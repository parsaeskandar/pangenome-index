#!/usr/bin/env python3
"""Long-lived middleware around `vg giraffe-server` framed-output mode.

Concurrency model (single loaded index copy, many concurrent callers):

  * `vg giraffe-server` already maps a batch across its `-t N` worker threads
    from ONE loaded set of indexes, and its framed output is keyed by read
    name (`READ<TAB>name<TAB>count`). So we do NOT need multiple engine copies
    to serve concurrent users — we need to keep the single subprocess busy and
    route its output back to the right caller.

  * A single **dispatcher** thread is the sole writer to the subprocess stdin.
    It drains a shared work queue, coalescing reads submitted by any number of
    concurrent callers into shared `PROCESS_BATCH` batches (so the N mapping
    threads work on everyone's reads at once), and writing surject blocks
    atomically.

  * A single **demuxer** thread reads the subprocess stdout, parses each
    `READ` frame, and hands it to the request that owns that name.

  * Callers (`map_reads`, `surject_with_anchors`) assign globally-unique
    internal names, enqueue their work, and block only on *their own* frames
    arriving — never on a global lock held across the whole round trip. Public
    signatures and return values are unchanged from the previous serial
    implementation, so callers (PangenomeMiddleware / the HTTP service) don't
    change; they can now just call these methods from multiple threads.
"""
from __future__ import annotations

import itertools
import os
import subprocess
import sys
import threading
import time
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# Reads are accepted in either form:
#   (name, sequence, quality)
#   (name, sequence, quality, surjection_target)
# An empty string in either trailing field disables that feature for that read.
FastqRead = Tuple[str, ...]


@dataclass
class GiraffeServerConfig:
    vg_binary: str
    gbz_path: str
    minimizer_path: str
    distance_path: str
    zipcode_path: str
    threads: int = 8
    # 0 means "use giraffe-server's built-in BLAT default (100)". Any value > 0
    # is forwarded as `-M N` and overrides the BLAT default.
    max_multimaps: int = 0
    batch_size: int = 256
    output_timeout_s: float = 60.0
    extra_args: Optional[Sequence[str]] = None
    # Haplotype path names to pre-index for surjection (forwarded as
    # --surject-target). If empty, per-read surjection requests are reported as
    # "not_indexed" by the server; the anchor pipeline doesn't need these.
    surject_target_paths: Sequence[str] = ()


@dataclass(eq=False)   # identity-based hash/eq so instances can live in a set
class _Request:
    """One in-flight caller waiting for a set of named frames."""
    outstanding: set                      # unique names not yet received
    frames: Dict[str, List[str]] = field(default_factory=dict)
    event: threading.Event = field(default_factory=threading.Event)
    error: Optional[str] = None


class GiraffeServerMiddleware:
    """Multiplexed client for a single `vg giraffe-server` subprocess."""

    def __init__(self, cfg: GiraffeServerConfig) -> None:
        self.cfg = cfg
        self._proc: Optional[subprocess.Popen[str]] = None
        self._start_lock = threading.Lock()

        # Unique-name generator for demux (never reused within a process).
        self._ids = itertools.count(1)

        # Routing table + in-flight set, guarded by _route_lock.
        self._route_lock = threading.Lock()
        self._sinks: Dict[str, _Request] = {}       # unique name -> request
        self._live: set = set()                     # requests awaiting frames
        self._closed_msg: Optional[str] = None      # set once the subprocess dies

        # Write queue drained by the dispatcher thread. Items are
        # ("read", line) or ("surject", block); lines/blocks already end in \n.
        self._write_q: deque = deque()
        self._write_cond = threading.Condition()
        self._shutdown = False

        # Threads.
        self._dispatch_thread: Optional[threading.Thread] = None
        self._stdout_thread: Optional[threading.Thread] = None
        self._stderr_thread: Optional[threading.Thread] = None

        # stderr rolling buffer (surfaced on errors).
        self._stderr_tail: deque = deque(maxlen=200)
        self._stderr_lock = threading.Lock()

    # ------------------------------------------------------------------ #
    # Lifecycle
    # ------------------------------------------------------------------ #

    def start(self) -> None:
        with self._start_lock:
            if self._proc is not None:
                return
            for p in (
                self.cfg.vg_binary, self.cfg.gbz_path, self.cfg.minimizer_path,
                self.cfg.distance_path, self.cfg.zipcode_path,
            ):
                if not Path(p).exists():
                    raise FileNotFoundError(f"Required path does not exist: {p}")

            cmd: List[str] = [
                self.cfg.vg_binary, "giraffe-server",
                "-Z", self.cfg.gbz_path,
                "-m", self.cfg.minimizer_path,
                "-d", self.cfg.distance_path,
                "-z", self.cfg.zipcode_path,
                "-t", str(self.cfg.threads),
                "-b", str(self.cfg.batch_size),
                "--framed-output",
            ]
            if self.cfg.max_multimaps > 0:
                cmd.extend(["-M", str(self.cfg.max_multimaps)])
            for target in self.cfg.surject_target_paths:
                cmd.extend(["--surject-target", target])
            if self.cfg.extra_args:
                cmd.extend(self.cfg.extra_args)

            self._proc = subprocess.Popen(
                cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True, bufsize=1,
            )
            self._shutdown = False
            self._closed_msg = None

            self._dispatch_thread = threading.Thread(
                target=self._dispatch_loop, daemon=True)
            self._dispatch_thread.start()
            if self._proc.stdout is not None:
                self._stdout_thread = threading.Thread(
                    target=self._demux_stdout, args=(self._proc.stdout,), daemon=True)
                self._stdout_thread.start()
            if self._proc.stderr is not None:
                self._stderr_thread = threading.Thread(
                    target=self._drain_stderr, args=(self._proc.stderr,), daemon=True)
                self._stderr_thread.start()

    def stop(self) -> None:
        proc = self._proc
        if proc is None:
            return
        self._proc = None
        with self._write_cond:
            self._shutdown = True
            self._write_cond.notify_all()
        self._fail_all("giraffe-server middleware stopped.")
        try:
            if proc.stdin:
                proc.stdin.close()
            proc.terminate()
            proc.wait(timeout=5)
        except Exception:
            proc.kill()
        finally:
            for stream in (proc.stdout, proc.stderr):
                if stream:
                    try:
                        stream.close()
                    except Exception:
                        pass

    @property
    def pid(self) -> Optional[int]:
        return None if self._proc is None else self._proc.pid

    def is_running(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    # ------------------------------------------------------------------ #
    # Public request API (safe to call from many threads concurrently)
    # ------------------------------------------------------------------ #

    def map_reads(
        self,
        reads: Iterable[FastqRead],
        surject_target: Optional[str] = None,
    ) -> List[List[str]]:
        """Map a batch of reads; return one list of GAF lines per input read,
        in input order. Each read is (name, seq, qual) or (name, seq, qual,
        target); `surject_target` supplies a default target for reads without
        one. Thread-safe: concurrent calls are coalesced onto the one
        subprocess and demultiplexed back by name."""
        if self._proc is None:
            self.start()

        batch = list(reads)
        if not batch:
            return []

        rid = next(self._ids)
        write_items: List[Tuple[str, str]] = []
        unique: List[str] = []
        for i, entry in enumerate(batch):
            if len(entry) == 3:
                name, seq, qual = entry
                target = surject_target or ""
            elif len(entry) == 4:
                name, seq, qual, target = entry
                if not target and surject_target:
                    target = surject_target
            else:
                raise ValueError(
                    "Each read must be (name, seq, qual) or (name, seq, qual, target)")
            if not name or not seq:
                raise ValueError("Each read needs non-empty name and sequence")
            if qual and len(seq) != len(qual):
                raise ValueError(f"Read '{name}' has mismatched sequence/quality lengths")

            u = f"m{rid}_{i}"               # unique framing name for demux
            if target:
                line = f"{u}\t{seq}\t{qual}\t{target}\n"
            elif qual:
                line = f"{u}\t{seq}\t{qual}\n"
            else:
                line = f"{u}\t{seq}\n"
            write_items.append(("read", line))
            unique.append(u)

        frames = self._submit(write_items, unique, self.cfg.output_timeout_s)
        return [frames[u] for u in unique]

    def map_sequences(
        self,
        sequences: Iterable[str],
        surject_target: Optional[str] = None,
    ) -> List[List[str]]:
        reads: List[FastqRead] = []
        for i, seq in enumerate(sequences):
            s = seq.strip().upper()
            if not s:
                continue
            reads.append((f"read_{i}", s, "I" * len(s)))
        return self.map_reads(reads, surject_target=surject_target)

    def surject_with_anchors(
        self,
        graph_alignment_gaf: str,
        anchors,
        target_haplotype: str,
        target_path_length: int = 0,
        read_name: Optional[str] = None,
    ) -> List[str]:
        """Surject a graph alignment onto `target_haplotype` using pre-computed
        anchors (SURJECT_WITH_ANCHORS). Returns the surjected GAF line(s).
        Thread-safe. The frame is keyed by a unique internal name, so the GAF
        content keeps its own read name; `read_name` is accepted for API
        compatibility but the wire framing uses a unique id."""
        if self._proc is None:
            self.start()

        anchor_list = list(anchors) if anchors is not None else []
        n_anchors = len(anchor_list)
        path_len = int(target_path_length) if target_path_length else 0

        rid = next(self._ids)
        u = f"s{rid}"                       # unique framing name for demux

        parts = [
            f"SURJECT_WITH_ANCHORS\t{u}\t{target_haplotype}\t{path_len}\t{n_anchors}\n"
        ]
        for a in anchor_list:
            parts.append(
                f"{a.gbwt_edge_begin_node}\t{a.gbwt_edge_begin_offset}\t"
                f"{a.gbwt_edge_end_node}\t{a.gbwt_edge_end_offset}\t"
                f"{a.path_offset_step_begin}\t{a.path_offset_step_end}\t"
                f"{a.read_begin_offset}\t{a.read_end_offset}\t"
                f"{a.source_mapping_begin}\t{a.source_mapping_end}\n"
            )
        gaf = graph_alignment_gaf if graph_alignment_gaf.endswith("\n") \
            else graph_alignment_gaf + "\n"
        parts.append(gaf)
        block = "".join(parts)

        frames = self._submit([("surject", block)], [u], self.cfg.output_timeout_s)
        return frames.get(u, [])

    def wait_until_ready(self) -> None:
        """Block until giraffe-server answers a probe read — the only reliable
        signal that all indexes have finished loading."""
        saved = self.cfg.output_timeout_s
        self.cfg.output_timeout_s = float("inf")     # wait forever for the probe
        try:
            self.map_reads([("__probe__", "ACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIII")])
        finally:
            self.cfg.output_timeout_s = saved

    # ------------------------------------------------------------------ #
    # Internals: submit / dispatch / demux
    # ------------------------------------------------------------------ #

    def _submit(
        self,
        write_items: List[Tuple[str, str]],
        unique_names: List[str],
        timeout_s: float,
    ) -> Dict[str, List[str]]:
        """Register the request, enqueue its work, and wait for all its frames."""
        req = _Request(outstanding=set(unique_names))
        with self._route_lock:
            if self._closed_msg is not None:
                raise RuntimeError(self._closed_msg + "\n" + self._format_stderr_tail())
            if self._proc is None or self._proc.poll() is not None:
                raise RuntimeError(
                    "giraffe-server is not running.\n" + self._format_stderr_tail())
            for u in unique_names:
                self._sinks[u] = req
            self._live.add(req)

        with self._write_cond:
            self._write_q.extend(write_items)
            self._write_cond.notify()

        timeout = None if timeout_s == float("inf") else timeout_s
        got = req.event.wait(timeout=timeout)

        with self._route_lock:
            self._live.discard(req)
            if not got:
                for u in unique_names:
                    self._sinks.pop(u, None)

        if req.error is not None:
            raise RuntimeError(req.error + "\n" + self._format_stderr_tail())
        if not got:
            raise RuntimeError(
                "Timed out waiting for giraffe-server output.\n" + self._format_stderr_tail())
        return req.frames

    def _dispatch_loop(self) -> None:
        while True:
            with self._write_cond:
                while not self._write_q and not self._shutdown:
                    self._write_cond.wait()
                if self._shutdown and not self._write_q:
                    return
                items = list(self._write_q)
                self._write_q.clear()

            proc = self._proc
            if proc is None or proc.stdin is None:
                self._fail_all("giraffe-server stdin unavailable.")
                return
            try:
                pending_reads = False
                for kind, payload in items:
                    if kind == "read":
                        proc.stdin.write(payload)
                        pending_reads = True
                    else:  # "surject": self-contained block; flush reads first
                        if pending_reads:
                            proc.stdin.write("PROCESS_BATCH\n")
                            pending_reads = False
                        proc.stdin.write(payload)
                if pending_reads:
                    proc.stdin.write("PROCESS_BATCH\n")
                proc.stdin.flush()
            except Exception as exc:  # broken pipe = subprocess died
                self._fail_all(f"giraffe-server stdin write failed: {exc}")
                return

    def _demux_stdout(self, stream) -> None:
        try:
            line = stream.readline()
            while line:
                header = line.rstrip("\n")
                if header.startswith("READ\t"):
                    parts = header.split("\t")
                    if len(parts) == 3:
                        name = parts[1]
                        try:
                            count = int(parts[2])
                        except ValueError:
                            count = 0
                        mapped: List[str] = []
                        broken = False
                        for _ in range(count):
                            body = stream.readline()
                            if not body:
                                broken = True
                                break
                            mapped.append(body.rstrip("\n"))
                        if broken:
                            break
                        self._deliver(name, mapped)
                    # malformed header: ignore and resync on the next READ line
                # non-READ stray line: ignore
                line = stream.readline()
        except Exception:
            pass
        finally:
            self._fail_all("giraffe-server closed stdout unexpectedly.")

    def _deliver(self, name: str, lines: List[str]) -> None:
        with self._route_lock:
            req = self._sinks.pop(name, None)
            if req is None:
                return                       # unknown/late frame — drop
            req.frames[name] = lines
            req.outstanding.discard(name)
            done = not req.outstanding
        if done:
            req.event.set()

    def _fail_all(self, msg: str) -> None:
        with self._route_lock:
            if self._closed_msg is None:
                self._closed_msg = msg
            reqs = list(self._live)
            self._live.clear()
            self._sinks.clear()
        for req in reqs:
            if req.error is None:
                req.error = msg
            req.event.set()

    def _drain_stderr(self, stderr_stream) -> None:
        # Rolling buffer surfaced on errors; also echo timing lines (and
        # everything when GIRAFFE_SERVER_STDERR is set) so server diagnostics
        # are visible.
        echo_all = os.environ.get("GIRAFFE_SERVER_STDERR") not in (None, "", "0")
        try:
            for line in stderr_stream:
                text = line.rstrip("\n")
                with self._stderr_lock:
                    self._stderr_tail.append(text)
                if echo_all or "[surject" in text:
                    print(text, file=sys.stderr, flush=True)
        except Exception:
            pass

    def _format_stderr_tail(self) -> str:
        with self._stderr_lock:
            if not self._stderr_tail:
                return "(no stderr captured from giraffe-server)"
            return "giraffe-server stderr (tail):\n" + "\n".join(self._stderr_tail)
