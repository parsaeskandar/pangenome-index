#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import threading
import time
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


# Reads are accepted in either form:
#   (name, sequence, quality)
#   (name, sequence, quality, surjection_target)
# An empty string in either of the trailing fields disables that feature for
# that read (e.g. ("r1", "ACGT", "", "HG002#1#chr1") = no quality, surject).
FastqRead = Tuple[str, ...]


@dataclass
class GiraffeServerConfig:
    vg_binary: str
    gbz_path: str
    minimizer_path: str
    distance_path: str
    zipcode_path: str
    threads: int = 8
    # 0 means "use giraffe-server's built-in BLAT default (100)". Any value
    # > 0 is forwarded as `-M N` and overrides the BLAT default.
    max_multimaps: int = 0
    batch_size: int = 256
    output_timeout_s: float = 60.0
    extra_args: Optional[Sequence[str]] = None
    # Haplotype path names to pre-index for surjection. Each name is forwarded
    # to giraffe-server as `--surject-target NAME`. If empty, per-read
    # surjection requests will be reported as "not_indexed" by the server.
    surject_target_paths: Sequence[str] = ()


class GiraffeServerMiddleware:
    """Long-lived middleware around `vg giraffe-server` framed output mode."""

    def __init__(self, cfg: GiraffeServerConfig) -> None:
        self.cfg = cfg
        self._proc: Optional[subprocess.Popen[str]] = None
        self._lock = threading.Lock()
        self._stderr_tail: deque[str] = deque(maxlen=200)
        self._stderr_cond = threading.Condition()  # notified on every new stderr line
        self._stderr_thread: Optional[threading.Thread] = None
        self._stdout_lines: deque[str] = deque()
        self._stdout_cond = threading.Condition()
        self._stdout_closed = False
        self._stdout_thread: Optional[threading.Thread] = None

    def start(self) -> None:
        if self._proc is not None:
            return
        for p in (
            self.cfg.vg_binary,
            self.cfg.gbz_path,
            self.cfg.minimizer_path,
            self.cfg.distance_path,
            self.cfg.zipcode_path,
        ):
            if not Path(p).exists():
                raise FileNotFoundError(f"Required path does not exist: {p}")

        cmd: List[str] = [
            self.cfg.vg_binary,
            "giraffe-server",
            "-Z",
            self.cfg.gbz_path,
            "-m",
            self.cfg.minimizer_path,
            "-d",
            self.cfg.distance_path,
            "-z",
            self.cfg.zipcode_path,
            "-t",
            str(self.cfg.threads),
            "-b",
            str(self.cfg.batch_size),
            "--framed-output",
        ]
        # Only forward -M when explicitly set; 0 means "use the engine's
        # BLAT default (100)".
        if self.cfg.max_multimaps > 0:
            cmd.extend(["-M", str(self.cfg.max_multimaps)])
        for target in self.cfg.surject_target_paths:
            cmd.extend(["--surject-target", target])
        if self.cfg.extra_args:
            cmd.extend(self.cfg.extra_args)

        self._proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        if self._proc.stdout is not None:
            self._stdout_closed = False
            self._stdout_thread = threading.Thread(
                target=self._drain_stdout,
                args=(self._proc.stdout,),
                daemon=True,
            )
            self._stdout_thread.start()
        if self._proc.stderr is not None:
            self._stderr_thread = threading.Thread(
                target=self._drain_stderr,
                args=(self._proc.stderr,),
                daemon=True,
            )
            self._stderr_thread.start()

    def stop(self) -> None:
        if self._proc is None:
            return
        proc = self._proc
        self._proc = None
        try:
            if proc.stdin:
                proc.stdin.close()
            proc.terminate()
            proc.wait(timeout=5)
        except Exception:
            proc.kill()
        finally:
            if proc.stdout:
                try:
                    proc.stdout.close()
                except Exception:
                    pass
            if proc.stderr:
                try:
                    proc.stderr.close()
                except Exception:
                    pass

    def map_reads(
        self,
        reads: Iterable[FastqRead],
        surject_target: Optional[str] = None,
    ) -> List[List[str]]:
        """Send a batch of reads to giraffe-server and return one list of GAF
        lines per input read (preserving order).

        Each read may be a 3-tuple (name, sequence, quality) or a 4-tuple
        (name, sequence, quality, surjection_target). The optional
        `surject_target` argument applies a default target to every read in
        the batch that doesn't have one explicitly. Empty string means "no
        surjection".
        """
        if self._proc is None:
            self.start()
        assert self._proc is not None
        proc = self._proc

        batch = list(reads)
        if not batch:
            return []

        with self._lock:
            if proc.stdin is None or proc.stdout is None:
                raise RuntimeError("giraffe-server process streams are unavailable")
            if proc.poll() is not None:
                raise RuntimeError(
                    "giraffe-server exited before request.\n" + self._format_stderr_tail()
                )

            names: List[str] = []
            for entry in batch:
                if len(entry) == 3:
                    name, seq, qual = entry
                    target = surject_target or ""
                elif len(entry) == 4:
                    name, seq, qual, target = entry
                    if not target and surject_target:
                        target = surject_target
                else:
                    raise ValueError(
                        "Each read must be (name, seq, qual) or (name, seq, qual, target)"
                    )
                if not name or not seq:
                    raise ValueError("Each read needs non-empty name and sequence")
                if qual and len(seq) != len(qual):
                    raise ValueError(f"Read '{name}' has mismatched sequence/quality lengths")

                if target:
                    # Always emit four fields when a surjection target is set;
                    # an empty quality column is fine on the server side.
                    proc.stdin.write(f"{name}\t{seq}\t{qual}\t{target}\n")
                elif qual:
                    proc.stdin.write(f"{name}\t{seq}\t{qual}\n")
                else:
                    proc.stdin.write(f"{name}\t{seq}\n")
                names.append(name)
            # Tell server to process current buffered reads immediately.
            proc.stdin.write("FLUSH_NOW\n")
            proc.stdin.flush()

            return self._read_framed_batch(proc, names)

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

    @property
    def pid(self) -> Optional[int]:
        if self._proc is None:
            return None
        return self._proc.pid

    def is_running(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    def _read_framed_batch(self, proc: subprocess.Popen[str], names: List[str]) -> List[List[str]]:
        deadline = time.monotonic() + self.cfg.output_timeout_s
        out: List[List[str]] = []
        for expected_name in names:
            header = self._readline_with_timeout(proc, deadline)
            if not header.startswith("READ\t"):
                raise RuntimeError(
                    f"Unexpected framed header: {header!r}\n" + self._format_stderr_tail()
                )
            parts = header.split("\t")
            if len(parts) != 3:
                raise RuntimeError(f"Malformed frame header: {header!r}")
            read_name = parts[1]
            if read_name != expected_name:
                raise RuntimeError(
                    f"Framed output out of order: expected {expected_name!r}, got {read_name!r}"
                )
            count = int(parts[2])
            mapped: List[str] = []
            for _ in range(count):
                mapped.append(self._readline_with_timeout(proc, deadline))
            out.append(mapped)
        return out

    def _readline_with_timeout(self, proc: subprocess.Popen[str], deadline: float) -> str:
        while True:
            with self._stdout_cond:
                if self._stdout_lines:
                    return self._stdout_lines.popleft()
                if self._stdout_closed:
                    raise RuntimeError(
                        "giraffe-server closed stdout unexpectedly.\n" + self._format_stderr_tail()
                    )
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise RuntimeError(
                        "Timed out waiting for giraffe-server output.\n" + self._format_stderr_tail()
                    )
                self._stdout_cond.wait(timeout=remaining)

    def _drain_stdout(self, stdout_stream) -> None:
        try:
            for line in stdout_stream:
                with self._stdout_cond:
                    self._stdout_lines.append(line.rstrip("\n"))
                    self._stdout_cond.notify()
        except Exception:
            pass
        finally:
            with self._stdout_cond:
                self._stdout_closed = True
                self._stdout_cond.notify_all()

    def wait_until_ready(self) -> None:
        """Block indefinitely until giraffe-server responds to a probe read.
        This is the only reliable signal that all indexes have finished loading."""
        saved_timeout = self.cfg.output_timeout_s
        # float("inf") as deadline means _readline_with_timeout waits forever.
        self.cfg.output_timeout_s = float("inf")
        try:
            self.map_reads([("__probe__", "ACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIII")])
        finally:
            self.cfg.output_timeout_s = saved_timeout

    def _drain_stderr(self, stderr_stream) -> None:
        try:
            for line in stderr_stream:
                with self._stderr_cond:
                    self._stderr_tail.append(line.rstrip("\n"))
                    self._stderr_cond.notify_all()
        except Exception:
            pass

    def _format_stderr_tail(self) -> str:
        with self._stderr_cond:
            if not self._stderr_tail:
                return "(no stderr captured from giraffe-server)"
            return "giraffe-server stderr (tail):\n" + "\n".join(self._stderr_tail)
