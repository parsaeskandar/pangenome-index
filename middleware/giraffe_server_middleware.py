#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import threading
import time
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


FastqRead = Tuple[str, str, str]


@dataclass
class GiraffeAlignment:
    """One GAF record and the graph path (node id, reverse on node, from_length per Mapping)."""

    gaf_line: str
    graph_path: List[Tuple[int, bool, int]]


@dataclass
class GiraffeServerConfig:
    vg_binary: str
    gbz_path: str
    minimizer_path: str
    distance_path: str
    zipcode_path: str
    threads: int = 8
    max_multimaps: int = 1
    batch_size: int = 256
    output_timeout_s: float = 60.0
    extra_args: Optional[Sequence[str]] = None
    emit_graph_path: bool = False


class GiraffeServerMiddleware:
    """Long-lived middleware around `vg giraffe-server` framed output mode."""

    def __init__(self, cfg: GiraffeServerConfig) -> None:
        self.cfg = cfg
        self._proc: Optional[subprocess.Popen[str]] = None
        self._lock = threading.Lock()
        self._stderr_tail: deque[str] = deque(maxlen=200)
        self._stderr_lock = threading.Lock()
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
            "-M",
            str(self.cfg.max_multimaps),
            "-b",
            str(self.cfg.batch_size),
            "--framed-output",
        ]
        if self.cfg.emit_graph_path:
            cmd.append("--emit-graph-path")
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

    def map_reads(self, reads: Iterable[FastqRead]) -> List[List[str]]:
        if self.cfg.emit_graph_path:
            raise ValueError(
                "map_reads cannot parse output when emit_graph_path is enabled; use map_reads_with_graph_paths"
            )
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

            for name, seq, qual in batch:
                if not name or not seq:
                    raise ValueError("Each read needs non-empty name and sequence")
                if qual and len(seq) != len(qual):
                    raise ValueError(f"Read '{name}' has mismatched sequence/quality lengths")
                if qual:
                    proc.stdin.write(f"{name}\t{seq}\t{qual}\n")
                else:
                    proc.stdin.write(f"{name}\t{seq}\n")
            # Tell server to process current buffered reads immediately.
            proc.stdin.write("@@FLUSH@@\n")
            proc.stdin.flush()

            return self._read_framed_batch(proc, [r[0] for r in batch])

    def map_reads_with_graph_paths(self, reads: Iterable[FastqRead]) -> List[List[GiraffeAlignment]]:
        if not self.cfg.emit_graph_path:
            raise ValueError("map_reads_with_graph_paths requires GiraffeServerConfig.emit_graph_path=True")
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

            for name, seq, qual in batch:
                if not name or not seq:
                    raise ValueError("Each read needs non-empty name and sequence")
                if qual and len(seq) != len(qual):
                    raise ValueError(f"Read '{name}' has mismatched sequence/quality lengths")
                if qual:
                    proc.stdin.write(f"{name}\t{seq}\t{qual}\n")
                else:
                    proc.stdin.write(f"{name}\t{seq}\n")
            proc.stdin.write("@@FLUSH@@\n")
            proc.stdin.flush()

            return self._read_framed_batch_with_graph(proc, [r[0] for r in batch])

    def map_sequences(self, sequences: Iterable[str]) -> List[List[str]]:
        reads: List[FastqRead] = []
        for i, seq in enumerate(sequences):
            s = seq.strip().upper()
            if not s:
                continue
            reads.append((f"read_{i}", s, "I" * len(s)))
        return self.map_reads(reads)

    def map_sequences_with_graph_paths(self, sequences: Iterable[str]) -> List[List[GiraffeAlignment]]:
        reads: List[FastqRead] = []
        for i, seq in enumerate(sequences):
            s = seq.strip().upper()
            if not s:
                continue
            reads.append((f"read_{i}", s, "I" * len(s)))
        return self.map_reads_with_graph_paths(reads)

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
            if not header.startswith("@READ\t"):
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

    def _read_framed_batch_with_graph(
        self, proc: subprocess.Popen[str], names: List[str]
    ) -> List[List[GiraffeAlignment]]:
        deadline = time.monotonic() + self.cfg.output_timeout_s
        out: List[List[GiraffeAlignment]] = []
        for expected_name in names:
            header = self._readline_with_timeout(proc, deadline)
            if not header.startswith("@READ\t"):
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
            mapped: List[GiraffeAlignment] = []
            for _ in range(count):
                gaf_line = self._readline_with_timeout(proc, deadline)
                graph_line = self._readline_with_timeout(proc, deadline)
                mapped.append(GiraffeAlignment(gaf_line=gaf_line, graph_path=_parse_graph_line(graph_line)))
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

    def _drain_stderr(self, stderr_stream) -> None:
        try:
            for line in stderr_stream:
                with self._stderr_lock:
                    self._stderr_tail.append(line.rstrip("\n"))
        except Exception:
            pass

    def _format_stderr_tail(self) -> str:
        with self._stderr_lock:
            if not self._stderr_tail:
                return "(no stderr captured from giraffe-server)"
            return "giraffe-server stderr (tail):\n" + "\n".join(self._stderr_tail)


def _parse_graph_line(line: str) -> List[Tuple[int, bool, int]]:
    if not line.startswith("@GRAPH"):
        raise RuntimeError(f"Expected @GRAPH line after GAF, got: {line!r}")
    fields = line.split("\t")[1:]
    if not fields:
        return []
    if len(fields) % 3 != 0:
        raise RuntimeError(f"Malformed @GRAPH line (expected triples): {line!r}")
    out: List[Tuple[int, bool, int]] = []
    for i in range(0, len(fields), 3):
        node_id = int(fields[i])
        is_rev = bool(int(fields[i + 1]))
        from_len = int(fields[i + 2])
        out.append((node_id, is_rev, from_len))
    return out
