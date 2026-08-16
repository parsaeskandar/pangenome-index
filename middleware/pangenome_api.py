#!/usr/bin/env python3
"""
Hardened HTTP service in front of PangenomeMiddleware — the private backend the
public hgPangenome CGI proxies to.

    internet -> hgPangenome CGI (public, trust boundary) -> this middleware (private)

Endpoints:
  POST /api/v1/map            {sequences:[{name,sequence}], options:{...}}
      -> 200 {job_id, status:"queued", n_sequences}
      (?sync=1 blocks and returns the finished envelope)
  GET  /api/v1/map/{job_id}   -> 200 {status, progress, results, error}
  POST /api/v1/liftover       {src, start, end, tgt}   (synchronous)
      src = full contig path; tgt = haplotype or contig, or an array of them
      -> 200 {intervals:[{haplotype, start, end, strand}, ...]}   (0..N per target)
  POST /api/v1/liftover/targets {src, start, end}   (synchronous, names only)
      -> 200 {haplotypes:[...]}   haplotypes this source interval CAN translate to
  GET  /api/v1/haplotypes     -> 200 {haplotypes:[...]}   (2-field names)
  GET  /healthz               -> 200 {status, ready, ...load/metrics}   (no auth)

Defense-in-depth (the CGI also rate-limits/validates, but we don't trust it):
  - concurrency cap + bounded queue with 429 load-shedding
  - per-job wall-clock timeout
  - independent input validation (count / length / total bytes / charset)
  - shared-secret auth (X-Pangenome-Token or Authorization: Bearer)
  - job TTL + eviction, capped job store
  - structured logs (no raw sequences) + metrics + /healthz

Stdlib only. Concurrency: the middleware multiplexes a single `vg
giraffe-server` subprocess (one loaded index copy) — concurrent jobs' reads are
coalesced into shared mapping batches and demultiplexed back by name, so
`--max-concurrent N` runs N worker threads that genuinely overlap on mapping and
surjection without loading the indexes N times. The coordinate-index
anchor-building step holds the GIL, so that stage serializes across workers (it
is cheap relative to mapping); releasing the GIL in the pybind11 bindings is a
later optimization. The queue + 429 still provide backpressure.
"""
from __future__ import annotations

import argparse
import hmac
import json
import os
import queue
import sys
import threading
import time
import uuid
from collections import deque
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Deque, Dict, List, Optional, Tuple
from urllib.parse import urlparse, parse_qs

from .pangenome_middleware import PangenomeMiddleware

_GAF_TAG_START = 12
VALID_BASES = set("ACGTNacgtn")

# Liftover limits. The span cap mirrors liftover_ext's own MAX_INTERVAL_LENGTH
# (validated again in C++); the target cap bounds fan-out on the multi-target
# form (one translate() call per target).
MAX_LIFTOVER_SPAN = 10_000_000
MAX_LIFTOVER_TARGETS = 1024


# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #

@dataclass
class ApiConfig:
    # §3 input limits (MUST match the CGI's)
    max_sequences: int = 50
    max_seq_len: int = 100_000          # 100 kb — supports long reads
    max_total_bytes: int = 10_000_000   # 10 MB request body cap
    # §1 concurrency / queue
    max_concurrent: int = 1             # N worker threads; the middleware multiplexes the single engine (see note above)
    max_queued: int = 32
    # §2 / §6 timeouts and lifecycle
    job_timeout_s: float = 120.0
    job_ttl_s: float = 1800.0           # 30 min
    max_jobs: int = 1000
    # §4 auth (None = disabled, with a warning)
    auth_token: Optional[str] = None
    # engine
    default_max_multimaps: int = 1


# --------------------------------------------------------------------------- #
# Metrics (§8)
# --------------------------------------------------------------------------- #

class Metrics:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self.counters: Dict[str, int] = {
            "submitted": 0, "completed": 0, "errored": 0, "timed_out": 0,
            "rejected_429": 0, "rejected_400": 0, "rejected_401": 0, "rejected_503": 0,
        }
        self._latencies: Deque[float] = deque(maxlen=500)   # run seconds
        self._completions: Deque[float] = deque(maxlen=2000)  # timestamps

    def inc(self, key: str, n: int = 1) -> None:
        with self._lock:
            self.counters[key] = self.counters.get(key, 0) + n

    def record_run(self, seconds: float) -> None:
        with self._lock:
            self._latencies.append(seconds)
            self._completions.append(time.time())

    def _pct(self, samples: List[float], p: float) -> Optional[float]:
        if not samples:
            return None
        s = sorted(samples)
        k = min(len(s) - 1, int(round((p / 100.0) * (len(s) - 1))))
        return round(s[k], 3)

    def snapshot(self) -> Dict[str, Any]:
        with self._lock:
            lat = list(self._latencies)
            now = time.time()
            jpm = sum(1 for t in self._completions if now - t <= 60.0)
            return {
                "counters": dict(self.counters),
                "jobs_per_min": jpm,
                "latency_p50_s": self._pct(lat, 50),
                "latency_p95_s": self._pct(lat, 95),
            }


# --------------------------------------------------------------------------- #
# GAF tag parsing -> JSON contract
# --------------------------------------------------------------------------- #

def _parse_tags(tag_fields: List[str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for f in tag_fields:
        parts = f.split(":", 2)
        if len(parts) != 3:
            continue
        tag, typ, val = parts
        if typ == "i":
            try:
                out[tag] = int(val)
            except ValueError:
                out[tag] = val
        elif typ == "f":
            try:
                out[tag] = float(val)
            except ValueError:
                out[tag] = val
        else:
            out[tag] = val
    return out


def _haplotypes_from_tags(tags: Dict[str, Any]) -> Dict[str, Any]:
    names = tags["hp"].split(",") if tags.get("hp") else []
    rep = tags.get("hb")
    if rep in ("*", ""):
        rep = None
    mosaic: List[Dict[str, Any]] = []
    hm = tags.get("hm")
    if hm:
        for seg in hm.split(";"):
            f = seg.split(":", 2)
            if len(f) == 3:
                try:
                    mosaic.append({
                        "covered_bp": int(f[0]),
                        "haplotype_count": int(f[1]),
                        "representative": None if f[2] in ("*", "") else f[2],
                    })
                except ValueError:
                    pass
    return {
        "count": tags.get("hn", len(names)),
        "names": names,
        "representative": rep,
        "representative_is_reference": tags.get("hf", 0) == 1,
        "fully_covered": tags.get("hq", 0) == 1,
        "coverage_percent": tags.get("hv"),
        "num_segments": tags.get("hs"),
        "mosaic": mosaic,
    }


def _surjection_from_tags(tags: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    status = tags.get("sj")
    if status is None:
        return None
    surj: Dict[str, Any] = {"status": status}
    if status == "ok":
        surj.update({
            "target": tags.get("sn"),
            "position": tags.get("sp"),
            "strand": "-" if tags.get("sr", 0) == 1 else "+",
            "score": tags.get("ss"),
            "mapping_quality": tags.get("sm"),
            "cigar": tags.get("sc"),
        })
    return surj


# --------------------------------------------------------------------------- #
# Jobs + Service
# --------------------------------------------------------------------------- #

class NotReady(Exception):
    pass


class Busy(Exception):
    pass


@dataclass
class Job:
    job_id: str
    sequences: List[Dict[str, str]]
    options: Dict[str, Any]
    status: str = "queued"                 # queued | running | done | error
    completed: int = 0
    results: List[Dict[str, Any]] = field(default_factory=list)
    error: Optional[str] = None
    created_at: float = field(default_factory=time.time)
    started_at: Optional[float] = None
    finished_at: Optional[float] = None
    total_bp: int = 0
    done_evt: threading.Event = field(default_factory=threading.Event)

    def envelope(self) -> Dict[str, Any]:
        return {
            "job_id": self.job_id,
            "status": self.status,
            "progress": {"completed": self.completed, "total": len(self.sequences)},
            "results": self.results,
            "error": self.error,
        }


class Service:
    def __init__(self, cfg: ApiConfig, mw=None, stub: bool = False, metrics: Optional[Metrics] = None) -> None:
        self._cfg = cfg
        self._mw = mw
        self._stub = stub
        self._ready = stub or (mw is not None)
        self._metrics = metrics or Metrics()
        self._jobs: Dict[str, Job] = {}
        self._jobs_lock = threading.Lock()
        self._queue: "queue.Queue[Job]" = queue.Queue(maxsize=max(1, cfg.max_queued))
        self._running = 0
        for _ in range(max(1, cfg.max_concurrent)):
            threading.Thread(target=self._worker, daemon=True).start()
        threading.Thread(target=self._reaper, daemon=True).start()

    # ---- lifecycle ----
    def set_engine(self, mw) -> None:
        self._mw = mw
        self._ready = True

    @property
    def metrics(self) -> Metrics:
        return self._metrics

    def health(self) -> Dict[str, Any]:
        return {
            "status": "ok" if self._ready else "starting",
            "ready": self._ready,
            "concurrency": self._running,
            "max_concurrent": self._cfg.max_concurrent,
            "queue_depth": self._queue.qsize(),
            "max_queued": self._cfg.max_queued,
            "jobs_retained": len(self._jobs),
            **self._metrics.snapshot(),
        }

    # ---- submission ----
    def submit(self, sequences: List[Dict[str, str]], options: Dict[str, Any]) -> Job:
        if not self._ready:
            raise NotReady()
        job = Job(job_id=uuid.uuid4().hex[:12], sequences=sequences, options=options)
        job.total_bp = sum(len(s["sequence"]) for s in sequences)
        with self._jobs_lock:
            self._jobs[job.job_id] = job
        try:
            self._queue.put_nowait(job)
        except queue.Full:
            with self._jobs_lock:
                self._jobs.pop(job.job_id, None)
            raise Busy()
        self._metrics.inc("submitted")
        return job

    def get(self, job_id: str) -> Optional[Job]:
        with self._jobs_lock:
            return self._jobs.get(job_id)

    # ---- synchronous coordinate translation (no queue) ----
    def liftover(self, src: str, start: int, end: int,
                 tgts: List[str]) -> List[Dict[str, Any]]:
        """Fold source→target correspondences into target intervals. Runs inline
        in the request thread (translation is fast); the middleware serializes
        coordinate-index access internally."""
        if not self._ready:
            raise NotReady()
        if self._stub:
            return []
        return self._mw.translate_intervals(src, start, end, tgts)

    def haplotypes(self) -> List[str]:
        if not self._ready:
            raise NotReady()
        if self._stub:
            return []
        return self._mw.get_haplotype_names()

    def liftover_targets(self, src: str, start: int, end: int,
                         scored: bool = False, min_coverage: float = 0.0,
                         max_nodes: int = 0):
        """The target haplotypes a source interval can translate to.

        `scored` returns [{haplotype, coverage, covered_bp}] instead of bare
        names, so a picker can rank destinations by how much of the region each
        haplotype actually shares."""
        if not self._ready:
            raise NotReady()
        if self._stub:
            return []
        if scored:
            return self._mw.translatable_haplotypes_scored(
                src, start, end, min_coverage, max_nodes)
        return self._mw.translatable_haplotypes(src, start, end)

    # ---- worker ----
    def _worker(self) -> None:
        while True:
            job = self._queue.get()
            self._running += 1
            try:
                self._process(job)
            except Exception as exc:
                job.status = "error"
                job.error = f"internal error: {exc}"
            finally:
                self._running -= 1
                if job.status not in ("error",):
                    job.status = "done"
                job.finished_at = time.time()
                job.done_evt.set()
                self._log_job(job)
                if job.status == "done":
                    self._metrics.inc("completed")
                    if job.started_at:
                        self._metrics.record_run(job.finished_at - job.started_at)
                elif job.error == "timed out":
                    self._metrics.inc("timed_out")
                else:
                    self._metrics.inc("errored")

    def _process(self, job: Job) -> None:
        job.started_at = time.time()
        job.status = "running"
        deadline = job.started_at + self._cfg.job_timeout_s
        if self._stub:
            self._process_stub(job)
            return

        surject = bool(job.options.get("surject", True))
        explicit_target = job.options.get("surject_target") or None
        # Per-haplotype coverage scores: on by default, with an optional cutoff
        # so a caller that only wants strong matches doesn't carry all ~464.
        want_coverage = bool(job.options.get("haplotype_coverage", True))
        try:
            min_coverage = float(job.options.get("min_haplotype_coverage") or 0.0)
        except (TypeError, ValueError):
            min_coverage = 0.0
        include_zero = bool(job.options.get("include_zero_coverage", False))

        reads = [(s["name"], s["sequence"], "I" * len(s["sequence"])) for s in job.sequences]
        # Bound the giraffe call to the job's remaining budget.
        self._mw.set_call_timeout(max(1.0, deadline - time.time()))
        mapped = self._mw.map_reads(reads)  # 1:1 with reads

        results: List[Dict[str, Any]] = []
        for i, alignments in enumerate(mapped):
            if time.time() > deadline:
                job.status = "error"
                job.error = "timed out"
                job.results = results
                return
            seq = job.sequences[i]
            qlen = len(seq["sequence"])
            if not alignments:
                results.append({"name": seq["name"], "status": "unmapped",
                                "error": None, "query_length": qlen, "alignments": []})
            else:
                alns = []
                for j, gaf in enumerate(alignments):
                    self._mw.set_call_timeout(max(1.0, deadline - time.time()))
                    alns.append(self._build_alignment(
                        gaf, surject, explicit_target, primary=(j == 0),
                        coverage=want_coverage, min_coverage=min_coverage,
                        include_zero=include_zero))
                results.append({"name": seq["name"], "status": "mapped",
                                "error": None, "query_length": qlen, "alignments": alns})
            job.completed = i + 1
            job.results = results
        job.results = results

    def _build_alignment(self, gaf: str, surject: bool,
                         explicit_target: Optional[str], primary: bool,
                         coverage: bool = True,
                         min_coverage: float = 0.0,
                         include_zero: bool = False) -> Dict[str, Any]:
        cols = gaf.split("\t")
        tags = _parse_tags(cols[_GAF_TAG_START:]) if len(cols) > _GAF_TAG_START else {}
        haplotypes = _haplotypes_from_tags(tags)

        # Graded score for EVERY haplotype, not just the exact-path carriers in
        # `haplotypes`. Never fail the alignment over this: it is supplementary.
        haplotype_coverage: List[Dict[str, Any]] = []
        if coverage:
            try:
                haplotype_coverage = self._mw.haplotype_coverage(
                    gaf, min_coverage, include_zero)
            except Exception:
                haplotype_coverage = []

        surjection: Optional[Dict[str, Any]] = None
        if surject:
            target = explicit_target or haplotypes["representative"]
            if target:
                try:
                    anchors, path_len, status = self._mw.build_surject_anchors(gaf, target)
                    if anchors:
                        lines = self._mw.surject_with_anchors(
                            gaf, anchors, target, target_path_length=path_len)
                        if lines:
                            scols = lines[0].split("\t")
                            surjection = _surjection_from_tags(
                                _parse_tags(scols[_GAF_TAG_START:]) if len(scols) > _GAF_TAG_START else {})
                        else:
                            surjection = {"status": "surjection_failed"}
                    else:
                        surjection = {"status": "surjection_failed", "detail": f"no_anchors ({status})"}
                except Exception as exc:
                    surjection = {"status": "surjection_failed", "detail": str(exc)}
            else:
                surjection = None

        return {
            "primary": primary,
            "score": tags.get("AS"),
            "mapping_quality": int(cols[11]) if len(cols) > 11 and cols[11].lstrip("-").isdigit() else None,
            "strand": cols[4] if len(cols) > 4 else None,
            "graph_path": cols[5] if len(cols) > 5 else None,
            "haplotypes": haplotypes,
            "haplotype_coverage": haplotype_coverage,
            "surjection": surjection,
        }

    def _process_stub(self, job: Job) -> None:
        results: List[Dict[str, Any]] = []
        for i, s in enumerate(job.sequences):
            results.append(_canned_result(s["name"], s["sequence"], i))
            job.completed = i + 1
            job.results = results
        job.results = results

    # ---- lifecycle: TTL eviction (§6) ----
    def _reaper(self) -> None:
        while True:
            time.sleep(60.0)
            now = time.time()
            with self._jobs_lock:
                # Drop finished jobs past their TTL.
                for jid in [j for j, job in self._jobs.items()
                            if job.finished_at and now - job.finished_at > self._cfg.job_ttl_s]:
                    self._jobs.pop(jid, None)
                # Cap total retained: evict oldest finished first.
                if len(self._jobs) > self._cfg.max_jobs:
                    finished = sorted(
                        (job for job in self._jobs.values() if job.finished_at),
                        key=lambda j: j.finished_at or 0.0)
                    for job in finished[: len(self._jobs) - self._cfg.max_jobs]:
                        self._jobs.pop(job.job_id, None)

    # ---- structured logging (§8; no raw sequences) ----
    def _log_job(self, job: Job) -> None:
        queue_wait = (job.started_at - job.created_at) if job.started_at else None
        run = (job.finished_at - job.started_at) if (job.finished_at and job.started_at) else None
        sys.stderr.write(json.dumps({
            "ev": "job",
            "job_id": job.job_id,
            "n_sequences": len(job.sequences),
            "total_bp": job.total_bp,
            "queue_wait_s": round(queue_wait, 3) if queue_wait is not None else None,
            "run_s": round(run, 3) if run is not None else None,
            "outcome": job.error or job.status,
        }) + "\n")
        sys.stderr.flush()


def _canned_result(name: str, seq: str, idx: int) -> Dict[str, Any]:
    qlen = len(seq)
    if idx % 3 == 2:
        return {"name": name, "status": "unmapped", "error": None,
                "query_length": qlen, "alignments": []}
    mosaic = (idx % 3 == 1)
    names = ["CHM13#0#chr10", "GRCh38#0#chr10",
             "HG00097#1#CM094066.1", "HG02300#1#CM086438.1"]
    seg1 = round(qlen * 0.72)
    haplotypes = {
        "count": 6 if mosaic else 12, "names": names,
        "representative": "CHM13#0#chr10", "representative_is_reference": True,
        "fully_covered": not mosaic, "coverage_percent": 72 if mosaic else 100,
        "num_segments": 2 if mosaic else 1,
        "mosaic": ([{"covered_bp": seg1, "haplotype_count": 6, "representative": "CHM13#0#chr10"},
                    {"covered_bp": qlen - seg1, "haplotype_count": 9, "representative": "HG02300#1#CM086438.1"}]
                   if mosaic else []),
    }
    surjection = {"status": "ok", "target": "CHM13#0#chr10",
                  "position": 4338779 + idx * 137, "strand": "+",
                  "cigar": f"{qlen}M", "score": qlen * 29, "mapping_quality": 60}
    return {"name": name, "status": "mapped", "error": None, "query_length": qlen,
            "alignments": [{"primary": True, "score": qlen * 29, "mapping_quality": 60,
                            "strand": "+", "graph_path": ">1>2>3",
                            "haplotypes": haplotypes, "surjection": surjection}]}


# --------------------------------------------------------------------------- #
# Validation (§3)
# --------------------------------------------------------------------------- #

def validate(sequences, cfg: ApiConfig) -> Tuple[bool, str]:
    if not isinstance(sequences, list) or not sequences:
        return False, "no sequences provided"
    if len(sequences) > cfg.max_sequences:
        return False, f"too many sequences ({len(sequences)}); max is {cfg.max_sequences}"
    for i, s in enumerate(sequences):
        if not isinstance(s, dict) or not s.get("sequence"):
            return False, f"sequence {i} missing 'sequence'"
        seq = s["sequence"].strip()
        if not seq:
            return False, f"sequence {i} is empty"
        if len(seq) > cfg.max_seq_len:
            return False, f"sequence {i} too long ({len(seq)}); max is {cfg.max_seq_len}"
        bad = set(seq) - VALID_BASES
        if bad:
            return False, f"sequence {i} has invalid bases: {''.join(sorted(bad))}"
        s["sequence"] = seq.upper()
        if not s.get("name"):
            s["name"] = f"seq_{i + 1}"
    return True, ""


# --------------------------------------------------------------------------- #
# HTTP handler
# --------------------------------------------------------------------------- #

def make_handler(service: Service, cfg: ApiConfig):
    class Handler(BaseHTTPRequestHandler):
        server_version = "PangenomeAPI/1.0"

        def _send_json(self, code: int, obj: Dict[str, Any]) -> None:
            body = json.dumps(obj).encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            try:
                self.wfile.write(body)
            except BrokenPipeError:
                pass

        def _error(self, code: int, message: str) -> None:
            self._send_json(code, {"status": "error", "error": message})

        def _handle_liftover(self, payload: Dict[str, Any]) -> None:
            src = payload.get("src")
            tgt = payload.get("tgt")
            start = payload.get("start")
            end = payload.get("end")

            if not isinstance(src, str) or not src:
                service.metrics.inc("rejected_400")
                self._error(400, "'src' must be a non-empty string (full contig path)")
                return
            # bool is a subclass of int — reject it explicitly.
            if (not isinstance(start, int) or isinstance(start, bool) or
                    not isinstance(end, int) or isinstance(end, bool)):
                service.metrics.inc("rejected_400")
                self._error(400, "'start' and 'end' must be integers")
                return
            if start < 0 or end < 0 or start > end:
                service.metrics.inc("rejected_400")
                self._error(400, f"invalid interval [{start}, {end})")
                return
            if end - start > MAX_LIFTOVER_SPAN:
                service.metrics.inc("rejected_400")
                self._error(400, f"interval length {end - start} exceeds maximum of {MAX_LIFTOVER_SPAN}")
                return
            if isinstance(tgt, str):
                tgts = [tgt] if tgt else []
            elif isinstance(tgt, list):
                tgts = tgt
            else:
                service.metrics.inc("rejected_400")
                self._error(400, "'tgt' must be a haplotype string or an array of strings")
                return
            if not tgts or not all(isinstance(t, str) and t for t in tgts):
                service.metrics.inc("rejected_400")
                self._error(400, "'tgt' must name at least one non-empty haplotype/contig")
                return
            if len(tgts) > MAX_LIFTOVER_TARGETS:
                service.metrics.inc("rejected_400")
                self._error(400, f"too many targets ({len(tgts)}); max is {MAX_LIFTOVER_TARGETS}")
                return

            try:
                intervals = service.liftover(src, start, end, tgts)
            except NotReady:
                service.metrics.inc("rejected_503")
                self._error(503, "service starting; indexes not loaded yet")
                return
            except ValueError as exc:
                # liftover_ext raises std::invalid_argument (→ ValueError) for an
                # unknown source haplotype or a bad interval.
                service.metrics.inc("rejected_400")
                self._error(400, str(exc))
                return
            except Exception as exc:
                service.metrics.inc("errored")
                self._error(500, f"liftover failed: {exc}")
                return

            self._send_json(200, {"intervals": intervals})

        def _handle_liftover_targets(self, payload: Dict[str, Any]) -> None:
            src = payload.get("src")
            start = payload.get("start")
            end = payload.get("end")

            if not isinstance(src, str) or not src:
                service.metrics.inc("rejected_400")
                self._error(400, "'src' must be a non-empty string (full contig path)")
                return
            if (not isinstance(start, int) or isinstance(start, bool) or
                    not isinstance(end, int) or isinstance(end, bool)):
                service.metrics.inc("rejected_400")
                self._error(400, "'start' and 'end' must be integers")
                return
            if start < 0 or end < 0 or start > end:
                service.metrics.inc("rejected_400")
                self._error(400, f"invalid interval [{start}, {end})")
                return
            if end - start > MAX_LIFTOVER_SPAN:
                service.metrics.inc("rejected_400")
                self._error(400, f"interval length {end - start} exceeds maximum of {MAX_LIFTOVER_SPAN}")
                return

            scored = bool(payload.get("scored", False))
            try:
                min_cov = float(payload.get("min_coverage") or 0.0)
            except (TypeError, ValueError):
                min_cov = 0.0
            try:
                max_nodes = int(payload.get("max_nodes") or 0)
            except (TypeError, ValueError):
                max_nodes = 0
            try:
                haplotypes = service.liftover_targets(src, start, end,
                                                      scored, min_cov, max_nodes)
            except NotReady:
                service.metrics.inc("rejected_503")
                self._error(503, "service starting; indexes not loaded yet")
                return
            except ValueError as exc:
                service.metrics.inc("rejected_400")
                self._error(400, str(exc))
                return
            except Exception as exc:
                service.metrics.inc("errored")
                self._error(500, f"liftover-targets failed: {exc}")
                return

            self._send_json(200, {"haplotypes": haplotypes})

        def log_message(self, fmt, *args):  # keep default request logging quiet
            return

        def _authed(self) -> bool:
            if cfg.auth_token is None:
                return True  # auth disabled (warned at startup)
            tok = self.headers.get("X-Pangenome-Token")
            if not tok:
                auth = self.headers.get("Authorization", "")
                if auth.startswith("Bearer "):
                    tok = auth[len("Bearer "):]
            return bool(tok) and hmac.compare_digest(tok, cfg.auth_token)

        def do_GET(self) -> None:
            parsed = urlparse(self.path)
            if parsed.path.rstrip("/") == "/healthz":
                self._send_json(200, service.health())   # no auth on health
                return
            if not self._authed():
                service.metrics.inc("rejected_401")
                self._error(401, "missing or invalid token")
                return
            if parsed.path.rstrip("/") == "/api/v1/haplotypes":
                try:
                    names = service.haplotypes()
                except NotReady:
                    self._error(503, "service starting; indexes not loaded yet")
                    return
                self._send_json(200, {"haplotypes": names})
                return
            prefix = "/api/v1/map/"
            if not parsed.path.startswith(prefix):
                self._error(404, f"unknown path: {parsed.path}")
                return
            job = service.get(parsed.path[len(prefix):].strip("/"))
            if job is None:
                self._error(404, "unknown or expired job")
                return
            self._send_json(200, job.envelope())

        def do_POST(self) -> None:
            parsed = urlparse(self.path)
            path = parsed.path.rstrip("/")
            if path not in ("/api/v1/map", "/api/v1/liftover", "/api/v1/liftover/targets"):
                self._error(404, f"unknown path: {parsed.path}")
                return
            if not self._authed():
                service.metrics.inc("rejected_401")
                self._error(401, "missing or invalid token")
                return

            length = int(self.headers.get("Content-Length", 0) or 0)
            if length > cfg.max_total_bytes:
                service.metrics.inc("rejected_400")
                self._error(400, f"payload too large ({length} bytes); max is {cfg.max_total_bytes}")
                return
            try:
                raw = self.rfile.read(length) if length else b"{}"
                payload = json.loads(raw or b"{}")
            except Exception as exc:
                service.metrics.inc("rejected_400")
                self._error(400, f"invalid JSON body: {exc}")
                return

            if path == "/api/v1/liftover":
                self._handle_liftover(payload)
                return
            if path == "/api/v1/liftover/targets":
                self._handle_liftover_targets(payload)
                return

            sequences = payload.get("sequences") or []
            options = payload.get("options") or {}
            if options.get("max_multimaps") is None:
                options["max_multimaps"] = cfg.default_max_multimaps

            ok, msg = validate(sequences, cfg)
            if not ok:
                service.metrics.inc("rejected_400")
                self._error(400, msg)
                return

            try:
                job = service.submit(sequences, options)
            except NotReady:
                service.metrics.inc("rejected_503")
                self._error(503, "service starting; indexes not loaded yet")
                return
            except Busy:
                service.metrics.inc("rejected_429")
                self._error(429, "server busy, please retry")
                return

            if parse_qs(parsed.query).get("sync", ["0"])[0] in ("1", "true", "yes"):
                job.done_evt.wait(timeout=cfg.job_timeout_s + 5.0)
                self._send_json(200, job.envelope())
            else:
                self._send_json(200, {"job_id": job.job_id, "status": "queued",
                                      "n_sequences": len(sequences)})

    return Handler


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pangenome mapping HTTP API")
    p.add_argument("--vg"); p.add_argument("--gbz")
    p.add_argument("--minimizer", help="LONG-READ minimizer index")
    p.add_argument("--dist"); p.add_argument("--zipcodes", help="LONG-READ zipcodes")
    p.add_argument("--ri"); p.add_argument("--tags"); p.add_argument("--gbwt-ri")
    p.add_argument("--t1")
    p.add_argument("--t2", default="",
                   help="Translation Table 2 (OPTIONAL). Omit to run table-free: "
                        "translation then uses only Table 1 + the GBWT/tag array.")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--max-multimaps", type=int, default=1)
    p.add_argument("--host", default="127.0.0.1", help="bind address (default localhost; use a tunnel)")
    p.add_argument("--port", type=int, default=8791)
    p.add_argument("--stub", action="store_true", help="canned results, no indexes")
    # hardening knobs
    p.add_argument("--max-sequences", type=int, default=50)
    p.add_argument("--max-seq-len", type=int, default=100_000)
    p.add_argument("--max-total-bytes", type=int, default=10_000_000)
    p.add_argument("--max-concurrent", type=int, default=1)
    p.add_argument("--max-queued", type=int, default=32)
    p.add_argument("--job-timeout", type=float, default=120.0)
    p.add_argument("--job-ttl", type=float, default=1800.0)
    p.add_argument("--max-jobs", type=int, default=1000)
    return p.parse_args()


def main() -> int:
    args = _parse_args()
    cfg = ApiConfig(
        max_sequences=args.max_sequences, max_seq_len=args.max_seq_len,
        max_total_bytes=args.max_total_bytes, max_concurrent=args.max_concurrent,
        max_queued=args.max_queued, job_timeout_s=args.job_timeout,
        job_ttl_s=args.job_ttl, max_jobs=args.max_jobs,
        auth_token=os.environ.get("PANGENOME_API_TOKEN") or None,
        default_max_multimaps=args.max_multimaps,
    )
    if cfg.auth_token is None:
        print("WARNING: PANGENOME_API_TOKEN not set — auth is DISABLED.", file=sys.stderr)
    if cfg.max_concurrent > 1:
        print(f"INFO: {cfg.max_concurrent} worker threads; the middleware multiplexes the "
              "single engine, so jobs overlap on one loaded index copy.", file=sys.stderr)

    metrics = Metrics()
    service = Service(cfg, mw=None, stub=args.stub, metrics=metrics)

    # Start listening immediately; load indexes in the background so /healthz and
    # 503-until-ready work (the CGI can surface "service starting").
    if not args.stub:
        # NOTE: t2 is deliberately absent — Table 2 is optional. Omitting it
        # selects the table-free translation path. Keeping it here made the
        # process exit before binding the port, so the proxy answered with its
        # own 503 instead of our JSON one.
        required = ("vg", "gbz", "minimizer", "dist", "zipcodes", "ri", "tags", "gbwt_ri", "t1")
        missing = [f for f in required if not getattr(args, f)]
        if missing:
            print("ERROR: missing required args (or pass --stub): "
                  + ", ".join("--" + m.replace("_", "-") for m in missing), file=sys.stderr)
            return 2

        def _load():
            print("Loading indexes (this can take a while)…", file=sys.stderr)
            mw = PangenomeMiddleware.from_paths(
                vg_binary=args.vg, gbz=args.gbz, minimizer=args.minimizer,
                distance=args.dist, zipcodes=args.zipcodes, ri=args.ri, tags=args.tags,
                gbwt_ri=args.gbwt_ri, t1=args.t1, t2=args.t2,
                threads=args.threads, max_multimaps=args.max_multimaps,
                output_timeout_s=args.job_timeout,
            )
            mw.wait_until_ready()
            service.set_engine(mw)
            print("Indexes ready — accepting jobs.", file=sys.stderr)
        threading.Thread(target=_load, daemon=True).start()

    httpd = ThreadingHTTPServer((args.host, args.port), make_handler(service, cfg))
    print(f"Pangenome API listening on http://{args.host}:{args.port} "
          f"(auth={'on' if cfg.auth_token else 'OFF'}, stub={args.stub})", file=sys.stderr)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
