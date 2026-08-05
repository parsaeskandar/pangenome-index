# Surject-with-Anchors — Implementation Plan

> **Status**: draft for review with Adam.
> **Owner**: Parsa.
> **Goal**: enable surjection onto *any* haplotype in a pangenome without pre-indexing every haplotype's positions at server startup.

---

## 1. Motivation

### 1.1 The problem we're solving

The current vg `Surjector` requires a `PathPositionHandleGraph` (typically built via `bdsg::ReferencePathOverlay`) that has been pre-indexed for every target haplotype the user might want to surject onto. For our pangenome use case (HPRC-v2 has ~94 K named paths), pre-indexing every haplotype at server startup is prohibitive:

- **RAM**: ~40–80 GB just for the overlay, on top of the GBZ (~20–40 GB) and other indexes.
- **Load time**: 30–90 minutes of single-threaded indexing at startup.

So the user has to pick a small set of pre-indexed targets up front. That's fine for fixed reference paths (GRCh38, CHM13) but rules out interactive "surject this alignment onto haplotype HG002#1#chr5" workflows.

### 1.2 The proposed solution

Use the **coordinate-translation pipeline** (already implemented in `pangenome-index`) to produce surjection **anchors** for any requested target haplotype on demand. The Surjector then runs against those anchors instead of querying a pre-built position graph.

Adam's framing: *"Modify Surjector so that it can be used with user-provided anchors and a PathHandleGraph. Then you would get the anchors from coordinate translation."*

### 1.3 Non-goals

- **No replacement of the existing surject pipeline.** `vg surject` CLI and the existing `Surjector::surject(...)` calls keep working bit-for-bit.
- **No spliced/multipath surjection** in this round. The new entry point handles linear (`realigning_surject`) only. Spliced can be added later if needed.
- **No removal of `HaplotypeSurjector`** from the engine in this round. It stays available for callers who pre-index a fixed set of targets (faster per-request).
- **No support for arbitrary `PathHandleGraph` implementations**. We assume `gbwtgraph::GBWTGraph` everywhere.

---

## 2. Architecture overview

Three long-lived processes, orchestrated by the Python middleware. Always-open indexes, fast per-request.

```
            ┌─────────────────────────────┐
            │     middleware (Python)      │
            │   orchestrates the pipeline  │
            └──────────────┬───────────────┘
                           │
        ┌──────────────────┼──────────────────┐
        │                  │                  │
        ▼                  │                  ▼
 ┌────────────────┐        │         ┌──────────────────────┐
 │ giraffe-server │        │         │ coord-trans-server   │
 │ (always open)  │        │         │ (always open)        │
 │ • GBZ          │        │         │ • GBZ                │
 │ • mapper       │        │         │ • FastLocate         │
 │ • surjector    │        │         │ • r-index            │
 │   ─ map()      │        │         │ • tag array          │
 │   ─ surject_   │        │         │ • translation tables │
 │     with_      │        │         │ • build_anchors()    │
 │     anchors()  │        │         │                      │
 └────────────────┘        │         └──────────────────────┘
        ▲                  │                  ▲
        │                  │                  │
        │       1. read    │   2. alignment   │
        │     ────────────►│ ────────────────►│
        │                  │                  │
        │       4. linear  │   3. anchors     │
        │     ◄────────────│ ◄────────────────│
        │     surjected    │                  │
        │     alignment    │                  │
        └─ 5. final ◄──────┘                  │
```

| Step | From → To | Payload |
|---|---|---|
| 1 | middleware → giraffe-server | read sequence (+ name, quality, target haplotype name) |
| 2 | giraffe-server → middleware | `Alignment` (graph path) |
| 3 | middleware → coord-trans-server | serialized `Alignment` + target haplotype name |
| 4 | coord-trans-server → middleware | serialized anchors + target_path_length + rev_strand |
| 5 | middleware → giraffe-server | serialized `Alignment` + anchors + target metadata |
| 6 | giraffe-server → middleware | surjected linear `Alignment` (GAF line) |

**Why three processes:**

1. Each process keeps only the indexes it needs (no FastLocate in vg-land, no MinimizerMapper in coord-trans-land).
2. Both servers stay open across many requests — startup cost amortized.
3. Restart/upgrade either side independently.
4. Per-request IPC overhead is ~1–3 ms, dwarfed by mapping (~10 ms) and surjection DP (~5–20 ms).

---

## 3. What an "anchor" is

The Surjector's anchor concept (`path_chunk_t` paired with a step range) tells the surjector: *"between read offsets `[r0, r1)`, the source graph alignment follows the target path through step range `[s_begin, s_end]`, which corresponds to base positions `[p_begin, p_end)` on the target path."*

Between consecutive anchors, the surjector does WFA / SW DP to fill in gaps where the source diverged from the target.

### 3.1 The new `PrecomputedAnchor` type

```cpp
namespace vg {
struct PrecomputedAnchor {
    // Read range this anchor pins.
    std::string::const_iterator read_begin;
    std::string::const_iterator read_end;

    // Graph alignment over that read range (same as path_chunk_t.second).
    path_t graph_path;

    // Step handles on the target path (begin/end of the chunk).
    step_handle_t step_begin;
    step_handle_t step_end;

    // Base positions on the target path for this chunk.
    size_t path_offset_begin;   // inclusive
    size_t path_offset_end;     // exclusive
};
}  // namespace vg
```

Free struct in `namespace vg`, declared in `surjector.hpp`. The coordinate translator produces these; the Surjector consumes them.

### 3.2 What anchor information replaces in surject

Walking through a 67-bp read and counting the original Surjector's position-graph calls (see `docs/surject_anchor_walkthrough.md` for the full trace), every `get_position_of_step`, `get_step_at_position`, and `get_path_length` call is answered by the anchor fields above. Total replacement: 9 position-graph calls per request → 0.

The only structural piece that still needs work is `SubpathOverlay`'s constructor typing (handled in §4.1).

---

## 4. Implementation plan, by component

### 4.1 Surjector side — vg / Giraffe_server

#### 4.1.1 New types and methods (`surjector.hpp`)

```cpp
// Free struct in namespace vg, before class Surjector.
struct PrecomputedAnchor { ... };  // as above

class Surjector : public AlignerClient {
public:
    // Existing:
    Surjector(const PathPositionHandleGraph* graph);
    vector<Alignment> surject(...);          // unchanged
    // ... existing methods ...

    /// NEW: surject against pre-computed anchors. Does NOT call
    /// extract_overlapping_paths, does NOT require PathPositionHandleGraph
    /// to have indexed the target path.
    ///
    /// Caller responsibilities:
    /// - target_path must exist in the underlying PathHandleGraph.
    /// - anchors must be sorted by read offset (begin ascending).
    /// - target_path_length is the total base length of the target path.
    /// - rev_strand: true if the source alignment runs along the target
    ///   path in reverse-complement (caller-determined).
    vector<Alignment> surject_with_anchors(
        const Alignment& source,
        path_handle_t target_path,
        bool rev_strand,
        size_t target_path_length,
        const std::vector<PrecomputedAnchor>& anchors,
        std::vector<std::tuple<std::string, int64_t, bool>>& positions_out,
        bool allow_negative_scores = false,
        bool preserve_deletions = false) const;

protected:
    // Existing realigning_surject / spliced_surject / etc., unchanged.

    /// NEW: sibling of realigning_surject that takes pre-computed per-chunk
    /// path intervals + path_length. Does not touch PathPositionHandleGraph.
    vector<pair<Alignment, pair<step_handle_t, step_handle_t>>>
    realigning_surject_anchored(
        const PathHandleGraph* graph,
        const Alignment& source,
        const path_handle_t& path_handle,
        bool rev_strand,
        size_t path_length,
        const std::vector<path_chunk_t>& path_chunks,
        const std::vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
        const std::vector<pair<size_t, size_t>>& chunk_path_intervals,
        bool allow_negative_scores,
        bool preserve_N_alignments = false) const;

    /// NEW: sibling of compute_disjoint_path_intervals, anchor-based.
    vector<tuple<size_t, size_t, vector<size_t>>>
    compute_disjoint_path_intervals_anchored(
        const Alignment& source,
        const std::vector<path_chunk_t>& path_chunks,
        const std::vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
        const std::vector<pair<size_t, size_t>>& chunk_path_intervals,
        size_t path_length,
        bool no_left_expansion,
        bool no_right_expansion,
        size_t max_gap) const;
};
```

#### 4.1.2 New implementation in `surjector.cpp`

- **`Surjector::surject_with_anchors`** — public entry. Converts `vector<PrecomputedAnchor>` → `(vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>, vector<pair<size_t, size_t>>)` then calls `realigning_surject_anchored`. Wraps result, populates `positions_out` using `anchors[0].path_offset_begin` directly (no `set_path_position` call needed).

- **`Surjector::realigning_surject_anchored`** — near-copy of `realigning_surject` (~150 lines) with:
  - `get_path_length(path_handle)` calls → use parameter.
  - `get_step_at_position(path_handle, X)` calls → leftmost/rightmost anchor's step handles.
  - `SubpathOverlay path_graph(path_position_graph, ...)` → `AnchoredSubpathOverlay path_graph(graph, ...)` (see below).

- **`Surjector::compute_disjoint_path_intervals_anchored`** — near-copy of `compute_disjoint_path_intervals` (~80 lines) with:
  - `get_position_of_step(ref_chunk.first/second)` calls → reads from `chunk_path_intervals[i].first/second`.
  - Rest of merging logic unchanged.

- **`AnchoredSubpathOverlay`** (new helper class, ~80 lines, in anonymous namespace of `surjector.cpp`). Implements `ExpandingOverlayGraph` over a step range on a `PathHandleGraph`. Mirrors `SubpathOverlay` API but constructor accepts `PathHandleGraph*` instead of `PathPositionHandleGraph*`. Internals are all `PathHandleGraph` operations: `get_next_step`, `get_handle_of_step`, `get_length`, `get_is_reverse`, etc.

#### 4.1.3 Files touched (vg side)

| File | Change |
|---|---|
| `src/surjector.hpp` | Add `PrecomputedAnchor`, declare 3 new methods on Surjector. |
| `src/surjector.cpp` | Define the 3 new methods + `AnchoredSubpathOverlay` helper class. |
| `src/subcommand/giraffe_server_main.cpp` | Add new stdin command `SURJECT_WITH_ANCHORS` (see §4.3). |
| `src/giraffe_engine.hpp` / `.cpp` | Add `engine.surject_with_anchors(...)` wrapper, parses the wire format. |

#### 4.1.4 Files NOT touched (regression-proofing)

- `extract_overlapping_paths` — unchanged.
- `realigning_surject` — unchanged.
- `spliced_surject` — unchanged.
- `compute_disjoint_path_intervals` / `compute_path_interval` — unchanged.
- `set_path_position` — unchanged.
- `surject_internal`, all public `surject(...)` overloads — unchanged.
- `SubpathOverlay` — unchanged (we add `AnchoredSubpathOverlay` next to it).
- `HaplotypeSurjector` in the engine — unchanged (still works for pre-indexed targets).

The existing call graph `surject() → surject_internal() → extract_overlapping_paths() → realigning_surject() → set_path_position()` keeps working bit-for-bit.

---

### 4.2 Coord-trans-server — pangenome-index

#### 4.2.1 New long-lived subcommand / binary

Currently `coordinate_translation` is a one-shot CLI binary that loads indexes per invocation. We add a long-lived server mode mirroring `giraffe-server`:

- New binary or subcommand: `coord-trans-server`.
- Loads at startup: GBZ, GBWT FastLocate (`.ri`), sampled tag array (`.tags`), translation tables (`.t1`, `.t2`) — same as the existing CLI but kept resident.
- Long-lived stdin loop with framed output (matches `giraffe-server`'s protocol).

#### 4.2.2 Stdin protocol

Two commands at launch:

```
BUILD_ANCHORS\t<read_name>\t<alignment_blob_base64>\t<target_haplotype_name>
TRANSLATE\t<src_haplotype>\t<start>\t<end>\t<tgt_haplotype>
PROCESS_BATCH         (alias: FLUSH_NOW)
```

`alignment_blob_base64` is the source `Alignment` proto, serialized via `vg::io::ProtobufIterator`-style encoding, then base64-encoded so it fits on one line.

#### 4.2.3 Response format (framed output)

For `BUILD_ANCHORS`:

```
ANCHORS\t<read_name>\t<status>\t<n_anchors>\t<target_path_length>\t<rev_strand>
<anchor_1_blob>
<anchor_2_blob>
...
<anchor_n_blob>
```

`<status>` is `ok` or one of `unknown_path` / `empty_alignment` / `no_common_nodes` for failure cases.

#### 4.2.4 Wire format for one anchor (per line)

Process-independent encoding. **No `step_handle_t` on the wire** — those are reconstructed by giraffe-server from path name + position.

```
read_begin\tread_end\tpath_offset_begin\tpath_offset_end\tn_mappings\t<mapping_1>;<mapping_2>;...;<mapping_n>
```

Each `<mapping_i>` is `node_id,is_reverse,node_offset,length[,edit_1,edit_2,...]` (matches encoded as empty edit list).

For BLAT-style queries this is small — typically a few hundred bytes per anchor.

Alternative: use protobuf for safety. Decision pending; ASCII works for v1, protobuf is the eventual right answer if the format stabilizes.

#### 4.2.5 New code in `pangenome-index`

| File | Status | Purpose |
|---|---|---|
| `include/pangenome_index/surject_anchor_builder.hpp` | new | declare `AnchorBuildResult`, `build_surject_anchors(...)` |
| `src/surject_anchor_builder.cpp` | new | implementation of anchor building |
| `src/coord_trans_server.cpp` | new | long-lived server wrapping the existing translation library + anchor builder |
| `include/pangenome_index/coord_trans_server.hpp` | new | server engine declaration |

#### 4.2.6 Anchor builder algorithm (`build_surject_anchors`)

Input:
- `gbwtgraph::GBZ& gbz`
- `gbwt::FastLocate& fast_locate`
- `const Alignment& source_alignment` (giraffe's graph path)
- `const std::string& target_haplotype_name`

Output:
```cpp
struct AnchorBuildResult {
    enum class Status { Ok, UnknownPath, EmptyAlignment, NoCommonNodes };
    Status status;
    std::vector<vg::PrecomputedAnchor> anchors;
    size_t target_path_length;
    vg::path_handle_t target_path_handle;
    bool target_rev_strand;
};
```

Algorithm (5 internal helpers):

**1. `resolve_target_path(gbz, target_haplotype_name) → (path_handle_t, gbwt_path_id)`**
Look up name in GBZ. Return `Status::UnknownPath` if not found.

**2. `compute_target_path_length(gbz, gbwt_path_id) → size_t`**
Walk the GBWT path summing `graph.get_length(handle)`. No position graph.

**3. `enumerate_source_visits(source_alignment) → vector<SourceVisit>`**
Walk `source.path().mapping(i)` for i=0..N-1. Each visit captures:
```
{ node_id, is_reverse, read_begin_offset, read_end_offset,
  mapping_pb_index, mapping (copy) }
```
Cumulative `mapping_to_length` produces read offsets.

**4. `find_target_visits_for_source(fast_locate, gbz, source_visits, gbwt_path_id) → vector<TargetVisit>`**

The work-horse. For each source visit, ask: does the target path also visit this `(node_id, orientation)`, and at what target base offset?

Sub-steps:
- Use `fast_locate.decompressSA(gbwt::Node::encode(node_id, is_reverse))` to find all GBWT sequences visiting this node in this orientation.
- Check if the target's forward seq_id (`2 * gbwt_path_id`) is among them; if not, try reverse seq_id (`2 * gbwt_path_id + 1`) → marks `target_rev_strand = true`.
- The target's GBWT offset on the node gives us a `gbwt::edge_type`. Convert to `step_handle_t` via `gbwtgraph::GBWTGraph::handle_to_path(...)` plus step-from-offset bridging.
- Get the cumulative base offset on the target by walking the target's GBWT path forward (or by exploiting `trace_coordinates_gbwt`'s machinery, which already does this).

This is essentially the same machinery `trace_coordinates_gbwt` already uses. Worth refactoring out a shared inner loop if maintainable.

Returns one `TargetVisit { source_mapping_idx, step_handle, target_base_offset_begin, target_base_offset_end }` per source visit that's also on the target. Source visits that the target doesn't traverse produce no entry.

**5. `group_visits_into_anchors(source_visits, target_visits) → vector<PrecomputedAnchor>`**

Walk source visits in order. Maintain a "current chunk in progress". On each step:

- If the source visit has a matching target visit AND the target step is `get_next_step(prev_target_step)` → **extend** current chunk: update its read_end, push the mapping into `graph_path`, update step_end and path_offset_end.
- If the source visit has a matching target visit but the target step is **not** the next step of the previous one (i.e., target took a detour we don't represent, or this is the first matching visit) → **close** the current chunk (if any), **open** a fresh chunk.
- If the source visit has no matching target visit → **close** the current chunk if open.

At chunk-close, build a `PrecomputedAnchor`:
- `read_begin/read_end` from first/last source visit's offsets.
- `graph_path` accumulated `path_t`.
- `step_begin/step_end` from first/last target visit's step handle.
- `path_offset_begin/path_offset_end` from first target visit's offset / last target visit's offset + last node length.

#### 4.2.7 What's reused vs. duplicated

- **Reused**: `find_first_and_last_common_nodes_gbwt`'s `decompressSA` + seq-id-check pattern.
- **Reused**: the inner GBWT-walk loop from `trace_coordinates_gbwt` (worth factoring out — see open question §6).
- **Duplicated**: the orientation-check logic from surjector's `extract_overlapping_paths` line 4135 (`is_reverse(handle) != is_reverse(handle_of_step)`).

---

### 4.3 Giraffe-server — new stdin command

Extend `giraffe-server`'s line-protocol with:

```
SURJECT_WITH_ANCHORS\t<read_name>\t<alignment_blob_base64>\t<target_path_name>\t<rev_strand>\t<target_path_length>\t<n_anchors>
<anchor_1_blob>
<anchor_2_blob>
...
<anchor_n_blob>
```

Dispatch:
1. Parse the source `Alignment` from base64.
2. Resolve `target_path_name` → `path_handle_t` via the GBZ.
3. Deserialize anchors. Reconstruct `step_handle_t` from `(target_path_name, path_offset_begin)` by walking the path. Cache the walk so multiple anchors in the same request are cheap.
4. Call `engine.surject_with_anchors(...)` which calls `surjector_->surject_with_anchors(...)`.
5. Emit the surjected GAF line (or framed `READ\t<name>\t<count>` + lines).

Add to help text. Add the option to giraffe-server's argv parsing if any startup config is needed (probably none — the server already has the GBZ loaded).

Files touched: `src/subcommand/giraffe_server_main.cpp`, plus a new public method `engine.surject_with_anchors(...)` on `GiraffeEngine` that takes the parsed inputs.

---

### 4.4 Middleware orchestrator — Python

#### 4.4.1 New `CoordinateTranslationMiddleware` class

Mirrors `GiraffeServerMiddleware`. Lives in `middleware/coord_trans_middleware.py` (new file).

```python
@dataclass
class CoordTransServerConfig:
    binary: str
    gbz_path: str
    gbwt_ri_path: str
    tags_path: str
    table1_path: str
    table2_path: str
    threads: int = 4

class CoordTransServerMiddleware:
    def start(self): ...
    def stop(self): ...
    def build_anchors(self, read_name, alignment_blob, target_haplotype) -> AnchorsResponse: ...
    def translate(self, src_hap, start, end, tgt_hap): ...
```

#### 4.4.2 New orchestrator: `PangenomeServerMiddleware`

Holds both `GiraffeServerMiddleware` and `CoordTransServerMiddleware`. Exposes the unified per-read API:

```python
class PangenomeServerMiddleware:
    def __init__(self, giraffe_cfg, coord_cfg):
        self.giraffe = GiraffeServerMiddleware(giraffe_cfg)
        self.coord   = CoordTransServerMiddleware(coord_cfg)

    def start(self):
        self.giraffe.start()
        self.coord.start()
        self._verify_gbz_match()   # fingerprint check (see §6.1)

    def map_and_surject(self, reads, target_haplotype) -> List[List[str]]:
        # 1. send reads to giraffe-server -> graph alignments
        graph_alns = self.giraffe.map_reads_get_alignments(reads)
        # 2. send each graph alignment to coord-trans-server -> anchors
        # 3. send (alignment, anchors, target metadata) back to giraffe-server -> linear surjections
        # 4. format and return
        ...
```

#### 4.4.3 Smoke test extension

Extend `smoke_test_middleware.py` with an `--orchestrate` mode that uses `PangenomeServerMiddleware`. Existing modes (giraffe-only, coord-only) keep working.

---

## 5. Phased implementation order

Each phase is independently testable and produces a working artifact.

### Phase 1 — Coord-trans-server skeleton

- Take the existing coordinate-translation library and wrap it in a long-lived `coord-trans-server` binary.
- Implement the `TRANSLATE` command (existing functionality, just over stdin instead of one-shot CLI).
- No anchor logic yet.

**Acceptance**: smoke test sends `TRANSLATE` commands and gets back the same results as the current CLI.

**Estimated effort**: 1 day.

### Phase 2 — Surjector::surject_with_anchors (vg-side, no coord-trans)

- Add `PrecomputedAnchor`, `surject_with_anchors`, `realigning_surject_anchored`, `compute_disjoint_path_intervals_anchored`, `AnchoredSubpathOverlay`.
- Test with hand-constructed anchors and a tiny test graph (the kind already used in `src/unittest/surject.cpp`).
- Verify against the existing `surject` output on identical inputs — outputs should be byte-identical when anchors are constructed from the same graph the old surject would produce.

**Acceptance**: unit test comparing `surject(...)` output to `surject_with_anchors(...)` output on a graph with pre-indexed paths produces identical alignments.

**Estimated effort**: 3–5 days, mostly because of the `AnchoredSubpathOverlay` correctness work.

### Phase 3 — Anchor builder (coord-trans-side)

- Add `build_surject_anchors(...)` to pangenome-index.
- Unit-test with constructed graph alignments against known target paths.
- Plug into coord-trans-server as `BUILD_ANCHORS` command.

**Acceptance**: given a graph alignment that passes through a known target path with a known set of detours, `build_surject_anchors` produces anchors matching the manually computed expected anchors.

**Estimated effort**: 3–4 days.

### Phase 4 — Wire format + giraffe-server SURJECT_WITH_ANCHORS

- Define and freeze the anchor wire format.
- Implement serialization in `coord-trans-server` and deserialization in `giraffe-server`.
- Add the `SURJECT_WITH_ANCHORS` stdin command.

**Acceptance**: manually constructed anchors fed via stdin produce correct surjections from giraffe-server.

**Estimated effort**: 2 days.

### Phase 5 — Middleware orchestrator

- Build `CoordTransServerMiddleware` Python class.
- Build `PangenomeServerMiddleware` that wires both servers.
- Extend smoke test.

**Acceptance**: `smoke_test_middleware.py --orchestrate --target-haplotype HG002#1#chr1` produces a surjected alignment for an arbitrary haplotype the server didn't pre-index.

**Estimated effort**: 2 days.

### Phase 6 — End-to-end testing + perf tuning

- Realistic queries against HPRC-v2 graph.
- Measure: total latency, per-step latency, IPC overhead, memory usage.
- Cache reconstructed `step_handle_t`s per request to amortize path-walk cost.

**Estimated effort**: 3–5 days.

**Total**: ~3 weeks of focused work.

---

## 6. Open questions for Adam

These are the things worth raising explicitly:

### 6.1 Cross-process GBZ consistency

Both servers MUST load the **same** GBZ file. If one reloads a fresher version while the other is still running, node IDs and path IDs may diverge, producing garbage with no error. The middleware should fingerprint (hash or mtime) the GBZ at startup of both processes and refuse to operate if they disagree.

**Question for Adam**: any existing convention in vg for fingerprinting graphs across processes? Or do we invent one?

### 6.2 `step_handle_t` reconstruction cost on giraffe-server

When giraffe-server receives anchors, it needs to convert `(target_path_name, path_offset_begin)` → `step_handle_t`. Today this requires `get_step_at_position(path_handle, path_offset)` which needs `PathPositionHandleGraph`. Three options to investigate:

1. Walk the path from start in `gbwtgraph::GBWTGraph` (PathHandleGraph) summing node lengths until we reach the target offset. O(steps) per request, but cacheable per (path, request).
2. Build a tiny on-demand position overlay covering just this one target path. ~hundreds of MB for a chromosome — comparable to caching one path's overlay anyway.
3. Have coord-trans-server include enough metadata in each anchor blob (e.g. cumulative steps from path start) that giraffe-server can reconstruct without walking.

Option 1 with per-request caching is probably fine for BLAT-style queries (read ~1 kb → maybe 10 anchors → 10 path walks of bounded length).

**Question for Adam**: any opinion on the right answer here? Does vg have helpers for step-handle reconstruction from base offsets that don't require PathPositionHandleGraph?

### 6.3 `SubpathOverlay` typing

`SubpathOverlay`'s constructor takes `PathPositionHandleGraph*` but its internals (per a quick read) use only `PathHandleGraph` operations during step iteration. Either:

(a) Loosen `SubpathOverlay` itself to accept `PathHandleGraph*` — vg-core change.
(b) Build `AnchoredSubpathOverlay` parallel to it — local change, no upstream impact.

**Question for Adam**: which would you prefer? If (a), I can do that as a separate PR.

### 6.4 Refactoring the GBWT-walk loop in coord translation

Both `trace_coordinates_gbwt` (existing) and `find_target_visits_for_source` (new) need to walk a target GBWT path forward producing per-node base offsets. Worth extracting that inner loop into a shared helper. Adam may have an opinion on where it lives (`src/coordinate_translation.cpp` is already 3k+ lines).

### 6.5 Spliced surject support

Should `surject_with_anchors` also support `preserve_deletions=true` (which dispatches to `spliced_surject`)? Spliced has its own DP machinery; supporting it would roughly double the implementation work. For BLAT-style use, spliced isn't needed.

**Question for Adam**: punt on spliced for v1? If yes, `surject_with_anchors` should reject `preserve_deletions=true` with a clear error.

### 6.6 Anchor wire format: ASCII vs. protobuf

ASCII is simpler for v1 — easy to debug, no schema dependencies. Protobuf is the eventual right answer for production. Question: ship ASCII and convert later, or pay the protobuf cost up front?

---

## 7. Risks and mitigations

| Risk | Severity | Mitigation |
|---|---|---|
| `AnchoredSubpathOverlay` behaves subtly differently from `SubpathOverlay` and produces wrong alignments. | High | Phase 2 unit tests: identical inputs through both paths → byte-identical outputs. Run on tiny graphs with known properties. |
| Orientation logic wrong in anchor builder → surjects to wrong strand silently. | High | Phase 3 unit tests on graphs with paths that share nodes in opposite orientations. Compare output strand to manual computation. |
| GBZ version mismatch between processes → silent garbage. | Medium | Fingerprint at middleware start (§6.1). Refuse to start if mismatch. |
| `step_handle_t` reconstruction is slow for long target paths. | Medium | Per-request cache. Profile early. |
| Cross-process IPC overhead dominates for short reads. | Low | Batch multiple reads per IPC round trip. |
| Coord-trans-server is a new long-lived process with its own lifecycle — leaks, crashes, hangs. | Medium | Reuse the supervision/restart pattern from `GiraffeServerMiddleware`. Stderr tail, output timeout, etc. |
| Anchor builder reports valid anchors but the alignment then has no overlap with the target → empty surject. | Low | Return `Status::NoCommonNodes` from `build_surject_anchors` cleanly; emit a sensible GAF "no surjection" record. |

---

## 8. Estimated counts of new code

| Component | Lines |
|---|---|
| `surjector.hpp` additions | ~40 |
| `surjector.cpp` additions (3 methods + `AnchoredSubpathOverlay`) | ~350 |
| `surject_anchor_builder.hpp` | ~50 |
| `surject_anchor_builder.cpp` | ~400 |
| `coord_trans_server.cpp` + header | ~300 |
| `giraffe_server_main.cpp` additions (new command) | ~80 |
| `GiraffeEngine` additions (anchor wrapper) | ~60 |
| Python middleware additions | ~250 |
| Unit tests | ~400 |
| Total | **~1.9 k lines**, ~3 weeks of focused work |

---

## 9. Backwards compatibility

- `vg surject` CLI: unchanged.
- `Surjector::surject(...)` public methods: unchanged.
- `GiraffeEngine::map_reads(...)` with `surjection_target_paths` pre-indexing: unchanged (fast path for known targets).
- Existing middleware `map_reads(...)`: unchanged.
- Existing coordinate-translation CLI and Python bindings: unchanged.

Everything in this plan is **additive**. The new path is opted into by either:
- Passing a per-read `surjection_target` that isn't in `surjection_target_paths` (engine falls back to anchor-based surjection if available), or
- Calling the new `PangenomeServerMiddleware.map_and_surject(...)` orchestrator.

---

## 10. Future extensions

- **Spliced surjection with anchors** (§6.5).
- **LRU eviction** of cached step-handle reconstructions and per-target temporary indexes.
- **Distributed deployment**: coord-trans-server and giraffe-server on different machines. Wire protocol already supports this since it's stdin/stdout-friendly and process-independent.
- **Protobuf wire format** if ASCII becomes a maintenance burden (§6.6).
- **vg core contribution**: if this anchor-based path proves robust, propose merging `PrecomputedAnchor` + `surject_with_anchors` upstream so other vg users benefit.

---

## Appendix A — Worked example for review

A 67-bp read going through `node1→node2→node3→node9→node5→node6→node7→node8` (variant path) is surjected onto target path `P_T = node1→node2→node3→node4→node5→node6→node7→node8` (node9 in source replaces node4 in target).

The coordinate translator produces:
- `target_path_length = 67`
- `target_rev_strand = false`
- 2 anchors:
  - **a0**: read[0, 23) ↔ target steps (S_T_1, S_T_3), path offsets [0, 23). Graph path: nodes 1→2→3.
  - **a1**: read[35, 67) ↔ target steps (S_T_5, S_T_8), path offsets [35, 67). Graph path: nodes 5→6→7→8.

The gap (read[23, 35)) is where the source went through node9. The surjector does DP between anchors against node4 on the target.

Full step-by-step trace through every Surjector function with concrete numbers: see `docs/surject_anchor_walkthrough.md` (companion document).

---

## Appendix B — Reading list for context

- [vg/src/surjector.hpp](../Giraffe_server/src/surjector.hpp) — Surjector public API.
- [vg/src/surjector.cpp](../Giraffe_server/src/surjector.cpp) — `extract_overlapping_paths` (line 4086), `realigning_surject` (line 3148), `compute_disjoint_path_intervals` (line 5142).
- [vg/src/subcommand/surject_main.cpp](../Giraffe_server/src/subcommand/surject_main.cpp) — how `vg surject` builds and uses `ReferencePathOverlay`.
- [vg/src/subcommand/giraffe_main.cpp:884-1072](../Giraffe_server/src/subcommand/giraffe_main.cpp) — preset infrastructure (for context on how blat etc. get applied).
- [pangenome-index/src/coordinate_translation.cpp](../src/coordinate_translation.cpp) — existing per-base coord translation; `trace_coordinates_gbwt` (line ~1531) is the closest existing analog to the anchor builder.
- [pangenome-index/src/coordinate_translation.cpp:1367-1373](../src/coordinate_translation.cpp) — `PathNode` struct (target-side per-node record that already carries everything an anchor needs).

---

*End of plan.*
