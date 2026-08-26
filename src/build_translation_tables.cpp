/**
 * build_translation_tables: build Translation Table 1 and Table 2 from a GBZ.
 *
 * Table 1 maps (named_path, global_interval) -> (path_id, local_interval).
 * Table 2 maps (source_path_id, target_haplotype) -> sorted segments, each
 * naming a TARGET PATH and the source range over which that path was observed.
 *
 * TABLE 2 IS ROUTING-ONLY, AND THAT IS DELIBERATE
 * ----------------------------------------------
 * Table 2 exists to answer one question for translate(): which target paths can
 * this source range reach? The coordinates are then produced by the GBWT trace.
 * It deliberately stores no target interval:
 *
 *   - Nothing in the query path reads one.
 *   - Producing one requires pairing a run's two end anchors onto a single
 *     target path, and that pairing is what caused misplaced blocks, five
 *     million abandoned runs and thirty-six million run splits when it was
 *     attempted. Routing needs no pairing and therefore cannot fail.
 *
 * The SOURCE RANGE per segment is what keeps queries fast, and it is not
 * optional. Give every target path of a run the run's full source range and a
 * query on a 150 Mb reference contig selects every target contig that run
 * touches anywhere -- around sixteen -- and translate() pays a full
 * find_tags_in_interval plus trace for each. Ranges derived from where each
 * path was actually probed keep that at roughly one candidate.
 *
 * ALGORITHM (no pairwise path comparison: that is O(paths^2 * length))
 *   Phase 1  Table 1, and each path's length, from one pass over every path.
 *   Phase 2a node -> bitmask of haplotypes visiting it (atomic OR, one pass).
 *   Phase 2b per source path: walk once, OR each node's mask into a coarse bin.
 *   Phase 2c per source path: find runs of covered bins per haplotype with an
 *            XOR delta between adjacent bins, merge runs closer than
 *            --merge-gap, then probe along each run with decompressSA. Record
 *            per target path WHICH probes saw it, as a bitmask, and emit one
 *            segment per contiguous stretch. The stretch matters: a repeat
 *            gives one isolated hit far from the real homology, and a plain
 *            min/max would stretch that path across the whole run.
 *
 * Usage:
 *   build_translation_tables <graph.gbz> <fastlocate.ri> <output.t1> <output.t2> [options]
 *
 * The r-index is required for Table 2 (decompressSA resolves target paths);
 * pass --only-table1 to build Table 1 alone, where it is unused.
 */

#include "pangenome_index/translation_tables.hpp"
#include <sdsl/simple_sds.hpp>
#include <gbwt/gbwt.h>
#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <tuple>
#include <atomic>
#include <cstdint>
#include <sstream>

using namespace std;
using namespace std::chrono;

// ============================================================
// Utilities
// ============================================================

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <graph.gbz> <fastlocate.ri> <output.t1> <output.t2> [options]" << endl;
    cerr << endl;
    cerr << "  graph.gbz       GBZ file (GBWT index + GBWTGraph)" << endl;
    cerr << "  fastlocate.ri   GBWT FastLocate r-index (needed for Table 2)" << endl;
    cerr << "  output.t1       Binary Translation Table 1 output" << endl;
    cerr << "  output.t2       Binary Translation Table 2 output" << endl;
    cerr << endl;
    cerr << "Options:" << endl;
    cerr << "  --dump              Print table contents to stderr after building" << endl;
    cerr << "  --debug             Enable verbose progress output" << endl;
    cerr << "  --only-table1       Build only Table 1 (r-index then unused)" << endl;
    cerr << "  --threads N         Threads (default: OMP_NUM_THREADS or max)" << endl;
    cerr << "  --bin-size N        Coverage granularity in bp (default 10000)" << endl;
    cerr << "  --merge-gap N       Merge runs separated by <= N bp (default 100000)" << endl;
    cerr << "  --min-run-bp N      Drop runs shorter than N bp (default 0 = keep all)" << endl;
    cerr << "  --probe-spacing N   One routing probe per ~N bp of run (default 1000000)." << endl;
    cerr << "                      Lower gives finer per-path source ranges and so" << endl;
    cerr << "                      fewer candidates per query, at more probes." << endl;
    cerr << "  --max-probes N      Cap probes per run (default 64, hard max 64:" << endl;
    cerr << "                      one bit per probe in the per-path bitmask)" << endl;
    cerr << "  --allow-same-haplotype  Also pair different contigs of one haplotype" << endl;
    cerr << "  --progress-every N  Progress line every N source paths (default 5000)" << endl;
    cerr << "  --help              Show this help" << endl;
}

/// Build the "sample#phase#contig" base name from GBWT metadata and a PathName.
static string build_base_name(const gbwt::Metadata& meta, const gbwt::PathName& pn) {
    string sample = (meta.samples() > 0 && pn.sample < meta.samples())
                    ? meta.sample(pn.sample) : to_string(pn.sample);
    string contig = (meta.contigs() > 0 && pn.contig < meta.contigs())
                    ? meta.contig(pn.contig) : to_string(pn.contig);
    return sample + "#" + to_string(pn.phase) + "#" + contig;
}

/// Build the "sample#phase" haplotype name (no contig) from a PathName.
static string build_haplotype_name(const gbwt::Metadata& meta, const gbwt::PathName& pn) {
    string sample = (meta.samples() > 0 && pn.sample < meta.samples())
                    ? meta.sample(pn.sample) : to_string(pn.sample);
    return sample + "#" + to_string(pn.phase);
}

// ============================================================
// Path metadata (lightweight; no node sequence)
// ============================================================

struct PathMetadata {
    size_t path_id        = 0;
    string base_name;
    string haplotype_name;
    string contig_name;
    size_t subpath_start  = 0;
    size_t length        = 0;  ///< length in bases (0 = skip)
};

// Full path data including node sequence; used only when loading a path for Phase 2.
struct PathData {
    size_t              path_id        = 0;
    string              base_name;
    string              haplotype_name;
    string              contig_name;
    size_t              subpath_start  = 0;
    size_t              length         = 0;
    gbwt::vector_type   nodes;
    vector<size_t>      node_offsets;
};

/// Load one path's node sequence and offsets from GBWT/graph (used in Phase 2 per-contig).
static PathData load_path_data(size_t path_id,
                               const gbwt::GBWT& gbwt_index,
                               const gbwtgraph::GBWTGraph& graph,
                               const gbwt::Metadata& meta,
                               const PathMetadata& pm)
{
    PathData pd;
    pd.path_id        = path_id;
    pd.base_name      = pm.base_name;
    pd.haplotype_name = pm.haplotype_name;
    pd.contig_name    = pm.contig_name;
    pd.subpath_start  = pm.subpath_start;
    pd.length         = pm.length;

    pd.nodes = gbwt_index.extract(gbwt::Path::encode(path_id, false));
    pd.node_offsets.reserve(pd.nodes.size() + 1);
    pd.node_offsets.push_back(0);
    size_t total = 0;
    for (gbwt::node_type node : pd.nodes) {
        if (node == gbwt::ENDMARKER) break;
        total += graph.get_length(
            graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node)));
        pd.node_offsets.push_back(total);
    }
    return pd;
}

// ============================================================
// Phase 1: Build Table 1 and collect path metadata only (no node sequences)
// ============================================================

static vector<PathMetadata> build_table1_and_collect_metadata(
        const gbwt::GBWT&        gbwt_index,
        const gbwtgraph::GBWTGraph& graph,
        const gbwt::Metadata&    meta,
        panindexer::TranslationTable1& table1,
        bool debug)
{
    size_t num_paths = meta.paths();
    vector<PathMetadata> all_meta(num_paths);

    // Parallel: extract each path only to compute length; do not store nodes.
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t path_id = 0; path_id < num_paths; ++path_id) {
        gbwt::PathName pn = meta.path(path_id);

        PathMetadata pm;
        pm.path_id        = path_id;
        pm.base_name      = build_base_name(meta, pn);
        pm.haplotype_name = build_haplotype_name(meta, pn);
        pm.contig_name    = (meta.contigs() > 0 && pn.contig < meta.contigs())
                            ? meta.contig(pn.contig) : to_string(pn.contig);
        pm.subpath_start  = static_cast<size_t>(pn.count);

        gbwt::vector_type nodes = gbwt_index.extract(gbwt::Path::encode(path_id, false));
        size_t total = 0;
        for (gbwt::node_type node : nodes) {
            if (node == gbwt::ENDMARKER) break;
            total += graph.get_length(
                graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node)));
        }
        pm.length = total;
        all_meta[path_id] = std::move(pm);
    }

    // Sequential: add to Table 1 in path_id order.
    size_t num_added = 0, num_skipped = 0;
    for (size_t path_id = 0; path_id < num_paths; ++path_id) {
        PathMetadata& pm = all_meta[path_id];
        if (pm.length == 0) {
            if (debug) {
                cerr << "  [skip] path_id=" << path_id
                     << " (" << pm.base_name << "[" << pm.subpath_start << "]) length=0" << endl;
            }
            ++num_skipped;
            continue;
        }
        table1.add_subpath(pm.base_name, path_id, pm.subpath_start, pm.length);
        ++num_added;
        if (debug) {
            cerr << "  [T1] path_id=" << path_id << " \"" << pm.base_name << "\""
                 << " offset=" << pm.subpath_start << " len=" << pm.length << endl;
        } else if (path_id % 500 == 0) {
            cerr << "  Phase1: " << path_id << "/" << num_paths << "\r" << flush;
        }
    }
    cerr << endl;
    cerr << "Table 1: " << num_added << " subpaths added, "
         << num_skipped << " skipped, "
         << table1.num_names() << " named paths." << endl;

    return all_meta;
}

// ============================================================
// Phase 2: Table 2 (routing) - which target paths a source range can reach
// ============================================================

/// One target path seen while probing a run, recorded as a BITMASK over probe
/// positions. A bitmask rather than a min/max range because a repeat produces
/// one isolated hit far from the real homology; min/max would stretch that path
/// across the entire run and make it a candidate for every query in it.
struct PathSpan { size_t pid; uint64_t mask; };

/// A target visit found by decompressSA: which path, and where along it.
struct Visit { size_t path_id; uint64_t sa_off; };

/// One emitted Table 2 segment, before it is written.
struct Row {
    size_t   src_path_id;
    uint32_t tgt_hap;      ///< index into hap_names
    size_t   src_start;
    size_t   src_end;
    size_t   tgt_path_id;
};

struct Table2Params {
    size_t bin_size       = 10000;
    size_t merge_gap      = 100000;
    size_t min_run_bp     = 0;
    size_t probe_spacing  = 1000000;
    size_t max_probes     = 64;     ///< hard max 64: one bit per probe
    size_t progress_every = 5000;
    bool   allow_same_hap = false;
};

static string commas(unsigned long long v) {
    string t = to_string(v), out;
    int c = 0;
    for (int i = (int)t.size() - 1; i >= 0; --i) {
        out.push_back(t[i]);
        if (++c % 3 == 0 && i > 0) out.push_back(',');
    }
    reverse(out.begin(), out.end());
    return out;
}

/**
 * Build Table 2 rows. Returns them unsorted; the caller sorts and writes.
 *
 * `hap_names` is filled with the haplotype name per index, and `path_hap` maps
 * each GBWT path id to its haplotype index.
 */
static vector<Row> build_table2_routing(
        const gbwt::GBWT& gbwt_index,
        const gbwtgraph::GBWTGraph& graph,
        const gbwt::Metadata& meta,
        const gbwt::FastLocate& rindex,
        const Table2Params& P,
        vector<string>& hap_names,
        vector<uint32_t>& path_hap)
{
    const size_t n_paths = meta.paths();
    path_hap.assign(n_paths, 0);
    hap_names.clear();
    {
        unordered_map<string, uint32_t> seen;
        for (size_t p = 0; p < n_paths; ++p) {
            gbwt::PathName pn = meta.path(p);
            string h = build_haplotype_name(meta, pn);
            auto it = seen.find(h);
            if (it == seen.end()) {
                uint32_t id = (uint32_t)hap_names.size();
                seen.emplace(h, id);
                hap_names.push_back(h);
                path_hap[p] = id;
            } else {
                path_hap[p] = it->second;
            }
        }
    }
    const size_t H = hap_names.size();
    const size_t WORDS = (H + 63) / 64;
    const size_t max_node = (size_t)graph.max_node_id();

    cerr << "  haplotypes: " << commas(H) << "   paths: " << commas(n_paths) << endl;
    if (n_paths > 10000000) {
        cerr << "  WARNING: " << commas(n_paths) << " paths. This looks like a"
             << " frequency-filtered graph; Table 2 scales with path count."
             << " Prefer the non-filtered graph." << endl;
    }

    // Node lengths: one dense array, reused by every walk. get_handle() +
    // get_length() per node was a large part of the per-node cost.
    cerr << "  node lengths ..." << endl;
    vector<uint32_t> node_len(max_node + 1, 0);
    #pragma omp parallel for schedule(static)
    for (size_t nid = 1; nid <= max_node; ++nid) {
        if (!graph.has_node((handlegraph::nid_t)nid)) continue;
        node_len[nid] = (uint32_t)graph.get_length(
            graph.get_handle((handlegraph::nid_t)nid, false));
    }

    // Phase 2a: node -> bitmask of haplotypes visiting it.
    cerr << "  node->haplotype masks ("
         << commas((max_node + 1) * WORDS * 8 / (1024 * 1024)) << " MB) ..." << endl;
    vector<uint64_t> nodemask;
    try {
        nodemask.assign((max_node + 1) * WORDS, 0);
    } catch (const std::bad_alloc&) {
        cerr << "Error: cannot allocate the node mask. Reduce haplotypes or add RAM." << endl;
        return {};
    }
    #pragma omp parallel for schedule(dynamic, 8)
    for (size_t p = 0; p < n_paths; ++p) {
        const uint32_t h = path_hap[p];
        const size_t word = h / 64;
        const uint64_t bit = 1ULL << (h % 64);
        gbwt::vector_type nodes = gbwt_index.extract(gbwt::Path::encode(p, false));
        for (gbwt::node_type node : nodes) {
            if (node == gbwt::ENDMARKER) break;
            const size_t nid = (size_t)gbwt::Node::id(node);
            if (nid > max_node) continue;
            uint64_t* slot = &nodemask[nid * WORDS + word];
            // Paths of different haplotypes share nodes, so this OR races.
            if ((__atomic_load_n(slot, __ATOMIC_RELAXED) & bit) == 0) {
                __atomic_fetch_or(slot, bit, __ATOMIC_RELAXED);
            }
        }
    }

    // Phases 2b/2c.
    vector<vector<Row>> per_thread(omp_get_max_threads());
    std::atomic<size_t> probes_done{0}, runs_dropped{0}, no_target{0}, done{0};
    const size_t merge_gap_bins = P.merge_gap / P.bin_size;
    auto t0 = high_resolution_clock::now();
    cerr << "  coverage bins and routing probes ..." << endl;

    #pragma omp parallel
    {
        vector<Row>& out = per_thread[omp_get_thread_num()];
        vector<uint64_t> bins;
        vector<uint32_t> run_start(H, 0), run_last(H, 0);
        vector<uint8_t>  run_open(H, 0);
        vector<uint32_t> src_nid, src_off;
        vector<Visit>    hits;
        vector<PathSpan> seen;
        vector<size_t>   probe_at;

        // Visits of `nid` by paths of haplotype h. Only EVEN sequence ids are
        // taken: a bidirectional GBWT stores each path twice and the reverse
        // copy measures offsets along the reversed sequence, so restricting to
        // the forward copy keeps one coordinate system. Probing both node
        // orientations still finds targets traversing it the other way.
        auto probe = [&](uint32_t nid, uint32_t h) {
            hits.clear();
            for (int orient = 0; orient < 2; ++orient) {
                gbwt::node_type gn = gbwt::Node::encode((gbwt::size_type)nid, orient == 1);
                vector<gbwt::size_type> sa = rindex.decompressSA(gn);
                ++probes_done;
                for (gbwt::size_type v : sa) {
                    const gbwt::size_type sid = rindex.seqId(v);
                    if (sid % 2 != 0) continue;
                    const size_t pid = (size_t)sid / 2;
                    if (pid < n_paths && path_hap[pid] == h) {
                        hits.push_back(Visit{pid, (uint64_t)rindex.seqOffset(v)});
                    }
                }
            }
        };

        #pragma omp for schedule(dynamic, 1)
        for (size_t sp = 0; sp < n_paths; ++sp) {
            // Phase 2b: one walk. Node list plus bin coverage.
            gbwt::vector_type nodes = gbwt_index.extract(gbwt::Path::encode(sp, false));
            src_nid.clear(); src_off.clear();
            size_t path_len = 0;
            for (gbwt::node_type node : nodes) {
                if (node == gbwt::ENDMARKER) break;
                const size_t nid = (size_t)gbwt::Node::id(node);
                if (nid > max_node) continue;
                src_nid.push_back((uint32_t)nid);
                src_off.push_back((uint32_t)path_len);
                path_len += node_len[nid];
            }
            if (src_nid.empty() || path_len == 0) { ++done; continue; }

            const size_t nbins = (path_len + P.bin_size - 1) / P.bin_size;
            bins.assign(nbins * WORDS, 0);
            for (size_t i = 0; i < src_nid.size(); ++i) {
                const uint64_t* m = &nodemask[(size_t)src_nid[i] * WORDS];
                // Mark EVERY bin the node spans. Marking only the bin holding
                // its start leaves the interior of any node longer than a bin
                // uncovered, which breaks runs apart and loses real homology.
                const size_t nlen = node_len[src_nid[i]];
                const size_t b_lo = src_off[i] / P.bin_size;
                const size_t b_hi = nlen > 0 ? (src_off[i] + nlen - 1) / P.bin_size : b_lo;
                for (size_t b = b_lo; b <= b_hi && b < nbins; ++b) {
                    uint64_t* dst = &bins[b * WORDS];
                    for (size_t w = 0; w < WORDS; ++w) dst[w] |= m[w];
                }
            }

            auto emit_run = [&](uint32_t h, size_t b0, size_t b1) {
                if (!P.allow_same_hap && h == path_hap[sp]) return;
                size_t src_start = b0 * P.bin_size;
                size_t src_end = min((b1 + 1) * P.bin_size, path_len);
                if (src_end <= src_start) return;
                // Filter on span, not bin count: a run can occupy two bins and
                // still be far shorter than the threshold.
                if (src_end - src_start < P.min_run_bp) { ++runs_dropped; return; }

                const size_t word = h / 64;
                const uint64_t bit = 1ULL << (h % 64);
                auto on_h = [&](size_t i) {
                    return (nodemask[(size_t)src_nid[i] * WORDS + word] & bit) != 0;
                };
                size_t lo = lower_bound(src_off.begin(), src_off.end(),
                                        (uint32_t)src_start) - src_off.begin();
                size_t hi = lower_bound(src_off.begin(), src_off.end(),
                                        (uint32_t)src_end) - src_off.begin();
                if (hi > src_nid.size()) hi = src_nid.size();
                // Bin coverage means "some node here is on h", so the run's
                // outermost nodes need not be; find ones that are.
                size_t i_first = hi, i_last = hi;
                for (size_t i = lo; i < hi; ++i) { if (on_h(i)) { i_first = i; break; } }
                if (i_first >= hi) { ++no_target; return; }
                for (size_t i = hi; i > i_first; --i) {
                    if (on_h(i - 1)) { i_last = i - 1; break; }
                }
                if (i_last >= hi) i_last = i_first;

                // Probe along the run. Count scales with LENGTH: a fixed handful
                // over a 15 Mb run leaves multi-megabase gaps, and a target
                // contig contributing only inside a gap is never recorded --
                // which a query landing there sees as "region does not exist".
                size_t n_probe = (src_end - src_start) / P.probe_spacing + 1;
                if (n_probe < 5) n_probe = 5;
                if (n_probe > P.max_probes) n_probe = P.max_probes;
                const size_t span_i = i_last - i_first;
                probe_at.assign(n_probe, src_end);
                seen.clear();
                for (size_t k = 0; k < n_probe; ++k) {
                    const size_t want_i =
                        i_first + (span_i * k) / (n_probe > 1 ? n_probe - 1 : 1);
                    for (size_t i = want_i; i <= i_last; ++i) {
                        if (!on_h(i)) continue;
                        probe(src_nid[i], h);
                        probe_at[k] = src_off[i];
                        const uint64_t kbit = 1ULL << k;
                        for (const Visit& v : hits) {
                            if (v.path_id == sp) continue;   // never map to itself
                            bool found = false;
                            for (PathSpan& ps : seen) {
                                if (ps.pid == v.path_id) { ps.mask |= kbit; found = true; break; }
                            }
                            if (!found) seen.push_back(PathSpan{v.path_id, kbit});
                        }
                        break;
                    }
                }
                if (seen.empty()) { ++no_target; return; }

                // One segment per maximal contiguous stretch of probes, padded
                // by a probe interval: a path seen at one probe may extend most
                // of the way to its neighbours, and the table must not
                // under-cover. An isolated repeat hit stays its own small
                // segment instead of widening the real one.
                const size_t pad = (src_end - src_start) / n_probe + 1;
                for (const PathSpan& ps : seen) {
                    for (size_t k = 0; k < n_probe; ) {
                        if (!(ps.mask & (1ULL << k))) { ++k; continue; }
                        size_t k0 = k;
                        while (k < n_probe && (ps.mask & (1ULL << k))) ++k;
                        const size_t p_lo = probe_at[k0], p_hi = probe_at[k - 1];
                        Row r;
                        r.src_path_id = sp;
                        r.tgt_hap     = h;
                        r.src_start   = (p_lo > src_start + pad) ? p_lo - pad : src_start;
                        r.src_end     = min(src_end, p_hi + pad);
                        if (r.src_end <= r.src_start) continue;
                        r.tgt_path_id = ps.pid;
                        out.push_back(r);
                    }
                }
            };

            // Phase 2c: runs via XOR delta between adjacent bins, gap-merged
            // inline so no intermediate run list is built.
            fill(run_open.begin(), run_open.end(), 0);
            for (size_t b = 0; b < nbins; ++b) {
                for (size_t w = 0; w < WORDS; ++w) {
                    const uint64_t cur  = bins[b * WORDS + w];
                    const uint64_t prev = (b == 0) ? 0ULL : bins[(b - 1) * WORDS + w];
                    uint64_t diff = cur ^ prev;
                    while (diff) {
                        const int t = __builtin_ctzll(diff);
                        diff &= diff - 1;
                        const uint32_t h = (uint32_t)(w * 64 + t);
                        if (h >= H) continue;
                        if (cur & (1ULL << t)) {
                            if (run_open[h] && b - run_last[h] - 1 <= merge_gap_bins) {
                                // close enough to the previous stretch: same run
                            } else {
                                if (run_open[h]) emit_run(h, run_start[h], run_last[h]);
                                run_start[h] = (uint32_t)b;
                            }
                            run_open[h] = 1;
                            run_last[h] = (uint32_t)b;
                        } else if (run_open[h]) {
                            run_last[h] = (uint32_t)(b - 1);
                        }
                    }
                }
            }
            for (size_t w = 0; w < WORDS; ++w) {
                uint64_t bits = bins[(nbins - 1) * WORDS + w];
                while (bits) {
                    const int t = __builtin_ctzll(bits);
                    bits &= bits - 1;
                    const uint32_t h = (uint32_t)(w * 64 + t);
                    if (h < H && run_open[h]) run_last[h] = (uint32_t)(nbins - 1);
                }
            }
            for (uint32_t h = 0; h < H; ++h) {
                if (run_open[h]) emit_run(h, run_start[h], run_last[h]);
            }

            const size_t d = ++done;
            if (P.progress_every && d % P.progress_every == 0) {
                const double el = duration_cast<milliseconds>(
                    high_resolution_clock::now() - t0).count() / 1000.0;
                const double rate = el > 0 ? d / el : 0;
                #pragma omp critical
                cerr << "    " << commas(d) << "/" << commas(n_paths) << " paths  "
                     << fixed << setprecision(0) << el << "s  eta "
                     << (rate > 0 ? (n_paths - d) / rate : 0) << "s  "
                     << commas(probes_done.load()) << " probes" << endl;
            }
        }
    }

    size_t total = 0;
    for (const auto& v : per_thread) total += v.size();
    vector<Row> rows;
    rows.reserve(total);
    for (auto& v : per_thread) {
        rows.insert(rows.end(), v.begin(), v.end());
        vector<Row>().swap(v);
    }
    cerr << "  probes: " << commas(probes_done.load())
         << "   runs with no target: " << commas(no_target.load())
         << "   runs dropped by --min-run-bp: " << commas(runs_dropped.load()) << endl;
    return rows;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        usage(argv[0]);
        return 1;
    }

    string gbz_file    = argv[1];
    string ri_file     = argv[2];
    string output_t1   = argv[3];
    string output_t2   = argv[4];
    bool dump          = false;
    bool debug         = false;
    bool only_table1   = false;
    int num_threads    = 0;  // 0 = use default (OMP_NUM_THREADS / omp_get_max_threads())
    Table2Params P;

    for (int i = 5; i < argc; ++i) {
        string arg = argv[i];
        auto need = [&](const char* what) -> bool {
            if (i + 1 >= argc) { cerr << "Error: " << what << " requires a value" << endl; return false; }
            return true;
        };
        if      (arg == "--dump")        dump = true;
        else if (arg == "--debug")       debug = true;
        else if (arg == "--only-table1") only_table1 = true;
        else if (arg == "--allow-same-haplotype") P.allow_same_hap = true;
        else if (arg == "--threads") {
            if (!need("--threads")) return 1;
            num_threads = atoi(argv[++i]);
            if (num_threads < 1) num_threads = 1;
        }
        else if (arg == "--bin-size")       { if (!need(arg.c_str())) return 1; P.bin_size = stoull(argv[++i]); }
        else if (arg == "--merge-gap")      { if (!need(arg.c_str())) return 1; P.merge_gap = stoull(argv[++i]); }
        else if (arg == "--min-run-bp")     { if (!need(arg.c_str())) return 1; P.min_run_bp = stoull(argv[++i]); }
        else if (arg == "--probe-spacing")  { if (!need(arg.c_str())) return 1; P.probe_spacing = stoull(argv[++i]); }
        else if (arg == "--max-probes")     { if (!need(arg.c_str())) return 1; P.max_probes = stoull(argv[++i]); }
        else if (arg == "--progress-every") { if (!need(arg.c_str())) return 1; P.progress_every = stoull(argv[++i]); }
        else if (arg == "--help" || arg == "-h") { usage(argv[0]); return 0; }
        else { cerr << "Unknown argument: " << arg << endl; usage(argv[0]); return 1; }
    }
    if (P.bin_size == 0)      { cerr << "Error: --bin-size must be > 0" << endl; return 1; }
    if (P.probe_spacing == 0) P.probe_spacing = 1;
    if (P.max_probes == 0)    P.max_probes = 1;
    // One bit per probe in the per-path bitmask.
    if (P.max_probes > 64)    P.max_probes = 64;

    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    cerr << "Using " << (num_threads > 0 ? num_threads : omp_get_max_threads()) << " thread(s)." << endl;

    auto t_start = high_resolution_clock::now();

    // ---- Load GBZ ----
    cerr << "Loading GBZ: " << gbz_file << " ..." << endl;
    gbwtgraph::GBZ gbz;
    try {
        sdsl::simple_sds::load_from(gbz, gbz_file);
    } catch (const exception& e) {
        cerr << "Error loading GBZ: " << e.what() << endl;
        return 1;
    }
    const gbwt::GBWT&          gbwt_index = gbz.index;
    const gbwtgraph::GBWTGraph& graph      = gbz.graph;

    auto t_loaded = high_resolution_clock::now();
    cerr << "GBZ loaded in "
         << duration_cast<milliseconds>(t_loaded - t_start).count()
         << " ms. Paths: " << gbwt_index.sequences() / 2
         << ", nodes: " << graph.get_node_count() << endl;

    if (!gbwt_index.hasMetadata()) {
        cerr << "Error: GBWT index has no metadata." << endl;
        return 1;
    }
    const gbwt::Metadata& meta = gbwt_index.metadata;
    cerr << "Metadata: " << meta.paths() << " paths, "
         << meta.samples() << " samples, "
         << meta.contigs() << " contigs" << endl;

    // ---- Phase 1: Table 1 (metadata only; no path node sequences stored) ----
    cerr << "\n=== Phase 1: Building Translation Table 1 ===" << endl;
    panindexer::TranslationTable1 table1;
    auto t_p1_start = high_resolution_clock::now();

    vector<PathMetadata> path_meta = build_table1_and_collect_metadata(
        gbwt_index, graph, meta, table1, debug);

    auto t_p1_end = high_resolution_clock::now();
    cerr << "Phase 1 done in "
         << duration_cast<milliseconds>(t_p1_end - t_p1_start).count() << " ms." << endl;

    // ---- Phase 2: Table 2 (routing) ----
    vector<Row> rows;
    vector<string> hap_names;
    vector<uint32_t> path_hap;
    if (!only_table1) {
        cerr << "\n=== Phase 2: Building Translation Table 2 (routing) ===" << endl;
        auto t_p2_start = high_resolution_clock::now();

        cerr << "  loading r-index: " << ri_file << " ..." << endl;
        gbwt::FastLocate rindex;
        {
            ifstream in(ri_file, ios::binary);
            if (!in) { cerr << "Error: cannot open r-index " << ri_file << endl; return 1; }
            rindex.load(in);
        }
        rindex.setGBWT(gbwt_index);

        rows = build_table2_routing(gbwt_index, graph, meta, rindex, P,
                                    hap_names, path_hap);

        // Key order must match TranslationTable2's: source path id, then
        // haplotype NAME. Interned ids follow first-seen order, not alphabetical.
        sort(rows.begin(), rows.end(), [&](const Row& a, const Row& b) {
            if (a.src_path_id != b.src_path_id) return a.src_path_id < b.src_path_id;
            if (a.tgt_hap != b.tgt_hap) return hap_names[a.tgt_hap] < hap_names[b.tgt_hap];
            if (a.src_start != b.src_start) return a.src_start < b.src_start;
            // Tiebreaker so the file cannot depend on which thread produced a
            // row: sort is not stable and rows arrive in thread order.
            return a.tgt_path_id < b.tgt_path_id;
        });

        auto t_p2_end = high_resolution_clock::now();
        cerr << "Phase 2 done in "
             << duration_cast<milliseconds>(t_p2_end - t_p2_start).count() << " ms."
             << "  segments: " << commas(rows.size()) << endl;
    }

    // ---- Dump (optional) ----
    if (dump) {
        cerr << "\n=== Table 1 contents ===" << endl;
        for (const string& name : table1.names()) {
            auto subpaths = table1.subpaths(name);
            cerr << name << "  (" << subpaths.size() << " subpath(s))" << endl;
            for (const auto& sp : subpaths) {
                cerr << "  path_id=" << sp.path_id
                     << "  offset=" << sp.subpath_start
                     << "  len=" << sp.length
                     << "  global=[" << sp.subpath_start << "," << sp.end() << ")" << endl;
            }
        }

        if (!only_table1) {
            cerr << "\n=== Table 2 contents ===" << endl;
            for (const Row& r : rows) {
                cerr << "src_path=" << r.src_path_id
                     << " tgt_hap=\"" << hap_names[r.tgt_hap] << "\""
                     << "  src=[" << r.src_start << "," << r.src_end << ")"
                     << " tgt_pid=" << r.tgt_path_id << endl;
            }
        }
    }

    // ---- Serialize Table 1 ----
    cerr << "\nWriting Table 1 to: " << output_t1 << " ..." << endl;
    {
        ofstream out(output_t1, ios::binary);
        if (!out) { cerr << "Error: cannot open " << output_t1 << endl; return 1; }
        table1.serialize(out);
    }

    // ---- Serialize Table 2 ----
    // Streamed, not via TranslationTable2: an all-pairs table has tens of
    // millions of keys and that map costs several GB of node, string and vector
    // overhead on top of the payload.
    if (!only_table1) {
        cerr << "Writing Table 2 to: " << output_t2 << " ..." << endl;
        ofstream out(output_t2, ios::binary);
        if (!out) { cerr << "Error: cannot open " << output_t2 << endl; return 1; }
        panindexer::TranslationTable2Writer writer(out, hap_names,
                                                   /*with_target_coords=*/false);
        size_t keys = 0;
        bool open_key = false;
        size_t cur_src = 0; uint32_t cur_hap = 0;
        for (const Row& r : rows) {
            if (!open_key || r.src_path_id != cur_src || r.tgt_hap != cur_hap) {
                if (!writer.begin_key(r.src_path_id, hap_names[r.tgt_hap])) {
                    cerr << "Error: key order violated at src_path " << r.src_path_id
                         << " hap " << hap_names[r.tgt_hap] << endl;
                    return 1;
                }
                open_key = true; cur_src = r.src_path_id; cur_hap = r.tgt_hap; ++keys;
            }
            panindexer::IntervalMapping m;
            m.src_start   = r.src_start;
            m.src_end     = r.src_end;
            m.tgt_path_id = r.tgt_path_id;
            writer.add_segment(m);
        }
        writer.finish();
        cerr << "  " << commas(keys) << " keys, " << commas(rows.size())
             << " segments." << endl;
    }

    auto t_end = high_resolution_clock::now();
    cerr << "\nDone. Total time: "
         << duration_cast<milliseconds>(t_end - t_start).count() << " ms." << endl;

    return 0;
}
