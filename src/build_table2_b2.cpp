/*
 * build_table2_b2.cpp
 *
 * Build an ALL-PAIRS Translation Table 2 in "B2" form: entries carry BOTH the
 * source and the target interval —
 *
 *     path X interval [a, b)  ->  path Y interval [c, d)   (+ orientation)
 *
 * as opposed to the stock table (build_translation_tables) and the coarse table
 * (build_table2_coarse), whose entries name only the target path and leave the
 * target coordinates to be recovered by a GBWT traversal at query time.
 *
 * WHY THE TARGET INTERVAL IS WORTH ITS BYTES
 * -----------------------------------------
 * Without [c, d) a Table 2 entry cannot be composed or refined: a lookup
 * returns "this source range touches path Y", so the caller still has to walk
 * path Y to find out where. On contig-level paths that walk has the whole
 * contig to search. Storing [c, d) turns the entry into a bound, so the
 * traversal starts already localized — and it makes two lookups chainable,
 * which a target-path-only entry never is.
 *
 * SIZE
 * ----
 * Entries are one per (source path, target path) syntenic block rather than one
 * per variant, so the count is governed by path count, not heterozygosity:
 *
 *   filtered graph      (223.8 M paths, ~6.4 kb each)   ~10^11 entries — do not
 *   non-filtered graph  ( 94.6 k paths, ~15 Mb each)     ~44 M entries, ~2.8 GB
 *
 * Measured on chrM (44 haplotypes): 1,892 segments = exactly 44 x 43, i.e. one
 * block per ordered haplotype pair, at 64.4 bytes each (44 payload + 20 key,
 * since a key here holds a single segment). The same graph at variant
 * resolution is ~93,000 segments, so this is ~49x smaller AND carries the
 * target coordinates the variant-resolution table does not.
 *
 * Use the NON-FILTERED graph. Frequency filtering chops each haplotype into
 * ~2,400 fragments per contig, and since entries scale with path count that
 * multiplies the table by the same factor for no extra information.
 *
 * SAFETY: block boundaries are only ever WIDENED (bin granularity, gap
 * merging), so the table can claim homology slightly too eagerly. A false
 * positive costs one traversal that then correctly returns nothing; a false
 * negative would silently lose a real translation. Nothing here narrows a
 * block, so the latter cannot happen.
 *
 * ALGORITHM (no pairwise path comparison — that would be O(paths^2 * length))
 * --------------------------------------------------------------------------
 *   Phase A  node -> bitmask of haplotypes visiting it. One pass over every
 *            path, atomic OR into a flat array. Cost: total path length.
 *   Phase B  per source path: walk it once, recording (node id, base offset)
 *            and OR-ing each node's mask into a coarse bin.
 *   Phase C  per source path: runs of covered bins per haplotype via an XOR
 *            delta between adjacent bins (presence rarely changes, so this
 *            costs bins*words, not bins*haplotypes); merge runs separated by
 *            less than --merge-gap. For each run, find the first and last
 *            source node actually on that haplotype and probe both with
 *            decompressSA to get the target path id and the two SA offsets.
 *   Phase D  per target path: walk it once to turn the SA offsets recorded in
 *            phase C into base offsets, giving [c, d) and the orientation.
 *   Phase E  sort by key and stream the file out.
 *
 * Phase D exists because decompressSA locates a visit in NODE units along the
 * path, not bases; converting needs that path's cumulative node lengths. Doing
 * it as a separate grouped pass costs one extra walk per target path — O(paths
 * * length) — instead of materializing prefix sums for every path at once,
 * which would not fit. Carrying the SA offset (not just the node id) is what
 * makes it exact: a node can occur many times on a path, and the SA value
 * identifies which occurrence, so no disambiguation heuristic is needed.
 *
 * Threading: phases A, B/C and D are each an independent parallel_for; B/C
 * accumulates into per-thread row vectors that are merged once at the end.
 *
 * Usage:
 *   build_table2_b2 <graph.gbz> <fastlocate.ri> <output.t2> [options]
 *
 * Table 1 is NOT needed and NOT rebuilt: source coordinates here are
 * path-local, exactly as Table 2 stores them, and path ids come from the GBWT.
 *
 * BUILD NOTE: compile and link against the SAME gbwtgraph (headers and
 * library). Mixing the vendored headers with a differently-built
 * libgbwtgraph.a silently corrupts the loaded GBZ and adjacent stack values
 * rather than failing to link — run `make gbwtgraph-lib` so the vendored
 * archive exists, or make sure $(LIB_DIR) matches the headers being included.
 */

#include <gbwt/fast_locate.h>
#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>
#include <sdsl/simple_sds.hpp>

#include "pangenome_index/translation_tables.hpp"

#include <algorithm>
#include <atomic>
#include <functional>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

/// How many times a run may be halved when no plausible anchor pairing exists.
/// 6 gives 64 pieces, far past any real number of target contigs per run.
constexpr int MAX_SPLIT_DEPTH = 6;

void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " <graph.gbz> <fastlocate.ri> <output.t1> <output.t2> [options]\n\n"
        << "Builds BOTH translation tables in one pass:\n"
        << "  Table 1  contig name + global interval -> (path id, local interval)\n"
        << "  Table 2  all-pairs, entries carrying the TARGET interval as well as\n"
        << "           the source interval (\"B2\" form, file version 3)\n\n"
        << "Use the NON-FILTERED graph: entries scale with path count, and\n"
        << "frequency filtering multiplies path count by ~2,400 for no gain.\n\n"
        << "  --threads N     worker threads (default: all)\n"
        << "  --bin-size N    coverage granularity in bp (default 10000)\n"
        << "  --merge-gap N   merge runs separated by <= N bp (default 100000)\n"
        << "  --min-run-bp N  drop runs shorter than N bp (default 0 = keep all).\n"
        << "                  Repeats make unrelated haplotypes share short\n"
        << "                  stretches; the report below shows how many runs sit\n"
        << "                  under each threshold so you can pick one.\n"
        << "  --max-paths N   process only the first N source paths (throughput probe)\n"
        << "  --progress-every N  progress line every N paths (default 1000)\n"
        << "  --progress      per-phase progress to stderr\n"
        << "  --no-t1         skip Table 1 (still requires the output.t1 argument)\n"
        << "  --target-intervals  also store the TARGET interval per segment\n"
        << "                  (file version 3). OFF by default: nothing in the\n"
        << "                  query path reads it, and producing it requires\n"
        << "                  pairing the run's two end anchors onto one target\n"
        << "                  path -- the step responsible for misplaced blocks,\n"
        << "                  abandoned runs and all the run splitting. Routing\n"
        << "                  alone needs no pairing and cannot fail.\n"
        << "  --allow-same-haplotype  also pair different contigs of the SAME\n"
        << "                  haplotype (off by default, matching the stock\n"
        << "                  builder). A path is never mapped to itself.\n";
}

std::string with_commas(unsigned long long v) {
    std::string s = std::to_string(v), out;
    int c = 0;
    for (int i = static_cast<int>(s.size()) - 1; i >= 0; --i) {
        out.push_back(s[i]);
        if (++c % 3 == 0 && i > 0) out.push_back(',');
    }
    std::reverse(out.begin(), out.end());
    return out;
}

std::string human_bytes(double bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    int u = 0;
    while (bytes >= 1024.0 && u < 4) { bytes /= 1024.0; ++u; }
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(bytes < 10 ? 2 : 1) << bytes << " " << units[u];
    return ss.str();
}

/**
 * One emitted row. The target interval is not known yet when phase C creates
 * this: `sa_first` / `sa_last` are the SA (seqOffset) positions of the run's
 * two anchor nodes on the target path, which phase D converts to base offsets.
 */
struct Row {
    size_t   src_path_id = 0;
    uint32_t tgt_hap     = 0;   ///< index into hap_names
    size_t   src_start   = 0;
    size_t   src_end     = 0;
    size_t   tgt_path_id = 0;
    uint64_t sa_first    = 0;   ///< seqOffset of the first anchor on the target
    uint64_t sa_last     = 0;   ///< seqOffset of the last anchor on the target
    size_t   tgt_start   = 0;   ///< filled by phase D
    size_t   tgt_end     = 0;   ///< filled by phase D
    bool     tgt_reverse = false;
    bool     resolved    = false;
};

/// A target visit found by decompressSA: which path, and where on it.
struct Visit {
    size_t   path_id;
    uint64_t sa_off;   ///< seqOffset; LARGER means EARLIER in the path
};

} // namespace

int main(int argc, char** argv) {
    std::string gbz_path, ri_path, t1_path, out_path;
    int threads = 0;
    size_t bin_size = 10000;
    size_t merge_gap = 100000;
    size_t min_run_bp = 0;
    size_t max_paths = 0;
    size_t progress_every = 1000;
    bool progress = false;
    bool allow_same_hap = false;
    bool build_t1 = true;
    bool want_tgt_intervals = false;

    {
        std::vector<std::string> pos;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
            else if (a == "--threads" && i + 1 < argc) threads = std::stoi(argv[++i]);
            else if (a == "--bin-size" && i + 1 < argc) bin_size = std::stoull(argv[++i]);
            else if (a == "--merge-gap" && i + 1 < argc) merge_gap = std::stoull(argv[++i]);
            else if (a == "--min-run-bp" && i + 1 < argc) min_run_bp = std::stoull(argv[++i]);
            else if (a == "--max-paths" && i + 1 < argc) max_paths = std::stoull(argv[++i]);
            else if (a == "--progress-every" && i + 1 < argc) progress_every = std::stoull(argv[++i]);
            else if (a == "--progress") progress = true;
            else if (a == "--allow-same-haplotype") allow_same_hap = true;
            else if (a == "--no-t1") build_t1 = false;
            else if (a == "--target-intervals") want_tgt_intervals = true;
            else if (!a.empty() && a[0] == '-') { usage(argv[0]); return 1; }
            else pos.push_back(a);
        }
        if (pos.size() != 4) { usage(argv[0]); return 1; }
        gbz_path = pos[0]; ri_path = pos[1]; t1_path = pos[2]; out_path = pos[3];
    }
    if (bin_size == 0) { std::cerr << "--bin-size must be > 0\n"; return 1; }
    if (progress_every == 0) progress_every = 1;

#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    std::cerr << "Threads: " << omp_get_max_threads() << std::endl;
#endif

    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
    auto secs = [](clk::time_point a, clk::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    // ---------------------------------------------------------------- load
    std::cerr << "Loading " << gbz_path << " ..." << std::endl;
    gbwtgraph::GBZ gbz;
    sdsl::simple_sds::load_from(gbz, gbz_path);
    const gbwt::GBWT& index = gbz.index;
    if (!index.hasMetadata() || !index.metadata.hasPathNames() ||
        !index.metadata.hasSampleNames()) {
        std::cerr << "Error: GBWT lacks path/sample metadata; cannot name haplotypes.\n";
        return 1;
    }
    const gbwt::Metadata& meta = index.metadata;

    std::cerr << "Loading " << ri_path << " ..." << std::endl;
    gbwt::FastLocate rindex;
    {
        std::ifstream in(ri_path, std::ios::binary);
        if (!in) { std::cerr << "Error: cannot open " << ri_path << "\n"; return 1; }
        rindex.load(in);
    }
    rindex.setGBWT(index);
    std::cerr << "  loaded in " << std::fixed << std::setprecision(1)
              << secs(t0, clk::now()) << "s" << std::endl;

    // ------------------------------------------------- haplotype numbering
    const size_t n_paths = meta.paths();
    std::vector<uint32_t> path_hap(n_paths, 0);
    std::vector<std::string> hap_names;
    {
        std::unordered_map<std::string, uint32_t> seen;
        for (size_t p = 0; p < n_paths; ++p) {
            gbwt::PathName pn = meta.path(p);
            std::string h = meta.sample(pn.sample) + "#" + std::to_string(pn.phase);
            auto it = seen.find(h);
            if (it == seen.end()) {
                uint32_t id = static_cast<uint32_t>(hap_names.size());
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
    std::cerr << "Haplotypes: " << with_commas(H)
              << "   paths: " << with_commas(n_paths)
              << "   bin: " << with_commas(bin_size) << " bp"
              << "   merge-gap: " << with_commas(merge_gap) << " bp" << std::endl;
    if (n_paths > 10000000) {
        std::cerr << "WARNING: " << with_commas(n_paths) << " paths. This looks like a"
                  << " frequency-filtered graph.\n         B2 entries scale with path"
                  << " count; expect ~10^11 entries. Use the non-filtered graph.\n";
    }

    const size_t max_node = static_cast<size_t>(gbz.graph.max_node_id());

    // Node length lookup: one dense array, reused by every path walk. Going
    // through get_handle()/get_length() per node was a large part of the cost.
    // Built BEFORE phase A so that pass can total each path's bases from it,
    // which is all Table 1 needs — no separate walk over every path.
    std::cerr << "Node lengths ("
              << human_bytes(static_cast<double>(max_node + 1) * 4.0) << ") ..." << std::endl;
    std::vector<uint32_t> node_len(max_node + 1, 0);
    {
        auto t_len = clk::now();
        #pragma omp parallel for schedule(static)
        for (size_t nid = 1; nid <= max_node; ++nid) {
            if (!gbz.graph.has_node(static_cast<handlegraph::nid_t>(nid))) continue;
            node_len[nid] = static_cast<uint32_t>(
                gbz.graph.get_length(gbz.graph.get_handle(static_cast<handlegraph::nid_t>(nid), false)));
        }
        std::cerr << "  node lengths in " << secs(t_len, clk::now()) << "s" << std::endl;
    }

    // ------------------------------------------- Phase A: node -> hap mask
    const size_t mask_entries = (max_node + 1) * WORDS;
    std::cerr << "Phase A: node->haplotype masks ("
              << human_bytes(static_cast<double>(mask_entries) * 8.0) << ") ..." << std::endl;
    std::vector<uint64_t> nodemask;
    try {
        nodemask.assign(mask_entries, 0);
    } catch (const std::bad_alloc&) {
        std::cerr << "Error: cannot allocate node mask ("
                  << human_bytes(static_cast<double>(mask_entries) * 8.0)
                  << "). Reduce haplotype count or add RAM.\n";
        return 1;
    }

    auto t_a = clk::now();
    // Total bases per path, accumulated in the same pass. This is the only
    // per-path quantity Table 1 needs, so building it here costs nothing.
    std::vector<size_t> path_bp(n_paths, 0);
    #pragma omp parallel for schedule(dynamic, 8)
    for (size_t p = 0; p < n_paths; ++p) {
        const uint32_t h = path_hap[p];
        const size_t word = h / 64;
        const uint64_t bit = 1ULL << (h % 64);
        gbwt::vector_type nodes = index.extract(gbwt::Path::encode(p, false));
        size_t total = 0;
        for (gbwt::node_type node : nodes) {
            if (node == gbwt::ENDMARKER) break;
            const size_t nid = static_cast<size_t>(gbwt::Node::id(node));
            if (nid > max_node) continue;
            total += node_len[nid];
            uint64_t* slot = &nodemask[nid * WORDS + word];
            // Paths of different haplotypes share nodes, so this OR races.
            if ((__atomic_load_n(slot, __ATOMIC_RELAXED) & bit) == 0) {
                __atomic_fetch_or(slot, bit, __ATOMIC_RELAXED);
            }
        }
        path_bp[p] = total;
    }
    std::cerr << "  phase A in " << secs(t_a, clk::now()) << "s" << std::endl;

    // ------------------------------------------------------ Table 1
    // Same construction as build_translation_tables: one subpath per GBWT path,
    // keyed by "sample#phase#contig", positioned at PathName::count. Added in
    // path_id order because add_subpath expects ascending subpath_start per name.
    if (build_t1) {
        auto t_t1 = clk::now();
        std::cerr << "Table 1: " << with_commas(n_paths) << " paths ..." << std::endl;
        panindexer::TranslationTable1 table1;
        size_t added = 0, skipped = 0;
        for (size_t p = 0; p < n_paths; ++p) {
            if (path_bp[p] == 0) { ++skipped; continue; }
            gbwt::PathName pn = meta.path(p);
            const std::string sample = (meta.samples() > 0 && pn.sample < meta.samples())
                ? meta.sample(pn.sample) : std::to_string(pn.sample);
            const std::string contig = (meta.contigs() > 0 && pn.contig < meta.contigs())
                ? meta.contig(pn.contig) : std::to_string(pn.contig);
            table1.add_subpath(sample + "#" + std::to_string(pn.phase) + "#" + contig,
                               p, static_cast<size_t>(pn.count), path_bp[p]);
            ++added;
        }
        {
            std::ofstream out(t1_path, std::ios::binary);
            if (!out) { std::cerr << "Error: cannot write " << t1_path << "\n"; return 1; }
            table1.serialize(out);
        }
        std::cerr << "  " << with_commas(added) << " subpaths, "
                  << with_commas(skipped) << " skipped (length 0), "
                  << with_commas(table1.num_names()) << " named paths, in "
                  << secs(t_t1, clk::now()) << "s -> " << t1_path << std::endl;
    } else {
        std::cerr << "Table 1: skipped (--no-t1)" << std::endl;
    }

    // ----------------------------- Phases B+C: coverage bins -> merged runs
    auto t_b = clk::now();
    std::cerr << "Phases B/C: coverage bins, run merging, target probes ..." << std::endl;

    const size_t merge_gap_bins = merge_gap / bin_size;
    const size_t n_src = (max_paths > 0 && max_paths < n_paths) ? max_paths : n_paths;

    std::vector<std::vector<Row>> per_thread;
#ifdef _OPENMP
    per_thread.resize(omp_get_max_threads());
#else
    per_thread.resize(1);
#endif

    std::atomic<size_t> probes_done{0};
    std::atomic<size_t> probe_failures{0};
    std::atomic<size_t> runs_dropped{0};
    std::atomic<size_t> widened_runs{0};
    std::atomic<size_t> runs_split{0};
    // Why runs were abandoned. Only the first is benign.
    std::atomic<size_t> fail_no_node_on_hap{0}, fail_first_unlocated{0},
                        fail_last_unlocated{0}, fail_no_shared_path{0},
                        fail_span_capped{0};
    std::atomic<size_t> src_done{0};

    #pragma omp parallel
    {
#ifdef _OPENMP
        std::vector<Row>& out = per_thread[omp_get_thread_num()];
#else
        std::vector<Row>& out = per_thread[0];
#endif
        // Per-source-path scratch, reused across paths to avoid reallocation.
        std::vector<uint64_t> bins;
        std::vector<uint32_t> run_start(H, 0), run_last(H, 0);
        std::vector<uint8_t>  run_open(H, 0);
        // The source path as walked: node id and base offset, ascending. Keeping
        // the whole walk (rather than one node per bin) is what lets a run's
        // anchors be the first and last node ACTUALLY on the target haplotype,
        // instead of whichever node happened to start a bin.
        std::vector<uint32_t> src_nid, src_off;

        // Collect this node's visits by paths belonging to haplotype `h`.
        // Only EVEN sequence ids are accepted: a bidirectional GBWT stores each
        // path twice, and on the reverse copy seqOffset is measured along the
        // reversed sequence. Restricting to the forward copy keeps every offset
        // in one coordinate system, and probing both node orientations still
        // finds targets that traverse the node the other way round.
        auto probe = [&](uint32_t nid, uint32_t h, std::vector<Visit>& hits) {
            hits.clear();
            for (int orient = 0; orient < 2; ++orient) {
                gbwt::node_type gn = gbwt::Node::encode(
                    static_cast<gbwt::size_type>(nid), orient == 1);
                std::vector<gbwt::size_type> sa = rindex.decompressSA(gn);
                ++probes_done;
                for (gbwt::size_type v : sa) {
                    const gbwt::size_type sid = rindex.seqId(v);
                    if (sid % 2 != 0) continue;
                    const size_t pid = static_cast<size_t>(sid) / 2;
                    if (pid < n_paths && path_hap[pid] == h) {
                        hits.push_back(Visit{pid, static_cast<uint64_t>(rindex.seqOffset(v))});
                    }
                }
            }
        };
        std::vector<Visit> hits_first, hits_last, hits_mid, hits_probe;
        std::vector<size_t> probe_paths[3];
        std::vector<size_t> seen_paths;

        #pragma omp for schedule(dynamic, 1)
        for (size_t sp = 0; sp < n_src; ++sp) {
            // ---- Phase B: one walk. Record the node list and bin coverage.
            gbwt::vector_type nodes = index.extract(gbwt::Path::encode(sp, false));
            src_nid.clear();
            src_off.clear();
            size_t path_len = 0;
            for (gbwt::node_type node : nodes) {
                if (node == gbwt::ENDMARKER) break;
                const size_t nid = static_cast<size_t>(gbwt::Node::id(node));
                if (nid > max_node) continue;
                src_nid.push_back(static_cast<uint32_t>(nid));
                src_off.push_back(static_cast<uint32_t>(path_len));
                path_len += node_len[nid];
            }
            if (src_nid.empty() || path_len == 0) { ++src_done; continue; }

            const size_t nbins = (path_len + bin_size - 1) / bin_size;
            bins.assign(nbins * WORDS, 0);
            for (size_t i = 0; i < src_nid.size(); ++i) {
                const uint64_t* m = &nodemask[static_cast<size_t>(src_nid[i]) * WORDS];
                // Mark EVERY bin the node spans, not just the one holding its
                // start. A node longer than bin_size would otherwise leave its
                // interior uncovered, which breaks runs apart and loses real
                // homology — the one error direction this table must not have.
                const size_t nlen = node_len[src_nid[i]];
                const size_t b_lo = src_off[i] / bin_size;
                const size_t b_hi = (nlen > 0)
                    ? (src_off[i] + nlen - 1) / bin_size : b_lo;
                for (size_t b = b_lo; b <= b_hi && b < nbins; ++b) {
                    uint64_t* dst = &bins[b * WORDS];
                    for (size_t w = 0; w < WORDS; ++w) dst[w] |= m[w];
                }
            }

            // Emit one row for the bin run [b0, b1] on haplotype h.
            // Recursive: a run with no plausible anchor pairing is usually a run
            // that spans TWO contigs of the target haplotype (assembly breaks are
            // common, and --merge-gap happily merges across them). One target
            // path id cannot describe it, so the run is split and each half
            // resolved on its own rather than forced onto a single path.
            std::function<void(uint32_t, size_t, size_t, int)> emit_run =
                [&](uint32_t h, size_t b0, size_t b1, int depth) {
                // The stock builder skips pairs within one haplotype; match it
                // by default so the two tables cover the same pair set.
                if (!allow_same_hap && h == path_hap[sp]) return;
                size_t src_start = b0 * bin_size;
                size_t src_end = (b1 + 1) * bin_size;
                if (src_end > path_len) src_end = path_len;
                if (src_end <= src_start) return;
                // Filter on actual span, not bin count: a run can occupy two
                // bins yet be far shorter than the threshold.
                if (src_end - src_start < min_run_bp) { ++runs_dropped; return; }

                // Anchors: first and last node in the run that is really on h.
                // Bin coverage means "some node in this bin is on h", so the
                // run's outermost nodes need not be.
                const size_t word = h / 64;
                const uint64_t bit = 1ULL << (h % 64);
                auto on_h = [&](size_t i) {
                    return (nodemask[static_cast<size_t>(src_nid[i]) * WORDS + word] & bit) != 0;
                };
                size_t lo = std::lower_bound(src_off.begin(), src_off.end(),
                                             static_cast<uint32_t>(src_start))
                            - src_off.begin();
                size_t hi = std::lower_bound(src_off.begin(), src_off.end(),
                                             static_cast<uint32_t>(src_end))
                            - src_off.begin();
                if (hi > src_nid.size()) hi = src_nid.size();
                size_t i_first = hi, i_last = hi;
                for (size_t i = lo; i < hi; ++i) { if (on_h(i)) { i_first = i; break; } }
                if (i_first >= hi) {
                    // Not a loss: no node in this sub-run is on the haplotype at
                    // all, so there is no homology here to record. Splitting
                    // produces these routinely when the on-h nodes all land in
                    // the other half.
                    ++fail_no_node_on_hap; ++probe_failures; return;
                }
                for (size_t i = hi; i > i_first; --i) {
                    if (on_h(i - 1)) { i_last = i - 1; break; }
                }
                if (i_last >= hi) i_last = i_first;

                if (!want_tgt_intervals) {
                    // ROUTING ONLY. Collect every distinct target path the run
                    // touches and emit one segment per path. There is no pairing
                    // and therefore no way to fail, to misplace a block, or to
                    // need a split: a run spanning two contigs simply yields two
                    // entries, each correct.
                    seen_paths.clear();
                    auto add_from = [&](const std::vector<Visit>& hits) {
                        for (const Visit& v : hits) {
                            if (v.path_id == sp) continue;   // never map to itself
                            if (std::find(seen_paths.begin(), seen_paths.end(),
                                          v.path_id) == seen_paths.end()) {
                                seen_paths.push_back(v.path_id);
                            }
                        }
                    };
                    // Sample across the run so every contig it crosses is seen.
                    const size_t span_i = i_last - i_first;
                    for (int frac = 0; frac <= 4; ++frac) {
                        const size_t want_i = i_first + (span_i * frac) / 4;
                        for (size_t i = want_i; i <= i_last; ++i) {
                            if (on_h(i)) { probe(src_nid[i], h, hits_probe);
                                           add_from(hits_probe); break; }
                        }
                    }
                    if (seen_paths.empty()) { ++fail_no_shared_path; ++probe_failures; return; }
                    for (size_t pid : seen_paths) {
                        Row r;
                        r.src_path_id = sp;
                        r.tgt_hap     = h;
                        r.src_start   = src_start;
                        r.src_end     = src_end;
                        r.tgt_path_id = pid;
                        r.resolved    = true;   // nothing for phase D to do
                        out.push_back(r);
                    }
                    return;
                }

                probe(src_nid[i_first], h, hits_first);
                if (hits_first.empty()) { ++fail_first_unlocated; ++probe_failures; return; }
                if (i_last != i_first) {
                    probe(src_nid[i_last], h, hits_last);
                } else {
                    hits_last = hits_first;
                }
                if (hits_last.empty()) { ++fail_last_unlocated; ++probe_failures; return; }

                // A third anchor from the MIDDLE of the run. Two anchors alone
                // cannot tell a correct pairing from a pair of unrelated repeat
                // copies that happen to share a path: both look like "same path,
                // some displacement". Requiring the midpoint to fall between them
                // on that same path rejects the pairing that lands the block in
                // the wrong place — the failure that produced blocks of the right
                // LENGTH but 90% of their shared nodes outside them.
                hits_mid.clear();
                for (int i = 0; i < 3; ++i) probe_paths[i].clear();
                if (i_last > i_first + 1) {
                    // Three interior probes, not one. A single midpoint can miss
                    // the contig boundary entirely (it may land in the half that
                    // agrees), leaving a two-contig run undetected. Sampling at
                    // 1/4, 1/2 and 3/4 catches a break anywhere in the middle.
                    const size_t span = i_last - i_first;
                    for (int frac = 1; frac <= 3; ++frac) {
                        const size_t want_i = i_first + (span * frac) / 4;
                        for (size_t i = want_i; i < i_last; ++i) {
                            if (on_h(i)) {
                                probe(src_nid[i], h, hits_probe);
                                // Keep each probe's paths SEPARATE. Pooling them
                                // defeats the purpose: across a contig boundary
                                // some probe lands on each side, so a pooled set
                                // "confirms" whichever pairing was picked.
                                probe_paths[frac - 1].clear();
                                for (const Visit& v : hits_probe)
                                    probe_paths[frac - 1].push_back(v.path_id);
                                hits_mid.insert(hits_mid.end(),
                                                hits_probe.begin(), hits_probe.end());
                                break;
                            }
                        }
                    }
                }

                // Pair the two anchors. Both must land on the SAME target path,
                // and among candidate pairings prefer the one whose node-count
                // separation best matches the source's — that is what discards
                // an unrelated paralogous copy when a node repeats.
                const long long want =
                    static_cast<long long>(i_last) - static_cast<long long>(i_first);
                // A pairing must span a comparable number of NODES to the source.
                // Without this floor any same-path pairing was accepted, however
                // absurd — a repeat copy giving got ~ 0 across a source run of
                // millions of nodes produced target intervals of a single base.
                const long long lo_ok = want / 4, hi_ok = want * 4;

                bool have = false, have_strict = false;
                size_t best_pid = 0;
                uint64_t best_a = 0, best_b = 0;
                long long best_cost = 0;
                // Fallback if nothing passes: the WIDEST same-path pairing.
                // Widening is the safe direction (a too-large block costs a
                // traversal that returns nothing); narrowing loses real hits.
                bool have_wide = false;
                size_t wide_pid = 0; uint64_t wide_a = 0, wide_b = 0;
                long long wide_span = -1;

                for (const Visit& a : hits_first) {
                    for (const Visit& b : hits_last) {
                        if (a.path_id != b.path_id) continue;
                        if (a.path_id == sp) continue;   // never map a path to itself
                        // seqOffset decreases along the path, so this signed
                        // difference is the target's node displacement.
                        const long long got = static_cast<long long>(a.sa_off)
                                            - static_cast<long long>(b.sa_off);
                        const long long mag = std::llabs(got);
                        if (mag > wide_span) {
                            have_wide = true; wide_span = mag;
                            wide_pid = a.path_id; wide_a = a.sa_off; wide_b = b.sa_off;
                        }
                        if (want > 0 && (mag < lo_ok || mag > hi_ok)) continue;

                        // Midpoint confirmation is PATH membership, not positional
                        // containment. Requiring the midpoint to fall strictly
                        // between the anchors rejected far too much — a repeat
                        // whose first located occurrence sits outside the span is
                        // still on the right contig — and the resulting splits ran
                        // away (36 M of them, and 5 M runs abandoned outright).
                        // What actually indicates a two-contig run is the midpoint
                        // being on a DIFFERENT path, which is checked below.
                        bool mid_ok = hits_mid.empty();
                        if (!mid_ok) {
                            for (const Visit& m : hits_mid) {
                                if (m.path_id == a.path_id) { mid_ok = true; break; }
                            }
                        }
                        const long long cost = std::llabs(got - want);
                        // Any midpoint-confirmed pairing beats every unconfirmed
                        // one, regardless of cost.
                        if (mid_ok) {
                            if (!have_strict || cost < best_cost) {
                                have = have_strict = true; best_cost = cost;
                                best_pid = a.path_id; best_a = a.sa_off; best_b = b.sa_off;
                            }
                        } else if (!have_strict && (!have || cost < best_cost)) {
                            have = true; best_cost = cost;
                            best_pid = a.path_id; best_a = a.sa_off; best_b = b.sa_off;
                        }
                    }
                }
                // Split only on POSITIVE evidence of a two-contig run: the
                // midpoint is on this haplotype, but on none of the paths any
                // candidate pairing used. Splitting merely because the midpoint
                // was unconfirmed fired on nearly every run.
                // Positive evidence of a multi-contig run: two interior probes
                // whose path sets are DISJOINT. One probe cannot show this, and a
                // probe merely agreeing with the chosen pairing cannot rule it
                // out — a run crossing a boundary has probes on both sides, so
                // some probe always agrees.
                bool split_wanted = false;
                for (int i = 0; i < 3 && !split_wanted; ++i) {
                    if (probe_paths[i].empty()) continue;
                    for (int j = i + 1; j < 3 && !split_wanted; ++j) {
                        if (probe_paths[j].empty()) continue;
                        bool overlap = false;
                        for (size_t a : probe_paths[i]) {
                            for (size_t b : probe_paths[j]) {
                                if (a == b) { overlap = true; break; }
                            }
                            if (overlap) break;
                        }
                        if (!overlap) split_wanted = true;
                    }
                }
                // Also split when the chosen pairing is on a path no interior
                // probe saw at all.
                if (have && !split_wanted && !hits_mid.empty()) {
                    bool mid_on_chosen = false;
                    for (const Visit& m : hits_mid) {
                        if (m.path_id == best_pid) { mid_on_chosen = true; break; }
                    }
                    split_wanted = !mid_on_chosen;
                }
                if (!have || split_wanted) {
                    if (depth < MAX_SPLIT_DEPTH && b1 > b0) {
                        const size_t mid = b0 + (b1 - b0) / 2;
                        ++runs_split;
                        emit_run(h, b0, mid, depth + 1);
                        emit_run(h, mid + 1, b1, depth + 1);
                        return;
                    }
                    // Never abandon a run that has ANY same-path pairing: a
                    // missing entry makes translate() report "this region does not
                    // exist", which is the false negative the whole table exists to
                    // avoid. An over-wide block only costs a traversal.
                    // Cap the fallback. An unbounded "widest pairing" produced
                    // blocks up to 36x their source span, which is safe but
                    // useless as a bound. Beyond the cap, prefer no entry over a
                    // meaningless one only if we truly cannot narrow it.
                    if (!have_wide) {
                        // Anchors located, but never together on one target path.
                        ++fail_no_shared_path; ++probe_failures; return;
                    }
                    if (want > 0 && wide_span > hi_ok) {
                        ++fail_span_capped; ++probe_failures; return;
                    }
                    ++widened_runs;
                    best_pid = wide_pid; best_a = wide_a; best_b = wide_b;
                }

                Row r;
                r.src_path_id = sp;
                r.tgt_hap     = h;
                r.src_start   = src_start;
                r.src_end     = src_end;
                r.tgt_path_id = best_pid;
                r.sa_first    = best_a;
                r.sa_last     = best_b;
                out.push_back(r);
            };

            // ---- Phase C: runs via XOR delta between adjacent bins, with gap
            // merging applied inline so no intermediate run list is built.
            std::fill(run_open.begin(), run_open.end(), 0);
            for (size_t b = 0; b < nbins; ++b) {
                for (size_t w = 0; w < WORDS; ++w) {
                    const uint64_t cur  = bins[b * WORDS + w];
                    const uint64_t prev = (b == 0) ? 0ULL : bins[(b - 1) * WORDS + w];
                    uint64_t diff = cur ^ prev;
                    while (diff) {
                        const int t = __builtin_ctzll(diff);
                        diff &= diff - 1;
                        const uint32_t h = static_cast<uint32_t>(w * 64 + t);
                        if (h >= H) continue;
                        if (cur & (1ULL << t)) {           // coverage starts
                            if (run_open[h] && b - run_last[h] - 1 <= merge_gap_bins) {
                                // close enough to the previous stretch: same run
                            } else {
                                if (run_open[h]) emit_run(h, run_start[h], run_last[h], 0);
                                run_start[h] = static_cast<uint32_t>(b);
                            }
                            run_open[h] = 1;
                            run_last[h] = static_cast<uint32_t>(b);
                        } else if (run_open[h]) {          // coverage stops
                            run_last[h] = static_cast<uint32_t>(b - 1);
                        }
                    }
                }
            }
            // Haplotypes still covered in the final bin run to the path end.
            for (size_t w = 0; w < WORDS; ++w) {
                uint64_t bits = bins[(nbins - 1) * WORDS + w];
                while (bits) {
                    const int t = __builtin_ctzll(bits);
                    bits &= bits - 1;
                    const uint32_t h = static_cast<uint32_t>(w * 64 + t);
                    if (h < H && run_open[h]) run_last[h] = static_cast<uint32_t>(nbins - 1);
                }
            }
            for (uint32_t h = 0; h < H; ++h) {
                if (run_open[h]) emit_run(h, run_start[h], run_last[h], 0);
            }

            if (progress) {
                const size_t d = ++src_done;
                if (d % progress_every == 0) {
                    const double el = secs(t_b, clk::now());
                    const double rate = (el > 0) ? static_cast<double>(d) / el : 0.0;
                    const double eta = (rate > 0) ? (n_src - d) / rate : 0.0;
                    #pragma omp critical
                    std::cerr << "  B/C " << with_commas(d) << "/" << with_commas(n_src)
                              << " paths  " << std::fixed << std::setprecision(0)
                              << el << "s elapsed  eta " << eta << "s  "
                              << with_commas(probes_done.load()) << " probes" << std::endl;
                }
            } else {
                ++src_done;
            }
        }
    }
    std::cerr << "  phases B/C in " << secs(t_b, clk::now()) << "s" << std::endl;

    // ------------------------------------------------------------ gather rows
    std::vector<Row> rows;
    {
        size_t total = 0;
        for (const auto& v : per_thread) total += v.size();
        rows.reserve(total);
        for (auto& v : per_thread) {
            rows.insert(rows.end(), v.begin(), v.end());
            std::vector<Row>().swap(v);   // release as we go
        }
    }
    std::cerr << "Rows: " << with_commas(rows.size()) << std::endl;

    // ------------------------------- Phase D: SA offsets -> target base coords
    auto t_d = clk::now();
    std::cerr << (want_tgt_intervals
                  ? "Phase D: resolving target coordinates ..."
                  : "Phase D: skipped (routing-only table)") << std::endl;
    // Group by target path so each one is walked exactly once.
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
        return a.tgt_path_id < b.tgt_path_id;
    });
    std::vector<size_t> group_begin;
    for (size_t i = 0; i < rows.size(); ) {
        group_begin.push_back(i);
        const size_t pid = rows[i].tgt_path_id;
        while (i < rows.size() && rows[i].tgt_path_id == pid) ++i;
    }
    group_begin.push_back(rows.size());

    std::atomic<size_t> unresolved{0};
    // OpenMP requires a plain relational test on the loop variable, so the
    // group count is hoisted rather than written as `g + 1 < size()`.
    const size_t n_groups = group_begin.empty() ? 0 : group_begin.size() - 1;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t g = 0; g < (want_tgt_intervals ? n_groups : 0); ++g) {
        const size_t lo = group_begin[g], hi = group_begin[g + 1];
        const size_t pid = rows[lo].tgt_path_id;

        // Base offset of every node position along the target path.
        gbwt::vector_type tn = index.extract(gbwt::Path::encode(pid, false));
        std::vector<uint32_t> off;
        off.reserve(tn.size());
        std::vector<uint32_t> len;
        len.reserve(tn.size());
        size_t acc = 0;
        for (gbwt::node_type node : tn) {
            if (node == gbwt::ENDMARKER) break;
            const size_t nid = static_cast<size_t>(gbwt::Node::id(node));
            const uint32_t l = (nid <= max_node) ? node_len[nid] : 0u;
            off.push_back(static_cast<uint32_t>(acc));
            len.push_back(l);
            acc += l;
        }
        const size_t L = off.size();
        if (L == 0) { unresolved += (hi - lo); continue; }

        // seqOffset counts from the END of the path (larger == earlier), so a
        // node position index is L - 1 - seqOffset.
        auto to_index = [&](uint64_t sa) -> size_t {
            if (sa >= L) return L;                       // out of range
            return L - 1 - static_cast<size_t>(sa);
        };

        for (size_t i = lo; i < hi; ++i) {
            Row& r = rows[i];
            const size_t ia = to_index(r.sa_first);
            const size_t ib = to_index(r.sa_last);
            if (ia >= L || ib >= L) { ++unresolved; continue; }
            const size_t a_start = off[ia], a_end = off[ia] + len[ia];
            const size_t b_start = off[ib], b_end = off[ib] + len[ib];
            // Store ascending regardless of direction; the flag records which
            // end of [tgt_start, tgt_end) the source's start corresponds to.
            r.tgt_start   = std::min(a_start, b_start);
            r.tgt_end     = std::max(a_end, b_end);
            r.tgt_reverse = (ib < ia);
            r.resolved    = true;
        }
    }
    std::cerr << "  phase D in " << secs(t_d, clk::now()) << "s" << std::endl;

    // ------------------------------------------------------ Phase E: write out
    auto t_w = clk::now();
    std::cerr << "Phase E: sorting and writing ..." << std::endl;
    // Key order must match TranslationTable2's: src_path_id, then haplotype
    // NAME (not the interned id, whose order is first-seen).
    std::sort(rows.begin(), rows.end(),
              [&](const Row& a, const Row& b) {
        if (a.src_path_id != b.src_path_id) return a.src_path_id < b.src_path_id;
        if (a.tgt_hap != b.tgt_hap) {
            return hap_names[a.tgt_hap] < hap_names[b.tgt_hap];
        }
        if (a.src_start != b.src_start) return a.src_start < b.src_start;
        // Tiebreaker so the result cannot depend on which thread produced a
        // row: std::sort is not stable and rows arrive in thread order.
        return a.tgt_path_id < b.tgt_path_id;
    });

    size_t written_keys = 0, written_segs = 0, skipped = 0;
    size_t under_100 = 0, under_1k = 0, under_10k = 0;
    {
        std::ofstream out(out_path, std::ios::binary);
        if (!out) { std::cerr << "Error: cannot write " << out_path << "\n"; return 1; }
        panindexer::TranslationTable2Writer writer(out, hap_names, want_tgt_intervals);
        bool key_open = false;
        size_t cur_src = 0;
        uint32_t cur_hap = 0;
        for (const Row& r : rows) {
            if (!r.resolved) { ++skipped; continue; }
            const size_t span = r.src_end - r.src_start;
            if (span < 100) ++under_100;
            if (span < 1000) ++under_1k;
            if (span < 10000) ++under_10k;
            if (!key_open || r.src_path_id != cur_src || r.tgt_hap != cur_hap) {
                if (!writer.begin_key(r.src_path_id, hap_names[r.tgt_hap])) {
                    std::cerr << "Error: key order violated at src_path "
                              << r.src_path_id << " hap " << hap_names[r.tgt_hap] << "\n";
                    return 1;
                }
                key_open = true;
                cur_src = r.src_path_id;
                cur_hap = r.tgt_hap;
                ++written_keys;
            }
            panindexer::IntervalMapping m;
            m.src_start   = r.src_start;
            m.src_end     = r.src_end;
            m.tgt_path_id = r.tgt_path_id;
            m.tgt_start   = r.tgt_start;
            m.tgt_end     = r.tgt_end;
            m.tgt_reverse = r.tgt_reverse;
            writer.add_segment(m);
            ++written_segs;
        }
        writer.finish();
    }
    std::cerr << "  written in " << secs(t_w, clk::now()) << "s" << std::endl;

    std::cout << "\n=============== B2 Table 2 ===============\n"
              << "  keys (src_path, tgt_haplotype): " << with_commas(written_keys) << "\n"
              << (want_tgt_intervals
                  ? "  segments (with target coords):  "
                  : "  segments (routing only):        ")
              << with_commas(written_segs) << "\n"
              << "  unresolved target coords:       " << with_commas(skipped) << "\n"
              << "  runs with no target anchor:     " << with_commas(probe_failures.load()) << "\n"
              << "  runs dropped by --min-run-bp:   " << with_commas(runs_dropped.load()) << "\n"
              << "  runs widened (no plausible pair):" << with_commas(widened_runs.load()) << "\n"
              << "  runs split (spanned 2 targets):  " << with_commas(runs_split.load()) << "\n"
              << "  abandoned, by cause:\n"
              << "    no node on that haplotype:    " << with_commas(fail_no_node_on_hap.load())
              << "   (benign: nothing to record)\n"
              << "    first anchor unlocated:       " << with_commas(fail_first_unlocated.load()) << "\n"
              << "    last anchor unlocated:        " << with_commas(fail_last_unlocated.load()) << "\n"
              << "    anchors never on one path:    " << with_commas(fail_no_shared_path.load())
              << "   <- LOST homology\n"
              << "    span exceeded 4x cap:         " << with_commas(fail_span_capped.load())
              << "   <- LOST homology\n"
              << "  decompressSA probes:            " << with_commas(probes_done.load()) << "\n"
              << "  segment source spans < 100 bp:  " << with_commas(under_100) << "\n"
              << "                       < 1 kb:    " << with_commas(under_1k) << "\n"
              << "                       < 10 kb:   " << with_commas(under_10k) << "\n"
              << "  table 1:                        "
              << (build_t1 ? t1_path : std::string("(skipped)")) << "\n"
              << "  table 2:                        " << out_path << "\n"
              << "  total time:                     " << std::fixed << std::setprecision(1)
              << secs(t0, clk::now()) << "s\n";
    std::cout.flush();
    std::cerr.flush();
    // Skip static destruction of the loaded GBZ/sdsl structures: it traps on
    // some toolchains after all work is complete and would turn a successful
    // build into a nonzero exit status.
    std::_Exit(0);
}
