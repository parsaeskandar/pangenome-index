/*
 * validate_table2.cpp
 *
 * Check that a Translation Table 2 tells the TRUTH about the graph.
 *
 * WHAT THIS TESTS THAT tags_check DOES NOT
 * ----------------------------------------
 * tags_check verifies Table 2 is self-consistent: that lookup() returns what
 * segments() says it should. That exercises the lookup code and says nothing
 * about whether the stored intervals correspond to anything real — a table full
 * of invented coordinates passes it.
 *
 * This validates the claim itself against the GBWT. For an entry
 *
 *     src path X interval [a, b)  ->  tgt path Y interval [c, d)
 *
 * the promise is: every graph node the source visits inside [a, b) that the
 * target also visits has one of its target occurrences inside [c, d). Block
 * boundaries are only ever WIDENED by the builder (bin granularity, gap
 * merging), so [c, d) must be a SUPERSET of the true correspondence. A node
 * landing outside is a false negative — the error class that silently loses a
 * real translation, and the one thing the table must never do.
 *
 * The converse is deliberately NOT an error: [c, d) may contain target
 * sequence with no source counterpart (an insertion in the target, or plain
 * over-widening). That costs a traversal that returns nothing, which is safe.
 *
 * WHAT IT REPORTS
 *   entries checked      how many (src, tgt) blocks were sampled
 *   clean entries        every shared node landed inside the target block
 *   shared nodes         nodes visited by both, inside the source block
 *   escaped nodes        shared nodes whose target occurrences ALL fell outside
 *                        [c, d) -- these are the real failures
 *   orientation mismatch entries whose tgt_reverse flag disagrees with the
 *                        direction the shared nodes actually run
 *   length ratio         |d - c| / |b - a|, summarized; wildly non-1 ratios on
 *                        a syntenic block suggest a mispaired paralog
 *
 * Exit status is nonzero if any escaped nodes were found, so it can gate a
 * build.
 *
 * Usage:
 *   validate_table2 <graph.gbz> <table.t2> [options]
 *
 *     --trials N     entries to sample (default 1000; 0 = every entry)
 *     --seed N       RNG seed (default 42) so runs are reproducible
 *     --threads N    worker threads (default: all)
 *     --verbose      print each failing entry
 */

#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>
#include <sdsl/simple_sds.hpp>

#include "pangenome_index/translation_tables.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " <graph.gbz> <table.t2> [options]\n\n"
        << "Validates that Table 2's stored intervals correspond to real shared\n"
        << "graph nodes, rather than merely being self-consistent.\n\n"
        << "  --trials N   entries to sample (default 1000; 0 = all)\n"
        << "  --seed N     RNG seed (default 42)\n"
        << "  --threads N  worker threads (default: all)\n"
        << "  --verbose    print each failing entry\n";
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

/// A path as base offsets: node id at position i covers [off[i], off[i]+len[i]).
struct PathWalk {
    std::vector<uint32_t> nid, off, len;
    /// node id -> every position index where the path visits it
    std::unordered_map<uint32_t, std::vector<uint32_t>> where;
};

PathWalk walk_path(const gbwt::GBWT& index, const gbwtgraph::GBWTGraph& graph,
                   size_t path_id) {
    PathWalk w;
    gbwt::vector_type nodes = index.extract(gbwt::Path::encode(path_id, false));
    size_t acc = 0;
    for (gbwt::node_type node : nodes) {
        if (node == gbwt::ENDMARKER) break;
        const auto id = gbwt::Node::id(node);
        const uint32_t l = static_cast<uint32_t>(
            graph.get_length(graph.get_handle(id, gbwt::Node::is_reverse(node))));
        w.where[static_cast<uint32_t>(id)].push_back(static_cast<uint32_t>(w.nid.size()));
        w.nid.push_back(static_cast<uint32_t>(id));
        w.off.push_back(static_cast<uint32_t>(acc));
        w.len.push_back(l);
        acc += l;
    }
    return w;
}

struct EntryResult {
    uint64_t shared = 0;
    uint64_t escaped = 0;
    bool orientation_ok = true;
    double  ratio = 0.0;
    bool    counted = false;
};

} // namespace

int main(int argc, char** argv) {
    std::string gbz_path, t2_path;
    size_t trials = 1000, seed = 42;
    int threads = 0;
    bool verbose = false;
    {
        std::vector<std::string> pos;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
            else if (a == "--trials" && i + 1 < argc) trials = std::stoull(argv[++i]);
            else if (a == "--seed" && i + 1 < argc) seed = std::stoull(argv[++i]);
            else if (a == "--threads" && i + 1 < argc) threads = std::stoi(argv[++i]);
            else if (a == "--verbose") verbose = true;
            else if (!a.empty() && a[0] == '-') { usage(argv[0]); return 1; }
            else pos.push_back(a);
        }
        if (pos.size() != 2) { usage(argv[0]); return 1; }
        gbz_path = pos[0]; t2_path = pos[1];
    }
#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
#endif

    std::cerr << "Loading " << gbz_path << " ..." << std::endl;
    gbwtgraph::GBZ gbz;
    sdsl::simple_sds::load_from(gbz, gbz_path);
    const gbwt::GBWT& index = gbz.index;

    std::cerr << "Loading " << t2_path << " ..." << std::endl;
    panindexer::TranslationTable2 t2;
    {
        std::ifstream in(t2_path, std::ios::binary);
        if (!in) { std::cerr << "Error: cannot open " << t2_path << "\n"; return 1; }
        t2.load(in);
    }
    std::cerr << "  " << with_commas(t2.num_entries()) << " keys, "
              << with_commas(t2.total_segments()) << " segments, target coords: "
              << (t2.has_target_coords() ? "yes" : "NO") << std::endl;

    if (!t2.has_target_coords()) {
        std::cerr << "\nThis table stores no target intervals, so there is no\n"
                  << "coordinate claim to validate. Only a B2 table (file version 3,\n"
                  << "built by build_table2_b2) can be checked here.\n";
        return 2;
    }

    // Flatten to (src_path, hap, segment) so sampling is uniform over segments
    // and so the work can be grouped by source path — each path is walked once.
    struct Item { size_t src_path; std::string hap; panindexer::IntervalMapping seg; };
    std::vector<Item> items;
    for (const auto& key : t2.keys())
        for (const auto& seg : t2.segments(key.first, key.second))
            items.push_back(Item{key.first, key.second, seg});
    if (items.empty()) { std::cerr << "Table is empty.\n"; return 2; }

    if (trials > 0 && trials < items.size()) {
        std::mt19937_64 rng(seed);
        std::shuffle(items.begin(), items.end(), rng);
        items.resize(trials);
    }
    // Group by source path: walking a whole-chromosome path per entry would
    // dominate the runtime, and adjacent entries usually share a source.
    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        if (a.src_path != b.src_path) return a.src_path < b.src_path;
        return a.seg.tgt_path_id < b.seg.tgt_path_id;
    });
    std::cerr << "Validating " << with_commas(items.size()) << " entries ..." << std::endl;

    std::vector<EntryResult> results(items.size());
    std::atomic<size_t> done{0};

    // One index over item positions grouped by source path.
    std::vector<size_t> group_start;
    for (size_t i = 0; i < items.size(); ) {
        group_start.push_back(i);
        const size_t sp = items[i].src_path;
        while (i < items.size() && items[i].src_path == sp) ++i;
    }
    const size_t n_groups = group_start.size();
    group_start.push_back(items.size());

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t g = 0; g < n_groups; ++g) {
        const size_t lo = group_start[g], hi = group_start[g + 1];
        PathWalk src = walk_path(index, gbz.graph, items[lo].src_path);
        if (src.nid.empty()) { done += (hi - lo); continue; }

        // Target walks are cached per group: consecutive entries on one source
        // often name the same target contig.
        std::unordered_map<size_t, PathWalk> tgt_cache;

        for (size_t i = lo; i < hi; ++i) {
            const auto& seg = items[i].seg;
            auto it = tgt_cache.find(seg.tgt_path_id);
            if (it == tgt_cache.end()) {
                if (tgt_cache.size() > 8) tgt_cache.clear();   // bound memory
                it = tgt_cache.emplace(seg.tgt_path_id,
                        walk_path(index, gbz.graph, seg.tgt_path_id)).first;
            }
            const PathWalk& tgt = it->second;
            if (tgt.nid.empty()) { ++done; continue; }

            EntryResult r;
            r.counted = true;
            r.ratio = (seg.src_end > seg.src_start)
                ? double(seg.tgt_end - seg.tgt_start) / double(seg.src_end - seg.src_start)
                : 0.0;

            // Walk the source's node positions inside [src_start, src_end).
            int64_t first_ti = -1, last_ti = -1;
            for (size_t k = 0; k < src.nid.size(); ++k) {
                const size_t s0 = src.off[k], s1 = src.off[k] + src.len[k];
                if (s1 <= seg.src_start) continue;
                if (s0 >= seg.src_end) break;

                auto wt = tgt.where.find(src.nid[k]);
                if (wt == tgt.where.end()) continue;   // not shared: nothing promised
                ++r.shared;

                // Does ANY occurrence on the target fall inside [tgt_start, tgt_end)?
                bool inside = false;
                for (uint32_t ti : wt->second) {
                    const size_t t0 = tgt.off[ti], t1 = tgt.off[ti] + tgt.len[ti];
                    if (t1 > seg.tgt_start && t0 < seg.tgt_end) {
                        inside = true;
                        if (first_ti < 0) first_ti = ti;
                        last_ti = ti;
                        break;
                    }
                }
                if (!inside) ++r.escaped;
            }
            // Orientation: if the shared nodes run backwards along the target,
            // tgt_reverse must say so.
            if (first_ti >= 0 && last_ti >= 0 && first_ti != last_ti) {
                const bool actually_reverse = (last_ti < first_ti);
                r.orientation_ok = (actually_reverse == seg.tgt_reverse);
            }
            results[i] = r;
            ++done;
        }
    }

    uint64_t entries = 0, clean = 0, shared = 0, escaped = 0, bad_orient = 0, no_shared = 0;
    double ratio_sum = 0.0; uint64_t ratio_n = 0;
    double ratio_min = 1e18, ratio_max = 0.0;
    for (size_t i = 0; i < results.size(); ++i) {
        const EntryResult& r = results[i];
        if (!r.counted) continue;
        ++entries;
        shared += r.shared;
        escaped += r.escaped;
        if (r.shared == 0) ++no_shared;
        if (r.escaped == 0) ++clean;
        if (!r.orientation_ok) ++bad_orient;
        if (r.ratio > 0) {
            ratio_sum += r.ratio; ++ratio_n;
            ratio_min = std::min(ratio_min, r.ratio);
            ratio_max = std::max(ratio_max, r.ratio);
        }
        if (verbose && (r.escaped > 0 || !r.orientation_ok)) {
            const auto& it = items[i];
            #pragma omp critical
            std::cerr << "  FAIL src_path=" << it.src_path << " hap=" << it.hap
                      << " src[" << it.seg.src_start << "," << it.seg.src_end << ")"
                      << " -> path " << it.seg.tgt_path_id
                      << " tgt[" << it.seg.tgt_start << "," << it.seg.tgt_end << ")"
                      << (it.seg.tgt_reverse ? " REV" : "")
                      << "  shared=" << r.shared << " escaped=" << r.escaped
                      << (r.orientation_ok ? "" : " ORIENTATION") << std::endl;
        }
    }

    std::cout << "\n=============== Table 2 validation ===============\n"
              << "  entries checked:        " << with_commas(entries) << "\n"
              << "  clean (no escapes):     " << with_commas(clean)
              << (entries ? "  (" + std::to_string(100.0 * clean / entries).substr(0, 5) + "%)" : "")
              << "\n"
              << "  entries with 0 shared:  " << with_commas(no_shared)
              << "   <- block claims homology where none exists\n"
              << "  shared nodes:           " << with_commas(shared) << "\n"
              << "  ESCAPED nodes:          " << with_commas(escaped)
              << "   <- shared node outside the target block (false negative)\n"
              << "  orientation mismatches: " << with_commas(bad_orient) << "\n";
    if (ratio_n) {
        std::cout << "  tgt/src length ratio:   mean "
                  << std::fixed << std::setprecision(3) << (ratio_sum / ratio_n)
                  << "  min " << ratio_min << "  max " << ratio_max << "\n";
    }
    std::cout << "  verdict:                "
              << ((escaped == 0 && bad_orient == 0) ? "PASS" : "FAIL") << "\n";
    std::cout.flush();
    std::cerr.flush();
    const int rc = (escaped == 0 && bad_orient == 0) ? 0 : 1;
    // Static destruction of the loaded GBZ traps on some toolchains after all
    // work is done; exit before it can turn a clean run into a crash.
    std::_Exit(rc);
}
