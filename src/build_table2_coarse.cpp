/*
 * build_table2_coarse.cpp
 *
 * Build an ALL-PAIRS Translation Table 2 at LOCUS resolution instead of
 * variant resolution.
 *
 * WHY
 * ---
 * Table 2's job inside Index::translate() is not coordinate math — it is
 * (1) routing: which target path(s) can this source range reach, and
 * (2) extent: how much of the source range to hand the GBWT trace.
 * The actual translation is done by find_tags_in_interval /
 * find_first_and_last_common_nodes_gbwt / trace_coordinates_gbwt.
 *
 * The stock builder (build_translation_tables) records a segment at every point
 * where colinearity breaks — i.e. roughly every SNP, ~1 segment per kb per
 * pair. Across n haplotypes that is O(n^2) segments and reaches tens of TB for
 * a full HPRC graph (see bin/measure_table2_size). Those breakpoints are noise
 * for routing: `[0,1000)->j` and `[1001,2000)->j` say the same thing as
 * `[0,2000)->j`. Coalescing them shrinks the table by orders of magnitude and
 * makes all-pairs affordable, which in turn lets sample haplotypes (not just
 * the chromosome-named references) be translated to and from.
 *
 * SAFETY: merging only ever WIDENS a source range, so the table can only claim
 * "homology may exist here" too eagerly. A false positive costs one traversal
 * that then correctly returns nothing. It can never produce a false negative,
 * which is the error that would silently lose a real translation.
 *
 * OUTPUT FORMAT IS UNCHANGED: same TranslationTable2 on disk, same
 * (src_path_id, tgt_haplotype) key, same IntervalMapping values. Query code
 * needs no changes — it just sees fewer, wider segments.
 *
 * ALGORITHM (no pairwise path comparisons; those are O(paths^2 * length))
 * ----------------------------------------------------------------------
 *   Phase A  node -> bitmask of haplotypes visiting it. One pass over every
 *            path; atomic OR into a flat array. Cost: total path length.
 *   Phase B  per source path: walk it once, OR each node's mask into a coarse
 *            bin (default 10 kb). A bin is "covered" for haplotype h iff any
 *            node in it is also visited by h.
 *   Phase C  per source path: find runs of covered bins per haplotype using an
 *            XOR delta between consecutive bins (presence rarely changes, so
 *            this costs bins*words, not bins*haplotypes), merge runs separated
 *            by less than --merge-gap, and resolve each run's target path id
 *            with a single GBWT probe.
 *
 * Orientation is deliberately ignored when marking coverage (node id only, not
 * node+strand): an inverted shared region is still homology worth tracing, and
 * over-approximating is the safe direction.
 *
 * Usage:
 *   build_table2_coarse <graph.gbz> <fastlocate.ri> <output.t2> [options]
 *
 *     --threads N       worker threads (default: all)
 *     --bin-size N      coverage granularity in bp (default 10000)
 *     --merge-gap N     merge runs separated by <= N bp (default 100000)
 *     --progress        report per-phase progress
 *
 * Table 1 is NOT needed and NOT rebuilt: source coordinates here are
 * path-local, exactly as Table 2 stores them, and path ids come from the GBWT.
 *
 * BUILD NOTE: compile and link against the SAME gbwtgraph (headers and
 * library). Mixing the vendored headers with a different installed
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

void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " <graph.gbz> <fastlocate.ri> <output.t2> [options]\n\n"
        << "Builds an all-pairs Translation Table 2 at locus resolution.\n\n"
        << "  --threads N     worker threads (default: all)\n"
        << "  --bin-size N    coverage granularity in bp (default 10000)\n"
        << "  --merge-gap N   merge runs separated by <= N bp (default 100000)\n"
        << "  --progress      per-phase progress to stderr\n";
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

/// One emitted row, before it goes into the table.
struct Row {
    size_t src_path_id;
    uint32_t tgt_hap;      ///< index into hap_names
    panindexer::IntervalMapping mapping;
};

} // namespace

int main(int argc, char** argv) {
    std::string gbz_path, ri_path, out_path;
    int threads = 0;
    size_t bin_size = 10000;
    size_t merge_gap = 100000;
    bool progress = false;

    {
        std::vector<std::string> pos;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
            else if (a == "--threads" && i + 1 < argc) threads = std::stoi(argv[++i]);
            else if (a == "--bin-size" && i + 1 < argc) bin_size = std::stoull(argv[++i]);
            else if (a == "--merge-gap" && i + 1 < argc) merge_gap = std::stoull(argv[++i]);
            else if (a == "--progress") progress = true;
            else if (!a.empty() && a[0] == '-') { usage(argv[0]); return 1; }
            else pos.push_back(a);
        }
        if (pos.size() != 3) { usage(argv[0]); return 1; }
        gbz_path = pos[0]; ri_path = pos[1]; out_path = pos[2];
    }
    if (bin_size == 0) { std::cerr << "--bin-size must be > 0\n"; return 1; }

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

    // ------------------------------------------- Phase A: node -> hap mask
    const size_t max_node = static_cast<size_t>(gbz.graph.max_node_id());
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
    std::atomic<size_t> paths_done{0};
    #pragma omp parallel for schedule(dynamic, 8)
    for (size_t p = 0; p < n_paths; ++p) {
        const uint32_t h = path_hap[p];
        const size_t word = h / 64;
        const uint64_t bit = 1ULL << (h % 64);
        gbwt::vector_type nodes = index.extract(gbwt::Path::encode(p, false));
        for (gbwt::node_type node : nodes) {
            if (node == gbwt::ENDMARKER) break;
            const size_t nid = static_cast<size_t>(gbwt::Node::id(node));
            if (nid > max_node) continue;
            uint64_t* slot = &nodemask[nid * WORDS + word];
            // Paths of different haplotypes share nodes, so this OR races.
            if ((__atomic_load_n(slot, __ATOMIC_RELAXED) & bit) == 0) {
                __atomic_fetch_or(slot, bit, __ATOMIC_RELAXED);
            }
        }
        if (progress) {
            size_t d = ++paths_done;
            if (d % 512 == 0) {
                #pragma omp critical
                std::cerr << "    A: " << with_commas(d) << "/" << with_commas(n_paths)
                          << " paths\r" << std::flush;
            }
        }
    }
    std::cerr << "  phase A in " << secs(t_a, clk::now()) << "s" << std::endl;

    // ----------------------------- Phases B+C: coverage bins -> merged runs
    auto t_b = clk::now();
    std::cerr << "Phases B/C: coverage bins and run merging ..." << std::endl;

    const size_t merge_gap_bins = merge_gap / bin_size;
    std::vector<std::vector<Row>> per_thread;
#ifdef _OPENMP
    per_thread.resize(omp_get_max_threads());
#else
    per_thread.resize(1);
#endif

    std::atomic<size_t> src_done{0};
    std::atomic<size_t> probe_failures{0};

    #pragma omp parallel
    {
#ifdef _OPENMP
        std::vector<Row>& out = per_thread[omp_get_thread_num()];
#else
        std::vector<Row>& out = per_thread[0];
#endif
        std::vector<uint64_t> bins;
        std::vector<size_t> first_node_of_bin;
        std::vector<gbwt::node_type> nodes;
        // Per-haplotype open run state for the XOR delta scan.
        std::vector<size_t> run_start(H, 0);
        std::vector<char> run_open(H, 0);
        std::vector<std::pair<size_t, size_t>> runs;   // (bin_begin, bin_end_exclusive)

        #pragma omp for schedule(dynamic, 4)
        for (size_t sp = 0; sp < n_paths; ++sp) {
            const uint32_t src_hap = path_hap[sp];

            gbwt::vector_type ext = index.extract(gbwt::Path::encode(sp, false));
            nodes.clear();
            nodes.reserve(ext.size());
            size_t path_len = 0;
            for (gbwt::node_type node : ext) {
                if (node == gbwt::ENDMARKER) break;
                nodes.push_back(node);
                path_len += gbz.graph.get_length(gbz.graph.get_handle(
                    gbwt::Node::id(node), gbwt::Node::is_reverse(node)));
            }
            if (nodes.empty() || path_len == 0) continue;

            const size_t nbins = path_len / bin_size + 1;
            bins.assign(nbins * WORDS, 0);
            first_node_of_bin.assign(nbins, static_cast<size_t>(-1));

            // Phase B: OR each node's haplotype mask into its bin.
            {
                size_t off = 0;
                for (size_t k = 0; k < nodes.size(); ++k) {
                    const size_t nid = static_cast<size_t>(gbwt::Node::id(nodes[k]));
                    const size_t b = off / bin_size;
                    if (b < nbins) {
                        if (first_node_of_bin[b] == static_cast<size_t>(-1)) {
                            first_node_of_bin[b] = k;
                        }
                        if (nid <= max_node) {
                            const uint64_t* src = &nodemask[nid * WORDS];
                            uint64_t* dst = &bins[b * WORDS];
                            for (size_t w = 0; w < WORDS; ++w) dst[w] |= src[w];
                        }
                    }
                    off += gbz.graph.get_length(gbz.graph.get_handle(
                        gbwt::Node::id(nodes[k]), gbwt::Node::is_reverse(nodes[k])));
                }
            }

            // Phase C: per-haplotype runs via XOR delta between adjacent bins.
            // Presence changes rarely, so this costs bins*WORDS rather than
            // bins*haplotypes.
            std::fill(run_open.begin(), run_open.end(), 0);
            std::unordered_map<uint32_t, std::vector<std::pair<size_t, size_t>>> raw;
            for (size_t b = 0; b < nbins; ++b) {
                for (size_t w = 0; w < WORDS; ++w) {
                    const uint64_t cur = bins[b * WORDS + w];
                    const uint64_t prev = (b == 0) ? 0ULL : bins[(b - 1) * WORDS + w];
                    uint64_t diff = cur ^ prev;
                    while (diff) {
                        const int t = __builtin_ctzll(diff);
                        diff &= diff - 1;
                        const uint32_t h = static_cast<uint32_t>(w * 64 + t);
                        if (h >= H) continue;
                        if (cur & (1ULL << t)) {
                            run_start[h] = b;
                            run_open[h] = 1;
                        } else if (run_open[h]) {
                            raw[h].emplace_back(run_start[h], b);
                            run_open[h] = 0;
                        }
                    }
                }
            }
            // Close runs still open at the end of the path.
            for (size_t w = 0; w < WORDS; ++w) {
                uint64_t cur = bins[(nbins - 1) * WORDS + w];
                while (cur) {
                    const int t = __builtin_ctzll(cur);
                    cur &= cur - 1;
                    const uint32_t h = static_cast<uint32_t>(w * 64 + t);
                    if (h < H && run_open[h]) {
                        raw[h].emplace_back(run_start[h], nbins);
                        run_open[h] = 0;
                    }
                }
            }

            // Merge runs separated by a small gap, then emit one row each.
            for (auto& kv : raw) {
                const uint32_t h = kv.first;
                if (h == src_hap) continue;          // T2 skips same-haplotype pairs
                auto& r = kv.second;
                if (r.empty()) continue;
                std::sort(r.begin(), r.end());
                runs.clear();
                runs.push_back(r[0]);
                for (size_t i = 1; i < r.size(); ++i) {
                    if (r[i].first <= runs.back().second + merge_gap_bins) {
                        runs.back().second = std::max(runs.back().second, r[i].second);
                    } else {
                        runs.push_back(r[i]);
                    }
                }

                for (const auto& run : runs) {
                    const size_t src_start = run.first * bin_size;
                    size_t src_end = run.second * bin_size;
                    if (src_end > path_len || run.second >= nbins) src_end = path_len;
                    if (src_end <= src_start) continue;

                    // Resolve which path of haplotype h this run reaches: find a
                    // node in the run that h actually visits, then read the
                    // sequence ids off that node. One probe per emitted run.
                    size_t tgt_path_id = 0;
                    bool resolved = false;
                    for (size_t b = run.first; b < run.second && !resolved; ++b) {
                        size_t k = (b < nbins) ? first_node_of_bin[b] : static_cast<size_t>(-1);
                        if (k == static_cast<size_t>(-1)) continue;
                        const size_t nid = static_cast<size_t>(gbwt::Node::id(nodes[k]));
                        if (nid > max_node) continue;
                        if ((nodemask[nid * WORDS + h / 64] & (1ULL << (h % 64))) == 0) {
                            continue;   // h not on this node; try the next bin
                        }
                        // h visits this node — in either orientation.
                        for (int orient = 0; orient < 2 && !resolved; ++orient) {
                            gbwt::node_type gn = gbwt::Node::encode(
                                static_cast<gbwt::size_type>(nid), orient == 1);
                            std::vector<gbwt::size_type> sa = rindex.decompressSA(gn);
                            for (gbwt::size_type v : sa) {
                                const size_t seq = static_cast<size_t>(rindex.seqId(v));
                                const size_t pid = seq / 2;
                                if (pid < n_paths && path_hap[pid] == h) {
                                    tgt_path_id = pid;
                                    resolved = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (!resolved) { ++probe_failures; continue; }

                    Row row;
                    row.src_path_id = sp;
                    row.tgt_hap = h;
                    row.mapping.src_start = src_start;
                    row.mapping.src_end = src_end;
                    row.mapping.tgt_path_id = tgt_path_id;
                    out.push_back(row);
                }
            }

            if (progress) {
                size_t d = ++src_done;
                if (d % 256 == 0) {
                    #pragma omp critical
                    std::cerr << "    B/C: " << with_commas(d) << "/" << with_commas(n_paths)
                              << " paths\r" << std::flush;
                }
            }
        }
    }
    std::cerr << "  phases B/C in " << secs(t_b, clk::now()) << "s" << std::endl;

    // ------------------------------------------------------------- assemble
    auto t_w = clk::now();
    size_t total_rows = 0;
    for (const auto& v : per_thread) total_rows += v.size();
    std::cerr << "Assembling " << with_commas(total_rows) << " segments ..." << std::endl;

    panindexer::TranslationTable2 table2;
    for (auto& v : per_thread) {
        for (const Row& r : v) {
            table2.add_mapping(r.src_path_id, hap_names[r.tgt_hap], r.mapping);
        }
        std::vector<Row>().swap(v);   // release as we go
    }
    table2.finalize();

    {
        std::ofstream out(out_path, std::ios::binary);
        if (!out) { std::cerr << "Error: cannot write " << out_path << "\n"; return 1; }
        table2.serialize(out);
    }
    std::cerr << "  written in " << secs(t_w, clk::now()) << "s" << std::endl;

    std::cout << "\n=============== coarse Table 2 ===============\n"
              << "  keys (src_path, tgt_haplotype): " << with_commas(table2.num_entries()) << "\n"
              << "  segments:                       " << with_commas(table2.total_segments()) << "\n"
              << "  unresolved runs (skipped):      " << with_commas(probe_failures.load()) << "\n"
              << "  output:                         " << out_path << "\n"
              << "  total time:                     " << std::fixed << std::setprecision(1)
              << secs(t0, clk::now()) << "s\n";
    std::cout.flush();
    std::cerr.flush();
    // Skip static destruction of the loaded GBZ/sdsl structures: it traps on
    // some toolchains after all work is complete and would turn a successful
    // build into a nonzero exit status.
    std::_Exit(0);
}
