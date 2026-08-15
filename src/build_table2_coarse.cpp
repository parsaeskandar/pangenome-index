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
        << "  --min-run-bp N  drop runs shorter than N bp (default 0 = keep all).\n"
        << "                  Repeats make distant haplotypes share short stretches;\n"
        << "                  raising this cuts spurious pairs and build time sharply.\n"
        << "  --max-paths N   process only the first N source paths (throughput probe)\n"
        << "  --progress-every N  progress line every N paths (default 1000000)\n"
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
    size_t min_run_bp = 0;
    size_t max_paths = 0;
    size_t progress_every = 1000000;
    bool progress = false;

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
            const size_t d = ++paths_done;
            if (d % progress_every == 0) {
                const double el = secs(t_a, clk::now());
                #pragma omp critical
                std::cerr << "  A " << with_commas(d) << "/" << with_commas(n_paths)
                          << " paths  " << std::fixed << std::setprecision(0) << el << "s  "
                          << with_commas(static_cast<unsigned long long>(d / (el > 0 ? el : 1)))
                          << "/s" << std::endl;
            }
        }
    }
    std::cerr << "  phase A in " << secs(t_a, clk::now()) << "s" << std::endl;

    // Node length lookup. Phase B needs a base offset for every node it walks;
    // going through get_handle()/get_length() did that twice per node and was a
    // large part of the per-node cost. One dense array replaces both calls, and
    // it is reused by every source path.
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
    std::atomic<size_t> probes_done{0};
    std::atomic<size_t> runs_dropped{0};
    const size_t n_src = (max_paths > 0 && max_paths < n_paths) ? max_paths : n_paths;
    if (n_src != n_paths) {
        std::cerr << "  (limited to the first " << with_commas(n_src)
                  << " source paths by --max-paths)" << std::endl;
    }

    #pragma omp parallel
    {
#ifdef _OPENMP
        std::vector<Row>& out = per_thread[omp_get_thread_num()];
#else
        std::vector<Row>& out = per_thread[0];
#endif
        std::vector<uint64_t> bins;
        // First GBWT node seen in each bin, used to resolve the target path id.
        // Storing the node itself (not an index) means the path's node list does
        // not have to be materialised at all.
        std::vector<gbwt::node_type> first_node_of_bin;

        // Per-haplotype run state. Flat arrays reused across paths: the previous
        // version built an unordered_map of vectors per path, which allocated
        // heavily because repeats make distant haplotypes flicker in and out of
        // coverage, producing very many short raw runs.
        std::vector<size_t> run_start(H, 0), run_last(H, 0);
        std::vector<char>   run_open(H, 0);

        // Target-path probe cache. decompressSA() is the single most expensive
        // call here, and it used to run once per emitted run. Runs of the same
        // haplotype that are close together on this source path resolve to the
        // same target contig, so one probe serves all of them; runs further
        // apart than the merge gap still probe again, so a source path spanning
        // several contigs of the same haplotype stays correct.
        std::vector<size_t>   probe_tgt(H, 0);
        std::vector<uint64_t> probe_stamp(H, 0);
        std::vector<size_t>   probe_last_bin(H, 0);
        uint64_t stamp = 0;

        #pragma omp for schedule(dynamic, 1)
        for (size_t sp = 0; sp < n_src; ++sp) {
            const uint32_t src_hap = path_hap[sp];
            ++stamp;

            // ---- Phase B: one pass. Bin coverage and accumulate base offsets
            // together, with no copy of the node list and one length lookup.
            gbwt::vector_type ext = index.extract(gbwt::Path::encode(sp, false));
            bins.clear();
            first_node_of_bin.clear();
            size_t off = 0;
            for (gbwt::node_type node : ext) {
                if (node == gbwt::ENDMARKER) break;
                const size_t nid = static_cast<size_t>(gbwt::Node::id(node));
                const size_t b = off / bin_size;
                if (b >= first_node_of_bin.size()) {
                    first_node_of_bin.resize(b + 1, gbwt::ENDMARKER);
                    bins.resize((b + 1) * WORDS, 0);
                }
                if (first_node_of_bin[b] == gbwt::ENDMARKER) first_node_of_bin[b] = node;
                if (nid <= max_node) {
                    const uint64_t* src = &nodemask[nid * WORDS];
                    uint64_t* dst = &bins[b * WORDS];
                    for (size_t w = 0; w < WORDS; ++w) dst[w] |= src[w];
                    off += node_len[nid];
                }
            }
            const size_t path_len = off;
            const size_t nbins = first_node_of_bin.size();
            if (nbins == 0 || path_len == 0) {
                if (progress) {
                    const size_t d = ++src_done;
                    if (d % progress_every == 0) {
                        const double el = secs(t_b, clk::now());
                        #pragma omp critical
                        std::cerr << "  B/C " << with_commas(d) << "/" << with_commas(n_src)
                                  << " paths  " << std::fixed << std::setprecision(0) << el << "s"
                                  << std::endl;
                    }
                }
                continue;
            }

            // Emit one merged run, resolving its target path id (cached).
            auto emit_run = [&](uint32_t h, size_t b0, size_t b1) {   // b1 inclusive
                if (h == src_hap) return;          // T2 skips same-haplotype pairs
                if (b1 < b0) return;
                const size_t src_start = b0 * bin_size;
                size_t src_end = (b1 + 1) * bin_size;
                if (src_end > path_len) src_end = path_len;
                if (src_end <= src_start) return;
                // Filter on actual span, not bin count: a run can occupy two bins
                // yet be far shorter than the threshold.
                if (src_end - src_start < min_run_bp) { ++runs_dropped; return; }

                size_t tgt = 0;
                bool have = false;
                if (probe_stamp[h] == stamp &&
                    b0 <= probe_last_bin[h] + merge_gap_bins + 1) {
                    tgt = probe_tgt[h];
                    have = true;
                } else {
                    for (size_t b = b0; b <= b1 && !have; ++b) {
                        const gbwt::node_type gn0 = first_node_of_bin[b];
                        if (gn0 == gbwt::ENDMARKER) continue;
                        const size_t nid = static_cast<size_t>(gbwt::Node::id(gn0));
                        if (nid > max_node) continue;
                        if ((nodemask[nid * WORDS + h / 64] & (1ULL << (h % 64))) == 0) {
                            continue;   // this bin's first node is not on h
                        }
                        for (int orient = 0; orient < 2 && !have; ++orient) {
                            gbwt::node_type gn = gbwt::Node::encode(
                                static_cast<gbwt::size_type>(nid), orient == 1);
                            std::vector<gbwt::size_type> sa = rindex.decompressSA(gn);
                            ++probes_done;
                            for (gbwt::size_type v : sa) {
                                const size_t pid =
                                    static_cast<size_t>(rindex.seqId(v)) / 2;
                                if (pid < n_paths && path_hap[pid] == h) {
                                    tgt = pid; have = true; break;
                                }
                            }
                        }
                    }
                    if (have) {
                        probe_stamp[h] = stamp;
                        probe_tgt[h] = tgt;
                    }
                }
                if (!have) { ++probe_failures; return; }
                probe_last_bin[h] = b1;

                Row row;
                row.src_path_id = sp;
                row.tgt_hap = h;
                row.mapping.src_start = src_start;
                row.mapping.src_end = src_end;
                row.mapping.tgt_path_id = tgt;
                out.push_back(row);
            };

            // ---- Phase C: runs found by XOR delta between adjacent bins, with
            // gap merging applied inline so no intermediate run list is built.
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
                                if (run_open[h]) emit_run(h, run_start[h], run_last[h]);
                                run_start[h] = b;
                            }
                            run_open[h] = 1;
                            run_last[h] = b;
                        } else if (run_open[h]) {          // coverage stops
                            run_last[h] = b - 1;
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
                    if (h < H && run_open[h]) run_last[h] = nbins - 1;
                }
            }
            for (uint32_t h = 0; h < H; ++h) {
                if (run_open[h]) emit_run(h, run_start[h], run_last[h]);
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
