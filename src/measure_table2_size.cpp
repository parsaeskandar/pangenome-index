/*
 * measure_table2_size.cpp
 *
 * Compute — exactly, and without building anything — how many IntervalMapping
 * segments a Translation Table 2 would contain if it covered ALL haplotype
 * pairs, plus the resulting file size. Answers "how big does this get?" in
 * O(total GBWT runs) instead of running a multi-day, multi-terabyte build.
 *
 * THE IDEA
 * --------
 * A T2 segment is a maximal colinear run: source path `a` and target path `b`
 * traverse the same consecutive nodes. So the number of segments for the pair
 * (a, b) is exactly the number of positions where their shared walk *starts* —
 * i.e. where they come together after having been apart.
 *
 * At any node v, the paths visiting v arrive from their respective predecessor
 * nodes. Two paths continue an existing run iff they arrived from the SAME
 * predecessor; if they arrived from different predecessors, a new run begins at
 * v. So if v's visits partition by predecessor into groups of sizes w_1..w_d
 * (W = sum w_i), the number of ordered path pairs starting a run at v is
 *
 *     W^2 - sum_i w_i^2          (= 2*w_1*w_2 for d = 2)
 *
 * and the total segment count is that summed over every node. Path starts are
 * handled for free: paths begin at the GBWT endmarker, so a first node's visits
 * form their own predecessor group.
 *
 * The w_i are just GBWT edge weights — the number of paths taking edge u->v —
 * which we read by walking each record's runs. No
 * path extraction, no pair enumeration, no homology computation.
 *
 * ORIENTATION: a GBZ's GBWT is bidirectional, storing every path twice (forward
 * and reverse). The reverse copies mirror the forward ones exactly, so the raw
 * sum double-counts; we halve it.
 *
 * WHAT THIS IS AND ISN'T
 * ----------------------
 * Counts ordered pairs of GBWT paths that share nodes — which is what T2 keys
 * (src_path_id, tgt_haplotype) enumerate. It does NOT subtract the small number
 * of same-haplotype pairs (different subpaths of one haplotype), which T2 skips,
 * so it is a slight over-estimate. It is exact about the thing that dominates:
 * the O(n^2) growth in shared-walk breakpoints.
 *
 * Usage:
 *   measure_table2_size <graph.gbz> [--threads N] [--ref-prefix P]...
 *
 * --ref-prefix (repeatable; defaults to GRCh38, CHM13, HG002) marks which
 * haplotypes are references, so the tool also reports the reference-mediated
 * ("star") size: every haplotype paired only with the references, which supports
 * any->any translation via one reference hop at a fraction of the size.
 */

#include <gbwt/gbwt.h>
#include <gbwt/internal.h>
#include <gbwt/support.h>
#include <gbwtgraph/gbz.h>
#include <sdsl/simple_sds.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <atomic>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

/// Bytes on disk per IntervalMapping: src_start, src_end, tgt_path_id as
/// uint64 each (TranslationTable2::serialize).
constexpr double BYTES_PER_SEGMENT = 24.0;

void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " <graph.gbz> [--threads N] [--ref-prefix PREFIX]...\n\n"
        << "Reports the exact segment count and file size an all-pairs Translation\n"
        << "Table 2 would have, without building it.\n\n"
        << "  --threads N      OpenMP threads (default: all available)\n"
        << "  --ref-prefix P   sample name treated as a reference (repeatable;\n"
        << "                   default: GRCh38, CHM13, HG002)\n";
}

std::string human_bytes(double bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB", "PB"};
    int u = 0;
    while (bytes >= 1024.0 && u < 5) { bytes /= 1024.0; ++u; }
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(bytes < 10 ? 2 : 1) << bytes << " " << units[u];
    return ss.str();
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

} // namespace

int main(int argc, char** argv) {
    std::string gbz_path;
    int threads = 0;
    std::vector<std::string> ref_prefixes;
    bool debug_records = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else if (a == "--threads" && i + 1 < argc) threads = std::stoi(argv[++i]);
        else if (a == "--ref-prefix" && i + 1 < argc) ref_prefixes.push_back(argv[++i]);
        else if (a == "--debug-records") debug_records = true;
        else if (a.size() && a[0] == '-') { usage(argv[0]); return 1; }
        else if (gbz_path.empty()) gbz_path = a;
        else { usage(argv[0]); return 1; }
    }
    if (gbz_path.empty()) { usage(argv[0]); return 1; }
    if (ref_prefixes.empty()) ref_prefixes = {"GRCh38", "CHM13", "HG002"};

#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    std::cerr << "Threads: " << omp_get_max_threads() << std::endl;
#endif

    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    auto secs = [](clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    std::cerr << "Loading " << gbz_path << " ..." << std::endl;
    gbwtgraph::GBZ gbz;
    sdsl::simple_sds::load_from(gbz, gbz_path);
    const gbwt::GBWT& index = gbz.index;
    auto t_load = clock::now();
    std::cerr << "  loaded in " << std::fixed << std::setprecision(1)
              << secs(t0, t_load) << "s" << std::endl;

    const size_t sigma = index.sigma();
    const size_t effective = index.effective();
    const size_t n_paths = index.sequences() / 2;   // bidirectional: 2 seqs/path

    // Incoming edge weights per node, accumulated as W = sum(w_i) and
    // SQ = sum(w_i^2). Storing only these two aggregates (not the weight lists)
    // is all the W^2 - sum(w_i^2) formula needs.
    std::vector<uint64_t> W(sigma, 0), SQ(sigma, 0);

    std::cerr << "Scanning " << with_commas(effective) << " GBWT records for edge weights..."
              << std::endl;

    uint64_t total_runs = 0, total_visits = 0, total_edges = 0;
    // Progress: a whole-genome graph has hundreds of millions of records, so
    // the scan must show it is alive.
    const bool report_progress = (effective > (1u << 24));
    std::atomic<size_t> records_done{0};

    #pragma omp parallel for schedule(dynamic, 1024) \
            reduction(+:total_runs, total_visits, total_edges)
    for (size_t comp = 0; comp < effective; ++comp) {
        gbwt::node_type from = index.toNode(comp);
        gbwt::CompressedRecord rec = index.record(from);

        // MUST come before touching the run decoder: CompressedRecordIterator
        // constructs Run(outdegree), which divides by outdegree — so an empty
        // record (outdegree 0, common in a large graph's alphabet) is a division
        // by zero, i.e. SIGFPE. Checking rec.size() here instead would be both
        // too late and expensive (CompressedRecord::size() decodes the record).
        const size_t outdeg = rec.outdegree();
        if (outdeg == 0) continue;

        // Sum run lengths per outgoing edge: run_type is (outrank, length), so
        // this is O(runs) with no allocation — much cheaper than decompressing
        // each record's body, which matters at hundreds of millions of records.
        uint64_t local[16];
        std::vector<uint64_t> spill;
        uint64_t* weights = local;
        if (outdeg > 16) { spill.assign(outdeg, 0); weights = spill.data(); }
        else { for (size_t r = 0; r < outdeg; ++r) local[r] = 0; }

        uint64_t body = 0;
        for (gbwt::CompressedRecordIterator iter(rec); !(iter.end()); ++iter) {
            weights[iter->first] += iter->second;
            body += iter->second;
            ++total_runs;
        }
        if (body == 0) continue;
        total_visits += body;

        if (report_progress) {
            const size_t done = ++records_done;
            if ((done & 0xFFFFFF) == 0) {   // every ~16.8M records
                #pragma omp critical
                std::cerr << "    " << with_commas(done) << " / "
                          << with_commas(effective) << " records\r" << std::flush;
            }
        }

        if (debug_records) {
            #pragma omp critical
            {
                std::cerr << "  rec node=" << from << " body=" << body
                          << " outdeg=" << outdeg << " :";
                for (size_t r = 0; r < outdeg; ++r) {
                    std::cerr << " [" << rec.successor(r) << "]=" << weights[r];
                }
                std::cerr << std::endl;
            }
        }

        for (size_t r = 0; r < outdeg; ++r) {
            const gbwt::node_type to = rec.successor(r);
            const uint64_t w = weights[r];
            if (w == 0) continue;
            if (to == gbwt::ENDMARKER) continue;   // path ends start no run
            ++total_edges;
            const size_t idx = to;                 // node ids index W/SQ directly
            if (idx >= sigma) continue;
            #pragma omp atomic
            W[idx] += w;
            #pragma omp atomic
            SQ[idx] += w * w;
        }
    }

    auto t_scan = clock::now();
    std::cerr << "  scanned in " << secs(t_load, t_scan) << "s" << std::endl;

    // Sum W^2 - sum(w_i^2) over nodes. Use long double: W^2 summed over a whole
    // pangenome can exceed 2^64.
    // Serial: one pass over the node arrays, cheap next to the record scan.
    long double raw_pairs = 0.0L;
    uint64_t nodes_with_visits = 0, max_W = 0;
    for (size_t v = 0; v < sigma; ++v) {
        const uint64_t w_total = W[v];
        if (w_total == 0) continue;
        ++nodes_with_visits;
        if (w_total > max_W) max_W = w_total;
        const long double t = static_cast<long double>(w_total);
        raw_pairs += t * t - static_cast<long double>(SQ[v]);
    }
    std::cerr << "  summed " << with_commas(nodes_with_visits) << " visited nodes" << std::endl;

    // Halve: the bidirectional GBWT stores each path forward and reverse, and the
    // reverse copies produce a mirrored (identical) count.
    const long double all_pairs_segments = raw_pairs / 2.0L;
    const long double all_pairs_bytes = all_pairs_segments * BYTES_PER_SEGMENT;

    // Reference-mediated ("star") estimate: keep only pairs where at least one
    // side is a reference. Fraction of ordered pairs = 1 - (n-r)(n-r-1)/(n(n-1)).
    size_t n_ref_paths = 0;
    if (index.hasMetadata() && index.metadata.hasPathNames() &&
        index.metadata.hasSampleNames()) {
        const gbwt::Metadata& meta = index.metadata;
        for (size_t p = 0; p < meta.paths(); ++p) {
            std::string sample = meta.sample(meta.path(p).sample);
            for (const std::string& pre : ref_prefixes) {
                if (sample.size() >= pre.size() &&
                    sample.compare(0, pre.size(), pre) == 0) { ++n_ref_paths; break; }
            }
        }
    }
    long double star_fraction = 0.0L;
    if (n_paths > 1 && n_ref_paths <= n_paths) {
        const long double n = static_cast<long double>(n_paths);
        const long double m = static_cast<long double>(n_paths - n_ref_paths);
        star_fraction = 1.0L - (m * (m - 1.0L)) / (n * (n - 1.0L));
        if (star_fraction < 0.0L) star_fraction = 0.0L;
    }

    auto t_end = clock::now();

    std::cout << "\n================ graph ================\n"
              << "  GBWT paths (haplotype subpaths): " << with_commas(n_paths) << "\n"
              << "  reference paths (by --ref-prefix): " << with_commas(n_ref_paths) << "\n"
              << "  GBWT nodes (records):            " << with_commas(effective) << "\n"
              << "  nodes with visits:               " << with_commas(nodes_with_visits) << "\n"
              << "  total path visits (both strands):" << with_commas(total_visits) << "\n"
              << "  total runs scanned:              " << with_commas(total_runs) << "\n"
              << "  GBWT edges counted:              " << with_commas(total_edges) << "\n"
              << "  max visits at one node:          " << with_commas(max_W) << "\n";

    std::cout << "\n========= ALL-PAIRS Table 2 (every haplotype -> every haplotype) =========\n"
              << "  segments: " << std::scientific << std::setprecision(3)
              << static_cast<double>(all_pairs_segments) << "\n"
              << "  size:     " << human_bytes(static_cast<double>(all_pairs_bytes))
              << "   (at " << BYTES_PER_SEGMENT << " B/segment)\n";

    if (star_fraction > 0.0L) {
        const long double star_seg = all_pairs_segments * star_fraction;
        std::cout << "\n========= REFERENCE-MEDIATED Table 2 (all haplotypes <-> references) =========\n"
                  << "  fraction of all-pairs: " << std::fixed << std::setprecision(3)
                  << static_cast<double>(star_fraction * 100.0L) << "%\n"
                  << "  segments: " << std::scientific << std::setprecision(3)
                  << static_cast<double>(star_seg) << "\n"
                  << "  size:     "
                  << human_bytes(static_cast<double>(star_seg * BYTES_PER_SEGMENT)) << "\n"
                  << "  (supports any->any translation via one reference hop)\n";
    }

    std::cout << "\nTotal time: " << std::fixed << std::setprecision(1)
              << secs(t0, t_end) << "s\n"
              << "\nNote: counts ordered GBWT path pairs sharing nodes; does not subtract\n"
              << "same-haplotype subpath pairs (T2 skips those), so it is a slight\n"
              << "over-estimate. Excludes per-key name/index overhead.\n";
    std::cout.flush();
    std::cerr.flush();
    // Skip static destruction: tearing down the loaded GBZ/sdsl structures traps
    // on some toolchains (an abort AFTER all work and output is complete, which
    // would otherwise give this read-only tool a nonzero exit status). Nothing
    // here owns unflushed state or external resources, so exiting directly is
    // safe and keeps the exit code meaningful for scripts.
    std::_Exit(0);
}
