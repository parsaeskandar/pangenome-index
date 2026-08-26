#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include "pangenome_index/translation_tables.hpp"
#include "pangenome_index/surject_anchor_builder.hpp"
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/simple_sds.hpp>
#include <gbwt/gbwt.h>
#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>
#include <handlegraph/util.hpp>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <chrono>
#include <sys/mman.h>

#include "pangenome_server.hpp"

using panindexer::FastLocate;
using panindexer::SampledTagArray;
using panindexer::TranslationTable1;
using panindexer::TranslationTable2;
using panindexer::PathInterval;
using panindexer::SubpathInfo;
using panindexer::TargetInterval;
using panindexer::IntervalMapping;

// ── Types matching coordinate_translation.cpp at file scope ────────────────
// These must be identical to the definitions in coordinate_translation.cpp
// so that extern function declarations resolve correctly at link time.

struct TagInfo {
    uint64_t tag_code;
    std::vector<size_t> source_offsets;
    std::vector<size_t> source_bwt_positions;
    std::vector<size_t> source_packed_positions;
};

struct TranslationResult {
    size_t source_offset;
    size_t target_offset;
    size_t target_seq_id;
    uint64_t tag_code;
};

struct FindTagsInIntervalTiming {
    size_t num_lf_before_phase1 = 0;
    double init_skip_lf_ms = 0;
    size_t phase1_lf_count = 0;
    double phase1_lf_ms = 0;
    double phase1_tag_lookup_ms = 0;
    double phase1_position_lookup_ms = 0;
    double phase1_tag_map_ms = 0;
    double phase1_fast_path_ms = 0;
    bool used_fast_path = false;
    bool found_last_common = false;
    size_t last_common_source_base = 0;
};

struct CommonNodes {
    bool first_is_unique = false;
    bool last_is_unique = false;
    size_t first_source_offset;
    size_t first_target_offset;
    size_t first_source_base;
    size_t first_target_base;
    uint64_t first_tag_code;
    size_t last_source_offset;
    size_t last_target_offset;
    size_t last_source_base;
    size_t last_target_base;
    uint64_t last_tag_code;
    bool found;
};

// Mirror of coordinate_translation.cpp's file-scope FindSeqStats (kept in sync
// by hand, per the extern-struct convention used for the structs above). The
// thread-local accumulator there records find_sequences_for_tag's LF cost;
// we reset it before a build and read it after to report per-query diagnostics.
struct FindSeqStats {
    size_t calls = 0;
    size_t runs = 0;
    size_t lf_steps = 0;
    size_t visits = 0;
    size_t last_run_nav_steps = 0;
    size_t last_run_length = 0;
    double total_ms = 0.0;
};
extern thread_local FindSeqStats g_find_seq_stats;

// ── Extern declarations for functions defined in coordinate_translation.cpp ─

extern std::vector<TagInfo> find_tags_in_interval(
    FastLocate& r_index, SampledTagArray& sampled,
    size_t source_seq_id, size_t seq_start, size_t seq_end,
    const gbwt::GBWT* gbwt_index,
    const gbwt::FastLocate* gbwt_fast_locate,
    const gbwtgraph::GBWTGraph* graph,
    size_t target_seq_id,
    FindTagsInIntervalTiming* out_timing);

// File-scope (not inside a namespace) so the extern matches the definition in
// coordinate_translation.cpp, per the same convention as the structs above.
struct NodeVisit {
    size_t seq_id;
    size_t offset;
    size_t bwt_pos;
    size_t packed_pos;
    uint64_t tag_code;
};

extern std::vector<NodeVisit> find_sequences_for_tag(
    FastLocate& r_index, SampledTagArray& sampled, uint64_t tag_code);

extern CommonNodes find_first_and_last_common_nodes_gbwt(
    const gbwt::FastLocate& gbwt_fast_locate,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const std::vector<TagInfo>& source_tags,
    size_t source_seq_id, size_t target_seq_id);

extern std::vector<TranslationResult> trace_coordinates_gbwt(
    const gbwt::GBWT& gbwt_index,
    const gbwt::FastLocate& gbwt_fast_locate,
    const gbwtgraph::GBWTGraph& graph,
    size_t source_seq_id, size_t source_start, size_t source_end,
    size_t target_seq_id, size_t anchor_source_offset, size_t anchor_target_offset,
    size_t anchor_source_base, size_t anchor_target_base,
    uint64_t anchor_tag_code,
    size_t last_common_source_base, size_t last_common_target_base,
    uint64_t last_common_tag_code);

// ── Local helpers ──────────────────────────────────────────────────────────

namespace {

std::unordered_map<size_t, std::pair<std::string, size_t>>
build_path_id_to_global(const TranslationTable1& t1) {
    std::unordered_map<size_t, std::pair<std::string, size_t>> out;
    std::vector<std::string> names = t1.names();
    for (const std::string& name : names) {
        std::vector<SubpathInfo> subpaths = t1.subpaths(name);
        for (const SubpathInfo& sp : subpaths) {
            out[sp.path_id] = {name, sp.subpath_start};
        }
    }
    return out;
}

std::vector<std::string>
path_names_for_haplotype(const TranslationTable1& t1,
                         const std::string& haplotype_prefix) {
    std::vector<std::string> names = t1.names();
    std::vector<std::string> result;

    // 1) Exact full-name match first: a full contig/path name like
    //    "CHM13#0#chr10" resolves to just that contig, so a query interval
    //    means that one locus (not the same offset on every contig).
    for (const std::string& name : names) {
        if (name == haplotype_prefix) {
            result.push_back(name);
        }
    }
    if (!result.empty()) {
        return result;
    }

    // 2) Otherwise treat the argument as a haplotype prefix ("CHM13#0") and
    //    match every contig under it. (Same behavior as before; note the
    //    interval is then looked up in each matched contig's coordinates.)
    std::string prefix = haplotype_prefix;
    if (prefix.empty() || prefix.back() != '#') {
        prefix += '#';
    }
    for (const std::string& name : names) {
        if (name.size() >= prefix.size() &&
            name.compare(0, prefix.size(), prefix) == 0) {
            result.push_back(name);
        }
    }
    return result;
}

} // anonymous namespace

// ── Index implementation ───────────────────────────────────────────────────

static constexpr int64_t MAX_INTERVAL_LENGTH = 10'000'000;

void Index::load(const std::string& gbz_path,
                 const std::string& ri_path,
                 const std::string& tags_path,
                 const std::string& gbwt_ri_path,
                 const std::string& table1_path,
                 const std::string& table2_path) {
    // Pinning pages in RAM (mlockall) avoids page-fault stalls during queries,
    // but mlockall(MCL_CURRENT | MCL_FUTURE) makes EVERY later allocation fail
    // when RLIMIT_MEMLOCK (`ulimit -l`) is too small to lock the working set —
    // which aborts index loading on most shared clusters (the default memlock
    // limit is often only 64 KB). So pinning is OPT-IN: set PANINDEX_MLOCK=1
    // only after raising `ulimit -l` (or with an admin-set unlimited limit).
    // PANINDEX_NO_MLOCK is still honored (and now redundant) for back-compat.
    if (std::getenv("PANINDEX_MLOCK") != nullptr &&
        std::getenv("PANINDEX_NO_MLOCK") == nullptr) {
        if (mlockall(MCL_CURRENT | MCL_FUTURE) != 0) {
            std::perror("[Index::load] mlockall failed; continuing without pinning");
        }
    }

    // Per-step progress logs so a bad_alloc identifies which artifact failed.
    // Silence by setting PANINDEX_QUIET_LOAD=1.
    const bool verbose = (std::getenv("PANINDEX_QUIET_LOAD") == nullptr);
    auto log_step = [&](const char* msg) {
        if (verbose) std::cerr << "[Index::load] " << msg << std::endl;
    };

    // 1. RLBWT r-index
    log_step(("[1/5] RLBWT r-index: " + ri_path).c_str());
    {
        std::ifstream rin(ri_path, std::ios::binary);
        if (!rin)
            throw std::runtime_error("Cannot open RLBWT r-index: " + ri_path);
        rlbwt_rindex_.load_encoded(rin);
        log_step("[1/5]   load_encoded done");
        rlbwt_rindex_.ensure_last_rank();
        log_step("[1/5]   ensure_last_rank done");
        rlbwt_rindex_.ensure_last_select();
        log_step("[1/5]   ensure_last_select done");
        // The third lazy support. Like the two above it is `mutable` + built
        // under std::call_once, so leaving it lazy is safe — but then the first
        // query to need it pays the construction cost and every concurrent query
        // blocks on that call_once. Building all three here keeps the first
        // query as cheap as the rest.
        rlbwt_rindex_.ensure_blocks_start_select();
        log_step("[1/5]   ensure_blocks_start_select done");
    }

    // 2. GBZ (GBWT + GBWTGraph)
    log_step(("[2/5] GBZ: " + gbz_path).c_str());
    gbz_ = std::make_unique<gbwtgraph::GBZ>();
    sdsl::simple_sds::load_from(*gbz_, gbz_path);
    log_step("[2/5]   loaded");

    // 3. GBWT FastLocate
    log_step(("[3/5] GBWT FastLocate: " + gbwt_ri_path).c_str());
    {
        std::ifstream gin(gbwt_ri_path, std::ios::binary);
        if (!gin)
            throw std::runtime_error("Cannot open GBWT FastLocate: " + gbwt_ri_path);
        gin.seekg(0, std::ios::end);
        if (static_cast<size_t>(gin.tellg()) == 0)
            throw std::runtime_error("GBWT FastLocate file is empty: " + gbwt_ri_path);
        gin.seekg(0, std::ios::beg);
        gbwt_rindex_ = std::make_unique<gbwt::FastLocate>();
        gbwt_rindex_->load(gin);
        gbwt_rindex_->setGBWT(gbz_->index);
        log_step("[3/5]   loaded");
    }

    // 4. Sampled tag array
    log_step(("[4/5] Sampled tag array: " + tags_path).c_str());
    {
        std::ifstream sin(tags_path, std::ios::binary);
        if (!sin)
            throw std::runtime_error("Cannot open sampled tags: " + tags_path);
        sampled_.load(sin);
        log_step("[4/5]   load done");
        sampled_.ensure_run_rank();
        log_step("[4/5]   ensure_run_rank done");
        sampled_.ensure_run_select();
        log_step("[4/5]   ensure_run_select done");
    }

    // 5. Translation tables
    log_step(("[5/5] Table 1: " + table1_path).c_str());
    {
        std::ifstream t1in(table1_path, std::ios::binary);
        if (!t1in)
            throw std::runtime_error("Cannot open Table 1: " + table1_path);
        table1_.load(t1in);
    }
    log_step("[5/5]   T1 loaded");
    // Table 2 is OPTIONAL. Pass an empty path to run without it: translation
    // then routes through the GBWT/tag array (translate_no_table2), which needs
    // no per-path-pair table and so is unaffected by how finely the graph
    // fragments haplotypes into GBWT paths.
    if (table2_path.empty()) {
        has_table2_ = false;
        log_step("[5/5] Table 2: (none) - using table-free translation");
    } else {
        log_step(("[5/5] Table 2: " + table2_path).c_str());
        std::ifstream t2in(table2_path, std::ios::binary);
        if (!t2in)
            throw std::runtime_error("Cannot open Table 2: " + table2_path);
        table2_.load(t2in);
        has_table2_ = true;
        // Report the form, since both are valid inputs and they behave very
        // differently: a "B2" table carries target intervals and one coarse
        // block per path pair, the older form only names the target path.
        log_step(table2_.has_target_coords()
                 ? "[5/5]   T2 loaded (B2 form: with target intervals)"
                 : "[5/5]   T2 loaded (target path ids only, no target intervals)");
    }

    path_to_global_ = build_path_id_to_global(table1_);

    // Cache the haplotype name list. Deriving it walks every GBWT path, which
    // is hundreds of millions of entries on a fragmented graph — far too slow
    // to repeat per query.
    log_step("[5/5] caching haplotype names ...");
    {
        const gbwt::GBWT& gi = gbz_->index;
        if (gi.hasMetadata() && gi.metadata.hasPathNames() &&
            gi.metadata.hasSampleNames()) {
            const gbwt::Metadata& meta = gi.metadata;
            std::unordered_set<std::string> seen;
            for (size_t p = 0; p < meta.paths(); ++p) {
                gbwt::PathName pn = meta.path(p);
                std::string h = meta.sample(pn.sample) + "#" + std::to_string(pn.phase);
                if (seen.insert(h).second) haplotype_names_.push_back(std::move(h));
            }
            std::sort(haplotype_names_.begin(), haplotype_names_.end());
        }
    }
    log_step(("[5/5]   " + std::to_string(haplotype_names_.size()) +
              " haplotypes cached").c_str());

    loaded_ = true;
    log_step("[done] all indexes loaded");
}

std::vector<HaplotypeCoverage>
Index::translatable_haplotypes_scored(const std::string& src_haplotype,
                                      int64_t start, int64_t end,
                                      double min_coverage,
                                      size_t max_nodes) const {
    if (!loaded_)
        throw std::runtime_error("Index::translatable_haplotypes_scored called before load()");
    if (end - start > MAX_INTERVAL_LENGTH)
        throw std::invalid_argument(
            "Interval length " + std::to_string(end - start) +
            " exceeds maximum of " + std::to_string(MAX_INTERVAL_LENGTH) + " bases");
    if (start < 0 || end < 0 || start > end)
        throw std::invalid_argument("Invalid interval [" +
            std::to_string(start) + ", " + std::to_string(end) + "]");

    std::vector<HaplotypeCoverage> out;

    // Inverse of SampledTagArray::encode_value(node_id, is_rev) =
    //   1 + (((node_id - 1) << 1) | is_rev)
    auto decode_tag = [](uint64_t code, bool& is_rev) -> int64_t {
        if (code == 0) return 0;             // 0 is reserved for gaps
        const uint64_t v = code - 1;
        is_rev = (v & 1ULL) != 0;
        return static_cast<int64_t>(v >> 1) + 1;
    };

    FastLocate& rindex = const_cast<FastLocate&>(rlbwt_rindex_);
    SampledTagArray& sampled = const_cast<SampledTagArray&>(sampled_);
    const gbwt::GBWT& gbwt_index = gbz_->index;
    const gbwtgraph::GBWTGraph& graph = gbz_->graph;
    const bool have_meta = gbwt_index.hasMetadata() &&
                           gbwt_index.metadata.hasPathNames() &&
                           gbwt_index.metadata.hasSampleNames();

    std::vector<std::string> source_path_names =
        path_names_for_haplotype(table1_, src_haplotype);
    if (source_path_names.empty())
        throw std::invalid_argument(
            "No paths found for source haplotype: " + src_haplotype);

    std::unordered_map<size_t, std::string> pid_to_hap;   // resolved lazily
    std::unordered_map<std::string, uint64_t> covered;
    std::unordered_set<std::string> here;
    uint64_t total_bp = 0;

    for (const std::string& name : source_path_names) {
        // lookup() end is EXCLUSIVE; the API is half-open, so pass `end` as-is.
        for (const PathInterval& pi : table1_.lookup(name,
                                                     static_cast<size_t>(start),
                                                     static_cast<size_t>(end))) {
            if (pi.end <= pi.start) continue;
            const size_t src_seq_id = 2 * pi.path_id;

            // Nodes the source visits here, unscoped by target.
            std::vector<TagInfo> tags = find_tags_in_interval(
                rindex, sampled, src_seq_id, pi.start, pi.end - 1,
                &gbwt_index, gbwt_rindex_.get(), &graph,
                std::numeric_limits<size_t>::max(), nullptr);
            if (tags.empty()) continue;

            // Optionally subsample: probing every node is exact but each probe
            // enumerates a node's whole pangenome usage, so wide intervals can
            // trade a little precision for a lot of speed.
            size_t stride = 1;
            if (max_nodes > 0 && tags.size() > max_nodes) {
                stride = (tags.size() + max_nodes - 1) / max_nodes;
            }

            for (size_t i = 0; i < tags.size(); i += stride) {
                const TagInfo& tag = tags[i];
                if (tag.source_offsets.empty()) continue;
                // source_offsets holds one entry per VISIT to this node (the tag
                // array is run-length encoded, so the walk steps run-by-run), not
                // one per base. Weight each visit by the node's length to get a
                // real base count — counting visits reports ~1/31 of the truth on
                // a graph whose nodes average ~31 bp.
                bool tag_rev = false;
                const int64_t nid = decode_tag(tag.tag_code, tag_rev);
                if (nid < 1) continue;
                const uint64_t node_bp = static_cast<uint64_t>(
                    graph.get_length(graph.get_handle(nid, tag_rev)));
                if (node_bp == 0) continue;
                uint64_t bp = node_bp * tag.source_offsets.size();
                bp *= stride;                 // a sampled node stands for its stride
                total_bp += bp;

                here.clear();
                for (const NodeVisit& v : find_sequences_for_tag(rindex, sampled, tag.tag_code)) {
                    const size_t pid = v.seq_id / 2;
                    auto it = pid_to_hap.find(pid);
                    if (it == pid_to_hap.end()) {
                        std::string hn;
                        if (have_meta && pid < gbwt_index.metadata.paths()) {
                            gbwt::PathName pn = gbwt_index.metadata.path(pid);
                            hn = gbwt_index.metadata.sample(pn.sample) + "#" +
                                 std::to_string(pn.phase);
                        } else {
                            hn = "path_" + std::to_string(pid);
                        }
                        it = pid_to_hap.emplace(pid, std::move(hn)).first;
                    }
                    here.insert(it->second);
                }
                for (const std::string& h : here) covered[h] += bp;
            }
        }
    }

    if (total_bp == 0) return out;
    out.reserve(covered.size());
    for (const auto& kv : covered) {
        double pct = 100.0 * static_cast<double>(kv.second) /
                             static_cast<double>(total_bp);
        if (pct > 100.0) pct = 100.0;        // stride rounding
        if (pct < min_coverage) continue;
        HaplotypeCoverage hc;
        hc.haplotype = kv.first;
        hc.covered_bp = kv.second;
        hc.coverage = pct;
        out.push_back(std::move(hc));
    }
    std::sort(out.begin(), out.end(),
              [](const HaplotypeCoverage& a, const HaplotypeCoverage& b) {
                  if (a.covered_bp != b.covered_bp) return a.covered_bp > b.covered_bp;
                  return a.haplotype < b.haplotype;
              });
    return out;
}

std::vector<TranslatedInterval>
Index::translate_no_table2(const std::string& src_haplotype,
                           int64_t start, int64_t end,
                           const std::string& tgt_haplotype,
                           double timeout_ms,
                           bool* timed_out,
                           TranslationDiagnostics* diag) const {
    if (timed_out) *timed_out = false;
    // Cap the per-fragment detail: a whole-chromosome request can split into
    // thousands of fragments and the aggregate counters answer most questions.
    constexpr size_t MAX_FRAGMENT_DETAIL = 500;
    if (!loaded_)
        throw std::runtime_error("Index::translate_no_table2 called before load()");
    if (end - start > MAX_INTERVAL_LENGTH)
        throw std::invalid_argument(
            "Interval length " + std::to_string(end - start) +
            " exceeds maximum of " + std::to_string(MAX_INTERVAL_LENGTH) + " bases");
    if (start < 0 || end < 0 || start > end)
        throw std::invalid_argument("Invalid interval [" +
            std::to_string(start) + ", " + std::to_string(end) + "]");

    // How many nodes to probe from each end of the source interval before
    // concluding the target is not present. A homologous target is normally hit
    // on the first probe; the cap only bounds the cost of a query whose target
    // genuinely shares nothing here.
    static const size_t probe_cap = []() -> size_t {
        const char* e = std::getenv("PANGENOME_TRANSLATE_PROBE_CAP");
        if (e) { long v = std::atol(e); if (v > 0) return static_cast<size_t>(v); }
        // Discovery probes now use the cheap GBWT gate rather than the RLBWT
        // enumeration, so a larger budget is affordable and buys completeness on
        // wide intervals.
        return 1024;
    }();

    // Target may be a haplotype ("HG002#1") or a full contig ("HG002#1#chr1").
    std::string tgt_key = tgt_haplotype, tgt_contig_filter;
    {
        size_t h1 = tgt_haplotype.find('#');
        if (h1 != std::string::npos) {
            size_t h2 = tgt_haplotype.find('#', h1 + 1);
            if (h2 != std::string::npos) {
                tgt_key = tgt_haplotype.substr(0, h2);
                tgt_contig_filter = tgt_haplotype;
            }
        }
    }
    // Split the two-field key so candidate paths can be tested against GBWT
    // metadata directly — with hundreds of millions of path fragments, listing
    // a haplotype's paths up front is not an option.
    std::string tgt_sample = tgt_key;
    unsigned tgt_phase = 0;
    {
        size_t h = tgt_key.rfind('#');
        if (h != std::string::npos) {
            tgt_sample = tgt_key.substr(0, h);
            tgt_phase = static_cast<unsigned>(std::atoi(tgt_key.c_str() + h + 1));
        }
    }
    const gbwt::GBWT& gbwt_index = gbz_->index;
    const bool have_meta = gbwt_index.hasMetadata() &&
                           gbwt_index.metadata.hasPathNames() &&
                           gbwt_index.metadata.hasSampleNames();
    if (!have_meta)
        throw std::runtime_error("translate_no_table2 requires GBWT path/sample metadata");

    auto path_is_target = [&](size_t pid) -> bool {
        if (pid >= gbwt_index.metadata.paths()) return false;
        gbwt::PathName pn = gbwt_index.metadata.path(pid);
        if (static_cast<unsigned>(pn.phase) != tgt_phase) return false;
        if (gbwt_index.metadata.sample(pn.sample) != tgt_sample) return false;
        if (!tgt_contig_filter.empty()) {
            auto itn = path_to_global_.find(pid);
            if (itn == path_to_global_.end() ||
                itn->second.first != tgt_contig_filter) return false;
        }
        return true;
    };

    FastLocate& rindex = const_cast<FastLocate&>(rlbwt_rindex_);
    SampledTagArray& sampled = const_cast<SampledTagArray&>(sampled_);
    const gbwt::GBWT* gbwt_index_ptr = &gbz_->index;
    const gbwtgraph::GBWTGraph& graph = gbz_->graph;

    std::vector<std::string> source_path_names =
        path_names_for_haplotype(table1_, src_haplotype);
    if (source_path_names.empty())
        throw std::invalid_argument(
            "No paths found for source haplotype: " + src_haplotype);

    std::vector<PathInterval> source_intervals;
    for (const std::string& name : source_path_names) {
        // lookup() end is EXCLUSIVE; the API is half-open, so pass `end` as-is.
        for (PathInterval& pi : table1_.lookup(name,
                                               static_cast<size_t>(start),
                                               static_cast<size_t>(end))) {
            source_intervals.push_back(pi);
        }
    }
    if (source_intervals.empty()) return {};

    struct HaplotypeTranslation {
        size_t source_haplotype_offset = 0;
        size_t target_haplotype_offset = 0;
        size_t target_path_id = 0;
    };
    std::vector<HaplotypeTranslation> all_raw;

    // Cooperative deadline. A C++ call holds the GIL for its whole duration, so
    // no caller-side timeout can interrupt it — the only way to bound one slow
    // target is to check the clock at safe points and stop ourselves.
    const auto deadline_t0 = std::chrono::steady_clock::now();
    bool hit_deadline = false;
    auto expired = [&]() -> bool {
        if (timeout_ms <= 0.0 || hit_deadline) return hit_deadline;
        if (std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - deadline_t0).count() > timeout_ms) {
            hit_deadline = true;
        }
        return hit_deadline;
    };

    for (const PathInterval& pi : source_intervals) {
        if (expired()) break;                       // between source fragments
        FragmentDiag fd;
        const size_t raw_before = all_raw.size();
        size_t frag_min_src = std::numeric_limits<size_t>::max(), frag_max_src = 0;
        if (diag) {
            diag->fragments++;
            fd.src_path_id = pi.path_id;
            fd.extent_start = pi.start;
            fd.extent_end = pi.end;
            fd.extent_bp = (pi.end > pi.start) ? (pi.end - pi.start) : 0;
            diag->requested_bp += fd.extent_bp;
        }
        auto close_fragment = [&]() {
            if (!diag) return;
            fd.points = all_raw.size() - raw_before;
            if (fd.points) {
                for (size_t k = raw_before; k < all_raw.size(); ++k) {
                    const size_t so = all_raw[k].source_haplotype_offset;
                    frag_min_src = std::min(frag_min_src, so);
                    frag_max_src = std::max(frag_max_src, so);
                }
                fd.mapped_span = frag_max_src - frag_min_src + 1;
                diag->productive++;
            }
            diag->mapped_bp += fd.points;
            if (fd.scoped_tags == 1) diag->single_anchor++;
            if (diag->detail.size() < MAX_FRAGMENT_DETAIL) diag->detail.push_back(fd);
        };
        const size_t src_path_id = pi.path_id;
        const size_t local_start = pi.start;
        const size_t local_end = pi.end;
        if (local_end <= local_start) continue;
        const size_t src_seq_id = 2 * src_path_id;

        auto it_src = path_to_global_.find(src_path_id);
        if (it_src == path_to_global_.end()) continue;
        const size_t src_subpath_start = it_src->second.second;

        const size_t seq_start_incl = local_start;
        const size_t seq_end_incl = local_end - 1;
        size_t probe_pad = 0;

        // 1) Nodes the source visits in this interval, unscoped by target.
        std::vector<TagInfo> all_tags = find_tags_in_interval(
            rindex, sampled, src_seq_id, seq_start_incl, seq_end_incl,
            gbwt_index_ptr, gbwt_rindex_.get(), &graph,
            std::numeric_limits<size_t>::max(), nullptr);
        if (diag) fd.unscoped_tags = all_tags.size();
        if (all_tags.empty()) { if (diag) { diag->no_tags++; close_fragment(); } continue; }

        // 2) Walk inward from both ends until a node the target also visits.
        //    find_tags_in_interval returns tags sorted by source offset, so the
        //    first hit going forward is the first common node and the first hit
        //    going backward is the last common node. Their target paths are the
        //    only candidates worth tracing — no table, no path enumeration.
        // candidate target path -> [min, max] source offset where it was observed.
        // Each target fragment covers only part of the interval; remembering where
        // lets each candidate be traced over its OWN window instead of the whole
        // request, which is the difference between linear and quadratic cost.
        std::unordered_map<size_t, std::pair<size_t, size_t>> candidates;

        // Discovery uses the GBWT decompressSA gate, NOT find_sequences_for_tag:
        // the latter enumerates a node's entire pangenome-wide visit list (the
        // ~1ms primitive that dominated the anchor builder), while discovery only
        // needs sequence ids. Nodes visited by an implausible number of sequences
        // are REPEATS: they belong to fragments from all over the haplotype, so
        // they invent candidates that share no synteny here and stretch a
        // candidate's window across the whole request. Skipping them keeps both
        // the candidate set and the windows honest.
        static const size_t repeat_visits = []() -> size_t {
            const char* e = std::getenv("PANGENOME_PROBE_REPEAT_LIMIT");
            if (e) { long v = std::atol(e); if (v > 0) return static_cast<size_t>(v); }
            return 1024;
        }();
        auto probe_tag = [&](const TagInfo& tag) -> bool {
            // Inverse of SampledTagArray::encode_value.
            if (tag.tag_code == 0) return false;
            const uint64_t code = tag.tag_code - 1;
            const bool tag_rev = (code & 1ULL) != 0;
            const int64_t nid = static_cast<int64_t>(code >> 1) + 1;

            size_t lo = std::numeric_limits<size_t>::max(), hi = 0;
            for (size_t off : tag.source_offsets) { lo = std::min(lo, off); hi = std::max(hi, off); }
            if (lo == std::numeric_limits<size_t>::max()) { lo = hi = 0; }

            bool hit = false;
            for (int flip = 0; flip < 2; ++flip) {
                gbwt::node_type node = gbwt::Node::encode(nid, flip ? !tag_rev : tag_rev);
                std::vector<gbwt::size_type> sa = gbwt_rindex_->decompressSA(node);
                if (sa.size() > repeat_visits) continue;      // repeat: unusable anchor
                for (gbwt::size_type v : sa) {
                    const size_t pid = static_cast<size_t>(gbwt_rindex_->seqId(v)) / 2;
                    if (!path_is_target(pid)) continue;
                    auto it = candidates.find(pid);
                    if (it == candidates.end()) candidates.emplace(pid, std::make_pair(lo, hi));
                    else {
                        it->second.first  = std::min(it->second.first, lo);
                        it->second.second = std::max(it->second.second, hi);
                    }
                    hit = true;
                }
            }
            return hit;
        };

        {
            const size_t n = all_tags.size();
            // Space probes by SOURCE DISTANCE, not by a fixed count. A fixed count
            // over-probes a small interval (10 kb needs ~1 fragment but paid 256
            // probes) and under-probes a large one (1 Mb strided ~3.9 kb, so any
            // target fragment shorter than that was missed and its bases lost).
            // Aim for one probe per `probe_spacing` bases of source, which is what
            // actually determines the shortest fragment we can still find.
            static const size_t probe_spacing = []() -> size_t {
                const char* e = std::getenv("PANGENOME_PROBE_SPACING_BP");
                if (e) { long v = std::atol(e); if (v > 0) return static_cast<size_t>(v); }
                return 2000;
            }();
            const size_t extent_bp = seq_end_incl - seq_start_incl + 1;
            size_t wanted = extent_bp / probe_spacing + 1;
            if (wanted < 4) wanted = 4;
            if (wanted > probe_cap) wanted = probe_cap;
            const size_t stride = (n > wanted) ? (n / wanted) : 1;
            probe_pad = (extent_bp * stride) / (n ? n : 1) + 1;
            for (size_t i = 0; i < n; i += stride) {
                if (expired()) break;
                probe_tag(all_tags[i]);             // collect ALL, do not stop
            }
            if (!hit_deadline && n) {
                probe_tag(all_tags[0]);
                probe_tag(all_tags[n - 1]);
            }
        }
        if (hit_deadline) break;

        if (diag) fd.candidates = static_cast<uint32_t>(candidates.size());
        if (candidates.empty()) {           // target shares nothing here
            if (diag) { diag->no_candidates++; close_fragment(); }
            continue;
        }
        if (diag) diag->traced++;

        // 3) Trace each candidate with the standard pipeline. The extent is the
        //    whole query interval: without Table 2 there is no precomputed
        //    homologous sub-range to narrow it to, and the common-node search
        //    below keeps only what source and target actually share.
        for (const auto& cand_entry : candidates) {
            if (expired()) break;                   // between candidate targets
            const size_t tgt_path_id = cand_entry.first;
            const size_t tgt_seq_id = 2 * tgt_path_id;

            // Trace only the window this target fragment was actually seen in,
            // padded by the probe resolution. Running every candidate over the
            // whole interval made the work candidates x interval; this makes the
            // total proportional to the interval, since the fragments partition it.
            const size_t win_lo = (cand_entry.second.first > seq_start_incl + probe_pad)
                                  ? cand_entry.second.first - probe_pad : seq_start_incl;
            const size_t win_hi = std::min(seq_end_incl, cand_entry.second.second + probe_pad);
            if (win_hi < win_lo) continue;

            std::vector<TagInfo> tags = find_tags_in_interval(
                rindex, sampled, src_seq_id, win_lo, win_hi,
                gbwt_index_ptr, gbwt_rindex_.get(), &graph, tgt_seq_id, nullptr);
            if (diag) fd.scoped_tags = std::max<uint64_t>(fd.scoped_tags, tags.size());
            if (tags.empty()) continue;

            CommonNodes common = find_first_and_last_common_nodes_gbwt(
                *gbwt_rindex_, rindex, sampled, tags, src_seq_id, tgt_seq_id);
            if (!common.found) continue;
            if (diag) {
                fd.first_source_base = common.first_source_base;
                fd.first_target_base = common.first_target_base;
                fd.last_source_base  = common.last_source_base;
                fd.last_target_base  = common.last_target_base;
                fd.first_unique = common.first_is_unique;
                fd.last_unique  = common.last_is_unique;
            }

            std::vector<TranslationResult> trans = trace_coordinates_gbwt(
                *gbwt_index_ptr, *gbwt_rindex_, graph,
                src_seq_id, win_lo, win_hi, tgt_seq_id,
                common.first_source_offset, common.first_target_offset,
                common.first_source_base, common.first_target_base,
                common.first_tag_code,
                common.last_source_base, common.last_target_base,
                common.last_tag_code);

            auto it_tgt = path_to_global_.find(tgt_path_id);
            if (it_tgt == path_to_global_.end()) continue;
            const size_t tgt_subpath_start = it_tgt->second.second;

            for (const TranslationResult& tr : trans) {
                if (tr.target_offset == 0) continue;
                HaplotypeTranslation ht;
                ht.source_haplotype_offset = src_subpath_start + tr.source_offset;
                ht.target_haplotype_offset = tgt_subpath_start + tr.target_offset;
                ht.target_path_id = tgt_path_id;
                all_raw.push_back(ht);
            }
        }
        if (diag) {
            if (all_raw.size() == raw_before) diag->empty_trace++;
            close_fragment();
        }
    }

    std::sort(all_raw.begin(), all_raw.end(),
              [](const HaplotypeTranslation& a, const HaplotypeTranslation& b) {
                  if (a.source_haplotype_offset != b.source_haplotype_offset)
                      return a.source_haplotype_offset < b.source_haplotype_offset;
                  return a.target_haplotype_offset < b.target_haplotype_offset;
              });

    if (timed_out) *timed_out = hit_deadline;

    std::vector<TranslatedInterval> results;
    results.reserve(all_raw.size());
    for (const HaplotypeTranslation& ht : all_raw) {
        TranslatedInterval ti;
        auto it_name = path_to_global_.find(ht.target_path_id);
        ti.haplotype = (it_name != path_to_global_.end()) ? it_name->second.first
                                                          : tgt_haplotype;
        ti.start = static_cast<int64_t>(ht.source_haplotype_offset);
        ti.end   = static_cast<int64_t>(ht.target_haplotype_offset);
        ti.strand = '+';
        results.push_back(ti);
    }
    return results;
}

DiagnosedTranslation
Index::translate_diagnosed(const std::string& src_haplotype,
                           int64_t start, int64_t end,
                           const std::string& tgt_haplotype) const {
    DiagnosedTranslation out;
    const auto t0 = std::chrono::steady_clock::now();
    // Every entry point routes through translate(), so a loaded Table 2 is used
    // here too. The per-fragment counters describe the table-free path's probing
    // and are not collected on the Table 2 path; diagnostics.table2_path says
    // which ran, so all-zero counters are not mistaken for "found nothing".
    // Set PANGENOME_TRANSLATE_NO_T2=1 to get the counters back.
    static const bool no_t2_env = (std::getenv("PANGENOME_TRANSLATE_NO_T2") != nullptr);
    if (no_t2_env || !has_table2_) {
        out.intervals = translate_no_table2(src_haplotype, start, end, tgt_haplotype,
                                            0.0, nullptr, &out.diagnostics);
    } else {
        out.diagnostics.table2_path = true;
        out.intervals = translate(src_haplotype, start, end, tgt_haplotype);
    }
    out.elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    return out;
}

TranslationRun
Index::translate_checked(const std::string& src_haplotype,
                         int64_t start, int64_t end,
                         const std::string& tgt_haplotype,
                         double timeout_ms) const {
    TranslationRun run;
    const auto t0 = std::chrono::steady_clock::now();
    bool timed_out = false;
    // Goes through translate(), so it uses Table 2 whenever one is loaded. It
    // previously called translate_no_table2 directly, which meant the serving
    // path — the API always passes a timeout, so it always lands here — could
    // never use a Table 2 no matter what was loaded.
    run.intervals = translate(src_haplotype, start, end, tgt_haplotype,
                              timeout_ms, &timed_out);
    run.timed_out = timed_out;
    run.elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    return run;
}

std::vector<TranslatedInterval>
Index::translate(const std::string& src_haplotype,
                 int64_t start, int64_t end,
                 const std::string& tgt_haplotype,
                 double timeout_ms,
                 bool* timed_out) const {
    if (!loaded_)
        throw std::runtime_error("Index::translate called before load()");
    if (timed_out) *timed_out = false;

    // Table 2 is the default whenever one is loaded. The table-free path (via
    // first/last common node through the GBWT/tag array) is the fallback when
    // there is no Table 2, and the escape hatch when PANGENOME_TRANSLATE_NO_T2
    // is set — useful for A/B comparison and if a table turns out to be wrong.
    static const bool no_t2_env = (std::getenv("PANGENOME_TRANSLATE_NO_T2") != nullptr);
    if (no_t2_env || !has_table2_) {
        return translate_no_table2(src_haplotype, start, end, tgt_haplotype,
                                   timeout_ms, timed_out);
    }

    // Cooperative deadline, same contract as translate_no_table2: checked
    // between source fragments and between candidate target paths, so a
    // pathological target is abandoned rather than stalling the request. It
    // cannot preempt a single long call, so the stop can overshoot slightly.
    const auto deadline_t0 = std::chrono::steady_clock::now();
    bool hit_deadline = false;
    auto past_deadline = [&]() {
        if (timeout_ms <= 0.0 || hit_deadline) return hit_deadline;
        if (std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - deadline_t0).count() > timeout_ms) {
            hit_deadline = true;
        }
        return hit_deadline;
    };
    if (end - start > MAX_INTERVAL_LENGTH)
        throw std::invalid_argument(
            "Interval length " + std::to_string(end - start) +
            " exceeds maximum of " + std::to_string(MAX_INTERVAL_LENGTH) + " bases");
    if (start < 0 || end < 0 || start > end)
        throw std::invalid_argument("Invalid interval [" +
            std::to_string(start) + ", " + std::to_string(end) + "]");

    size_t global_start = static_cast<size_t>(start);
    size_t global_end   = static_cast<size_t>(end);

    // The target may be a haplotype ("GRCh38#0") or a full contig
    // ("GRCh38#0#chr10"). T2 is keyed by the 2-field haplotype name, so we
    // always look up with the first two '#'-fields; when a contig is named,
    // we additionally keep only target paths on that contig.
    std::string tgt_key = tgt_haplotype;
    std::string tgt_contig_filter;
    {
        size_t h1 = tgt_haplotype.find('#');
        if (h1 != std::string::npos) {
            size_t h2 = tgt_haplotype.find('#', h1 + 1);
            if (h2 != std::string::npos) {
                tgt_key = tgt_haplotype.substr(0, h2);   // e.g. "GRCh38#0"
                tgt_contig_filter = tgt_haplotype;        // e.g. "GRCh38#0#chr10"
            }
        }
    }

    // The underlying query functions take non-const references (they read, not write,
    // but were not declared const).  const_cast is safe here because those functions
    // do not mutate the objects after construction.
    FastLocate& rindex   = const_cast<FastLocate&>(rlbwt_rindex_);
    SampledTagArray& sampled = const_cast<SampledTagArray&>(sampled_);

    std::vector<std::string> source_path_names =
        path_names_for_haplotype(table1_, src_haplotype);
    if (source_path_names.empty())
        throw std::invalid_argument(
            "No paths found for source haplotype: " + src_haplotype);

    std::vector<PathInterval> source_intervals;
    for (const std::string& name : source_path_names) {
        std::vector<PathInterval> pis =
            // lookup() end is EXCLUSIVE and the API is half-open: pass as-is.
            table1_.lookup(name, global_start, global_end);
        for (PathInterval& pi : pis)
            source_intervals.push_back(pi);
    }
    if (source_intervals.empty())
        return {};

    const gbwt::GBWT* gbwt_index_ptr = &gbz_->index;
    const gbwtgraph::GBWTGraph& graph = gbz_->graph;

    struct HaplotypeTranslation {
        size_t source_haplotype_offset = 0;
        size_t target_haplotype_offset = 0;
        size_t target_path_id = 0;
    };
    std::vector<HaplotypeTranslation> all_raw;

    for (const PathInterval& pi : source_intervals) {
        if (past_deadline()) break;
        size_t src_path_id = pi.path_id;
        size_t local_start = pi.start;
        size_t local_end   = pi.end;
        if (local_end <= local_start) continue;

        std::vector<TargetInterval> tgt_intervals =
            table2_.lookup(src_path_id, tgt_key, local_start, local_end);
        if (tgt_intervals.empty()) continue;

        size_t src_seq_id = 2 * src_path_id;
        auto it_src = path_to_global_.find(src_path_id);
        if (it_src == path_to_global_.end()) continue;
        size_t src_subpath_start = it_src->second.second;

        std::unordered_set<size_t> distinct_tgt_paths;
        for (const TargetInterval& ti : tgt_intervals)
            distinct_tgt_paths.insert(ti.tgt_path_id);

        // If a specific target contig was named, drop target paths on any
        // other contig (e.g. paralogous hits on a different chromosome).
        if (!tgt_contig_filter.empty()) {
            std::unordered_set<size_t> filtered;
            for (size_t pid : distinct_tgt_paths) {
                auto itn = path_to_global_.find(pid);
                if (itn != path_to_global_.end() &&
                    itn->second.first == tgt_contig_filter) {
                    filtered.insert(pid);
                }
            }
            distinct_tgt_paths.swap(filtered);
        }
        if (distinct_tgt_paths.empty()) continue;

        std::vector<IntervalMapping> segs =
            table2_.segments(src_path_id, tgt_key);

        for (size_t tgt_path_id : distinct_tgt_paths) {
            if (past_deadline()) break;
            size_t extent_start = local_end, extent_end = local_start;
            for (const IntervalMapping& m : segs) {
                if (m.tgt_path_id != tgt_path_id) continue;
                size_t overlap_start = std::max(local_start, m.src_start);
                size_t overlap_end   = std::min(local_end,   m.src_end);
                if (overlap_start >= overlap_end) continue;
                extent_start = std::min(extent_start, overlap_start);
                extent_end   = std::max(extent_end,   overlap_end);
            }
            if (extent_start >= extent_end) continue;

            size_t extent_start_incl = extent_start;
            size_t extent_end_incl   = extent_end - 1;
            size_t tgt_seq_id = 2 * tgt_path_id;

            std::vector<TagInfo> tags = find_tags_in_interval(
                rindex, sampled, src_seq_id,
                extent_start_incl, extent_end_incl,
                gbwt_index_ptr, gbwt_rindex_.get(), &graph,
                tgt_seq_id, nullptr);
            if (tags.empty()) continue;

            CommonNodes common = find_first_and_last_common_nodes_gbwt(
                *gbwt_rindex_, rindex, sampled, tags,
                src_seq_id, tgt_seq_id);
            if (!common.found) continue;

            std::vector<TranslationResult> trans = trace_coordinates_gbwt(
                *gbwt_index_ptr, *gbwt_rindex_, graph,
                src_seq_id, extent_start_incl, extent_end_incl,
                tgt_seq_id,
                common.first_source_offset, common.first_target_offset,
                common.first_source_base,   common.first_target_base,
                common.first_tag_code,
                common.last_source_base,    common.last_target_base,
                common.last_tag_code);

            auto it_tgt = path_to_global_.find(tgt_path_id);
            if (it_tgt == path_to_global_.end()) continue;
            size_t tgt_subpath_start = it_tgt->second.second;

            for (const TranslationResult& tr : trans) {
                if (tr.target_offset == 0) continue;
                HaplotypeTranslation ht;
                ht.source_haplotype_offset = src_subpath_start + tr.source_offset;
                ht.target_haplotype_offset = tgt_subpath_start + tr.target_offset;
                ht.target_path_id = tgt_path_id;
                all_raw.push_back(ht);
            }
        }
    }

    std::sort(all_raw.begin(), all_raw.end(),
              [](const HaplotypeTranslation& a, const HaplotypeTranslation& b) {
                  if (a.source_haplotype_offset != b.source_haplotype_offset)
                      return a.source_haplotype_offset < b.source_haplotype_offset;
                  return a.target_haplotype_offset < b.target_haplotype_offset;
              });

    std::vector<TranslatedInterval> results;
    results.reserve(all_raw.size());
    for (const HaplotypeTranslation& ht : all_raw) {
        TranslatedInterval ti;
        // Report the actual target CONTIG the position resolved to (e.g.
        // "GRCh38#0#chr10"), not just the queried haplotype ("GRCh38#0").
        // The contig is known from the resolved target path id; fall back to
        // the queried haplotype name if it isn't in the path→name map.
        auto it_name = path_to_global_.find(ht.target_path_id);
        ti.haplotype = (it_name != path_to_global_.end())
                       ? it_name->second.first
                       : tgt_haplotype;
        ti.start     = static_cast<int64_t>(ht.source_haplotype_offset);
        ti.end       = static_cast<int64_t>(ht.target_haplotype_offset);
        ti.strand    = '+';
        results.push_back(ti);
    }
    if (timed_out) *timed_out = hit_deadline;
    return results;
}

std::vector<std::string>
Index::translatable_haplotypes(const std::string& src_haplotype,
                               int64_t start, int64_t end) const {
    if (!loaded_)
        throw std::runtime_error("Index::translatable_haplotypes called before load()");

    // Without Table 2 there is nothing to look up: derive the same answer from
    // the graph instead, and return just the names for API compatibility.
    if (!has_table2_) {
        std::vector<std::string> names;
        for (const HaplotypeCoverage& hc :
                 translatable_haplotypes_scored(src_haplotype, start, end)) {
            names.push_back(hc.haplotype);
        }
        std::sort(names.begin(), names.end());
        return names;
    }
    if (end - start > MAX_INTERVAL_LENGTH)
        throw std::invalid_argument(
            "Interval length " + std::to_string(end - start) +
            " exceeds maximum of " + std::to_string(MAX_INTERVAL_LENGTH) + " bases");
    if (start < 0 || end < 0 || start > end)
        throw std::invalid_argument("Invalid interval [" +
            std::to_string(start) + ", " + std::to_string(end) + "]");

    size_t s = static_cast<size_t>(start);
    size_t e = static_cast<size_t>(end);

    // Resolve the source contig/haplotype to its GBWT path(s) and map the query
    // interval into path-local coordinates (same front half as translate()).
    std::vector<std::string> source_path_names =
        path_names_for_haplotype(table1_, src_haplotype);
    if (source_path_names.empty())
        throw std::invalid_argument(
            "No paths found for source haplotype: " + src_haplotype);

    std::vector<PathInterval> source_intervals;
    for (const std::string& name : source_path_names) {
        // lookup() end is EXCLUSIVE and the API is half-open: pass as-is.
        std::vector<PathInterval> pis = table1_.lookup(name, s, e);
        for (PathInterval& pi : pis)
            source_intervals.push_back(pi);
    }
    if (source_intervals.empty())
        return {};

    // Candidate target haplotypes per source path come straight from T2's keys.
    // T2 is sparse — it only stores (src_path_id, tgt_haplotype) pairs that share
    // at least one graph node — so this is the homology-pruned candidate set, not
    // every haplotype in the graph.
    // Candidates for THIS source path only. Building a map of every key in the
    // table (what this used to do via keys()) copies one string per key — tens
    // of millions on an all-pairs table — on every single request, which made
    // this endpoint far slower than the translation it precedes.
    std::unordered_set<std::string> found;
    for (const PathInterval& pi : source_intervals) {
        if (pi.end <= pi.start) continue;
        for (const std::string& tgt : table2_.target_haplotypes_for(pi.path_id)) {
            if (found.count(tgt)) continue;  // already confirmed via another segment
            std::vector<TargetInterval> hits =
                table2_.lookup(pi.path_id, tgt, pi.start, pi.end);
            if (!hits.empty()) found.insert(tgt);
        }
    }

    std::vector<std::string> out(found.begin(), found.end());
    std::sort(out.begin(), out.end());
    return out;
}

// ── GAF parsing helpers ──────────────────────────────────────────────────

namespace {

/// Tab-split helper that preserves empty fields.
std::vector<std::string> split_tabs_strict(const std::string& s) {
    std::vector<std::string> out;
    size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == '\t') {
            out.push_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
    return out;
}

/// Parse the GAF path field (e.g. ">1>2<3>4") into (node_id, is_reverse) pairs.
/// Throws on malformed input. Returns empty vector for "*" (missing path).
std::vector<std::pair<int64_t, bool>>
parse_gaf_path_field(const std::string& path_str) {
    std::vector<std::pair<int64_t, bool>> result;
    if (path_str == "*" || path_str.empty()) return result;
    size_t i = 0;
    while (i < path_str.size()) {
        char c = path_str[i];
        if (c != '>' && c != '<') {
            // Stable-path form ("chr1:100-200" etc.) — not what giraffe emits.
            // Return empty so caller falls back to "no nodes" handling.
            return {};
        }
        bool is_rev = (c == '<');
        ++i;
        size_t j = i;
        while (j < path_str.size() &&
               path_str[j] != '>' && path_str[j] != '<') {
            ++j;
        }
        if (j == i) {
            throw std::runtime_error("Empty node id in GAF path: " + path_str);
        }
        std::string token = path_str.substr(i, j - i);
        // Reject stable-path step format ("name:start-end") — we only handle
        // segment IDs (positive integers) at this layer.
        if (token.find(':') != std::string::npos) return {};
        int64_t nid = std::stoll(token);
        result.emplace_back(nid, is_rev);
        i = j;
    }
    return result;
}

/// Parse a cg:Z: CIGAR string ("5M2D3I4M") into (op, len) pairs.
std::vector<std::pair<char, size_t>> parse_cg_cigar(const std::string& cg) {
    std::vector<std::pair<char, size_t>> ops;
    size_t i = 0;
    while (i < cg.size()) {
        size_t j = i;
        while (j < cg.size() && std::isdigit(static_cast<unsigned char>(cg[j]))) ++j;
        if (j == i || j >= cg.size()) break;
        size_t len = std::stoull(cg.substr(i, j - i));
        char op = cg[j];
        ops.emplace_back(op, len);
        i = j + 1;
    }
    return ops;
}

/// Parse a cs:Z: difference string into (op, len) pairs in CIGAR form, using
/// '=' for matches and 'X' for mismatches so the existing walk consumes it
/// unchanged. cs is what vg actually emits (alignment_to_gaf defaults
/// cs_cigar=true); it never writes cg:Z:.
///
///   :42     42 matching bases        -> ('=', 42)
///   *ac     mismatch, ref a, query c -> ('X', 1)
///   +ACGT   inserted in the query    -> ('I', 4)
///   -ACGT   deleted from the query   -> ('D', 4)
std::vector<std::pair<char, size_t>> parse_cs_cigar(const std::string& cs) {
    std::vector<std::pair<char, size_t>> ops;
    size_t i = 0;
    while (i < cs.size()) {
        const char tok = cs[i++];
        if (tok == ':') {
            size_t j = i;
            while (j < cs.size() && std::isdigit(static_cast<unsigned char>(cs[j]))) ++j;
            if (j == i) break;
            ops.emplace_back('=', std::stoull(cs.substr(i, j - i)));
            i = j;
        } else if (tok == '*') {
            // Exactly two bases: reference then query.
            if (i + 2 > cs.size()) break;
            ops.emplace_back('X', 1);
            i += 2;
        } else if (tok == '+' || tok == '-') {
            size_t j = i;
            while (j < cs.size() && std::isalpha(static_cast<unsigned char>(cs[j]))) ++j;
            if (j == i) break;
            ops.emplace_back(tok == '+' ? 'I' : 'D', j - i);
            i = j;
        } else {
            break;   // unrecognised: stop rather than misattribute
        }
    }
    return ops;
}

/// Convert a GAF line into per-node SourceMappings using the path field
/// and (preferred) the cg:Z: CIGAR. Falls back to proportional distribution
/// of the query interval across node lengths when CIGAR is missing.
///
/// On any parse failure returns an empty vector; caller can detect this and
/// return status="parse_error".
std::vector<panindexer::SourceMapping>
gaf_to_source_mappings(const std::string& gaf_str,
                       const gbwtgraph::GBWTGraph& graph) {
    auto fields = split_tabs_strict(gaf_str);
    if (fields.size() < 12) return {};

    // Required numeric fields (return empty on malformed input).
    size_t query_start, query_end, path_start, path_end_field;
    try {
        query_start     = std::stoull(fields[2]);
        query_end       = std::stoull(fields[3]);
        path_start      = std::stoull(fields[7]);
        path_end_field  = std::stoull(fields[8]);
    } catch (const std::exception&) {
        return {};
    }
    (void) query_end;
    (void) path_end_field;

    const std::string& path_str = fields[5];
    std::vector<std::pair<int64_t, bool>> path_nodes;
    try {
        path_nodes = parse_gaf_path_field(path_str);
    } catch (const std::exception&) {
        return {};
    }
    if (path_nodes.empty()) return {};

    // Look up node lengths in the graph (returns 0 if node is missing).
    std::vector<size_t> node_lengths;
    node_lengths.reserve(path_nodes.size());
    for (const auto& [nid, is_rev] : path_nodes) {
        if (!graph.has_node(static_cast<handlegraph::nid_t>(nid))) {
            return {};
        }
        auto handle = graph.get_handle(
            static_cast<handlegraph::nid_t>(nid), is_rev);
        node_lengths.push_back(graph.get_length(handle));
    }

    // Find the alignment string among the optional fields. vg emits cs:Z:, so
    // that is checked FIRST; the cg:Z: branch below predates this and never
    // fired against giraffe output, which silently sent every alignment down
    // the proportional-estimate fallback.
    std::string cg_str, cs_str;
    for (size_t k = 12; k < fields.size(); ++k) {
        if (cs_str.empty() && fields[k].size() > 5 &&
            fields[k].compare(0, 5, "cs:Z:") == 0) {
            cs_str = fields[k].substr(5);
        } else if (cg_str.empty() && fields[k].size() > 5 &&
                   fields[k].compare(0, 5, "cg:Z:") == 0) {
            cg_str = fields[k].substr(5);
        }
    }

    std::vector<panindexer::SourceMapping> result;
    result.reserve(path_nodes.size());

    if (!cs_str.empty() || !cg_str.empty()) {
        // ── CIGAR-driven walk ─────────────────────────────────────────────
        // cs carries '=' / 'X' so matches are countable; cg's 'M' is not.
        const bool have_matches = !cs_str.empty();
        auto cigar_ops = have_matches ? parse_cs_cigar(cs_str)
                                      : parse_cg_cigar(cg_str);

        size_t query_pos = query_start;
        size_t path_pos  = path_start;            // absolute path position
        size_t cigar_i   = 0;
        size_t cigar_rem = cigar_ops.empty() ? 0 : cigar_ops[0].second;
        size_t path_node_start = 0;               // path offset at start of current node

        for (size_t ni = 0; ni < path_nodes.size(); ++ni) {
            const auto& [nid, is_rev] = path_nodes[ni];
            size_t nlen = node_lengths[ni];
            size_t path_node_end = path_node_start + nlen;

            // Skip nodes the alignment hasn't reached yet (rare with proper GAF).
            if (path_node_end <= path_pos) {
                path_node_start = path_node_end;
                continue;
            }

            size_t node_offset = (path_pos >= path_node_start)
                                 ? (path_pos - path_node_start) : 0;
            size_t read_begin = query_pos;
            size_t path_at_node_start = path_pos;
            size_t node_matched = 0;

            // Walk CIGAR until we exhaust this node or run out of ops.
            while (path_pos < path_node_end && cigar_i < cigar_ops.size()) {
                char op = cigar_ops[cigar_i].first;
                bool consumes_path  = (op == 'M' || op == '=' || op == 'X'
                                       || op == 'D' || op == 'N');
                bool consumes_query = (op == 'M' || op == '=' || op == 'X'
                                       || op == 'I' || op == 'S');

                size_t steps;
                if (consumes_path) {
                    size_t path_left_in_node = path_node_end - path_pos;
                    steps = std::min(cigar_rem, path_left_in_node);
                } else {
                    // I / S / H / P — consume the rest of the op without
                    // advancing path. Attribute to this node.
                    steps = cigar_rem;
                }

                if (op == '=') node_matched += steps;
                if (consumes_path)  path_pos  += steps;
                if (consumes_query) query_pos += steps;

                cigar_rem -= steps;
                if (cigar_rem == 0) {
                    ++cigar_i;
                    if (cigar_i < cigar_ops.size()) {
                        cigar_rem = cigar_ops[cigar_i].second;
                    }
                }
            }

            panindexer::SourceMapping sm;
            sm.node_id = nid;
            sm.is_reverse = is_rev;
            sm.read_begin_offset = read_begin;
            sm.read_end_offset = query_pos;
            sm.node_offset_in_node = node_offset;
            sm.mapping_from_length = path_pos - path_at_node_start;
            sm.matched_bases = have_matches ? node_matched : 0;
            result.push_back(sm);

            path_node_start = path_node_end;
        }
    } else {
        // ── No CIGAR: distribute the query range proportionally to path coverage.
        size_t path_node_start = 0;
        size_t total_path_len = 0;
        for (size_t l : node_lengths) total_path_len += l;
        if (path_end_field == 0) path_end_field = total_path_len;
        size_t covered_total = (path_end_field > path_start)
                               ? (path_end_field - path_start) : 0;
        size_t query_total = (query_end > query_start)
                             ? (query_end - query_start) : 0;

        for (size_t ni = 0; ni < path_nodes.size(); ++ni) {
            const auto& [nid, is_rev] = path_nodes[ni];
            size_t nlen = node_lengths[ni];
            size_t node_end = path_node_start + nlen;

            if (node_end <= path_start || path_node_start >= path_end_field) {
                path_node_start = node_end;
                continue;
            }
            size_t cov_start = std::max(path_node_start, path_start);
            size_t cov_end   = std::min(node_end, path_end_field);
            size_t cov_len   = cov_end - cov_start;

            size_t read_b, read_e;
            if (covered_total > 0) {
                read_b = query_start + (cov_start - path_start) * query_total / covered_total;
                read_e = query_start + (cov_end   - path_start) * query_total / covered_total;
            } else {
                read_b = query_start;
                read_e = query_start;
            }

            panindexer::SourceMapping sm;
            sm.node_id = nid;
            sm.is_reverse = is_rev;
            sm.read_begin_offset = read_b;
            sm.read_end_offset = read_e;
            sm.node_offset_in_node = cov_start - path_node_start;
            sm.mapping_from_length = cov_len;
            result.push_back(sm);

            path_node_start = node_end;
        }
    }

    return result;
}

/// Translate panindexer::PrecomputedAnchor records into the portable
/// AnchorRecord struct used at the Python boundary.
AnchorRecord to_anchor_record(const panindexer::PrecomputedAnchor& a) {
    AnchorRecord r;
    r.source_mapping_begin = a.source_mapping_begin;
    r.source_mapping_end   = a.source_mapping_end;
    r.read_begin_offset    = a.read_begin_offset;
    r.read_end_offset      = a.read_end_offset;
    r.path_offset_step_begin = a.path_offset_step_begin;
    r.path_offset_step_end   = a.path_offset_step_end;
    r.gbwt_edge_begin_node   = a.gbwt_edge_begin.first;
    r.gbwt_edge_begin_offset = a.gbwt_edge_begin.second;
    r.gbwt_edge_end_node     = a.gbwt_edge_end.first;
    r.gbwt_edge_end_offset   = a.gbwt_edge_end.second;
    return r;
}

const char* status_token(panindexer::AnchorBuildResult::Status s) {
    using S = panindexer::AnchorBuildResult::Status;
    switch (s) {
        case S::Ok:             return "ok";
        case S::EmptyAlignment: return "empty_alignment";
        case S::UnknownPath:    return "unknown_path";
        case S::NoCommonNodes:  return "no_common_nodes";
    }
    return "unknown";
}

} // anonymous namespace

std::vector<HaplotypeCoverage>
Index::haplotype_coverage(const std::string& gaf_str, double min_coverage,
                          bool include_zero) const {
    if (!loaded_) {
        throw std::runtime_error("Index::haplotype_coverage called before load()");
    }

    std::vector<HaplotypeCoverage> out;
    auto mappings = gaf_to_source_mappings(gaf_str, gbz_->graph);
    if (mappings.empty()) return out;

    const gbwt::GBWT& gbwt_index = gbz_->index;
    const bool have_names = gbwt_index.hasMetadata() &&
                            gbwt_index.metadata.hasPathNames() &&
                            gbwt_index.metadata.hasSampleNames();

    // GBWT path id -> two-field haplotype name, resolved lazily: an alignment
    // touches only a small number of paths, so this avoids walking all of them.
    std::unordered_map<size_t, std::string> pid_to_hap;
    std::unordered_map<std::string, uint64_t> covered;
    std::unordered_map<std::string, uint64_t> matched;
    // Matches are knowable only when the GAF carried cs:Z:. One mapping with a
    // non-zero count proves the string was there; without it, identity is
    // unknown rather than zero.
    bool any_matches = false;
    std::unordered_set<std::string> here;
    uint64_t total_bp = 0;

    for (const panindexer::SourceMapping& m : mappings) {
        const uint64_t bp = (m.read_end_offset > m.read_begin_offset)
                            ? (m.read_end_offset - m.read_begin_offset) : 0;
        if (bp == 0) continue;   // pure insertion: no graph node to credit
        total_bp += bp;
        if (m.matched_bases > 0) any_matches = true;

        // Which haplotypes visit this node? Count a haplotype once per node even
        // if it visits repeatedly, and accept either orientation so that
        // reverse-strand alignments (and inversions) still score.
        here.clear();
        for (int flip = 0; flip < 2; ++flip) {
            const bool is_rev = (flip == 0) ? m.is_reverse : !m.is_reverse;
            gbwt::node_type node = gbwt::Node::encode(m.node_id, is_rev);
            std::vector<gbwt::size_type> sa = gbwt_rindex_->decompressSA(node);
            for (gbwt::size_type v : sa) {
                const size_t pid = static_cast<size_t>(gbwt_rindex_->seqId(v)) / 2;
                auto it = pid_to_hap.find(pid);
                if (it == pid_to_hap.end()) {
                    std::string name;
                    if (have_names && pid < gbwt_index.metadata.paths()) {
                        gbwt::PathName pn = gbwt_index.metadata.path(pid);
                        name = gbwt_index.metadata.sample(pn.sample) + "#" +
                               std::to_string(pn.phase);
                    } else {
                        name = "path_" + std::to_string(pid);
                    }
                    it = pid_to_hap.emplace(pid, std::move(name)).first;
                }
                here.insert(it->second);
            }
        }
        for (const std::string& h : here) {
            covered[h] += bp;
            // Same loop, same decompressSA result: identity for every haplotype
            // costs nothing beyond the coverage pass already being made.
            matched[h] += m.matched_bases;
        }
    }

    if (total_bp == 0) return out;

    // Optionally round the list out to every haplotype in the graph, so callers
    // that want a complete ranked table see non-overlapping ones as an explicit
    // 0 rather than as a silent absence.
    if (include_zero && min_coverage <= 0.0) {
        for (const std::string& name : get_haplotype_names()) {
            covered.emplace(name, 0);   // no-op where already scored
        }
    }

    out.reserve(covered.size());
    for (const auto& kv : covered) {
        const double pct = 100.0 * static_cast<double>(kv.second) /
                                   static_cast<double>(total_bp);
        if (pct < min_coverage) continue;
        HaplotypeCoverage hc;
        hc.haplotype = kv.first;
        hc.covered_bp = kv.second;
        hc.coverage = pct;
        if (any_matches) {
            auto mit = matched.find(kv.first);
            hc.matched_bp = (mit == matched.end()) ? 0 : mit->second;
            hc.identity = 100.0 * static_cast<double>(hc.matched_bp) /
                                  static_cast<double>(total_bp);
            hc.has_identity = true;
        }
        out.push_back(std::move(hc));
    }
    // Rank by identity when it is available: it is the number that answers
    // "how well does this haplotype match", where coverage only answers
    // "how much of the read is present at all".
    std::sort(out.begin(), out.end(),
              [](const HaplotypeCoverage& a, const HaplotypeCoverage& b) {
                  if (a.has_identity && b.has_identity && a.matched_bp != b.matched_bp)
                      return a.matched_bp > b.matched_bp;
                  if (a.covered_bp != b.covered_bp) return a.covered_bp > b.covered_bp;
                  return a.haplotype < b.haplotype;   // stable, deterministic order
              });
    return out;
}

AnchorBuildPyResult Index::build_surject_anchors(
    const std::string& gaf_str,
    const std::string& target_haplotype) const
{
    AnchorBuildPyResult out;

    if (!loaded_) {
        throw std::runtime_error("Index::build_surject_anchors called before load()");
    }

    auto source_mappings = gaf_to_source_mappings(gaf_str, gbz_->graph);
    if (source_mappings.empty()) {
        // Differentiate "couldn't parse anything" vs "alignment is empty":
        // peek at the GAF query name to choose. Either way we report empty.
        out.status = "parse_error";
        return out;
    }

    FastLocate& rindex = const_cast<FastLocate&>(rlbwt_rindex_);
    SampledTagArray& sampled = const_cast<SampledTagArray&>(sampled_);

    // Resolve `target_haplotype` to one or more T1 subpath ids.
    //   1) Exact-name match against T1.names() — covers the surjection case
    //      where the user passes a full path name (e.g. "GRCh38#0#chrM").
    //   2) Fall back to haplotype-prefix match (target + "#") — covers the
    //      translate-style case where the user passes a haplotype prefix
    //      (e.g. "HG002#1" → matches "HG002#1#chr1", "HG002#1#chr2", …).
    //
    // panindexer::build_surject_anchors only does (2), so calling it with a
    // full path name returns UnknownPath. We do the resolution here so both
    // forms work transparently.
    // Carry (path_id, path_length) together: SubpathInfo.length is the target
    // path's base length, precomputed at index-build time by the same
    // extract-and-sum loop build_surject_anchors_for_path would otherwise run.
    // Passing it avoids a full target-path extraction per query.
    std::vector<std::pair<size_t, size_t>> target_path_ids;
    {
        std::vector<std::string> t1_names = table1_.names();
        std::vector<std::string> matched_names;
        for (const auto& nm : t1_names) {
            if (nm == target_haplotype) {
                matched_names.push_back(nm);
            }
        }
        if (matched_names.empty()) {
            std::string prefix = target_haplotype;
            if (prefix.empty() || prefix.back() != '#') prefix += '#';
            for (const auto& nm : t1_names) {
                if (nm.size() >= prefix.size() &&
                    nm.compare(0, prefix.size(), prefix) == 0) {
                    matched_names.push_back(nm);
                }
            }
        }
        for (const auto& nm : matched_names) {
            auto sps = table1_.subpaths(nm);
            for (const auto& sp : sps) {
                target_path_ids.emplace_back(sp.path_id, sp.length);
            }
        }
    }
    if (target_path_ids.empty()) {
        out.status = "unknown_path";
        return out;
    }

    // Reset the diagnostics accumulators so they reflect only this query
    // (all subpath attempts): the find_sequences_for_tag LF cost and the
    // target-path walk LF cost.
    g_find_seq_stats = FindSeqStats{};
    panindexer::g_anchor_walk_stats = panindexer::AnchorWalkStats{};

    // ── Restrict to the subpaths the read actually touches ──────────────────
    // The previous code ran the full per-path build against EVERY resolved
    // subpath. For a chromosome target that's thousands of subpaths, and each
    // runs a detect_strand scan (~2 × n_source_mappings decompressSA) before
    // discovering the read isn't on it — an O(subpaths × mappings) decompressSA
    // blow-up (measured: >1M decompressSA calls and ~300 s for one read).
    //
    // A subpath the read does not visit produces zero anchors and is never the
    // "best" result, so skipping it CANNOT change the output. We find the
    // touched subpaths in a single pass: decompressSA each read node once and
    // keep the resolved subpaths whose GBWT path id (seqId/2, either
    // orientation) appears. That is O(n_source_mappings) decompressSA total,
    // independent of how many subpaths the name resolved to.
    std::unordered_set<size_t> target_pid_set;
    std::unordered_map<size_t, size_t> len_by_pid;
    target_pid_set.reserve(target_path_ids.size() * 2);
    len_by_pid.reserve(target_path_ids.size() * 2);
    for (const auto& [pid, plen] : target_path_ids) {
        target_pid_set.insert(pid);
        len_by_pid[pid] = plen;
    }

    std::vector<std::pair<size_t, size_t>> touched;  // (path_id, length), first-seen order
    std::unordered_set<size_t> touched_seen;
    for (const panindexer::SourceMapping& sm : source_mappings) {
        gbwt::node_type node = gbwt::Node::encode(
            static_cast<gbwt::node_type>(sm.node_id), sm.is_reverse);
        const auto _ds_t0 = std::chrono::high_resolution_clock::now();
        std::vector<gbwt::size_type> sa = gbwt_rindex_->decompressSA(node);
        panindexer::g_anchor_walk_stats.decompress_sa_ms +=
            std::chrono::duration<double, std::milli>(
                std::chrono::high_resolution_clock::now() - _ds_t0).count();
        panindexer::g_anchor_walk_stats.decompress_sa_calls++;
        panindexer::g_anchor_walk_stats.decompress_sa_entries += sa.size();
        for (gbwt::size_type v : sa) {
            size_t pid = static_cast<size_t>(gbwt_rindex_->seqId(v)) / 2;
            if (target_pid_set.count(pid) && touched_seen.insert(pid).second) {
                touched.emplace_back(pid, len_by_pid[pid]);
            }
        }
    }

    // Build anchors only for the touched subpaths, keeping the same "pick the
    // result with the most anchors" semantics as before.
    std::vector<panindexer::AnchorBuildResult> results;
    results.reserve(touched.size());
    for (const auto& [pid, plen] : touched) {
        results.push_back(panindexer::build_surject_anchors_for_path(
            *gbz_, rindex, sampled, *gbwt_rindex_, source_mappings, pid, plen));
    }

    if (results.empty()) {
        // Name resolved, but the read shares no node with any of its subpaths.
        out.status = "no_common_nodes";
    } else {
        // Pick the subpath the read overlaps most. We count DISTINCT covered
        // source mappings (read nodes), not raw anchors: the multiple-candidate
        // builder emits one anchor per target occurrence of a node, so a
        // repeat-heavy subpath would otherwise be over-credited by raw anchor
        // count and could beat the subpath the read truly aligns to.
        auto covered_source_mappings = [](const panindexer::AnchorBuildResult& r) {
            std::unordered_set<size_t> covered;
            for (const auto& a : r.anchors) {
                for (size_t i = a.source_mapping_begin; i < a.source_mapping_end; ++i) {
                    covered.insert(i);
                }
            }
            return covered.size();
        };
        const panindexer::AnchorBuildResult* best = &results.front();
        size_t best_cov = covered_source_mappings(*best);
        for (const auto& r : results) {
            size_t cov = covered_source_mappings(r);
            if (cov > best_cov) {
                best = &r;
                best_cov = cov;
            }
        }
        out.status = status_token(best->status);
        out.target_path_length = best->target_path_length;
        out.target_rev_strand = best->target_rev_strand;
        out.anchors.reserve(best->anchors.size());
        for (const auto& a : best->anchors) {
            out.anchors.push_back(to_anchor_record(a));
        }
    }

    // Report the find_sequences_for_tag LF cost accumulated over this query.
    out.find_seq_calls     = g_find_seq_stats.calls;
    out.find_seq_runs      = g_find_seq_stats.runs;
    out.find_seq_lf_steps  = g_find_seq_stats.lf_steps;
    out.find_seq_visits    = g_find_seq_stats.visits;
    out.last_run_nav_steps = g_find_seq_stats.last_run_nav_steps;
    out.last_run_length    = g_find_seq_stats.last_run_length;
    // Report the target-path walk LF cost.
    out.walk_lf_steps     = panindexer::g_anchor_walk_stats.walk_lf_steps;
    out.walk_span         = panindexer::g_anchor_walk_stats.walk_span;
    out.first_anchor_base = panindexer::g_anchor_walk_stats.first_anchor_base;
    out.last_anchor_base  = panindexer::g_anchor_walk_stats.last_anchor_base;
    // Wall-clock attribution of the build time.
    out.find_seq_ms            = g_find_seq_stats.total_ms;
    out.decompress_sa_ms       = panindexer::g_anchor_walk_stats.decompress_sa_ms;
    out.walk_ms                = panindexer::g_anchor_walk_stats.walk_ms;
    out.decompress_sa_calls    = panindexer::g_anchor_walk_stats.decompress_sa_calls;
    out.decompress_sa_entries  = panindexer::g_anchor_walk_stats.decompress_sa_entries;
    // Confirm the call explosion = subpath loop × per-subpath strand scan.
    out.n_target_subpaths = target_path_ids.size();
    out.n_source_mappings = source_mappings.size();
    return out;
}

std::vector<std::string> Index::get_haplotype_names() const {
    if (!loaded_)
        throw std::runtime_error("Index::get_haplotype_names called before load()");

    // Built once during load(); recomputing means walking every GBWT path.
    if (!haplotype_names_.empty()) return haplotype_names_;

    const gbwt::GBWT& gbwt_index = gbz_->index;
    if (!gbwt_index.hasMetadata())
        return {};

    const gbwt::Metadata& metadata = gbwt_index.metadata;
    std::unordered_set<std::string> seen;
    std::vector<std::string> result;

    for (size_t path_id = 0; path_id < metadata.paths(); ++path_id) {
        gbwt::PathName pn = metadata.path(path_id);
        std::string sample_name = metadata.sample(pn.sample);
        std::string hap_name = sample_name + "#" + std::to_string(pn.phase);
        if (seen.insert(hap_name).second) {
            result.push_back(hap_name);
        }
    }

    std::sort(result.begin(), result.end());
    return result;
}
