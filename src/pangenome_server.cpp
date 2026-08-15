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
    log_step(("[5/5] Table 2: " + table2_path).c_str());
    {
        std::ifstream t2in(table2_path, std::ios::binary);
        if (!t2in)
            throw std::runtime_error("Cannot open Table 2: " + table2_path);
        table2_.load(t2in);
    }
    log_step("[5/5]   T2 loaded");

    path_to_global_ = build_path_id_to_global(table1_);
    loaded_ = true;
    log_step("[done] all indexes loaded");
}

std::vector<TranslatedInterval>
Index::translate(const std::string& src_haplotype,
                 int64_t start, int64_t end,
                 const std::string& tgt_haplotype) const {
    if (!loaded_)
        throw std::runtime_error("Index::translate called before load()");
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
            table1_.lookup(name, global_start, global_end + 1);
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
    return results;
}

std::vector<std::string>
Index::translatable_haplotypes(const std::string& src_haplotype,
                               int64_t start, int64_t end) const {
    if (!loaded_)
        throw std::runtime_error("Index::translatable_haplotypes called before load()");
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
        std::vector<PathInterval> pis = table1_.lookup(name, s, e + 1);
        for (PathInterval& pi : pis)
            source_intervals.push_back(pi);
    }
    if (source_intervals.empty())
        return {};

    // Candidate target haplotypes per source path come straight from T2's keys.
    // T2 is sparse — it only stores (src_path_id, tgt_haplotype) pairs that share
    // at least one graph node — so this is the homology-pruned candidate set, not
    // every haplotype in the graph.
    std::unordered_map<size_t, std::vector<std::string>> tgts_by_src;
    for (const std::pair<size_t, std::string>& key : table2_.keys())
        tgts_by_src[key.first].push_back(key.second);

    // Confirm each candidate actually overlaps THIS interval (a binary-searched
    // segment lookup — still no coordinate trace).
    std::unordered_set<std::string> found;
    for (const PathInterval& pi : source_intervals) {
        if (pi.end <= pi.start) continue;
        auto it = tgts_by_src.find(pi.path_id);
        if (it == tgts_by_src.end()) continue;
        for (const std::string& tgt : it->second) {
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

    // Find cg:Z: CIGAR among optional fields.
    std::string cg_str;
    for (size_t k = 12; k < fields.size(); ++k) {
        if (fields[k].size() > 5 && fields[k].compare(0, 5, "cg:Z:") == 0) {
            cg_str = fields[k].substr(5);
            break;
        }
    }

    std::vector<panindexer::SourceMapping> result;
    result.reserve(path_nodes.size());

    if (!cg_str.empty()) {
        // ── CIGAR-driven walk ─────────────────────────────────────────────
        auto cigar_ops = parse_cg_cigar(cg_str);

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
    std::unordered_set<std::string> here;
    uint64_t total_bp = 0;

    for (const panindexer::SourceMapping& m : mappings) {
        const uint64_t bp = (m.read_end_offset > m.read_begin_offset)
                            ? (m.read_end_offset - m.read_begin_offset) : 0;
        if (bp == 0) continue;   // pure insertion: no graph node to credit
        total_bp += bp;

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
        for (const std::string& h : here) covered[h] += bp;
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
        out.push_back(std::move(hc));
    }
    std::sort(out.begin(), out.end(),
              [](const HaplotypeCoverage& a, const HaplotypeCoverage& b) {
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
