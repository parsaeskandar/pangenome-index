#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include "pangenome_index/translation_tables.hpp"
#include <sdsl/wavelet_trees.hpp>
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
    mlockall(MCL_CURRENT | MCL_FUTURE);

    // 1. RLBWT r-index
    {
        std::ifstream rin(ri_path, std::ios::binary);
        if (!rin)
            throw std::runtime_error("Cannot open RLBWT r-index: " + ri_path);
        rlbwt_rindex_.load_encoded(rin);
        rlbwt_rindex_.ensure_last_rank();
        rlbwt_rindex_.ensure_last_select();
    }

    // 2. GBZ (GBWT + GBWTGraph)
    gbz_ = std::make_unique<gbwtgraph::GBZ>();
    sdsl::simple_sds::load_from(*gbz_, gbz_path);

    // 3. GBWT FastLocate
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
    }

    // 4. Sampled tag array
    {
        std::ifstream sin(tags_path, std::ios::binary);
        if (!sin)
            throw std::runtime_error("Cannot open sampled tags: " + tags_path);
        sampled_.load(sin);
        sampled_.ensure_run_rank();
        sampled_.ensure_run_select();
    }

    // 5. Translation tables
    {
        std::ifstream t1in(table1_path, std::ios::binary);
        if (!t1in)
            throw std::runtime_error("Cannot open Table 1: " + table1_path);
        table1_.load(t1in);
    }
    {
        std::ifstream t2in(table2_path, std::ios::binary);
        if (!t2in)
            throw std::runtime_error("Cannot open Table 2: " + table2_path);
        table2_.load(t2in);
    }

    path_to_global_ = build_path_id_to_global(table1_);
    loaded_ = true;
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
            table2_.lookup(src_path_id, tgt_haplotype, local_start, local_end);
        if (tgt_intervals.empty()) continue;

        size_t src_seq_id = 2 * src_path_id;
        auto it_src = path_to_global_.find(src_path_id);
        if (it_src == path_to_global_.end()) continue;
        size_t src_subpath_start = it_src->second.second;

        std::unordered_set<size_t> distinct_tgt_paths;
        for (const TargetInterval& ti : tgt_intervals)
            distinct_tgt_paths.insert(ti.tgt_path_id);

        std::vector<IntervalMapping> segs =
            table2_.segments(src_path_id, tgt_haplotype);

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
        ti.haplotype = tgt_haplotype;
        ti.start     = static_cast<int64_t>(ht.source_haplotype_offset);
        ti.end       = static_cast<int64_t>(ht.target_haplotype_offset);
        ti.strand    = '+';
        results.push_back(ti);
    }
    return results;
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
