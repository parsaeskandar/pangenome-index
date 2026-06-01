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

    std::vector<panindexer::AnchorBuildResult> results =
        panindexer::build_surject_anchors(
            *gbz_, rindex, sampled, *gbwt_rindex_,
            table1_, source_mappings, target_haplotype);

    if (results.empty()) {
        out.status = "unknown_path";
        return out;
    }

    // Pick the subpath result with the most anchors. If none have anchors,
    // fall back to the first (so status/target_path_length get reported).
    const panindexer::AnchorBuildResult* best = &results.front();
    for (const auto& r : results) {
        if (r.anchors.size() > best->anchors.size()) {
            best = &r;
        }
    }

    out.status = status_token(best->status);
    out.target_path_length = best->target_path_length;
    out.target_rev_strand = best->target_rev_strand;
    out.anchors.reserve(best->anchors.size());
    for (const auto& a : best->anchors) {
        out.anchors.push_back(to_anchor_record(a));
    }
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
