#include "pangenome_index/surject_anchor_builder.hpp"

#include <gbwt/utils.h>
#include <handlegraph/util.hpp>

#include <algorithm>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// ── File-scope (NOT inside namespace panindexer) so the extern matches
//    coordinate_translation.cpp's definition, which is also at file scope.
//    Mirrors the pattern used in pangenome_server.cpp:78-106.
struct NodeVisit {
    size_t seq_id;
    size_t offset;
    size_t bwt_pos;
    size_t packed_pos;
    uint64_t tag_code;
};

extern std::vector<NodeVisit> find_sequences_for_tag(
    panindexer::FastLocate& r_index,
    panindexer::SampledTagArray& sampled,
    uint64_t tag_code);

namespace panindexer {

namespace {

// ── Helpers ───────────────────────────────────────────────────────────────

/// Build a per-node lookup from source mappings: node → ordered list of
/// source-mapping indices that visit that node (in read order).
std::unordered_map<gbwt::node_type, std::vector<size_t>>
build_source_visits(const std::vector<SourceMapping>& source_mappings) {
    std::unordered_map<gbwt::node_type, std::vector<size_t>> out;
    out.reserve(source_mappings.size() * 2);
    for (size_t i = 0; i < source_mappings.size(); ++i) {
        const SourceMapping& m = source_mappings[i];
        gbwt::node_type node = gbwt::Node::encode(m.node_id, m.is_reverse);
        out[node].push_back(i);
    }
    return out;
}

/// Walk the target's GBWT path once, summing node lengths.
size_t compute_target_path_length(const gbwtgraph::GBZ& gbz,
                                  size_t target_gbwt_path_id) {
    gbwt::vector_type target_nodes = gbz.index.extract(
        gbwt::Path::encode(target_gbwt_path_id, false));
    size_t total = 0;
    for (gbwt::node_type node : target_nodes) {
        if (node == gbwt::ENDMARKER) break;
        total += gbz.graph.get_length(
            gbz.graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node)));
    }
    return total;
}

/// Locate a single (gbwt::edge_type, target_base_offset) anchor for a given
/// source mapping by pairing decompressSA's target visits with
/// find_sequences_for_tag's target base offsets.
///
/// `prefer_last` selects which visit to return when the node is visited
/// multiple times by the target:
///   - false: smallest base offset (largest GBWT seqOffset) — used for first anchor.
///   - true:  largest base offset (smallest GBWT seqOffset) — used for last anchor.
///
/// Returns true if at least one target visit exists; output params are then
/// populated. False otherwise.
bool locate_target_visit(
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const SourceMapping& mapping,
    size_t target_seq_id_fwd,
    bool prefer_last,
    gbwt::edge_type& out_edge,
    size_t& out_target_base_offset)
{
    gbwt::node_type node = gbwt::Node::encode(mapping.node_id, mapping.is_reverse);
    uint64_t tag_code = SampledTagArray::encode_value(mapping.node_id, mapping.is_reverse);

    // 1) RLBWT: all haplotype visits to this node with base offsets.
    std::vector<NodeVisit> rlbwt_visits = find_sequences_for_tag(rlbwt_rindex, sampled, tag_code);
    std::vector<size_t> target_base_offsets;
    target_base_offsets.reserve(rlbwt_visits.size());
    for (const NodeVisit& v : rlbwt_visits) {
        if (v.seq_id == target_seq_id_fwd) {
            target_base_offsets.push_back(v.offset);
        }
    }
    if (target_base_offsets.empty()) return false;
    std::sort(target_base_offsets.begin(), target_base_offsets.end());

    // 2) GBWT: all visits to this node from FastLocate, with their per-visit
    //    GBWT seqOffsets. The INDEX into decompressSA's output is the
    //    offset_in_record component of gbwt::edge_type (see
    //    coordinate_translation.cpp:663-675).
    std::vector<gbwt::size_type> sa_values = gbwt_fast_locate.decompressSA(node);
    // (seqOffset, sa_index_in_decompressSA) for target visits only.
    std::vector<std::pair<size_t, size_t>> target_gbwt_visits;
    target_gbwt_visits.reserve(sa_values.size());
    for (size_t i = 0; i < sa_values.size(); ++i) {
        if (gbwt_fast_locate.seqId(sa_values[i]) == target_seq_id_fwd) {
            target_gbwt_visits.emplace_back(
                static_cast<size_t>(gbwt_fast_locate.seqOffset(sa_values[i])),
                i);
        }
    }
    if (target_gbwt_visits.empty()) {
        // Inconsistent: RLBWT says there's a visit but GBWT doesn't agree.
        // This shouldn't happen if both indexes are built from the same GBZ.
        return false;
    }
    // Sort by seqOffset DESCENDING. Larger GBWT seqOffset == earlier in path
    // (see coordinate_translation.cpp:1406 comment), so the i-th entry here
    // pairs with the i-th smallest base offset in target_base_offsets.
    std::sort(target_gbwt_visits.begin(), target_gbwt_visits.end(),
              [](const std::pair<size_t, size_t>& a,
                 const std::pair<size_t, size_t>& b) {
                  return a.first > b.first;
              });

    // 3) Pair the chosen base offset with the matching GBWT visit.
    const size_t visit_count = std::min(target_base_offsets.size(),
                                        target_gbwt_visits.size());
    if (visit_count == 0) return false;
    const size_t pick = prefer_last ? (visit_count - 1) : 0;
    out_target_base_offset = target_base_offsets[pick];
    size_t offset_in_record = target_gbwt_visits[pick].second;
    out_edge = gbwt::edge_type(node, offset_in_record);
    return true;
}

/// Construct a gbwtgraph::GBWTGraph step_handle from a gbwt::edge_type using
/// the documented encoding (gbwtgraph.cpp:911-1024).
handlegraph::step_handle_t edge_to_step_handle(const gbwt::edge_type& edge) {
    handlegraph::step_handle_t step;
    handlegraph::as_integers(step)[0] = edge.first;
    handlegraph::as_integers(step)[1] = edge.second;
    return step;
}

/// Information about a single matched step during the LF walk.
struct WalkMatch {
    size_t source_mapping_index = 0;
    gbwt::edge_type edge{gbwt::ENDMARKER, 0};
    size_t target_base_offset = 0;
};

/// Walk the target GBWT path forward from `start_edge` (= first anchor),
/// accumulating base offsets. At each step, if the node is in source_visits,
/// pair with the next unmatched source mapping and emit a WalkMatch.
///
/// Stops when either:
///   - cursor edge equals `end_edge` (= last anchor) AND we've processed it, OR
///   - cursor base offset > `end_base_offset` (safety bound), OR
///   - LF returns ENDMARKER (path ended).
std::vector<WalkMatch> walk_target_collecting_matches(
    const gbwtgraph::GBZ& gbz,
    gbwt::edge_type start_edge,
    size_t start_base_offset,
    gbwt::edge_type end_edge,
    size_t end_base_offset,
    const std::unordered_map<gbwt::node_type, std::vector<size_t>>& source_visits)
{
    std::vector<WalkMatch> matches;

    // Per-node consumption queue: which source-mapping indices for this node
    // haven't been paired yet. Initialised from source_visits on first encounter.
    std::unordered_map<gbwt::node_type, std::deque<size_t>> pending;
    pending.reserve(source_visits.size());

    auto try_match = [&](const gbwt::edge_type& cursor, size_t base_off) {
        auto it = source_visits.find(cursor.first);
        if (it == source_visits.end()) return;
        auto& queue = pending.try_emplace(cursor.first,
                                          std::deque<size_t>(it->second.begin(),
                                                             it->second.end()))
                              .first->second;
        if (queue.empty()) return;
        size_t src_idx = queue.front();
        queue.pop_front();
        WalkMatch m;
        m.source_mapping_index = src_idx;
        m.edge = cursor;
        m.target_base_offset = base_off;
        matches.push_back(m);
    };

    gbwt::edge_type cursor = start_edge;
    size_t cursor_base = start_base_offset;

    // Process the starting step.
    try_match(cursor, cursor_base);

    // If start == end (one-step anchor range), we're done.
    if (cursor == end_edge) return matches;

    // Walk forward.
    while (true) {
        // Advance.
        handlegraph::handle_t handle = gbz.graph.get_handle(
            gbwt::Node::id(cursor.first), gbwt::Node::is_reverse(cursor.first));
        size_t node_len = gbz.graph.get_length(handle);
        gbwt::edge_type next = gbz.index.LF(cursor);
        if (next.first == gbwt::ENDMARKER) break;
        cursor = next;
        cursor_base += node_len;

        try_match(cursor, cursor_base);

        // Stop conditions.
        if (cursor == end_edge) break;
        // Safety net: don't run away past the last anchor's expected position.
        // Allow a small slack in case end_base_offset was the last anchor's
        // first base and we walk one more node before catching up.
        if (cursor_base > end_base_offset) break;
    }

    return matches;
}

/// Group `matches` (in target path order) into PrecomputedAnchor chunks.
/// New chunk starts when source_mapping_index isn't the previous + 1 (gap
/// in the source's read coverage) or when target steps aren't adjacent
/// (gap in target path coverage).
std::vector<PrecomputedAnchor> group_matches_into_chunks(
    const std::vector<WalkMatch>& matches,
    const std::vector<SourceMapping>& source_mappings)
{
    std::vector<PrecomputedAnchor> anchors;
    if (matches.empty()) return anchors;

    auto begin_chunk = [&](const WalkMatch& m) {
        PrecomputedAnchor a;
        a.source_mapping_begin = m.source_mapping_index;
        a.source_mapping_end = m.source_mapping_index + 1;
        a.read_begin_offset = source_mappings[m.source_mapping_index].read_begin_offset;
        a.read_end_offset   = source_mappings[m.source_mapping_index].read_end_offset;
        a.path_offset_step_begin = m.target_base_offset;
        a.path_offset_step_end   = m.target_base_offset;
        a.gbwt_edge_begin = m.edge;
        a.gbwt_edge_end   = m.edge;
        a.step_begin = edge_to_step_handle(m.edge);
        a.step_end   = edge_to_step_handle(m.edge);
        anchors.push_back(a);
    };
    auto extend_chunk = [&](const WalkMatch& m) {
        PrecomputedAnchor& a = anchors.back();
        a.source_mapping_end = m.source_mapping_index + 1;
        a.read_end_offset = source_mappings[m.source_mapping_index].read_end_offset;
        a.path_offset_step_end = m.target_base_offset;
        a.gbwt_edge_end = m.edge;
        a.step_end = edge_to_step_handle(m.edge);
    };

    begin_chunk(matches.front());
    for (size_t k = 1; k < matches.size(); ++k) {
        const WalkMatch& prev = matches[k - 1];
        const WalkMatch& curr = matches[k];
        // Same chunk iff source mappings are read-adjacent AND target moved
        // forward (not staying / not jumping). Target adjacency in path-step
        // terms is implicit — the walk only advances by one step per LF call,
        // and we record matches in walk order.
        bool source_adjacent = (curr.source_mapping_index == prev.source_mapping_index + 1);
        // Target base must strictly advance — walking via LF guarantees this
        // because each step has positive node length.
        if (source_adjacent) {
            extend_chunk(curr);
        } else {
            begin_chunk(curr);
        }
    }
    return anchors;
}

// path_names_for_haplotype is declared at file scope in pangenome_server.cpp
// for the same TranslationTable1 lookup. Reimplemented here so this file
// doesn't depend on pangenome_server.cpp.
std::vector<std::string> path_names_for_haplotype_local(
    const TranslationTable1& t1, const std::string& haplotype_prefix)
{
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

} // namespace

// ── Public entry points ────────────────────────────────────────────────────

AnchorBuildResult build_surject_anchors_for_path(
    const gbwtgraph::GBZ& gbz,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const std::vector<SourceMapping>& source_mappings,
    size_t target_gbwt_path_id)
{
    AnchorBuildResult result;

    // Validate target path.
    if (target_gbwt_path_id >= gbz.index.sequences() / 2) {
        result.status = AnchorBuildResult::Status::UnknownPath;
        return result;
    }
    // Translate to a libhandlegraph path_handle for the caller's convenience.
    result.target_path_handle = gbz.graph.path_to_handle(target_gbwt_path_id);

    if (source_mappings.empty()) {
        result.status = AnchorBuildResult::Status::EmptyAlignment;
        return result;
    }

    const size_t target_seq_id_fwd = 2 * target_gbwt_path_id;

    // Source-side hashmap for the walk lookup.
    auto source_visits = build_source_visits(source_mappings);

    // Target path length — used by AnchorBackedPositionGraph::get_path_length()
    // and as a safety guard for the LF walk.
    result.target_path_length = compute_target_path_length(gbz, target_gbwt_path_id);

    // 1. Find the FIRST anchor: walk source mappings forward, take the first
    //    one whose node is visited by the target.
    gbwt::edge_type first_edge{gbwt::ENDMARKER, 0};
    size_t first_target_base = 0;
    bool first_found = false;
    for (size_t i = 0; i < source_mappings.size(); ++i) {
        if (locate_target_visit(rlbwt_rindex, sampled, gbwt_fast_locate,
                                source_mappings[i], target_seq_id_fwd,
                                /*prefer_last=*/false,
                                first_edge, first_target_base)) {
            first_found = true;
            break;
        }
    }
    if (!first_found) {
        result.status = AnchorBuildResult::Status::NoCommonNodes;
        return result;
    }

    // 2. Find the LAST anchor: walk source mappings backward, take the first
    //    one (from the end) with a target match. Picks the largest target base
    //    offset for that mapping (prefer_last=true) so the walk terminates at
    //    the right path position even when the target visits the node
    //    multiple times.
    gbwt::edge_type last_edge{gbwt::ENDMARKER, 0};
    size_t last_target_base = 0;
    bool last_found = false;
    for (size_t i = source_mappings.size(); i-- > 0; ) {
        if (locate_target_visit(rlbwt_rindex, sampled, gbwt_fast_locate,
                                source_mappings[i], target_seq_id_fwd,
                                /*prefer_last=*/true,
                                last_edge, last_target_base)) {
            last_found = true;
            break;
        }
    }
    if (!last_found) {
        // Should not happen if first_found succeeded, but handle defensively.
        result.status = AnchorBuildResult::Status::NoCommonNodes;
        return result;
    }

    // Edge case: if first and last anchors are the same mapping, the walk is
    // a single step.
    if (first_target_base > last_target_base) {
        std::swap(first_target_base, last_target_base);
        std::swap(first_edge, last_edge);
    }

    // 3. Walk target forward from first to last, collecting matches.
    std::vector<WalkMatch> matches = walk_target_collecting_matches(
        gbz, first_edge, first_target_base, last_edge, last_target_base,
        source_visits);

    if (matches.empty()) {
        // Defensive: at least the first anchor's match should have been emitted.
        result.status = AnchorBuildResult::Status::NoCommonNodes;
        return result;
    }

    // 4. Group into chunks.
    result.anchors = group_matches_into_chunks(matches, source_mappings);

    // The walk emits matches in target path order. For colinear alignments
    // that equals read order, but when the target visits source nodes in a
    // different order than the read (e.g. path y in test0_target_loop loops
    // back through node 1, visiting node 4 before node 2), the chunks end
    // up in target order. Surject's downstream pipeline expects chunks in
    // read order, so sort by read_begin_offset here.
    std::sort(result.anchors.begin(), result.anchors.end(),
              [](const PrecomputedAnchor& a, const PrecomputedAnchor& b) {
                  return a.read_begin_offset < b.read_begin_offset;
              });

    result.status = AnchorBuildResult::Status::Ok;
    return result;
}

std::vector<AnchorBuildResult> build_surject_anchors(
    const gbwtgraph::GBZ& gbz,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const TranslationTable1& table1,
    const std::vector<SourceMapping>& source_mappings,
    const std::string& target_haplotype_name)
{
    std::vector<AnchorBuildResult> results;

    std::vector<std::string> subpath_names =
        path_names_for_haplotype_local(table1, target_haplotype_name);
    if (subpath_names.empty()) {
        AnchorBuildResult r;
        r.status = AnchorBuildResult::Status::UnknownPath;
        results.push_back(r);
        return results;
    }

    // Each subpath name in T1 corresponds to a set of GBWT path ids
    // (one SubpathInfo per fragment). Build anchors per fragment.
    for (const std::string& name : subpath_names) {
        std::vector<SubpathInfo> subpaths = table1.subpaths(name);
        for (const SubpathInfo& sp : subpaths) {
            results.push_back(build_surject_anchors_for_path(
                gbz, rlbwt_rindex, sampled, gbwt_fast_locate,
                source_mappings, sp.path_id));
        }
    }
    return results;
}

} // namespace panindexer
