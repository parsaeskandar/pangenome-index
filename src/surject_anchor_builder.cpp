#include "pangenome_index/surject_anchor_builder.hpp"

#include <gbwt/utils.h>
#include <handlegraph/util.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
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

// Diagnostics accumulator for the target-path LF walk (see header). Defined
// here (external linkage); read from pangenome_server.cpp.
thread_local AnchorWalkStats g_anchor_walk_stats;

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
///
/// PERFORMANCE — mirrors coordinate_translation.cpp's check_common_node():
/// the cheap GBWT decompressSA() probe runs FIRST and the function returns
/// early when the target does not visit this node, so the expensive RLBWT
/// find_sequences_for_tag() (which enumerates every haplotype's visits to the
/// node via locateNext) runs ONLY for nodes the target actually visits. When
/// callers scan source mappings until the first/last target-visited node, this
/// means find_sequences_for_tag() is invoked ~once per anchor (i.e. ~twice per
/// read) instead of once per candidate mapping. The original order
/// (find_sequences_for_tag first) made it the per-candidate cost and was the
/// whole-genome bottleneck. The returned values are identical to that order:
/// this is purely a reordering of two independent probes plus an early exit.
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

    // 1) GBWT decompressSA FIRST — the cheap gate (check_common_node:1399-1433).
    //    Collect this node's target visits as (seqOffset, sa_index_in_decompressSA).
    //    The INDEX into decompressSA's output is the offset_in_record component
    //    of gbwt::edge_type (see coordinate_translation.cpp:663-675).
    const auto _ds_t0 = std::chrono::high_resolution_clock::now();
    std::vector<gbwt::size_type> sa_values = gbwt_fast_locate.decompressSA(node);
    g_anchor_walk_stats.decompress_sa_ms +=
        std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - _ds_t0).count();
    g_anchor_walk_stats.decompress_sa_calls++;
    g_anchor_walk_stats.decompress_sa_entries += sa_values.size();
    std::vector<std::pair<size_t, size_t>> target_gbwt_visits;
    target_gbwt_visits.reserve(sa_values.size());
    for (size_t i = 0; i < sa_values.size(); ++i) {
        if (gbwt_fast_locate.seqId(sa_values[i]) == target_seq_id_fwd) {
            target_gbwt_visits.emplace_back(
                static_cast<size_t>(gbwt_fast_locate.seqOffset(sa_values[i])),
                i);
        }
    }
    // Target does not visit this node → cheap reject, WITHOUT touching the RLBWT.
    if (target_gbwt_visits.empty()) return false;
    // Sort by seqOffset DESCENDING. Larger GBWT seqOffset == earlier in path
    // (see coordinate_translation.cpp:1406 comment), so the i-th entry here
    // pairs with the i-th smallest base offset in target_base_offsets.
    std::sort(target_gbwt_visits.begin(), target_gbwt_visits.end(),
              [](const std::pair<size_t, size_t>& a,
                 const std::pair<size_t, size_t>& b) {
                  return a.first > b.first;
              });

    // 2) RLBWT find_sequences_for_tag ONLY now that the target is confirmed to
    //    visit this node — to read off the target base offset(s). This is the
    //    single expensive call check_common_node makes for a confirmed common
    //    node (coordinate_translation.cpp:1470).
    uint64_t tag_code = SampledTagArray::encode_value(mapping.node_id, mapping.is_reverse);
    std::vector<NodeVisit> rlbwt_visits = find_sequences_for_tag(rlbwt_rindex, sampled, tag_code);
    std::vector<size_t> target_base_offsets;
    target_base_offsets.reserve(rlbwt_visits.size());
    for (const NodeVisit& v : rlbwt_visits) {
        if (v.seq_id == target_seq_id_fwd) {
            target_base_offsets.push_back(v.offset);
        }
    }
    if (target_base_offsets.empty()) {
        // Inconsistent: GBWT says there's a target visit but RLBWT doesn't agree.
        // This shouldn't happen if both indexes are built from the same GBZ.
        return false;
    }
    std::sort(target_base_offsets.begin(), target_base_offsets.end());

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

/// Enumerate ALL of the target sequence's visits to `mapping`'s node, as
/// (target_base_offset, gbwt::edge_type) pairs sorted ascending by base offset.
///
/// This is the multiple-candidate sibling of locate_target_visit(): instead of
/// collapsing a repeated node to a single chosen occurrence (first/last), it
/// returns every occurrence so the caller can hand them all to the Surjector,
/// whose colinear chunk-chaining picks the occurrence the read actually came
/// from (see extract_overlapping_paths in surjector.cpp: it starts a distinct
/// ref chunk per step of a node). Same cheap-gate-first ordering as
/// locate_target_visit — the GBWT decompressSA probe runs FIRST and the RLBWT
/// find_sequences_for_tag() runs only for nodes the target actually visits.
///
/// Returns true (and fills out_visits) iff the target visits this node.
bool locate_all_target_visits(
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const SourceMapping& mapping,
    size_t target_seq_id_fwd,
    std::vector<std::pair<size_t, gbwt::edge_type>>& out_visits)
{
    out_visits.clear();
    gbwt::node_type node = gbwt::Node::encode(mapping.node_id, mapping.is_reverse);

    // 1) GBWT decompressSA gate FIRST (cheap). Collect this node's target visits
    //    as (seqOffset, index-in-decompressSA); the index is the offset_in_record
    //    component of gbwt::edge_type.
    const auto _ds_t0 = std::chrono::high_resolution_clock::now();
    std::vector<gbwt::size_type> sa_values = gbwt_fast_locate.decompressSA(node);
    g_anchor_walk_stats.decompress_sa_ms +=
        std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - _ds_t0).count();
    g_anchor_walk_stats.decompress_sa_calls++;
    g_anchor_walk_stats.decompress_sa_entries += sa_values.size();
    std::vector<std::pair<size_t, size_t>> target_gbwt_visits;
    target_gbwt_visits.reserve(sa_values.size());
    for (size_t i = 0; i < sa_values.size(); ++i) {
        if (gbwt_fast_locate.seqId(sa_values[i]) == target_seq_id_fwd) {
            target_gbwt_visits.emplace_back(
                static_cast<size_t>(gbwt_fast_locate.seqOffset(sa_values[i])),
                i);
        }
    }
    if (target_gbwt_visits.empty()) return false;  // off-target node: cheap reject.
    // Sort by seqOffset DESCENDING: larger GBWT seqOffset == earlier in path, so
    // the j-th entry here pairs with the j-th SMALLEST base offset below.
    std::sort(target_gbwt_visits.begin(), target_gbwt_visits.end(),
              [](const std::pair<size_t, size_t>& a,
                 const std::pair<size_t, size_t>& b) {
                  return a.first > b.first;
              });

    // 2) RLBWT find_sequences_for_tag only now that a target visit is confirmed.
    uint64_t tag_code = SampledTagArray::encode_value(mapping.node_id, mapping.is_reverse);
    std::vector<NodeVisit> rlbwt_visits = find_sequences_for_tag(rlbwt_rindex, sampled, tag_code);
    std::vector<size_t> target_base_offsets;
    target_base_offsets.reserve(rlbwt_visits.size());
    for (const NodeVisit& v : rlbwt_visits) {
        if (v.seq_id == target_seq_id_fwd) {
            target_base_offsets.push_back(v.offset);
        }
    }
    if (target_base_offsets.empty()) return false;  // index disagreement; skip.
    std::sort(target_base_offsets.begin(), target_base_offsets.end());

    // 3) Pair each base offset (ascending) with its matching GBWT visit
    //    (seqOffset descending). Emit EVERY pair as a candidate.
    const size_t visit_count = std::min(target_base_offsets.size(),
                                        target_gbwt_visits.size());
    out_visits.reserve(visit_count);
    for (size_t j = 0; j < visit_count; ++j) {
        gbwt::edge_type edge(node, target_gbwt_visits[j].second);
        out_visits.emplace_back(target_base_offsets[j], edge);
    }
    return !out_visits.empty();
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
        g_anchor_walk_stats.walk_lf_steps++;   // diagnostics: this is the real hot path
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

// ── Reverse-strand support ──────────────────────────────────────────────────

/// Cheap probe: does the target's forward sequence (target_seq_id_fwd) visit
/// `node_id` in the given orientation? Uses only the GBWT decompressSA (the
/// same cheap gate locate_target_visit uses); no RLBWT work.
bool target_visits_node(const gbwt::FastLocate& gbwt_fast_locate,
                        int64_t node_id, bool is_reverse,
                        size_t target_seq_id_fwd) {
    gbwt::node_type node = gbwt::Node::encode(node_id, is_reverse);
    const auto _ds_t0 = std::chrono::high_resolution_clock::now();
    std::vector<gbwt::size_type> sa_values = gbwt_fast_locate.decompressSA(node);
    g_anchor_walk_stats.decompress_sa_ms +=
        std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - _ds_t0).count();
    g_anchor_walk_stats.decompress_sa_calls++;
    g_anchor_walk_stats.decompress_sa_entries += sa_values.size();
    for (gbwt::size_type sa : sa_values) {
        if (gbwt_fast_locate.seqId(sa) == target_seq_id_fwd) return true;
    }
    return false;
}

/// Determine which strand of the read aligns to the target, from the first
/// source mapping the target visits in EITHER orientation:
///   - returns 0 (forward): target visits the node in the read's orientation,
///   - returns 1 (reverse): target visits the node in the FLIPPED orientation
///     (the read is the reverse complement of this target region),
///   - returns -1: no source node is on the target in either orientation.
/// Probing both orientations per candidate stops at the first shared node, so
/// a reverse-strand read is detected in O(1) probes instead of scanning all
/// mappings in the wrong orientation (which is what made revcomp reads slow).
/// For a colinear read sampled from the target this hits the very first mapping.
int detect_strand(const gbwt::FastLocate& gbwt_fast_locate,
                  const std::vector<SourceMapping>& mappings,
                  size_t target_seq_id_fwd) {
    for (const SourceMapping& m : mappings) {
        if (target_visits_node(gbwt_fast_locate, m.node_id, m.is_reverse, target_seq_id_fwd)) {
            return 0;  // forward strand
        }
        if (target_visits_node(gbwt_fast_locate, m.node_id, !m.is_reverse, target_seq_id_fwd)) {
            return 1;  // reverse strand
        }
    }
    return -1;  // no common node either strand
}

/// Build strand-normalized source mappings for the reverse-strand case.
///
/// A reverse-complement read traverses the target's nodes in the OPPOSITE
/// orientation and in reverse order (see the orientation analysis in the
/// header / docs). This rewrites the list as if the read had been on the
/// target's forward strand:
///   - flip each node's orientation, so encode(node_id, is_reverse) now matches
///     the orientation the target path traverses the node in,
///   - reverse the order (read order → target-forward order), and
///   - flip read offsets about the read length R (read_begin' = R - read_end,
///     read_end' = R - read_begin) so they increase in target-forward order.
/// Running find_anchors_for_strand() on this list yields anchors whose
/// step_handles and base offsets are on the target's FORWARD strand — exactly
/// what AnchorBackedPositionGraph needs. The Surjector then derives the '-'
/// strand itself from the (reverse-oriented) graph alignment; these anchors
/// only describe the target path's geometry, which is strand-independent.
std::vector<SourceMapping> make_reverse_strand_mappings(
    const std::vector<SourceMapping>& mappings) {
    size_t read_len = 0;
    for (const SourceMapping& m : mappings) {
        read_len = std::max(read_len, m.read_end_offset);
    }
    std::vector<SourceMapping> out;
    out.reserve(mappings.size());
    for (size_t k = mappings.size(); k-- > 0; ) {
        SourceMapping s = mappings[k];
        s.is_reverse = !s.is_reverse;
        const size_t rb = mappings[k].read_begin_offset;
        const size_t re = mappings[k].read_end_offset;
        s.read_begin_offset = (read_len >= re) ? (read_len - re) : 0;
        s.read_end_offset   = (read_len >= rb) ? (read_len - rb) : 0;
        out.push_back(s);
    }
    return out;
}

/// Multiple-candidate anchor search (DEFAULT). For each source mapping the
/// target visits, emit ONE singleton anchor per target visit — i.e. hand the
/// Surjector every occurrence of every read node on the target path and let its
/// own colinear chunk-chaining pick the occurrence the read came from.
///
/// This replaces the "pick one occurrence (first/last) of the boundary nodes,
/// then LF-walk the span between them" strategy. That strategy had two failure
/// modes on repeated nodes: (1) if the read came from a LATER occurrence of the
/// first boundary node, the chosen (earliest) occurrence sat chromosome-scale
/// away from the true locus, so the LF walk between anchors traversed the whole
/// intervening span (slow); and (2) the FIFO per-node pairing during that walk
/// collapsed each node to a single occurrence, so the wrong copy could be the
/// only one seeded. Enumerating every occurrence removes both: there is no walk
/// (each visit's base offset comes straight from the RLBWT), and the Surjector
/// sees all copies. The AnchorBackedPositionGraph records every anchor's step
/// in node_target_steps_ (no dedup), so for_each_step_on_handle returns them all
/// and Surjector::extract_overlapping_paths seeds a ref chunk per copy.
///
/// Returns anchors sorted by (read_begin_offset, target base offset); empty if
/// the target shares no node with the read.
std::vector<PrecomputedAnchor> find_anchors_all_candidates(
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const std::vector<SourceMapping>& mappings,
    size_t target_seq_id_fwd) {
    std::vector<PrecomputedAnchor> anchors;
    std::vector<std::pair<size_t, gbwt::edge_type>> visits;
    for (size_t i = 0; i < mappings.size(); ++i) {
        // off-target nodes (read insertions relative to the target) are rejected
        // by the cheap decompressSA gate inside locate_all_target_visits.
        if (!locate_all_target_visits(rlbwt_rindex, sampled, gbwt_fast_locate,
                                      mappings[i], target_seq_id_fwd, visits)) {
            continue;
        }
        for (const auto& v : visits) {
            const size_t base_off = v.first;
            const gbwt::edge_type& edge = v.second;
            PrecomputedAnchor a;
            a.source_mapping_begin = i;
            a.source_mapping_end   = i + 1;
            a.read_begin_offset = mappings[i].read_begin_offset;
            a.read_end_offset   = mappings[i].read_end_offset;
            a.path_offset_step_begin = base_off;
            a.path_offset_step_end   = base_off;
            a.gbwt_edge_begin = edge;
            a.gbwt_edge_end   = edge;
            a.step_begin = edge_to_step_handle(edge);
            a.step_end   = edge_to_step_handle(edge);
            anchors.push_back(a);
        }
    }
    // Read order (ties broken by target position) — deterministic and what the
    // tests / consumer expect. Multiple anchors MAY share a read range: they are
    // alternative placements of the same read node, which is the whole point.
    std::sort(anchors.begin(), anchors.end(),
              [](const PrecomputedAnchor& a, const PrecomputedAnchor& b) {
                  if (a.read_begin_offset != b.read_begin_offset)
                      return a.read_begin_offset < b.read_begin_offset;
                  return a.path_offset_step_begin < b.path_offset_step_begin;
              });
    return anchors;
}

/// Core anchor search for ONE strand. `mappings` are oriented so each node
/// matches the target path's forward-strand traversal (the read's own mappings
/// for a forward-strand read; make_reverse_strand_mappings()'s output for a
/// reverse-strand read). Returns anchors sorted by read_begin_offset, or an
/// empty vector if no common node with the target was found.
///
/// Dispatches between the default multiple-candidate search and the legacy
/// single-occurrence-plus-LF-walk search. The walk path is retained behind the
/// PANGENOME_SURJECT_ANCHOR_WALK env var so the two can be A/B compared on the
/// cluster (correctness + speed) without a rebuild; unset (the default) uses the
/// multiple-candidate search.
std::vector<PrecomputedAnchor> find_anchors_walk(
    const gbwtgraph::GBZ& gbz,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const std::vector<SourceMapping>& mappings,
    size_t target_seq_id_fwd);

std::vector<PrecomputedAnchor> find_anchors_for_strand(
    const gbwtgraph::GBZ& gbz,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const std::vector<SourceMapping>& mappings,
    size_t target_seq_id_fwd) {
    static const bool use_walk =
        (std::getenv("PANGENOME_SURJECT_ANCHOR_WALK") != nullptr);
    if (use_walk) {
        return find_anchors_walk(gbz, rlbwt_rindex, sampled, gbwt_fast_locate,
                                 mappings, target_seq_id_fwd);
    }
    return find_anchors_all_candidates(rlbwt_rindex, sampled, gbwt_fast_locate,
                                       mappings, target_seq_id_fwd);
}

/// Legacy single-occurrence + LF-walk search (see find_anchors_for_strand).
/// This is the former body of find_anchors_for_strand, unchanged.
std::vector<PrecomputedAnchor> find_anchors_walk(
    const gbwtgraph::GBZ& gbz,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const std::vector<SourceMapping>& mappings,
    size_t target_seq_id_fwd) {
    auto source_visits = build_source_visits(mappings);

    // 1. FIRST anchor: first source mapping (read order) the target visits.
    gbwt::edge_type first_edge{gbwt::ENDMARKER, 0};
    size_t first_target_base = 0;
    bool first_found = false;
    for (size_t i = 0; i < mappings.size(); ++i) {
        if (locate_target_visit(rlbwt_rindex, sampled, gbwt_fast_locate,
                                mappings[i], target_seq_id_fwd,
                                /*prefer_last=*/false, first_edge, first_target_base)) {
            first_found = true;
            break;
        }
    }
    if (!first_found) return {};

    // 2. LAST anchor: last source mapping (read order) the target visits.
    gbwt::edge_type last_edge{gbwt::ENDMARKER, 0};
    size_t last_target_base = 0;
    bool last_found = false;
    for (size_t i = mappings.size(); i-- > 0; ) {
        if (locate_target_visit(rlbwt_rindex, sampled, gbwt_fast_locate,
                                mappings[i], target_seq_id_fwd,
                                /*prefer_last=*/true, last_edge, last_target_base)) {
            last_found = true;
            break;
        }
    }
    if (!last_found) return {};

    // Edge case: keep first/last ordered by target base offset.
    if (first_target_base > last_target_base) {
        std::swap(first_target_base, last_target_base);
        std::swap(first_edge, last_edge);
    }

    // Diagnostics: how far apart the chosen anchors are on the target. The walk
    // below traverses this whole span via LF, so a large span (repeated boundary
    // nodes whose chosen occurrences are chromosome-scale apart) = a slow build.
    g_anchor_walk_stats.walks++;
    g_anchor_walk_stats.first_anchor_base = first_target_base;
    g_anchor_walk_stats.last_anchor_base = last_target_base;
    g_anchor_walk_stats.walk_span = last_target_base - first_target_base;

    // 3. Walk target forward from first to last, collecting matches.
    const auto _walk_t0 = std::chrono::high_resolution_clock::now();
    std::vector<WalkMatch> matches = walk_target_collecting_matches(
        gbz, first_edge, first_target_base, last_edge, last_target_base, source_visits);
    g_anchor_walk_stats.walk_ms +=
        std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - _walk_t0).count();
    if (matches.empty()) return {};

    // 4. Group into chunks; sort into read order (walk emits target order).
    std::vector<PrecomputedAnchor> anchors = group_matches_into_chunks(matches, mappings);
    std::sort(anchors.begin(), anchors.end(),
              [](const PrecomputedAnchor& a, const PrecomputedAnchor& b) {
                  return a.read_begin_offset < b.read_begin_offset;
              });
    return anchors;
}

} // namespace

// ── Public entry points ────────────────────────────────────────────────────

AnchorBuildResult build_surject_anchors_for_path(
    const gbwtgraph::GBZ& gbz,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const std::vector<SourceMapping>& source_mappings,
    size_t target_gbwt_path_id,
    size_t precomputed_target_path_length)
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

    // Target path length — used by AnchorBackedPositionGraph::get_path_length().
    // Prefer the caller's precomputed value (e.g. TranslationTable1's
    // SubpathInfo.length, computed by the identical extract-and-sum loop at
    // index-build time) to avoid extracting the whole target path here, which
    // is O(path length) — a full chromosome for chromosome-scale targets.
    // Coordinate translation never extracts the full path; this matches that.
    result.target_path_length = (precomputed_target_path_length > 0)
        ? precomputed_target_path_length
        : compute_target_path_length(gbz, target_gbwt_path_id);

    // Determine which strand of the read aligns to the target. A revcomp read
    // traverses the target's nodes in the opposite orientation, so it shares
    // no node with the target's forward sequence in the read's own orientation
    // — detect that up front (O(1) probes for a colinear read) instead of
    // scanning every mapping in the wrong orientation.
    int strand = detect_strand(gbwt_fast_locate, source_mappings, target_seq_id_fwd);
    if (strand < 0) {
        result.status = AnchorBuildResult::Status::NoCommonNodes;
        return result;
    }

    std::vector<PrecomputedAnchor> anchors;
    if (strand == 0) {
        // Forward strand: the read's own mappings. Byte-for-byte the old path.
        anchors = find_anchors_for_strand(gbz, rlbwt_rindex, sampled,
                                          gbwt_fast_locate, source_mappings,
                                          target_seq_id_fwd);
    } else {
        // Reverse strand: normalize to the target's forward orientation/order,
        // then run the identical search. The resulting step_handles and base
        // offsets are on the target's forward strand (what the position graph
        // needs); the Surjector emits the '-' strand from the graph alignment.
        std::vector<SourceMapping> rev_mappings =
            make_reverse_strand_mappings(source_mappings);
        anchors = find_anchors_for_strand(gbz, rlbwt_rindex, sampled,
                                          gbwt_fast_locate, rev_mappings,
                                          target_seq_id_fwd);
        result.target_rev_strand = true;
    }

    if (anchors.empty()) {
        result.status = AnchorBuildResult::Status::NoCommonNodes;
        return result;
    }

    result.anchors = std::move(anchors);
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
    // (one SubpathInfo per fragment). Build anchors per fragment, passing
    // sp.length so the builder doesn't re-extract the full target path.
    for (const std::string& name : subpath_names) {
        std::vector<SubpathInfo> subpaths = table1.subpaths(name);
        for (const SubpathInfo& sp : subpaths) {
            results.push_back(build_surject_anchors_for_path(
                gbz, rlbwt_rindex, sampled, gbwt_fast_locate,
                source_mappings, sp.path_id, sp.length));
        }
    }
    return results;
}

} // namespace panindexer
