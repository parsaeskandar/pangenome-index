#ifndef PANGENOME_INDEX_SURJECT_ANCHOR_BUILDER_HPP
#define PANGENOME_INDEX_SURJECT_ANCHOR_BUILDER_HPP

/**
 * surject_anchor_builder.hpp
 *
 * Build pre-computed surjection anchors between a graph alignment (Giraffe
 * output) and a target haplotype path.
 *
 * The output anchors carry the data Surjector's downstream pipeline needs
 * (step_handles + base positions on the target path) without requiring a
 * built bdsg::ReferencePathOverlay over the target. The consumer wraps
 * these in an AnchorBackedPositionGraph and hands that to Surjector.
 *
 * Algorithm (Shape A — mirrors trace_coordinates_gbwt's pattern):
 *   1. Find one anchor by scanning source mappings forward → first one
 *      whose graph node is also visited by the target haplotype.
 *   2. Find another anchor by scanning source mappings backward → last one
 *      with a target match.
 *   3. Walk the target's GBWT path from first to last anchor via
 *      gbwt::GBWT::LF, accumulating base offsets. At each step, check
 *      whether the node is in the source's mapping set.
 *   4. Group matches into chunks (contiguous on both read and target path).
 *
 * Cost per query: 2× find_sequences_for_tag + a few decompressSA + O(walk
 * distance) GBWT LF calls. No full path enumeration.
 */

#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include "pangenome_index/translation_tables.hpp"

#include <gbwt/fast_locate.h>
#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>
#include <handlegraph/types.hpp>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace panindexer {

/// Minimal per-mapping payload the anchor builder needs from a Giraffe-style
/// graph alignment. The caller (engine code that has vg::Alignment) converts
/// each Mapping into one of these so the liftover library doesn't have to
/// depend on vg::Alignment / vg.pb.h.
struct SourceMapping {
    /// Graph node id (1-based, as in GBWT/gbwtgraph).
    int64_t node_id = 0;
    /// True if the source visits the node on its reverse strand.
    bool is_reverse = false;
    /// Cumulative read offset where this mapping begins (0-based, on the
    /// forward strand of the read sequence).
    size_t read_begin_offset = 0;
    /// Cumulative read offset where this mapping ends (exclusive).
    size_t read_end_offset = 0;
    /// Offset within the node where this mapping starts.
    size_t node_offset_in_node = 0;
    /// Number of bases the mapping consumes on the graph side (= path_from_length
    /// of the mapping's edits combined).
    size_t mapping_from_length = 0;
};

/// One anchor: a contiguous run of source mappings whose graph nodes also
/// appear as consecutive steps on the target haplotype path.
///
/// Fields mirror Giraffe_server/src/anchor_backed_position_graph.hpp's
/// PrecomputedAnchor (which is the consumer of these records). Keeping the
/// two definitions in sync is a manual coordination today; the wire protocol
/// will eventually serialize this same payload.
struct PrecomputedAnchor {
    /// First step on the target path for this anchor.
    handlegraph::step_handle_t step_begin;
    /// Last step (inclusive) on the target path for this anchor.
    handlegraph::step_handle_t step_end;
    /// Base position of step_begin on the target path.
    size_t path_offset_step_begin = 0;
    /// Base position of step_end on the target path.
    size_t path_offset_step_end = 0;

    /// First source-mapping index covered by this anchor.
    size_t source_mapping_begin = 0;
    /// One past the last source-mapping index covered by this anchor.
    size_t source_mapping_end = 0;
    /// Read range covered (from source mappings).
    size_t read_begin_offset = 0;
    size_t read_end_offset = 0;

    /// GBWT search states (= gbwtgraph step_handle_t encoded as edge_type
    /// per gbwtgraph.cpp:911-1024). Carried for the future surject_with_anchors
    /// API; AnchorBackedPositionGraph only reads the step_handles above.
    gbwt::edge_type gbwt_edge_begin{gbwt::ENDMARKER, 0};
    gbwt::edge_type gbwt_edge_end{gbwt::ENDMARKER, 0};
};

/// Result of build_surject_anchors_for_path.
struct AnchorBuildResult {
    enum class Status {
        Ok,                  ///< anchors populated; caller can proceed
        EmptyAlignment,      ///< source had no graph mappings
        UnknownPath,         ///< the requested target GBWT path id is out of range
        NoCommonNodes,       ///< source and target share no graph node
    };

    Status status = Status::EmptyAlignment;
    std::vector<PrecomputedAnchor> anchors;

    /// Total length in bases of the target path. Needed by
    /// AnchorBackedPositionGraph::get_path_length().
    size_t target_path_length = 0;

    /// libhandlegraph handle for the target path. Convenience for the caller.
    handlegraph::path_handle_t target_path_handle{};

    /// Forward-strand match assumed (v1). When reverse-strand surjection is
    /// added, this flips to true and the anchor positions are measured from
    /// the target path's forward strand regardless.
    bool target_rev_strand = false;
};

/**
 * Build surjection anchors for a single target subpath identified by its
 * GBWT path id (= forward seq_id / 2).
 *
 * Inputs:
 *   - `gbz`           : graph + GBWT (Giraffe_server's gbwtgraph::GBZ).
 *   - `rlbwt_rindex`  : RLBWT r-index used by find_sequences_for_tag.
 *   - `sampled`       : sampled tag array used by find_sequences_for_tag.
 *   - `gbwt_fast_locate`: GBWT FastLocate for decompressSA queries.
 *   - `source_mappings`: graph alignment from Giraffe, in read order.
 *   - `target_gbwt_path_id`: the GBWT path id of the target subpath.
 *   - `precomputed_target_path_length`: total base length of the target path,
 *     if the caller already knows it (e.g. from TranslationTable1's
 *     SubpathInfo.length, which is computed by the identical extract-and-sum
 *     loop at index-build time). Pass 0 to have this function compute it by
 *     extracting the full target path — O(path length), a whole chromosome for
 *     chromosome-scale targets. Passing the precomputed value avoids that
 *     extraction entirely (coordinate translation never extracts the full path).
 *
 * Outputs (via the returned struct):
 *   - `anchors`: in read order. May be empty if NoCommonNodes.
 *   - `target_path_length`, `target_path_handle`.
 *   - Status code distinguishing failure modes.
 *
 * Notes for v1:
 *   - Forward strand only. If the source visits the target's reverse strand,
 *     it appears as "no common nodes" (a different gbwt::node_type key).
 *   - Tandem repeats: where the target visits a node multiple times, the i-th
 *     matched source mapping (in read order) is paired with the i-th target
 *     visit (in target path order). Same convention as Surjector.
 */
AnchorBuildResult build_surject_anchors_for_path(
    const gbwtgraph::GBZ& gbz,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const std::vector<SourceMapping>& source_mappings,
    size_t target_gbwt_path_id,
    size_t precomputed_target_path_length = 0);

/**
 * Convenience overload: resolve a target haplotype name (e.g.
 * "GRCh38#0#chr1") via TranslationTable1, then build anchors for each
 * matching subpath. Returns one AnchorBuildResult per subpath, in T1 order.
 */
std::vector<AnchorBuildResult> build_surject_anchors(
    const gbwtgraph::GBZ& gbz,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const gbwt::FastLocate& gbwt_fast_locate,
    const TranslationTable1& table1,
    const std::vector<SourceMapping>& source_mappings,
    const std::string& target_haplotype_name);

} // namespace panindexer

#endif // PANGENOME_INDEX_SURJECT_ANCHOR_BUILDER_HPP
