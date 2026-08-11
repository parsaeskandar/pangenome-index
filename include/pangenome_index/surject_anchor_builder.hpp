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
 * Algorithm (multiple-candidate — DEFAULT):
 *   For each source mapping the target haplotype visits, emit ONE anchor per
 *   target occurrence of that node (a "candidate"). Where the target visits a
 *   node several times (tandem repeat / circular contig), every occurrence is
 *   emitted. The consumer records all of them, so the Surjector's own colinear
 *   chunk-chaining picks the occurrence the read actually came from rather than
 *   the builder pre-guessing one. Each occurrence's target base offset comes
 *   straight from the RLBWT tag array, so there is NO LF walk between anchors.
 *
 *   Cost per query: ≤1 find_sequences_for_tag + 1 decompressSA per target-
 *   visited source mapping. No full path enumeration, no cross-locus walk.
 *
 * Legacy algorithm (single-occurrence + LF-walk — opt-in):
 *   Set the env var PANGENOME_SURJECT_ANCHOR_WALK to use the older strategy:
 *   pick the first/last target occurrence of the boundary nodes, then LF-walk
 *   the target path between them pairing nodes FIFO. Retained only for A/B
 *   comparison; on repeated boundary nodes it can walk a chromosome-scale span
 *   and seed the wrong occurrence.
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

/// Diagnostics for the target-path GBWT LF walk (walk_target_collecting_matches).
/// The walk's length = number of gbz.index.LF steps between the first and last
/// anchor on the target path; when boundary nodes recur, the chosen anchor
/// occurrences can be chromosome-scale apart and dominate build time. Reset
/// before a measured region and read afterwards. Thread-local; external linkage.
struct AnchorWalkStats {
    size_t walk_lf_steps = 0;     ///< total gbz.index.LF calls during the walk(s)
    size_t walks = 0;             ///< number of walk_target_collecting_matches calls
    size_t first_anchor_base = 0; ///< base offset of the first anchor (last successful walk)
    size_t last_anchor_base = 0;  ///< base offset of the last anchor (last successful walk)
    size_t walk_span = 0;         ///< last_anchor_base - first_anchor_base (bases spanned)
    double walk_ms = 0.0;         ///< wall-clock time in walk_target_collecting_matches
    size_t decompress_sa_calls = 0;   ///< number of gbwt FastLocate decompressSA() calls
    size_t decompress_sa_entries = 0; ///< total SA values returned across those calls
    double decompress_sa_ms = 0.0;    ///< wall-clock time in decompressSA()
};
extern thread_local AnchorWalkStats g_anchor_walk_stats;

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
