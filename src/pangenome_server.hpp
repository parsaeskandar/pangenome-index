#ifndef PANGENOME_SERVER_HPP
#define PANGENOME_SERVER_HPP

#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include "pangenome_index/surject_anchor_builder.hpp"
#include "pangenome_index/translation_tables.hpp"
#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbz.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

/// How much of one alignment a single haplotype accounts for.
///
/// Unlike the engine's "carried by" list — which reports only haplotypes that
/// thread the read's EXACT allele path and therefore drops a haplotype that
/// differs at a single variant — this is a graded score, so a haplotype that
/// matches everywhere except one small site scores near 100 instead of being
/// omitted.
struct HaplotypeCoverage {
    std::string haplotype;    ///< two-field haplotype name, e.g. "HG00097#1"
    uint64_t covered_bp = 0;  ///< aligned read bases on nodes this haplotype visits
    double coverage = 0.0;    ///< covered_bp as a percentage of aligned bases, 0..100
};

struct TranslatedInterval {
    std::string haplotype;
    int64_t start;
    int64_t end;
    char strand;  // '+' or '-'
};

/// Per-source-fragment accounting for one translation, so it is possible to see
/// WHERE bases are lost instead of only that the total came up short. A large
/// interval is split by Table 1 into many GBWT path fragments, each translated
/// independently; any fragment that fails contributes nothing and is otherwise
/// invisible.
struct FragmentDiag {
    uint64_t src_path_id = 0;
    uint64_t extent_start = 0;     ///< path-local extent handed to the trace
    uint64_t extent_end = 0;       ///< exclusive
    uint64_t extent_bp = 0;        ///< bases this fragment was asked to cover
    uint64_t unscoped_tags = 0;    ///< nodes the source visits here
    uint64_t scoped_tags = 0;      ///< nodes shared with the target (1 == the
                                   ///< "extended search" case: a single anchor
                                   ///< outside the interval, which maps ~nothing)
    uint32_t candidates = 0;       ///< candidate target paths found by probing
    uint64_t points = 0;           ///< per-base correspondences produced
    uint64_t mapped_span = 0;      ///< last mapped source base - first + 1
    // Anchor geometry: the two numbers that decide whether the anchors are
    // colinear. If (last_target_base - first_target_base) is nothing like
    // (last_source_base - first_source_base), a paralogous copy was anchored.
    uint64_t first_source_base = 0;
    uint64_t first_target_base = 0;
    uint64_t last_source_base = 0;
    uint64_t last_target_base = 0;
    bool first_unique = false;     ///< anchor sat on a node unique to both
    bool last_unique = false;
    uint32_t diag_version = 2;     ///< bump on change; 0/absent => stale build
};

/// Aggregate view over all fragments of one translation.
struct TranslationDiagnostics {
    uint64_t fragments = 0;        ///< fragments Table 1 split the request into
    uint64_t no_tags = 0;          ///< no source tags at all
    uint64_t no_candidates = 0;    ///< probing found no target path
    uint64_t traced = 0;           ///< a trace was attempted
    uint64_t productive = 0;       ///< produced at least one point
    uint64_t empty_trace = 0;      ///< traced but produced nothing
    uint64_t single_anchor = 0;    ///< scoped_tags == 1 (extended-search case)
    uint64_t requested_bp = 0;     ///< sum of fragment extents
    uint64_t mapped_bp = 0;        ///< total points produced
    std::vector<FragmentDiag> detail;   ///< capped, see MAX_FRAGMENT_DETAIL
};

/// translate() plus the per-fragment accounting above.
struct DiagnosedTranslation {
    std::vector<TranslatedInterval> intervals;
    TranslationDiagnostics diagnostics;
    double elapsed_ms = 0.0;
};

/// Outcome of a translation that may be cut short by a deadline.
struct TranslationRun {
    std::vector<TranslatedInterval> intervals;
    /// True if the deadline fired before the query finished. `intervals` then
    /// holds whatever was produced up to that point (possibly empty), so a
    /// caller can surface a partial answer plus a warning rather than nothing.
    bool timed_out = false;
    double elapsed_ms = 0.0;
};

/// Plain-data version of a PrecomputedAnchor suitable for Python bindings.
/// Same fields as panindexer::PrecomputedAnchor, with step_handle_t replaced
/// by a portable (node, offset) pair (the gbwtgraph step encoding —
/// gbwtgraph.cpp:911-1024 stores edge.first/edge.second into as_integers(step)).
struct AnchorRecord {
    /// Source-mapping range covered by this anchor (read order).
    uint64_t source_mapping_begin = 0;
    uint64_t source_mapping_end = 0;
    /// Read range covered.
    uint64_t read_begin_offset = 0;
    uint64_t read_end_offset = 0;
    /// Base positions on the target path.
    uint64_t path_offset_step_begin = 0;
    uint64_t path_offset_step_end = 0;
    /// GBWT search states (= gbwt::edge_type). Pairs of (node, offset_in_record).
    /// step_handle reconstruction: as_integers(step)[0]=node, [1]=offset.
    uint64_t gbwt_edge_begin_node = 0;
    uint64_t gbwt_edge_begin_offset = 0;
    uint64_t gbwt_edge_end_node = 0;
    uint64_t gbwt_edge_end_offset = 0;
};

struct AnchorBuildPyResult {
    /// One of: "ok", "empty_alignment", "unknown_path", "no_common_nodes".
    std::string status;
    std::vector<AnchorRecord> anchors;
    uint64_t target_path_length = 0;
    bool target_rev_strand = false;

    /// Diagnostics for the dominant cost: find_sequences_for_tag's RLBWT
    /// enumeration. Measured over the whole build (all subpath attempts).
    ///   find_seq_calls       — number of find_sequences_for_tag invocations
    ///   find_seq_runs        — total tag runs ("vectors") iterated
    ///   find_seq_lf_steps    — total locateNext (LF) calls (navigation + walk)
    ///   find_seq_visits      — total node visits enumerated (node's pangenome usage)
    ///   last_run_nav_steps   — LF steps to navigate to the LAST run's start
    ///   last_run_length      — number of positions in the LAST run iterated
    uint64_t find_seq_calls = 0;
    uint64_t find_seq_runs = 0;
    uint64_t find_seq_lf_steps = 0;
    uint64_t find_seq_visits = 0;
    uint64_t last_run_nav_steps = 0;
    uint64_t last_run_length = 0;

    /// Diagnostics for the target-path GBWT LF walk (the other candidate
    /// bottleneck). walk_lf_steps is the number of gbz.index.LF calls; walk_span
    /// is how far apart (in target bases) the chosen first/last anchors are —
    /// a large span means the walk traverses a big chunk of the target path.
    uint64_t walk_lf_steps = 0;
    uint64_t walk_span = 0;
    uint64_t first_anchor_base = 0;
    uint64_t last_anchor_base = 0;

    /// Wall-clock attribution of the build (ms): which phase actually eats the
    /// time when LF counts are small (points at per-call cost / cold memory).
    double find_seq_ms = 0.0;       ///< time in find_sequences_for_tag
    double decompress_sa_ms = 0.0;  ///< time in gbwt decompressSA
    double walk_ms = 0.0;           ///< time in the target-path walk
    uint64_t decompress_sa_calls = 0;
    uint64_t decompress_sa_entries = 0;

    /// How many target subpaths the name resolved to (we currently run the full
    /// per-path build against EACH), and how many source mappings the GAF
    /// produced. decompress_sa_calls should track ~n_target_subpaths × 2 ×
    /// n_source_mappings — these confirm the call explosion is the subpath loop.
    uint64_t n_target_subpaths = 0;
    uint64_t n_source_mappings = 0;
};

class Index {
public:
    Index() = default;

    /// Load all indexes into memory.
    ///
    /// @param gbz_path        Path to the GBZ file (GBWT + GBWTGraph).
    /// @param ri_path         Path to the RLBWT r-index (.ri, encoded format).
    /// @param tags_path       Path to the sampled tag array (.tags).
    /// @param gbwt_ri_path    Path to the GBWT FastLocate (.ri).
    /// @param table1_path     Path to Translation Table 1 (.t1).
    /// @param table2_path     Path to Translation Table 2 (.t2).
    /// `table2_path` may be empty: translation then uses the table-free path
    /// (translate_no_table2), which needs only Table 1 plus the GBWT/tag array.
    void load(const std::string& gbz_path,
              const std::string& ri_path,
              const std::string& tags_path,
              const std::string& gbwt_ri_path,
              const std::string& table1_path,
              const std::string& table2_path = "");

    /// True if a Table 2 was loaded. When false, translate() and
    /// translatable_haplotypes() automatically use their table-free forms.
    bool has_table2() const { return has_table2_; }

    /// Run coordinate translation from source haplotype interval to target.
    /// Throws std::invalid_argument if (end - start) > 10 000 000.
    std::vector<TranslatedInterval>
    translate(const std::string& src_haplotype,
              int64_t start, int64_t end,
              const std::string& tgt_haplotype) const;

    /// Discovery query: for a source contig interval, return the set of target
    /// haplotype names (2-field) that have a homologous region overlapping it —
    /// i.e. every haplotype this interval CAN be translated to. This is a pure
    /// Table-2 overlap check (no coordinate trace), so it is much cheaper than
    /// translating to each. Result is sorted and de-duplicated; may include the
    /// source's own haplotype (identity homology).
    /// Throws std::invalid_argument on a bad/oversized interval or unknown source.
    std::vector<std::string>
    translatable_haplotypes(const std::string& src_haplotype,
                            int64_t start, int64_t end) const;

    /// Score every haplotype by how much of a graph alignment it accounts for.
    ///
    /// For each node the alignment visits, the tag array / GBWT reports which
    /// haplotypes also visit that node; a haplotype is credited with the read
    /// bases aligned there. The score is those bases as a percentage of all
    /// aligned bases, so 100 means the haplotype visits every node the read
    /// does. Both node orientations count, so reverse-strand alignments and
    /// inverted haplotypes are scored correctly.
    ///
    /// Results are sorted by descending coverage. `min_coverage` (0..100) drops
    /// haplotypes below the threshold; pass 0 (the default) to report every
    /// haplotype that shares any node, however partial the match.
    ///
    /// `include_zero` additionally lists haplotypes that share NO node with the
    /// alignment, scored 0, so the result covers every haplotype in the graph
    /// rather than only those with some overlap.
    std::vector<HaplotypeCoverage>
    haplotype_coverage(const std::string& gaf_str, double min_coverage = 0.0,
                       bool include_zero = false) const;

    /// Which haplotypes a source interval can reach, each with a coverage score
    /// — the table-free counterpart of translatable_haplotypes().
    ///
    /// Walks the nodes the source visits in the interval and asks the tag array
    /// which haplotypes also visit them, crediting each haplotype the source
    /// bases on the nodes it shares. `coverage` is those bases as a percentage
    /// of the interval, so 100 means the haplotype covers the whole region and
    /// a small value means it shares only a fragment (often a repeat).
    ///
    /// `min_coverage` (0..100) filters the result; `max_nodes` caps how many
    /// nodes are probed (0 = all), trading precision for speed on wide
    /// intervals. Sorted by descending coverage.
    std::vector<HaplotypeCoverage>
    translatable_haplotypes_scored(const std::string& src_haplotype,
                                   int64_t start, int64_t end,
                                   double min_coverage = 0.0,
                                   size_t max_nodes = 0) const;

    /// Coordinate translation WITHOUT Table 2.
    ///
    /// Table 2 exists only to route a source interval to candidate target paths
    /// and to narrow the trace extent. This finds the candidates directly: it
    /// walks inward from both ends of the source interval until it meets a node
    /// the target haplotype also visits (the first and last common nodes), then
    /// runs the same GBWT trace. Nothing is precomputed per path pair, so it is
    /// unaffected by how finely the graph fragments haplotypes into GBWT paths.
    ///
    /// Selected at runtime by setting PANGENOME_TRANSLATE_NO_T2=1, which makes
    /// translate() delegate here; PANGENOME_TRANSLATE_PROBE_CAP bounds how many
    /// nodes are probed from each end before giving up (default 256).
    /// `timeout_ms > 0` bounds the work: the deadline is checked between source
    /// fragments, between candidate target paths, and between node probes, so a
    /// pathological target is abandoned instead of stalling the whole request.
    /// It cannot preempt a single long call inside those steps, so the actual
    /// stop can overshoot slightly. `*timed_out` is set when it fires.
    std::vector<TranslatedInterval>
    translate_no_table2(const std::string& src_haplotype,
                        int64_t start, int64_t end,
                        const std::string& tgt_haplotype,
                        double timeout_ms = 0.0,
                        bool* timed_out = nullptr,
                        TranslationDiagnostics* diag = nullptr) const;

    /// Translate and report per-fragment accounting: how many fragments the
    /// request was split into, how many produced anything, and how many bases
    /// each covered. Use this to locate where a short result lost its bases.
    DiagnosedTranslation
    translate_diagnosed(const std::string& src_haplotype,
                        int64_t start, int64_t end,
                        const std::string& tgt_haplotype) const;

    /// translate() with a per-query deadline, reporting whether it fired.
    /// Use this when one slow haplotype must not hold up a multi-target query.
    TranslationRun
    translate_checked(const std::string& src_haplotype,
                      int64_t start, int64_t end,
                      const std::string& tgt_haplotype,
                      double timeout_ms) const;

    /// Return all valid haplotype names present in the loaded index.
    std::vector<std::string> get_haplotype_names() const;

    /// Build surjection anchors from a graph GAF alignment against the named
    /// target haplotype. Parses the GAF path + CIGAR into per-node source
    /// mappings, then runs panindexer::build_surject_anchors.
    ///
    /// Returns one AnchorBuildPyResult; if the haplotype maps to multiple
    /// subpaths, the result with the longest non-empty anchor list wins.
    /// status is one of: "ok", "empty_alignment", "unknown_path",
    /// "no_common_nodes", "parse_error".
    AnchorBuildPyResult build_surject_anchors(const std::string& gaf_str,
                                              const std::string& target_haplotype) const;

private:
    bool loaded_ = false;

    panindexer::FastLocate rlbwt_rindex_;
    panindexer::SampledTagArray sampled_;
    std::unique_ptr<gbwtgraph::GBZ> gbz_;
    std::unique_ptr<gbwt::FastLocate> gbwt_rindex_;
    panindexer::TranslationTable1 table1_;
    panindexer::TranslationTable2 table2_;
    bool has_table2_ = false;
    /// Haplotype names, cached at load(): deriving them walks every GBWT path.
    std::vector<std::string> haplotype_names_;
    std::unordered_map<size_t, std::pair<std::string, size_t>> path_to_global_;
};

#endif // PANGENOME_SERVER_HPP
