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

struct TranslatedInterval {
    std::string haplotype;
    int64_t start;
    int64_t end;
    char strand;  // '+' or '-'
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
    void load(const std::string& gbz_path,
              const std::string& ri_path,
              const std::string& tags_path,
              const std::string& gbwt_ri_path,
              const std::string& table1_path,
              const std::string& table2_path);

    /// Run coordinate translation from source haplotype interval to target.
    /// Throws std::invalid_argument if (end - start) > 10 000 000.
    std::vector<TranslatedInterval>
    translate(const std::string& src_haplotype,
              int64_t start, int64_t end,
              const std::string& tgt_haplotype) const;

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
    std::unordered_map<size_t, std::pair<std::string, size_t>> path_to_global_;
};

#endif // PANGENOME_SERVER_HPP
