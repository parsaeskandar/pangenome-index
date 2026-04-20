#ifndef PANGENOME_SERVER_HPP
#define PANGENOME_SERVER_HPP

#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include "pangenome_index/translation_tables.hpp"
#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbz.h>
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

/// One step along a GBWT-graph alignment path (e.g. from Giraffe). Used for downstream haplotype resolution.
struct GraphPathVisit {
    int64_t node_id = 0;
    bool is_reverse = false;
    int32_t from_length = 0;
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

    /// Receive graph coordinates from an upstream mapper. Stub until haplotype walk is implemented.
    void accept_graph_path(const std::vector<GraphPathVisit>& path);

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
