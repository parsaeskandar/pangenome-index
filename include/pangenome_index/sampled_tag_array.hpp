#ifndef PANGENOME_INDEX_SAMPLED_TAG_ARRAY_HPP
#define PANGENOME_INDEX_SAMPLED_TAG_ARRAY_HPP

#include <cstdint>
#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <gbwtgraph/utils.h>

namespace panindexer {

    class SampledTagArray {
    public:
        SampledTagArray() = default;

        // Build from a stream of runs: for each input run (pos_t, length),
        // emit value = encode(node_id,is_rev) if offset==0, else GAP_CODE (0), merging consecutive runs with equal value.
        void build_from_runs(const std::vector<std::pair<handlegraph::pos_t, uint64_t>>& runs, size_t bwt_size);

        // Build from a callback enumerator (yields many runs without materializing all of them)
        void build_from_enumerator(const std::function<void(const std::function<void(handlegraph::pos_t,uint64_t)>&)>& enumerator,
                                   size_t bwt_size);

        // Serialization
        void serialize(std::ostream& out) const;
        void load(std::istream& in);

        // Accessors
        inline const sdsl::wm_int<>& values() const { return sampled_values; }
        inline const sdsl::sd_vector<>& run_starts() const { return bwt_intervals; }

        // Encode (node_id,is_rev) into integer code; 0 reserved for gaps
        static inline uint64_t encode_value(int64_t node_id, bool is_rev) {
            // node_id >= 1; map to 1.. via shift and set bit0 for strand
            return 1 + ((static_cast<uint64_t>(node_id - 1) << 1) | static_cast<uint64_t>(is_rev));
        }

    private:
        sdsl::wm_int<> sampled_values; // one value per run (0 for gap, otherwise encode(node_id,is_rev))
        sdsl::sd_vector<> bwt_intervals; // 1 at BWT positions where a run starts
    };

}

#endif // PANGENOME_INDEX_SAMPLED_TAG_ARRAY_HPP


