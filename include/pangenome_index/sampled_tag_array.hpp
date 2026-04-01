#ifndef PANGENOME_INDEX_SAMPLED_TAG_ARRAY_HPP
#define PANGENOME_INDEX_SAMPLED_TAG_ARRAY_HPP

#include <cstdint>
#include <cassert>
#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_gmr.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/util.hpp>
#include <ostream>
#include <gbwtgraph/utils.h>

namespace panindexer {

    class SampledTagArray {
    public:
        SampledTagArray();
        SampledTagArray(const SampledTagArray& source);
        SampledTagArray& operator=(const SampledTagArray& source);
        SampledTagArray(SampledTagArray&& source) noexcept;
        SampledTagArray& operator=(SampledTagArray&& source) noexcept;

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
        inline const sdsl::wt_gmr<sdsl::int_vector<>, sdsl::inv_multi_perm_support<4, sdsl::int_vector<>>>& values() const { return sampled_values; }
        inline const sdsl::sd_vector<>& run_starts() const { return bwt_intervals; }
        inline bool is_first_run_gap() const { return first_run_is_gap; }

        // Compatibility no-ops: supports are eagerly initialized.
        inline void ensure_run_rank() const {}
        inline void ensure_run_select() const {}

        // Debug: print is_first_run_gap, then bwt_intervals bit and rank for positions [0, limit)
        void print_bwt_intervals_and_rank(size_t limit, std::ostream& out = std::cerr) const;

        // Helpers for queries
        inline size_t total_runs() const { 
            // Total runs includes both gap and non-gap runs
            // bwt_intervals has one bit per run (including gaps)
            return run_rank_support(bwt_intervals.size());
        }

        // Return run id (0-based) that contains BWT position pos.
        // Uses rank: run containing pos is the one whose start is the last 1-bit at or before pos.
        // SDSL rank_1(i) = number of 1s in [0..i-1], so rank_1(pos+1) = number of 1s in [0..pos].
        // So run_id = rank_1(pos+1) - 1 (0-based); when rank==0 we return 0 (pos before first run start).
        inline size_t run_id_at(size_t pos) const {
            if (pos >= bwt_intervals.size()) {
                pos = bwt_intervals.size() - 1;
            }
            size_t r = run_rank_support(pos + 1);
            return (r == 0) ? 0 : (r - 1);
        }

        // Return [start,end] BWT span for run_id
        inline std::pair<size_t,size_t> run_span(size_t run_id) const {
            size_t start = run_select_support(run_id + 1);
            size_t end;
            if (run_id + 1 < total_runs()) {
                end = run_select_support(run_id + 2) - 1;
            } else {
                end = bwt_intervals.size() - 1;
            }
            return { start, end };
        }

        // Return encoded tag value for a run
        // If run_id corresponds to a gap run, returns 0
        // Otherwise maps run_id to correct index in wt_gmr based on gap flag
        inline uint64_t run_value(size_t run_id) const {
            // Check if this run_id is a gap run
            if ((run_id + static_cast<size_t>(first_run_is_gap)) % 2 == 1) {
                // This run_id is a gap; return 0
                return 0;
            } else {
                // This run_id is a normal tag; access wt_gmr at (run_id - first_run_is_gap)/2
                return sampled_values[(run_id - static_cast<size_t>(first_run_is_gap)) / 2];
            }
        }

        // Encode (node_id,is_rev) into integer code; 0 reserved for gaps
        static inline uint64_t encode_value(int64_t node_id, bool is_rev) {
            // node_id >= 1; map to 1.. via shift and set bit0 for strand
            return 1 + ((static_cast<uint64_t>(node_id - 1) << 1) | static_cast<uint64_t>(is_rev));
        }

    private:
        void init_supports();

        sdsl::wt_gmr<sdsl::int_vector<>, sdsl::inv_multi_perm_support<4, sdsl::int_vector<>>> sampled_values; // only non-gap values (gaps not stored)
        sdsl::sd_vector<> bwt_intervals; // 1 at BWT positions where a run starts (including zero-length gap runs)
        bool first_run_is_gap = true; // 1 if first run is gap, 0 if first run is normal tag

        // Eager supports for run_starts
        sdsl::sd_vector<>::rank_1_type run_rank_support;
        sdsl::sd_vector<>::select_1_type run_select_support;
    };

}

#endif // PANGENOME_INDEX_SAMPLED_TAG_ARRAY_HPP


