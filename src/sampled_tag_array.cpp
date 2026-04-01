#include "pangenome_index/sampled_tag_array.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include <sdsl/util.hpp>

namespace panindexer {

    using handlegraph::pos_t;

    static inline uint64_t encode_val_from_pos(pos_t p) {
        if (gbwtgraph::offset(p) != 0) return 0; // gap within node
        if (gbwtgraph::id(p) == 0) return 0;     // treat special zero-tag as gap
        return SampledTagArray::encode_value(gbwtgraph::id(p), gbwtgraph::is_rev(p));
    }

    SampledTagArray::SampledTagArray() {
        this->init_supports();
    }

    SampledTagArray::SampledTagArray(const SampledTagArray& source) :
            sampled_values(source.sampled_values),
            bwt_intervals(source.bwt_intervals),
            first_run_is_gap(source.first_run_is_gap) {
        this->init_supports();
    }

    SampledTagArray& SampledTagArray::operator=(const SampledTagArray& source) {
        if (this != &source) {
            this->sampled_values = source.sampled_values;
            this->bwt_intervals = source.bwt_intervals;
            this->first_run_is_gap = source.first_run_is_gap;
            this->init_supports();
        }
        return *this;
    }

    SampledTagArray::SampledTagArray(SampledTagArray&& source) noexcept :
            sampled_values(std::move(source.sampled_values)),
            bwt_intervals(std::move(source.bwt_intervals)),
            first_run_is_gap(source.first_run_is_gap) {
        this->init_supports();
        source.init_supports();
    }

    SampledTagArray& SampledTagArray::operator=(SampledTagArray&& source) noexcept {
        if (this != &source) {
            this->sampled_values = std::move(source.sampled_values);
            this->bwt_intervals = std::move(source.bwt_intervals);
            this->first_run_is_gap = source.first_run_is_gap;
            this->init_supports();
            source.init_supports();
        }
        return *this;
    }

    void SampledTagArray::init_supports() {
        sdsl::util::init_support(run_rank_support, &bwt_intervals);
        sdsl::util::init_support(run_select_support, &bwt_intervals);
    }

    // Helper function to construct wt_gmr from values
    // wt_gmr uses grammar-based compression which is more memory-efficient and stable on macOS
    static void construct_wt_gmr_from_values(sdsl::wt_gmr<sdsl::int_vector<>, sdsl::inv_multi_perm_support<4, sdsl::int_vector<>>>& target, const std::vector<uint64_t>& values) {
        if (values.empty()) {
            target = sdsl::wt_gmr<sdsl::int_vector<>, sdsl::inv_multi_perm_support<4, sdsl::int_vector<>>>();
            return;
        }
        
        uint64_t maxv = 0;
        for (auto v : values) { if (v > maxv) maxv = v; }
        size_t width = std::max<size_t>(1, sdsl::bits::hi(maxv) + 1);
        
        // Create int_vector and populate it
        sdsl::int_vector<> iv(values.size(), 0, width);
        for (size_t i = 0; i < values.size(); ++i) {
            iv[i] = values[i];
        }
        
        // Construct wt_gmr directly - it should handle memory more reliably than wm_int
        sdsl::construct_im(target, iv);
        
        // iv will be destroyed when it goes out of scope
        // target should have its own independent copy of the data
    }

    void SampledTagArray::build_from_runs(const std::vector<std::pair<pos_t, uint64_t>>& runs, size_t bwt_size) {
        // First, transform runs into (start position markers) using cumulative lengths, and values merged.
        // We'll collect run_starts positions and run values in order.
        std::vector<uint64_t> values;
        values.reserve(runs.size());

        std::vector<uint64_t> starts;
        starts.reserve(runs.size());

        uint64_t cum = 0;
        uint64_t last_val = std::numeric_limits<uint64_t>::max();

        for (const auto& pr : runs) {
            pos_t p = pr.first;
            uint64_t len = pr.second;
            uint64_t val = encode_val_from_pos(p);
            if (values.empty() || val != last_val) {
                starts.push_back(cum);
                values.push_back(val);
                last_val = val;
            }
            cum += len;
        }

        // Determine if first run is gap (false if gap, true if normal tag)
        first_run_is_gap = (values.empty() || values[0] == 0);

        // Build bwt_intervals with zero-length gap runs inserted between adjacent normal tags
        // Only store non-gap values in wt_gmr
        std::vector<uint64_t> bwt_starts; // All run starts including gaps
        std::vector<uint64_t> non_gap_values; // Only non-gap values for wt_gmr
        
        for (size_t i = 0; i < values.size(); ++i) {
            uint64_t val = values[i];
            uint64_t start = starts[i];
            
            // Add start position for this run
            bwt_starts.push_back(start);
            
            // If this is a non-gap value, store it in wt_gmr
            if (val != 0) {
                non_gap_values.push_back(val);
                
                // Insert zero-length gap run between this and next run if next run is also non-gap
                if (i + 1 < values.size() && values[i + 1] != 0) {
                    // Insert gap run at the start position of next run (zero-length gap)
                    bwt_starts.push_back(starts[i + 1]);
                }
            }
        }

        // Build sd_vector for run starts with multiset=true
        {
            sdsl::sd_vector_builder builder(bwt_size + 1, bwt_starts.size(), true); // multiset=true
            for (uint64_t s : bwt_starts) builder.set(s);
            bwt_intervals = sdsl::sd_vector<>(builder);
        }

        // Build wt_gmr from only non-gap values
        construct_wt_gmr_from_values(sampled_values, non_gap_values);
        this->init_supports();
    }

    void SampledTagArray::build_from_enumerator(const std::function<void(const std::function<void(pos_t,uint64_t)>&)>& enumerator,
                                                size_t bwt_size) {
        std::vector<uint64_t> values;
        std::vector<uint64_t> starts;
        uint64_t cum = 0;
        uint64_t last_val = std::numeric_limits<uint64_t>::max();

        auto sink = [&](pos_t p, uint64_t len) {
            uint64_t val = encode_val_from_pos(p);
            if (values.empty() || val != last_val) {
                starts.push_back(cum);
                values.push_back(val);
                last_val = val;
            }
            cum += len;
        };

        enumerator(sink);

        if (cum != bwt_size) {
            std::cerr << "WARNING: Sampled tag array build: streamed run lengths sum to " << cum
                      << " but bwt_size is " << bwt_size << ". Alignment may be wrong." << std::endl;
        }

        // Determine if first run is gap (with bounds check for empty input)
        first_run_is_gap = (values.empty() || values[0] == 0);

        // Build bwt_intervals with zero-length gap runs inserted between adjacent normal tags
        // Only store non-gap values in wt_gmr
        std::vector<uint64_t> bwt_starts; // All run starts including gaps
        std::vector<uint64_t> non_gap_values; // Only non-gap values for wt_gmr
        
        for (size_t i = 0; i < values.size(); ++i) {
            uint64_t val = values[i];
            uint64_t start = starts[i];
            
            // Add start position for this run
            bwt_starts.push_back(start);
            
            // If this is a non-gap value, store it in wt_gmr
            if (val != 0) {
                non_gap_values.push_back(val);
                
                // Insert zero-length gap run between this and next run if next run is also non-gap
                if (i + 1 < values.size() && values[i + 1] != 0) {
                    // Insert gap run at the start position of next run (zero-length gap)
                    bwt_starts.push_back(starts[i + 1]);
                }
            }
        }
        
        // Build sd_vector for run starts with multiset=true
        {
            sdsl::sd_vector_builder builder(bwt_size + 1, bwt_starts.size(), true); // multiset=true
            for (uint64_t s : bwt_starts) builder.set(s);
            bwt_intervals = sdsl::sd_vector<>(builder);
        }

        std::cerr << "The size of the bwt is: " << bwt_size << std::endl;
        std::cerr << "The size of the bwt_starts vector is: " << bwt_starts.size() << std::endl;
        std::cerr << "The size of the bwt_intervals vector is: " << bwt_intervals.size() << std::endl;

        std::cerr << "The size of the non_gap_values vector is: " << non_gap_values.size() << std::endl;
        // Build wt_gmr from only non-gap values
        construct_wt_gmr_from_values(sampled_values, non_gap_values);
        this->init_supports();
        
        std::cerr << "Finished building sampled_tag_array" << std::endl;
    }

    void SampledTagArray::print_bwt_intervals_and_rank(size_t limit, std::ostream& out) const {
        // Same rank as run_id_at: rank_1(pos+1) = number of 1s in [0..pos]; run_id = rank - 1 (so run_id_at is correct).
        size_t n = std::min(limit, static_cast<size_t>(bwt_intervals.size()));
        size_t nruns = total_runs();
        out << "is_first_run_gap=" << first_run_is_gap << "\n";
        out << "bwt_intervals (pos, bit, rank) up to " << n << " positions (rank = rank_1(pos+1), run_id = rank-1):\n";
        for (size_t i = 0; i < n; ++i) {
            size_t r = run_rank_support(i + 1);
            out << "  pos=" << i << " bit=" << static_cast<int>(bwt_intervals[i]) << " rank=" << r << "\n";
        }
        out << "select_1(1..10) [position of j-th run start]:\n";
        for (size_t j = 1; j <= 10 && j <= nruns; ++j) {
            size_t sel = run_select_support(j);
            out << "  select(" << j << ")=" << sel << "\n";
        }
    }

    void SampledTagArray::serialize(std::ostream& out) const {
        // On macOS, ensure supports are not accessed during serialization
        // Serialize the structures directly without accessing lazy supports
        sdsl::serialize(sampled_values, out);
        sdsl::serialize(bwt_intervals, out);
        sdsl::write_member(first_run_is_gap, out);
    }

    void SampledTagArray::load(std::istream& in) {
        sampled_values.load(in);
        bwt_intervals.load(in);
        sdsl::read_member(first_run_is_gap, in);
        this->init_supports();
    }

}


