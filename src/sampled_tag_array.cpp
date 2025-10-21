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

        // Build sd_vector for run starts
        sdsl::sd_vector_builder builder(bwt_size + 1, starts.size());
        for (uint64_t s : starts) builder.set(s);
        bwt_intervals = sdsl::sd_vector<>(builder);

        // Build wm_int from values (guard for empty and set width explicitly)
        if (values.empty()) {
            sampled_values = sdsl::wm_int<>();
        } else {
            uint64_t maxv = 0;
            for (auto v : values) { if (v > maxv) maxv = v; }
            size_t width = std::max<size_t>(1, sdsl::bits::hi(maxv) + 1);
            sdsl::int_vector<> iv(values.size(), 0, width);
            for (size_t i = 0; i < values.size(); ++i) iv[i] = values[i];
            sdsl::construct_im(sampled_values, iv);
        }
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
        std::cerr << "Enumerated runs" << std::endl;

        // print the last 10 values of the starts vector and the last 10 values of the values vector
        std::cerr << "Last 10 values of starts: ";
        for (size_t i = 0; i < 10; ++i) {
            std::cerr << starts[starts.size() - 10 + i] << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Last 10 values of values: ";
        for (size_t i = 0; i < 10; ++i) {
            // Decode the value before printing (assume pos_t is int64_t and direction encoded in bit 0, id in high bits)
            uint64_t encoded_val = values[values.size() - 10 + i];
            if (encoded_val == 0) {
                std::cerr << "[GAP] ";
            } else {
                int64_t node_id = static_cast<int64_t>((encoded_val - 1) >> 1) + 1;
                bool is_rev = (encoded_val - 1) & 1;
                std::cerr << "(node_id=" << node_id << ",is_rev=" << is_rev << ") ";
            }
        }
        std::cerr << std::endl;
        // print value size and start size  
        std::cerr << "Value size: " << values.size() << std::endl;
        std::cerr << "Start size: " << starts.size() << std::endl;
        // print all values
        std::cerr << "All values: ";
        for (size_t i = 0; i < values.size(); ++i) {
            std::cerr << values[i] << " ";
        }
        std::cerr << std::endl;
        sdsl::sd_vector_builder builder(bwt_size + 1, starts.size());
        for (uint64_t s : starts) builder.set(s);
        bwt_intervals = sdsl::sd_vector<>(builder);
        std::cerr << "Built bwt_intervals" << std::endl;
        if (values.empty()) {
            sampled_values = sdsl::wm_int<>();
        } else {
            uint64_t maxv = 0;
            for (auto v : values) { if (v > maxv) maxv = v; }
            size_t width = std::max<size_t>(1, sdsl::bits::hi(maxv) + 1);
            sdsl::int_vector<> iv(values.size(), 0, width);
            std::cerr << "Created int_vector of size: " << values.size() << " width=" << width << std::endl;
            // print bwt size
            std::cerr << "Bwt size: " << bwt_size << std::endl;
            for (size_t i = 0; i < values.size(); ++i) iv[i] = values[i];
            std::cerr << "Filled int_vector" << std::endl;
            sdsl::construct_im(sampled_values, iv);
            std::cerr << "Constructed sampled_values" << std::endl;
        }
        std::cerr << "Finished building sampled_tag_array" << std::endl;
    }

    void SampledTagArray::serialize(std::ostream& out) const {
        sdsl::serialize(sampled_values, out);
        sdsl::serialize(bwt_intervals, out);
    }

    void SampledTagArray::load(std::istream& in) {
        sampled_values.load(in);
        bwt_intervals.load(in);
    }

}


