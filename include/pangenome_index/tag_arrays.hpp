
#ifndef VG_TAG_ARRAYS_HPP
#define VG_TAG_ARRAYS_HPP


#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/sd_vector.hpp>
#include <gbwt/utils.h>
#include <gbwt/internal.h>
#include <iostream>
#include <cstdint>
#include <bitset>
#include <stdexcept>
#include "bplus_tree.hpp"
#include <gbwtgraph/utils.h>
//#include <vg/io/vpkg.hpp>


namespace panindexer {

    using namespace std;
    using namespace gbwtgraph;


    class TagArray {
    public:
        TagArray();

//    void TagArray::add_run(size_t offset, bool is_reverse, int length, size_t node_id);

        void load_bptree(BplusTree<Run> &bptree, size_t bwt_size);
        void serialize(std::ostream& out);
        void load(std::istream& in);
        void query(size_t start, size_t end);

        void serialize_run_by_run(const std::vector<std::pair<pos_t, uint8_t>>& tag_runs, const std::string& file_name);
        void deserialize_run_by_run(const std::string& file_name);

        void store_blocks(std::ostream& out);
        static std::pair<pos_t, uint8_t> load_block_at(std::istream& in, size_t &next_block_start);

        vector<pair<pos_t, uint8_t>> get_tag_runs(){
            return tag_runs;
        };

    private:

        //
        std::vector<gbwt::byte_type> encoded_runs;


        sdsl::bit_vector encoded_runs_starts;


        sdsl::sd_vector<> bwt_intervals;

        sdsl::bit_vector::select_1_type encoded_runs_starts_select;
        sdsl::sd_vector<>::rank_1_type bwt_intervals_rank;

        vector<pair<pos_t, uint8_t>> tag_runs;


    };

}

#endif //VG_TAG_ARRAYS_HPP