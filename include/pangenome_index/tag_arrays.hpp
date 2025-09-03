
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

//    using namespace std;
    using namespace gbwtgraph;
    using handlegraph::pos_t;


    class TagArray {
    public:
//        TagArray(int length_bits);
        TagArray();

//    void TagArray::add_run(size_t offset, bool is_reverse, int length, size_t node_id);

        void load_bptree(BplusTree<Run> &bptree, size_t bwt_size);
        void load_bptree_lite(BplusTree <Run> &bptree);
        void serialize(std::ostream& out);
        void load(std::istream& in);
        void query(size_t start, size_t end);

        void serialize_run_by_run(std::ofstream& out, const std::vector<std::pair<pos_t, uint16_t>>& tag_runs);
        void serialize_run_by_run_batch(sdsl::int_vector_buffer<8>& out, const std::vector<std::pair<gbwtgraph::Position, uint16_t>>& tag_runs);
        void deserialize_run_by_run(const std::string& file_name);
        void serialize_bptree_lite(std::string filename, BplusTree <Run> &bptree);


        void store_blocks(std::ostream& out);
        static std::pair<pos_t, uint16_t> load_block_at(std::istream& in, size_t &next_block_start);

        void store_blocks_sdsl(std::string filename);

        static std::pair<pos_t, uint16_t> decode_run(gbwt::size_type decc);
        gbwt::size_type encode_run_length(size_t offset, bool is_rev, uint16_t length, int64_t node_id);

        vector<pair<pos_t, uint16_t>> get_tag_runs(){
            return tag_runs;
        };

        void compressed_serialize(std::ostream &main_out, std::ostream &encoded_starts_file, std::ostream &bwt_intervals_file, std::vector<std::pair<pos_t, uint16_t>> &tag_runs);
        void merge_compressed_files(const std::string filename, const std::string encoded_starts_file, const std::string bwt_intervals_file);
        void load_compressed_tags(std::istream &in);
        void query_compressed(size_t start, size_t end, size_t &number_of_runs);

        // Statistics helpers for compressed format
        size_t number_of_runs_compressed() const;
        size_t bytes_encoded_runs() const;
        size_t bytes_encoded_runs_starts_sd() const;
        size_t bytes_bwt_intervals() const;
        size_t bytes_total_compressed() const;

    private:

        //
        std::vector<gbwt::byte_type> encoded_runs;


        sdsl::bit_vector encoded_runs_starts;


        sdsl::sd_vector<> encoded_runs_starts_sd;
        sdsl::sd_vector<>::select_1_type encoded_runs_sd_starts_select;

        sdsl::sd_vector<> bwt_intervals;

        sdsl::bit_vector::select_1_type encoded_runs_starts_select;
        sdsl::sd_vector<>::rank_1_type bwt_intervals_rank;

        vector<pair<pos_t, uint16_t>> tag_runs;



        // Variables for encoding/decoding
        const static int length_bits = 9;


        // Variables for the serialization
        size_t encoded_start_every_k_run = 10;
        size_t remaining_run_to_write_start = 0;
        size_t cumulative_starts = 0;
        size_t encoded_start_ones = 0;
        size_t start_pos = 0;
        size_t cumulative_run_bwt_position = 0;


    };

}

#endif //VG_TAG_ARRAYS_HPP