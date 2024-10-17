//
// Created by seeskand on 9/13/24.
//

#ifndef PANGENOME_INDEX_TAG_ARRAYS_HPP
#define PANGENOME_INDEX_TAG_ARRAYS_HPP

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
//#include <vg/io/vpkg.hpp>


using namespace std;
using handlegraph::pos_t;

namespace panindexer {


    class TagArray {
    public:
        TagArray();

//    void TagArray::add_run(size_t offset, bool is_reverse, int length, size_t node_id);

        void load_bptree(BplusTree<Run> &bptree, size_t bwt_size);

        void serialize(std::ostream &out);

        void load(std::istream &in);

        void query(size_t start, size_t end);

    private:

        std::vector <gbwt::byte_type> encoded_runs;
        sdsl::bit_vector encoded_runs_starts;
        sdsl::sd_vector<> bwt_intervals;

        sdsl::bit_vector::select_1_type encoded_runs_starts_select;
        sdsl::sd_vector<>::rank_1_type bwt_intervals_rank;
    };

}
#endif //PANGENOME_INDEX_TAG_ARRAYS_HPP
