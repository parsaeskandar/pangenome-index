//
// Created by seeskand on 12/9/24.
//

#include "pangenome_index/r-index.hpp"
#include <iostream>


using namespace panindexer;

int main(int argc, char **argv) {

    std::string rlbwt_file = std::string(argv[1]);


    int threads = 8;

    FastLocate idx(rlbwt_file);
    // idx.serialize(std::cout);
    idx.serialize_encoded(std::cout);

}