//
// Created by seeskand on 11/14/24.
//

/*
 * This file is used to merge the tags from different chromosomes
 * This file use the whole genome r-index and merge the tags from different chromosomes into one whole genome tag array index
 * The inputs are the whole-genome r-index and the folder containing the chri.rl_bwt
 */


#include "pangenome_index/r-index.hpp"



// if not define Time define it
#ifndef TIME
#define TIME 1
#endif



using namespace gbwtgraph;
using namespace panindexer;



int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "usage: ... " << std::endl;
        exit(0);
    }


    std::string r_index_file = std::string(argv[1]);
    std::string tag_array_index_dir = std::string(argv[2]);



    cerr << "Reading the whole genome r-index file" << endl;
    FastLocate idx;
    std::ifstream in(r_index_file);
    idx.load(in);


    // get the list of files in the directory
    std::vector<std::string> files = get_files_in_dir(tag_array_index_dir);


    std::vector<>







}
