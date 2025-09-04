//
// Created by seeskand on 9/18/24.
//
#include "pangenome_index/bplus_tree.hpp"
#include "pangenome_index/algorithm.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include "pangenome_index/unique_kmer.hpp"
#include <getopt.h>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <queue>
#include <mutex>
#include <condition_variable>

#ifndef TIME
#define TIME 1
#endif
using namespace std;
using namespace gbwtgraph;
using namespace panindexer;


std::vector<std::string> readSequencesFromFile(const std::string &filename) {
    std::ifstream file(filename);
    std::vector<std::string> sequences;
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            sequences.push_back(line);
        }
    }
    return sequences;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage: ./bin/query_tags <r_index.ri> <compressed_tags.tags> <reads.txt>" << std::endl;
        return 1;
    }
    std::string r_index_file = std::string(argv[1]);
    std::string tag_array_index = std::string(argv[2]);
    std::string reads_file = std::string(argv[3]);

    cerr << "Reading the rindex file" << endl;

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
#endif

    FastLocate r_index;
    if (!sdsl::load_from_file(r_index, r_index_file)) {
        std::cerr << "Cannot load the r-index from " << r_index_file << std::endl;
        std::exit(EXIT_FAILURE);
    }

#if TIME
    auto time2 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = time2 - time1;
    std::cerr << "Loading r-index into memory took " << duration1.count() << " seconds" << std::endl;
#endif


    cerr << "Reading the tag array index" << endl;
    TagArray tag_array;
    std::ifstream in_ds(tag_array_index, std::ios::binary);
    tag_array.load_compressed_tags_compact(in_ds);

#if TIME
    auto time3 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = time3 - time2;
    std::cerr << "Loading tag arrays took " << duration2.count() << " seconds" << std::endl;
#endif

    // Read reads from file
    std::vector<std::string> reads = readSequencesFromFile(reads_file);
    if (reads.empty()) {
        std::cerr << "No reads found in file." << std::endl;
        return 1;
    }

    // For each read, query r-index and tag arrays, then output summary
    for (size_t i = 0; i < reads.size(); ++i) {
        std::string &read = reads[i];
        auto range = r_index.count(read);
        size_t number_of_runs = 0;
        if (range.first > range.second) {
            std::cerr << "Read " << i << " has no matches" << std::endl;
            continue;
        }
        tag_array.query_compressed_compact(range.first, range.second, number_of_runs);
        std::cout << "read_index=" << i
                  << "\tlen=" << read.size()
                  << "\tbwt_start=" << range.first
                  << "\tbwt_end=" << range.second
                  << "\truns=" << number_of_runs
                  << std::endl;
    }

}





