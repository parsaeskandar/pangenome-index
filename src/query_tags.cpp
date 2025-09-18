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
    {
        std::ifstream rin(r_index_file, std::ios::binary);
        if (!rin) { std::cerr << "Cannot open r-index: " << r_index_file << std::endl; return 1; }
        r_index.load_encoded(rin);
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
        gbwt::range_type range;
        if (r_index.is_encoded()) {
            range = {0, r_index.bwt_size() - 1};
            std::cerr << "[encoded count] read_index=" << i << " len=" << read.size() << std::endl;
            for (size_t pos = read.size(); pos > 0; --pos) {
                char c = read[pos - 1];
                gbwt::range_type prev = range;
                range = r_index.count_encoded(read);
                // range = r_index.LF_encoded(range, static_cast<size_t>(c));
                // std::cerr << "  step=" << (read.size() - pos + 1)
                //           << " char='" << c << "' prev=[" << prev.first << "," << prev.second
                //           << "] new=[" << range.first << "," << range.second << "]" << std::endl;
                if (range.first > range.second) break;
            }
            std::cerr << "  final range=[" << range.first << "," << range.second << "]\n";
        } else {
            range = r_index.count(read);
        }
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





