//
// Created by seeskand on 9/18/24.
//

//
// Created by seeskand on 9/13/24.
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
//    if (argc != 3) {
//        std::cerr << "usage: ... " << std::endl;
//        exit(0);
//    }
    std::string r_index_file = std::string(argv[1]);
    std::string tag_array_index = std::string(argv[2]);
    int k = std::stoi(argv[3]);
    int threads = 8;

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
    std::ifstream in_ds(tag_array_index);
    tag_array.load_compressed_tags(in_ds);

#if TIME
    auto time3 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = time3 - time2;
    std::cerr << "Loading tag arrays took " << duration2.count() << " seconds" << std::endl;
#endif


    std::vector<std::string> sequences = readSequencesFromFile("/Users/seeskand/CLionProjects/pangenome-index/test_data/two_contig_graph/contigs_XY.");
    if (sequences.empty()) {
        std::cerr << "No sequences found in file." << std::endl;
        return 1;
    }
//
    cerr << "Finished reading the test sequences" << endl;
//
    std::vector<std::string> kmers;
    std::default_random_engine generator(static_cast<long unsigned int>(std::time(0)));
    std::uniform_int_distribution<int> sequenceDist(0, sequences.size() - 1);

    for (int i = 0; i < 100000; ++i) {
        const std::string &sequence = sequences[sequenceDist(generator)];
        if (sequence.size() < k) {
            continue; // Skip this sequence if it's shorter than k
        }
        std::uniform_int_distribution<int> positionDist(0, sequence.size() - k);
        int pos = positionDist(generator);
        kmers.push_back(sequence.substr(pos, k));
    }


#if TIME
    auto time10 = chrono::high_resolution_clock::now();
#endif

    cerr << "sample kmers created" << endl;

    for (string kmer: kmers) {



        //find the interval in the bwt
//        auto range = idx.count(kmer);
        auto range = r_index.count(kmer);
//        cout << range.first << " " << range.second << endl;
//        cerr << kmer << " " << range.first << " " << range.second << endl;
        tag_array.query_compressed(range.first, range.second);



//        if (range.first <= range.second) {
////            cerr << "The kmer is: " << kmer << endl;
////            cout << "bwt range is " << range.first << " " << range.second << endl;
//            tag_array.query(range.first, range.second);
//        }

    }

    cerr << "All queries done" << endl;
#if TIME
    auto time11 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration10 = time11 - time10;
    std::cerr << "Querying " << kmers.size() << " kmers with size " << k << " took " << duration10.count() << " seconds" << std::endl;
#endif







}





