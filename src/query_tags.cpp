//
// Created by seeskand on 9/18/24.
//

//
// Created by seeskand on 9/13/24.
//

#include "../r-index/internal/r_index.hpp"
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


#define TIME 1

using namespace std;
using namespace ri;
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
    if (argc != 3) {
        std::cerr << "usage: ... " << std::endl;
        exit(0);
    }
    std::string r_index = std::string(argv[1]);
    std::string tag_array_index = std::string(argv[2]);
    int threads = 8;
    size_t k = 31;


    cerr << "Reading the rindex file" << endl;
    FastLocate idx(r_index);
//    std::ifstream in(index_file);
//    bool fast;
//    //fast or small index?
//    in.read((char *) &fast, sizeof(fast));
//    r_index<> idx;
//    idx.load(in);



    cerr << "Reading the tag array index" << endl;
    TagArray tag_array;
    std::ifstream in_ds(tag_array_index);

//    tag_array.load(in_ds);
//    tag_array.query(3, 40);
    size_t next_block_start = 0;
    for (int i = 0; i < 4; ++i) {
        tag_array.load_block_at(in_ds, next_block_start);
    }
//    tag_array.load_block_at(in_ds, next_block_start);

//    cerr << "Finished loading the tag array index" << endl;



//    std::vector<std::string> sequences = readSequencesFromFile("/Users/seeskand/CLionProjects/pangenome-index/test_data/x.giraffe");
////        std::vector<std::string> sequences = readSequencesFromFile("/Users/seeskand/Documents/pangenome-index/test_data/1mb.chr20");
//    if (sequences.empty()) {
//        std::cerr << "No sequences found in file." << std::endl;
//        return 1;
//    }
//
//    cerr << "Finished reading the test sequences" << endl;
//
//    std::vector<std::string> kmers;
//    std::default_random_engine generator(static_cast<long unsigned int>(std::time(0)));
//    std::uniform_int_distribution<int> sequenceDist(0, sequences.size() - 1);
//
//    for (int i = 0; i < 1000; ++i) {
//        const std::string &sequence = sequences[sequenceDist(generator)];
//        if (sequence.size() < k) {
//            continue; // Skip this sequence if it's shorter than k
//        }
//        std::uniform_int_distribution<int> positionDist(0, sequence.size() - k);
//        int pos = positionDist(generator);
//        kmers.push_back(sequence.substr(pos, k));
//    }
//
//
//#if TIME
//    auto time10 = chrono::high_resolution_clock::now();
//#endif
//
//    cerr << "sample kmers created" << endl;
//
//    for (string kmer: kmers) {
//
//
//
//        //find the interval in the bwt
//        auto range = idx.count(kmer);
////        cout << range.first << " " << range.second << endl;
//            cerr << kmer << " " << range.first << " " << range.second << endl;
//
//
//
//        if (range.first <= range.second) {
////            cerr << "The kmer is: " << kmer << endl;
////            cout << "bwt range is " << range.first << " " << range.second << endl;
//            tag_array.query(range.first, range.second);
//        }
//
//    }
//
//    cerr << "All queries done" << endl;
//#if TIME
//    auto time11 = chrono::high_resolution_clock::now();
//    std::chrono::duration<double> duration10 = time11 - time10;
//    std::cerr << "Querying " << kmers.size() << " kmers with size " << k << " took " << duration10.count() << " seconds" << std::endl;
//#endif







}





