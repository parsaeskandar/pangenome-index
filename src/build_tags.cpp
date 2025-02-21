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
#include "../deps/grlBWT/scripts/fm_index.h"
#include "../deps/grlBWT/include/bwt_io.h"
#include "pangenome_index/r-index.hpp"

#ifndef TIME
#define TIME 1
#endif

using namespace std;
using namespace ri;
using namespace panindexer;


int main(int argc, char **argv) {
//    if (argc != 3) {
//        std::cout << "usage: ... " << std::endl;
//        exit(0);
//    }
    std::string graph_file = std::string(argv[1]);
//    std::string index_file = std::string(argv[2]);
    std::string rlbwt_file = std::string(argv[2]);
    std::string output_file = std::string(argv[3]);

    int threads = 8;
//    omp_set_num_threads(threads);
    size_t k = 31;


//    cerr << "Reading the rindex file" << endl;
//    std::ifstream in(index_file);
//    bool fast;
//    //fast or small index?
//    in.read((char *) &fast, sizeof(fast));
//    r_index<> idx;
//    idx.load(in);


    FastLocate idx(rlbwt_file);
//    idx.serialize(std::cout);
//    exit(0);
//    FastLocate idx;
//    std::ifstream in(rlbwt_file);
//    idx.load(in);

    // print all blocks
//    for (size_t i = 0; i < idx.blocks.size(); i++) {
//        idx.blocks[i].print();
//    }

//    cerr << "Finished reading the rindex file" << endl;




//    gbwt::range_type run(0, 13);
//    size_t run_id = 0;
//    bool starts_with_to = false;
//    size_t first_run = 0;
////    size_t symb = 'G';
//
////    auto a1 = idx.LF({0, 13}, 'A');
//    cerr << '1' << endl;
//    gbwt::range_type a = idx.LF({0, 13}, 'A');
//    gbwt::range_type b = idx.LF(a, 'T');
//
//
//    std::cerr << "LFFFF" << a << std::endl;
//    std::cerr << "LFFFF" << b << std::endl;
//
//
//    for (size_t i = 0; i < 14; i++) {
//        std::cerr << "rankat " << idx.rankAt(i, 'A') << std::endl;
//    }
//
//    auto test = idx.locate({10, 50});
//    std::cerr << "LOCATE " << std::endl;
//    for(auto i: test){
//        std::cerr << i << std::endl;
//    }

//    idx.serialize(std::cout);



//    exit(0);
    GBZ gbz;
    cerr << "Loading the graph file" << endl;
    sdsl::simple_sds::load_from(gbz, graph_file);

//    auto input = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, gbwtgraph::GBWTGraph, HandleGraph>(graph_file);
//    gbz = std::move(get<0>(input));


    typedef gbwtgraph::Key64::value_type kmer_type;

    hash_map<kmer_type, gbwtgraph::Position> index;
    cerr << "Computing the unique kmers in the graph" << endl;

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
#endif

    unique_kmers_parallel<gbwtgraph::Key64>(gbz.graph, index, k);

#if TIME
    auto time2 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = time2 - time1;
    std::cerr << "Indexing unique kmers took " << duration1.count() << " seconds" << std::endl;
#endif


    BplusTree<panindexer::Run> bptree(15); // TODO: determine the BPlusTree degree

//
//
//
    cerr << "Adding the kmers to the BPlusTree" << endl;

    gbwt::range_type whole_bwt(0, idx.bwt_size() - 1);
    parallel_kmers_to_bplustree(idx, bptree, index, k, whole_bwt);

//

#if TIME
    auto time4 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration3 = time4 - time2;
    std::cerr << "Adding kmers to bptree took " << duration3.count() << " seconds" << std::endl;
#endif


    // computing some statistics
    auto unique_kmers_size = index.size();
    auto bptree_items = bptree.get_bpt_size();
    auto bwt_size = idx.bwt_size();
    size_t tag_arrays_covered = 0;

    cerr << "The number of unique kmers in the index is: " << unique_kmers_size << endl;
    cerr << "The number of items in the BPlusTree is: " << bptree_items << endl;
    cerr << "The size of the BWT is: " << bwt_size << endl;


    cerr << "calculating the fraction of the tag arrays covered" << endl;
    // calculating the fraction of the tag arrays covered
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        panindexer::Run current_item = *it;
        if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
            auto next_it = it;
            ++next_it; // Move to the next element
            if (next_it != bptree.end()) { // Check if the next element is not the end
                panindexer::Run next_item = *next_it; // Get the next item
                tag_arrays_covered += (next_item.start_position - current_item.start_position);
            }
        }
    }


    cerr << "The fraction of the tag arrays covered by unique kmers is: " << tag_arrays_covered << " / " << bwt_size
         << " = "
         << (double) tag_arrays_covered / bwt_size << endl;


#if TIME
    auto time5 = chrono::high_resolution_clock::now();
#endif
    cerr << "Extending the kmers on the graph" << endl;
    auto extension_candidates = extend_kmers_bfs_parallel(gbz.graph, idx, bptree, 1024);

#if TIME
    auto time6 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration4 = time6 - time5;
    std::cerr << "Extending kmers took " << duration4.count() << " seconds" << std::endl;
#endif


//    cerr << "The number of items in the BPlusTree is: " << bptree.get_bpt_size() << endl;
    tag_arrays_covered = 0;

    cerr << "calculating the fraction of the tag arrays covered after one extension" << endl;
    // calculating the fraction of the tag arrays covered
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        panindexer::Run current_item = *it;
        if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
            auto next_it = it;
            ++next_it; // Move to the next element
            if (next_it != bptree.end()) { // Check if the next element is not the end
                panindexer::Run next_item = *next_it; // Get the next item
                tag_arrays_covered += (next_item.start_position - current_item.start_position);
//                cerr << "The current item is: " << current_item << " NEXT " << next_item << endl;
//                if (next_item.start_position - current_item.start_position > 500){
//                    cerr << "The current item is: " << current_item << " The next item is: " << next_item << endl;
//                }
            }
        }
    }

    cerr << "The fraction of the tag arrays covered after extending the kmers is: " << tag_arrays_covered << " / "
         << bwt_size << " = "
         << (double) tag_arrays_covered / bwt_size << endl;


    string pattern = "$";


    // Want to find the end of each sequence in the text
//    auto OCC = idx.ISA(pattern);
    auto OCC = idx.OCC();

    //print all elements of the OCC array
//    for (size_t i = 0; i < OCC.size(); i++){
//        cerr << "OCC[" << i << "] = " << OCC[i] << endl;
//    }


    // using the sort_end_of_seq function
    // the first item in the nth element of this is the end of the nth sequence in the bwt
    auto end_of_seq = sort_end_of_seq(OCC);

#if TIME
    auto time7 = chrono::high_resolution_clock::now();
#endif

    traverse_sequences_parallel(gbz, bptree, idx, end_of_seq);
#if TIME
    auto time8 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration5 = time8 - time7;
    std::cerr << "Traversing all paths and fill all the gaps took " << duration5.count() << " seconds" << std::endl;
#endif


    cerr << "The final number of items in the BPlusTree is: " << bptree.get_bpt_size() << endl;
    tag_arrays_covered = 0;

    cerr << "calculating the fraction of the tag arrays covered " << endl;
    // calculating the fraction of the tag arrays covered
    int count = 0;
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        panindexer::Run current_item = *it;
        if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
            auto next_it = it;
            ++next_it; // Move to the next element
            if (next_it != bptree.end()) { // Check if the next element is not the end
                panindexer::Run next_item = *next_it; // Get the next item
                tag_arrays_covered += (next_item.start_position - current_item.start_position);
                count++;
                // print the decoded items of current_item and next item
                if (count < 5) {

                    pos_t t1 = current_item.graph_position.decode();
                    cerr << "The current item is: " << current_item << endl;
                    cerr << "the graph position is " << t1 << endl;
                    cerr << "node id: " << id(t1) << " offset " << offset(t1) << " rev? " << is_rev(t1) << endl;
                    cerr << "The current start position is " << current_item.start_position << endl;
                }
//                if (next_item.start_position - current_item.start_position > 500){
//                    cerr << "The current item is: " << current_item << " The next item is: " << next_item << endl;
//                }
            }
        } else {
            std::cerr << "================ERROR================ " << count << std::endl;
            // this case should not happen other than when we are at the end
            pos_t t1 = current_item.graph_position.decode();
            cerr << "The current item is: " << current_item << endl;
            cerr << "the graph position is " << t1 << endl;
            cerr << "node id: " << id(t1) << " offset " << offset(t1) << " rev? " << is_rev(t1) << endl;
            cerr << "The current start position is " << current_item.start_position << endl;


            auto next_it = it;
            ++next_it;

            if (next_it != bptree.end()){

                panindexer::Run next_item = *next_it;

                pos_t t2 = next_item.graph_position.decode();
                cerr << "The current item is: " << next_item << endl;
                cerr << "the graph position is " << t2 << endl;
                cerr << "node id: " << id(t2) << " offset " << offset(t2) << " rev? " << is_rev(t2) << endl;
                cerr << "The current start position is " << next_item.start_position << endl;

                std::cerr << "==========================================================" << std::endl;

            }





        }
    }

    cerr << "The fraction of the tag arrays covered after filling the gaps is: " << tag_arrays_covered << " / "
         << bwt_size << " = "
         << (double) tag_arrays_covered / bwt_size << endl;

    panindexer::TagArray tag_array;
    tag_array.load_bptree(bptree, idx.bwt_size());

    cerr << "The bptree loaded into Tag array DS correctly " << endl;

//    auto node_to_comp = node_to_component(gbz);

    

    // query the tag array and print the first 10 runs
//    tag_array.query(0, 100);

//    tag_array.serialize(std::cout);
//    tag_array.store_blocks(std::cout);
    tag_array.store_blocks_sdsl(output_file);




}





