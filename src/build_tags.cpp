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


#define TIME 1

using namespace std;
using namespace ri;
using namespace panindexer;


int main(int argc, char **argv) {
//    if (argc != 3) {
//        std::cout << "usage: ... " << std::endl;
//        exit(0);
//    }
    std::string graph_file = std::string(argv[1]);
    std::string index_file = std::string(argv[2]);
    std::string rlbwt_file = std::string(argv[3]);

    int threads = 8;
    size_t k = 31;


    cerr << "Reading the rindex file" << endl;
    std::ifstream in(index_file);
    bool fast;
    //fast or small index?
    in.read((char *) &fast, sizeof(fast));
    r_index<> idx;
    idx.load(in);


//    std::cerr<<"Computing the FM-index from the input BCR BWT"<<std::endl;
//    fm_index fmi(rlbwt_file);

    bwt_buff_reader bwt_buff(rlbwt_file);
    cerr << "Runs " << bwt_buff.size() << endl;

    fm_index fmi(rlbwt_file);


    cerr << fmi.lf(0).first << " " << fmi.lf(0).second << endl;
    cerr << fmi.lf(1).first << " " << fmi.lf(1).second << endl;
    cerr << fmi.lf(2).first << " " << fmi.lf(2).second << endl;
    cerr << fmi.lf(3).first << " " << fmi.lf(3).second << endl;

    cerr << fmi.C[0] << " " << fmi.C[1] << " " << fmi.C[2] << " " << fmi.C[3] << endl;
//    cerr << fmi.sym_map['A'] << " " << fmi.sym_map['C'] << " " << fmi.sym_map['G'] << " " << fmi.sym_map['T'] << endl;
    cerr << fmi.bwt[0] << " " << fmi.bwt[1] << " " << fmi.bwt[2] << " " << fmi.bwt[3] << endl;

    // print stuff in fmi.sym_map
    for (size_t i = 0; i < 256; i++){
        if (fmi.sym_map[i] != 0){
            cerr << i << " " << fmi.sym_map[i] << endl;
        }
    }


    size_t symb = 'A';
    cerr << fmi.bwt << endl;


    FastLocate r_index(rlbwt_file);



    // checking the backward navigation psi
    cerr << "PSI " << r_index.psi(0) << endl;
    cerr << "PSI " << r_index.psi(1) << endl;
    cerr << "PSI " << r_index.psi(2) << endl;
    cerr << "PSI " << r_index.psi(3) << endl;
    cerr << "PSI " << r_index.psi(4) << endl;



    cerr << " 5 " << r_index.locateNext(5) << endl;
    cerr << " 6 " << r_index.locateNext(6) << endl;
    cerr << " 7 " << r_index.locateNext(7) << endl;
    cerr << " 8 " << r_index.locateNext(8) << endl;
    cerr << " 9 " << r_index.locateNext(9) << endl;
    cerr << " 10 " << r_index.locateNext(10) << endl;
    cerr << " 11 " << r_index.locateNext(11) << endl;
    cerr << " 12 " << r_index.locateNext(12) << endl;
    cerr << " 13 " << r_index.locateNext(13) << endl;
    cerr << " 14 " << r_index.locateNext(14) << endl;




    gbwt::range_type run(5, 10);
    size_t run_id = 0;
    bool starts_with_to = false;
    size_t first_run = 0;
//    size_t symb = 'G';



    auto temp = r_index.LF(run, symb, starts_with_to, first_run);
    cerr << "TEMP " << temp.first << " " << temp.second << endl;


    std::vector <size_type> res = r_index.locate(run);

    // print all the res
    for (size_t i = 0; i < res.size(); i++){
        cerr << res[i] << " ";
    }


    auto x = r_index.decompressDA();
    cerr << "DA " << x.size() << endl;
    for (size_t i = 0; i < x.size(); i++){
        cerr << x[i] << " ";
    }

    cerr << r_index.tot_strings() << endl;

    exit(0);


    GBZ gbz;
    cerr << "Loading the graph file" << endl;
    sdsl::simple_sds::load_from(gbz, graph_file);

//    auto input = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, gbwtgraph::GBWTGraph, HandleGraph>(graph_file);
//    gbz = std::move(get<0>(input));


    bool progress = true;
    auto threshold = 0;
    auto space_efficient_counting = false;

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
////    kmers_to_bplustree(idx, bptree, index, k, {0, idx.bwt_size() - 1}, "");
    parallel_kmers_to_bplustree(idx, bptree, index, k, {0, idx.bwt_size() - 1});
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


    // get the ISA values for the end of each sequence by searching for pattern # and locating the results
    auto OCC = idx.ISA(pattern);

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
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        panindexer::Run current_item = *it;
        if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
            auto next_it = it;
            ++next_it; // Move to the next element
            if (next_it != bptree.end()) { // Check if the next element is not the end
                panindexer::Run next_item = *next_it; // Get the next item
                tag_arrays_covered += (next_item.start_position - current_item.start_position);
                // print the decoded items of current_item and next item
                pos_t t1 = current_item.graph_position.decode();
//                cerr << "The current item is: " << current_item << endl;
//                cerr << "node id: " << vg::id(t1) << " offset " << vg::offset(t1) << " rev? " << vg::is_rev(t1) << endl;

//                if (next_item.start_position - current_item.start_position > 500){
//                    cerr << "The current item is: " << current_item << " The next item is: " << next_item << endl;
//                }
            }
        }
    }

    cerr << "The fraction of the tag arrays covered after filling the gaps is: " << tag_arrays_covered << " / "
         << bwt_size << " = "
         << (double) tag_arrays_covered / bwt_size << endl;

    TagArray tag_array;
    tag_array.load_bptree(bptree, idx.bwt_size());

    tag_array.serialize(std::cout);





}





