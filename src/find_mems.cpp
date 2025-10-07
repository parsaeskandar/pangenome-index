#include "pangenome_index/algorithm.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include <chrono>

using namespace std;
using namespace gbwtgraph;
using namespace panindexer;

#ifndef TIME
#define TIME 1
#endif

int main(int argc, char **argv) {
    string r_index_file = argv[1];
    string tag_array_index = argv[2];
    string reads_file = argv[3];
    size_t mem_length = std::stoi(argv[4]);
    size_t min_occ = std::stoi(argv[5]);

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
    double total_mem_time = 0.0;
    double total_tag_time = 0.0;
#endif

    cerr << "Reading the rindex file (encoded)" << endl;
    FastLocate r_index;
    {
        std::ifstream rin(r_index_file, std::ios::binary);
        if (!rin) { std::cerr << "Cannot open r-index: " << r_index_file << std::endl; std::exit(EXIT_FAILURE); }
        r_index.load_encoded(rin);
    }



//    std::cerr << "sym map " << (int) r_index.sym_map[NENDMARKER] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['A'] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['C'] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['G'] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['T'] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['N'] << endl;
//
//    // r_index.initialize_complement_table();
//
//     FastLocate::bi_interval bint = {0, 0, r_index.bwt_size()};
//     auto forw = r_index.forward_extend(bint, 'A');
//     forw = r_index.forward_extend(forw, 'C');
//     forw = r_index.forward_extend(forw, 'A');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//
//
//     // print forw
//     std::cerr << "forward: " << forw.forward << " size: " << forw.size << " reverse: " << forw.reverse << std::endl;
//     auto back = r_index.backward_extend(bint, 'G');
//
//     back = r_index.backward_extend(back, 'A');
//     back = r_index.backward_extend(back, 'C');
//     back = r_index.backward_extend(back, 'A');
//     std::cerr << "backward: " << back.forward << " size: " << back.size << " reverse: " << back.reverse << std::endl;
//     exit(0);


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

    // Open reads file
    std::ifstream reads(reads_file);
    if (!reads) {
        std::cerr << "Cannot open reads file: " << reads_file << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string read;
    int i = 0;
    while (std::getline(reads, read)) {
        if (read.empty()) continue;
        i++;

#if TIME
        auto time4 = chrono::high_resolution_clock::now();
#endif

        // Find MEMs for this read
        auto mems = find_all_mems(read, mem_length, min_occ, r_index);

#if TIME
        auto time5 = chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration3 = time5 - time4;
        total_mem_time += duration3.count();
//        std::cerr << "Finding MEMs took " << duration3.count() << " seconds" << std::endl;
#endif

        // Print results for this read
        std::cout << "Seq: " << i << std::endl;
        for (const auto& mem : mems) {
            // std::cout << "MEM START: " << mem.start << ", MEM END: " << mem.end + 1 << " OCC " << mem.bi_interval.size << std::endl;
            std::cout << "MEM START: " << mem.start << ", MEM END: " << mem.end << " BWT START: " << mem.bwt_start << " SIZE: " << mem.size << std::endl;
//            std::cout << "BWT interval start: " << mem.bi_interval.forward << ", size: " << mem.bi_interval.size << std::endl;

            size_t tag_nums = 0;

#if TIME
            auto time6 = chrono::high_resolution_clock::now();
#endif


            // Query the tag array for this MEM (compact)
            tag_array.query_compressed_compact(mem.bwt_start, mem.bwt_start + mem.size - 1, tag_nums);

#if TIME
            auto time7 = chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration4 = time7 - time6;
            total_tag_time += duration4.count();
//            std::cerr << "Querying tag array took " << duration4.count() << " seconds" << std::endl;
#endif
        }
        std::cout << std::endl;
    }

    reads.close();

#if TIME
    std::cout << "\nTotal time for finding all MEMs: " << total_mem_time << " seconds" << std::endl;
    std::cout << "Total time for all tag queries: " << total_tag_time << " seconds" << std::endl;
#endif
}

