#include "pangenome_index/memory_mapped_index.hpp"
#include "pangenome_index/algorithm.hpp"
#include <chrono>
#include <fstream>

using namespace std;
using namespace gbwtgraph;
using namespace panindexer;

#ifndef TIME
#define TIME 1
#endif

int main(int argc, char **argv) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <r_index_file> <tag_array_file> <reads_file> <mem_length> <min_occ>" << std::endl;
        std::cerr << "Example: " << argv[0] << " whole.ri whole.tags reads.txt 20 5" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    string r_index_file = argv[1];
    string tag_array_file = argv[2];
    string reads_file = argv[3];
    size_t mem_length = std::stoi(argv[4]);
    size_t min_occ = std::stoi(argv[5]);

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
    double total_mem_time = 0.0;
    double total_tag_time = 0.0;
#endif

    std::cerr << "Opening memory-mapped indices..." << std::endl;
    
    // Use memory-mapped index manager instead of loading into RAM
    MemoryMappedIndexManager index_manager;
    
    if (!index_manager.open(r_index_file, tag_array_file)) {
        std::cerr << "Failed to open memory-mapped indices" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Get references to the memory-mapped data structures
    FastLocate& r_index = index_manager.get_r_index();
    TagArray& tag_array = index_manager.get_tag_array();

#if TIME
    auto time2 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = time2 - time1;
    std::cerr << "Opening memory-mapped indices took " << duration1.count() << " seconds" << std::endl;
    std::cerr << "Total memory-mapped size: " << (index_manager.total_size() / (1024.0 * 1024.0 * 1024.0)) << " GB" << std::endl;
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

        // Find MEMs for this read using the memory-mapped r-index
        auto mems = find_all_mems(read, mem_length, min_occ, r_index);

#if TIME
        auto time5 = chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration3 = time5 - time4;
        total_mem_time += duration3.count();
#endif

        // Print results for this read
        std::cout << "Seq: " << i << std::endl;
        for (const auto& mem : mems) {
            std::cout << "MEM START: " << mem.start << ", MEM END: " << mem.end 
                      << " BWT START: " << mem.bwt_start << " SIZE: " << mem.size << std::endl;

            size_t tag_nums = 0;

#if TIME
            auto time6 = chrono::high_resolution_clock::now();
#endif

            // Query the memory-mapped tag array for this MEM
            tag_array.query_compressed(mem.bwt_start, mem.bwt_start + mem.size - 1, tag_nums);

#if TIME
            auto time7 = chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration4 = time7 - time6;
            total_tag_time += duration4.count();
#endif
        }
        std::cout << std::endl;
    }

    reads.close();

#if TIME
    std::cout << "\nTotal time for finding all MEMs: " << total_mem_time << " seconds" << std::endl;
    std::cout << "Total time for all tag queries: " << total_tag_time << " seconds" << std::endl;
    std::cout << "Memory-mapped indices used: " << (index_manager.total_size() / (1024.0 * 1024.0 * 1024.0)) << " GB" << std::endl;
#endif

    // Clean up (memory mapping is automatically cleaned up by destructor)
    std::cerr << "Cleaning up memory-mapped indices..." << std::endl;
}
