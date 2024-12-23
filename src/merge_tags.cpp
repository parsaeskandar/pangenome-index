//
// Created by seeskand on 11/14/24.
//

/*
 * This file is used to merge the tags from different chromosomes
 * This file use the whole genome r-index and merge the tags from different chromosomes into one whole genome tag array index
 * The inputs are the whole-genome r-index and the folder containing the chri.rl_bwt
 */


#include "pangenome_index/r-index.hpp"
#include "pangenome_index/algorithm.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include <gbwtgraph/utils.h>
#include <iostream>
#include <filesystem>
#include <omp.h>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <utility>
#include <stdexcept>


// if not define Time define it
#ifndef TIME
#define TIME 1
#endif



//using namespace gbwtgraph;
using namespace panindexer;
using namespace std;
using namespace gbwtgraph;

namespace fs = std::filesystem;



class FileReader {
public:
    FileReader(const std::vector<std::string>& files, size_t n_threads)
            : files(files), n_threads(n_threads), current_thread_id(0) {
        if (files.empty()) {
            throw std::invalid_argument("File list cannot be empty.");
        }

        // Initialize per-file mutexes
        for (size_t i = 0; i < files.size(); ++i) {
            mutexes.emplace_back(std::make_unique<std::mutex>());
        }

        // Initialize file data
        file_positions.resize(files.size(), 0);
        current_nodes.resize(files.size());
        remaining_lengths.resize(files.size(), 0);

        initializeFiles();
    }

    pos_t getNextBlock(size_t thread_id, size_t fileIndex) {
        waitForTurn(thread_id); // Ensure the thread waits for its turn

        pos_t currentTag;

        {
            // Lock file access using the correct mutex for the file
            std::lock_guard<std::mutex> lock(*mutexes[fileIndex]);

            // Check if there are remaining tags in the current block
            if (remaining_lengths[fileIndex] > 0) {
                currentTag = current_nodes[fileIndex];
                remaining_lengths[fileIndex]--;

                // If the block is exhausted, load the next block
                if (remaining_lengths[fileIndex] == 0) {
                    std::ifstream in(files[fileIndex]);
                    in.seekg(file_positions[fileIndex]);

                    auto nextBlock = panindexer::TagArray::load_block_at(in, file_positions[fileIndex]);

                    if (nextBlock.second == 0) {
                        std::cerr << "Empty block encountered in file:" << files[fileIndex] << std::endl;
                    }

                    current_nodes[fileIndex] = nextBlock.first;
                    remaining_lengths[fileIndex] = nextBlock.second;
                }
            } else {
                // Load the next block if no remaining length
                std::ifstream in(files[fileIndex]);
                in.seekg(file_positions[fileIndex]);

                auto nextBlock = panindexer::TagArray::load_block_at(in, file_positions[fileIndex]);

                if (nextBlock.second == 0) {
                    std::cerr << "Empty block encountered in file:" << files[fileIndex] << std::endl;
                }

                current_nodes[fileIndex] = nextBlock.first;
                remaining_lengths[fileIndex] = nextBlock.second - 1;
                currentTag = current_nodes[fileIndex];
            }
        }

        {
            // Output current processing details
            std::lock_guard<std::mutex> lock(output_mutex);
//            std::cerr << "Thread " << thread_id << " File " << fileIndex
//                      << " Block " << currentTag
//                      << " Remaining Length " << int(remaining_lengths[fileIndex])
//                      << std::endl;
        }

        notifyNext(thread_id); // Notify the next thread

        return currentTag;
    }

private:
    void initializeFiles() {
        for (size_t i = 0; i < files.size(); i++) {
            std::ifstream in(files[i]);
            if (!in.is_open()) {
                throw std::runtime_error("Cannot open file: " + files[i]);
            }
            auto firstBlock = panindexer::TagArray::load_block_at(in, file_positions[i]);
            current_nodes[i] = firstBlock.first;
            remaining_lengths[i] = firstBlock.second;
        }
    }

    void waitForTurn(size_t thread_id) {
        std::unique_lock<std::mutex> lock(thread_mutex);
        thread_cv.wait(lock, [this, thread_id] {
            return thread_id == current_thread_id;
        });
    }

    void notifyNext(size_t thread_id) {
        std::lock_guard<std::mutex> lock(thread_mutex);

        // Advance to the next thread
        current_thread_id++;
        if (current_thread_id == n_threads) {
            current_thread_id = 0; // Reset to thread 0
        }

        thread_cv.notify_all(); // Notify waiting threads
    }

    // Data Members
    std::vector<std::string> files;          // List of file paths
    size_t n_threads;                        // Number of threads
    size_t current_thread_id;                // Current thread ID to execute
    std::vector<std::unique_ptr<std::mutex>> mutexes; // Mutex for each file
    std::mutex thread_mutex;                 // Mutex for thread coordination
    std::condition_variable thread_cv;       // Condition variable for thread coordination
    std::mutex output_mutex;                 // Mutex for synchronized output

    // File Handling
    std::vector<size_t> file_positions;      // Current file positions
    std::vector<pos_t> current_nodes;        // Current nodes being processed
    std::vector<uint8_t> remaining_lengths;  // Remaining lengths of current blocks
};





void extract_tag(const FastLocate &r_index, FileReader& reader, size_t thread_id, int& current, bool first, std::vector<int> comp_to_file, unordered_map<size_t, size_t> seq_id_to_comp_id, pos_t& buffer, int thread_num){
    // current is the previous position which we have to run locate_next_nth on, this is when current is none zero
    // current equal to zero means this is the first time we are calling this thread_id
    if (first){

        current = r_index.locateFirst();
        for (size_t i = 0; i < r_index.tot_strings() + thread_id; i++){
            current = r_index.locateNext(current);
        }

    } else {
        // get the number of threads we set using omp_set_num_threads(threads);
//        int threads = omp_get_num_threads();
//        std::cerr << thread_num << std::endl;
        current = r_index.locate_next_nth(current, thread_num);
//        if (current == 0) std::cerr << "current" << current << std::endl;
    }


    auto seq_id = r_index.seqId(current);

    // want to get the file number that is associated with the seq id
    auto current_file = comp_to_file[seq_id_to_comp_id[seq_id]];


    auto block = reader.getNextBlock(thread_id, current_file);
    buffer = block;
    // print the block
//    std::cerr << "Block " << block << "from thread " << thread_id << std::endl;
}




std::vector<std::string> get_files_in_dir(const std::string& directoryPath) {
    std::vector<std::string> files;

    // Check if the path is a valid directory
    if (!fs::is_directory(directoryPath)) {
        std::cerr << "Path is not a directory: " << directoryPath << std::endl;
        return files;
    }

    // Iterate over entries in the directory
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (fs::is_regular_file(entry.status())) { // Only include regular files
            files.push_back(entry.path().string()); // Add file path as a string
        }
    }

    return files;
}



int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "usage: ... " << std::endl;
        exit(0);
    }

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
#endif





    int threads = 8;
    omp_set_num_threads(threads);

    std::string gbz_graph = std::string(argv[1]);
    std::string r_index_file = std::string(argv[2]);
    std::string tag_array_index_dir = std::string(argv[3]);




    GBZ gbz;
    cerr << "Loading the graph file" << endl;
    sdsl::simple_sds::load_from(gbz, gbz_graph);

    cerr << "Reading the whole genome r-index file" << endl;
    FastLocate r_index;
//    std::ifstream r_index(r_index_file);
//    idx.load(in);
    if(!sdsl::load_from_file(r_index, r_index_file))
    {
        std::cerr << "Cannot load the r-index from " << r_index_file << std::endl;
        std::exit(EXIT_FAILURE);
    }




    // get the list of files in the directory
    std::vector<std::string> files = get_files_in_dir(tag_array_index_dir);

    int number_of_file = files.size();


    // finding the components
    auto node_to_comp_map = node_to_component(gbz);



    // print all node_to_comp_map
//    for (auto it = node_to_comp_map.begin(); it != node_to_comp_map.end(); ++it){
//        std::cerr << it->first << " " << it->second << std::endl;
//    }


    std::vector<int> file_to_comp(number_of_file);
    std::vector<int> comp_to_file(number_of_file);

    std::vector<pos_t> current_node_of_files(number_of_file);
    std::vector<uint8_t> current_remaining_length_of_files(number_of_file);
    std::vector<size_t> current_file_positions(number_of_file, 0);



    // for each tag block files, we read the first block and read the first node
    for (auto i=0; i<number_of_file; i++){
        std::ifstream in(files[i]);
        std::pair<pos_t, uint8_t> first_block = panindexer::TagArray::load_block_at(in, current_file_positions[i]);
        current_node_of_files[i] = first_block.first;
        current_remaining_length_of_files[i] = first_block.second;

        // get the component of the node of the first block
        size_t comp = node_to_comp_map[id(current_node_of_files[i])];
        file_to_comp[i] = comp;
        comp_to_file[comp] = i;
    }


    // print all file_to_comp
//    for (size_t i = 0; i < number_of_file; i++){
//        std::cerr << "comp " << i << " file " << comp_to_file[i] << std::endl;
//    }

    unordered_map<size_t, size_t> seq_id_to_comp_id;

    // get the first node of each path and get the component id of the node
    for (size_t i = 0; i < r_index.tot_strings(); i++){
        auto seq_graph_nodes = gbz.index.extract(i * 2);

        // extracting the first node of the sequence component id
        seq_id_to_comp_id[i] = node_to_comp_map[gbwt::Node::id(seq_graph_nodes[0])];
    }


    // reading the r-index file and traversing the r-index and getting the sequence number of each BWT position

    auto current = r_index.locateFirst();
    // have to bypass the first #seqs because there are no tags for them
    for (size_t i = 0; i < r_index.tot_strings(); i++){
//        cerr << "The first node is " << current << endl;
        current = r_index.locateNext(current);
    }

//    cerr << "The first nodes " << current << endl;
    pos_t current_tag;


    // note that the first #num_seq tags are correspond to the ENDMARKERs
    auto total_tags_count = r_index.get_sequence_size() - r_index.tot_strings();
    vector<pos_t> tags(total_tags_count);
    vector<pair<pos_t, uint8_t>> tag_runs;
    size_t current_run = 0;


    // Single threaded version of extracting the tags for each position of the r-index
//    for (size_t i = r_index.tot_strings() ; i < r_index.get_sequence_size(); i++){
//        // seq_id is the same as GBWT path identifier
//        auto seq_id = r_index.seqId(current);
//
//        // want to get the file number that is associated with the seq id
//        auto current_file = comp_to_file[seq_id_to_comp_id[seq_id]];
//
//        // reading the next tag from the current file
//        if (current_remaining_length_of_files[current_file] == 0){
//            std::cerr << "Error: reaching the end of file too soon " << std::endl;
//        } else {
//            // Just getting the node and decrementing the remaining length
//            current_tag = current_node_of_files[current_file];
//            current_remaining_length_of_files[current_file]--;
//
//            if (current_remaining_length_of_files[current_file] == 0){
//                // read the next block
//                std::ifstream in(files[current_file]);
//                std::pair<pos_t, uint8_t> next_block = panindexer::TagArray::load_block_at(in, current_file_positions[current_file]);
//
//                current_node_of_files[current_file] = next_block.first;
//                current_remaining_length_of_files[current_file] = next_block.second;
//
//                if (next_block.second == 0){
//
//                    std::cerr << "Error: the next block is empty " << seq_id << std::endl;
//                    current_remaining_length_of_files[current_file] = 0;
//                    current_node_of_files[current_file] = next_block.first;
//                }
//
//            }
//
//
//            if (current_run == 0){
//                tag_runs.push_back(make_pair(current_tag, 1));
//                current_run++;
//            } else {
//                if (current_tag == tag_runs[current_run - 1].first){
//                    tag_runs[current_run - 1].second++;
//                } else {
//                    tag_runs.push_back(make_pair(current_tag, 1));
//                    current_run++;
//                }
//            }
//            tags[ i - r_index.tot_strings()] = current_tag;
//        }
//
//        if (i != r_index.get_sequence_size() - 1){
//            current = r_index.locateNext(current);
//        }
//    }






    // print the tags
//    cerr << "tags size " << tags.size() << " r index size " << r_index.get_sequence_size() << endl;

    // print the tag_runs
//    for (size_t i = 0; i < tag_runs.size(); i++){
//        cerr << "tag " << tag_runs[i].first << " length " << int(tag_runs[i].second) << endl;
//    }

    TagArray tag_array;

#if TIME
    auto time2 = chrono::high_resolution_clock::now();
#endif
//    tag_array.serialize_run_by_run(tag_runs, "whole_genome_tag_array.tags");

    FileReader reader(files, threads);
    std::vector<std::thread> threads_list; threads_list.reserve(threads);
    std::vector<pos_t> thread_buffers(threads);
    std::vector<int> current_position(threads, 0);
    int batch_size = 1024;
    vector<pos_t> tags_batch(batch_size);
    vector<pair<pos_t, uint8_t>> tag_runs2;
    vector<pos_t> tags2(total_tags_count);
    current_run = 0;
    vector<pair<pos_t, uint8_t>> tag_runs_check;



    for(size_t to_read = 0; to_read < threads && to_read < total_tags_count; to_read++){
        threads_list.emplace_back(extract_tag, std::ref(r_index), std::ref(reader), to_read, std::ref(current_position[to_read]), true, comp_to_file, seq_id_to_comp_id, std::ref(thread_buffers[to_read]), threads);
    }
//    extract_tag(r_index, reader, 0, 0, comp_to_file, seq_id_to_comp_id, thread_buffers[0]);


    for(size_t to_write = 0; to_write < total_tags_count; to_write++){
        size_t thread_id = to_write % threads;
        threads_list[thread_id].join();
        pos_t current_tag = thread_buffers[thread_id];
        // print the block
//        std::cerr << "Block " << current_tag << "from thread " << thread_id << std::endl;
        tags2[to_write] = current_tag;
        tags_batch[to_write % batch_size] = current_tag;

        if (current_run == 0){
            tag_runs2.push_back(make_pair(current_tag, 1));
            tag_runs_check.push_back(make_pair(current_tag, 1));
            current_run++;
        } else {
            if (current_tag == tag_runs2[current_run - 1].first){
                tag_runs2[current_run - 1].second++;
                tag_runs_check[tag_runs_check.size() - 1].second++;
            } else {
                tag_runs2.push_back(make_pair(current_tag, 1));
                tag_runs_check.push_back(make_pair(current_tag, 1));
                current_run++;
            }
        }

        // whenever the batch is full, we serialize the batch

        if (tag_runs2.size() == batch_size || to_write == total_tags_count - 1){
            // pop the last element
            auto temp = tag_runs2.back();
            if (to_write < total_tags_count - 1){
                tag_runs2.pop_back();
            }
            tag_array.serialize_run_by_run(tag_runs2, "whole_genome_tag_array_parallel.tags");
            tag_runs2.clear();
            current_run = 1;
            tag_runs2.push_back(temp);
        }



        if(to_write + threads < total_tags_count)
        {
            threads_list[thread_id] = std::thread(extract_tag, std::ref(r_index), std::ref(reader), to_write % threads, std::ref(current_position[thread_id]), false, comp_to_file, seq_id_to_comp_id, std::ref(thread_buffers[thread_id]), threads);
        }



    }
    for (auto& thread : threads_list) {
        if (thread.joinable()) {
            thread.join();
        }
    }

#if TIME
    auto time3 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = time3 - time2;
    std::cerr << "Converting tags using multiple threads took " << duration2.count() << " seconds" << std::endl;
#endif

    // see of tags and tags2 are the same
//    for (size_t i = 0; i < total_tags_count; i++){
//        if (tags[i] != tags2[i]){
//            std::cerr << "Error: tags and tags2 are not the same " << i << " " << tags[i] << " " << tags2[i] << std::endl;
//        }
//    }
//
//    // see if the tag_runs and tag_runs_check are the same
//    for (size_t i = 0; i < tag_runs.size(); i++){
//        if (tag_runs[i].first != tag_runs_check[i].first || tag_runs[i].second != tag_runs_check[i].second){
//            std::cerr << "Error: tag_runs and tag_runs_check are not the same " << i << " " << tag_runs[i].first << " " << tag_runs_check[i].first << " " << tag_runs[i].second << " " << tag_runs_check[i].second << std::endl;
//        }
//    }

//    tag_array.serialize_run_by_run(tag_runs_check, "whole_genome_tag_array_parallel.tags");
//    tag_array.deserialize_run_by_run("whole_genome_tag_array_parallel.tags");
//    TagArray tagArray2;
//    tagArray2.deserialize_run_by_run("whole_genome_tag_array.tags");
//
//    //check if the two tag arrays are the same by comparing their encoded runs
//    auto en1 = tag_array.get_tag_runs();
//    auto en2 = tagArray2.get_tag_runs();
//
//    std::cerr << en1.size() << " " << en2.size() << std::endl;
//    // print the two encoded runs
//    for (size_t i = 0; i < en1.size(); i++){
//        std::cerr << int(en1[i].second) << " " << int(en2[i].second) << std::endl;
//        if (en1[i].first != en2[i].first){
//            std::cerr << "Error: tag arrays are not the same " << i << " " << int(en1[i].second) << " " << int(en2[i].second) << std::endl;
//        }
//    }




    return 0;





}

