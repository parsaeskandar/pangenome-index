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
using handlegraph::pos_t;


namespace fs = std::filesystem;


class FileReader {
public:


    FileReader(const std::vector <std::string> &files, size_t n_threads, size_t batch_size)
            : files(files), n_threads(n_threads), current_thread_id(0), batch_size(batch_size) {
        if (files.empty()) {
            throw std::invalid_argument("File list cannot be empty.");
        }

        // Initialize per-file mutexes
        for (size_t i = 0; i < files.size(); ++i) {
            mutexes.emplace_back(std::make_unique<std::mutex>());
        }

        // Initialize file data
        file_positions.resize(files.size(), 0);
        file_end_reached.resize(files.size(), false);
        file_buffers.resize(files.size());
        tag_batches.resize(files.size());
        current_index_tag_batch.resize(files.size(), 0);
        current_length_tag_batch.resize(files.size(), 0);
        for (size_t i = 0; i < files.size(); ++i) {
            tag_batches[i].resize(batch_size);
        }

        initializeFiles();

    }

    int get_file_num(){
        return files.size();
    }


    pos_t get_first_tag(size_t fileIndex) {
        return tag_batches[fileIndex][0].first;
    }



    void closeAllFiles() {
        std::cerr << "Closing all files" << std::endl;
        for (size_t i = 0; i < file_buffers.size(); i++) {
            if (file_buffers[i].is_open()) {
                file_buffers[i].close(); // Explicitly close the file buffer
                std::cerr << "Closed file: " << files[i] << std::endl;
            }
        }
        file_buffers.clear(); // Optionally clear the vector to release all resources
        file_positions.clear(); // Clear positions as files are closed
        file_end_reached.clear(); // Clear end-of-file flags
        tag_batches.clear(); // Clear tag buffersStarting creating the request list
    }



    std::vector<std::vector<std::pair<pos_t, uint8_t>>> extract_requested_tags(
            size_t thread_id, const std::vector<size_t>& requests) {

        waitForTurn(thread_id);

        std::vector<std::vector<std::pair<pos_t, uint8_t>>> extracted_tags(files.size());

        for (size_t i = 0; i < files.size(); i++) {
            size_t current_extracted = 0;


            // TODO: the idea is to use a circular array
            while (current_extracted < requests[i]){
                if (current_length_tag_batch[i] == 0) {
                    std::cerr << "No more tags is possible to read from file: " << files[i] << std::endl;
                    break;
                }
                if (current_extracted + tag_batches[i][current_index_tag_batch[i]].second <= requests[i]){
                    // add the whole tag run to the extracted tags and delete it from the tag_batches
                    extracted_tags[i].push_back(tag_batches[i][current_index_tag_batch[i]]);
                    current_extracted += tag_batches[i][current_index_tag_batch[i]].second;
                    current_index_tag_batch[i]= (current_index_tag_batch[i] + 1) % batch_size;
                    current_length_tag_batch[i]--;

                    // remove the first element of tag_batches[i]
//                    tag_batches[i].erase(tag_batches[i].begin());
                } else {
                    // add the first part of the tag run to the extracted tags and update the tag run in the tag_batches
                    extracted_tags[i].push_back(std::make_pair(tag_batches[i][current_index_tag_batch[i]].first, requests[i] - current_extracted));
                    tag_batches[i][current_index_tag_batch[i]].second -= (requests[i] - current_extracted);
                    current_extracted = requests[i];
                }

            }



//            while (current_extracted < requests[i]) {
//                if (tag_batches[i].empty() || tag_batches[i][0].second == 0) {
//                    // TODO: FIX THIS
////                    refill_tags();
//                    if (tag_batches[i].empty() || tag_batches[i][0].second == 0) {
//                        std::cerr << "No more tags to read from file: " << files[i] << std::endl;
//                        break;
//                    }
//                }
//                if (current_extracted + tag_batches[i][0].second <= requests[i]) {
//                    // add the whole tag run to the extracted tags and delete it from the tag_batches
//                    extracted_tags[i].push_back(tag_batches[i][0]);
//                    current_extracted += tag_batches[i][0].second;
//                    // remove the first element of tag_batches[i]
//                    tag_batches[i].erase(tag_batches[i].begin());
//
//                } else {
//                    // add the first part of the tag run to the extracted tags and update the tag run in the tag_batches
//                    extracted_tags[i].push_back(std::make_pair(tag_batches[i][0].first, requests[i] - current_extracted));
//                    tag_batches[i][0].second -= (requests[i] - current_extracted);
//                    current_extracted = requests[i];
//                }
//
//            }

        }


        refill_tags();

        notifyNext(thread_id);
        return extracted_tags;
    }


    pos_t getNextBlock(size_t thread_id, size_t fileIndex) {
        waitForTurn(thread_id); // Ensure the thread waits for its turn

        pos_t currentTag;

//        {
//            // Lock file access using the correct mutex for the file
//            std::lock_guard <std::mutex> lock(*mutexes[fileIndex]);
//
//            // Check if there are remaining tags in the current block
//            if (remaining_lengths[fileIndex] > 0) {
//                currentTag = current_nodes[fileIndex];
//                remaining_lengths[fileIndex]--;
//
//                // If the block is exhausted, load the next block
//                if (remaining_lengths[fileIndex] == 0) {
//                    if (file_positions[fileIndex] < file_buffers[fileIndex].size()) {
//                        auto nextBlock = panindexer::TagArray::decode_run(
//                                gbwt::ByteCode::read(file_buffers[fileIndex], file_positions[fileIndex]));
//
//
////                    current_nodes[i] = block.first;
////                    remaining_lengths[i] = block.second;
//
////                    std::ifstream in(files[fileIndex]);
////                    in.seekg(file_positions[fileIndex]);
////
////                    auto nextBlock = panindexer::TagArray::load_block_at(in, file_positions[fileIndex]);
//
//                        if (nextBlock.second == 0) {
//                            std::cerr << "Emptyyy block encountered in file: " << files[fileIndex] << std::endl;
//                        }
//
//                        current_nodes[fileIndex] = nextBlock.first;
//                        remaining_lengths[fileIndex] = nextBlock.second;
//                    } else {
//                        std::cerr << "Empty block encountered in file: " << files[fileIndex] << std::endl;
//                    }
//                }
//            } else {
//                // Load the next block if no remaining length
////                std::ifstream in(files[fileIndex]);
////                in.seekg(file_positions[fileIndex]);
//
////                auto nextBlock = panindexer::TagArray::load_block_at(in, file_positions[fileIndex]);
//
//                if (file_positions[fileIndex] < file_buffers[fileIndex].size()) {
//                    auto nextBlock = panindexer::TagArray::decode_run(
//                            gbwt::ByteCode::read(file_buffers[fileIndex], file_positions[fileIndex]));
//
//
//                    if (nextBlock.second == 0) {
//                        std::cerr << "Empty block encountered in file:" << files[fileIndex] << std::endl;
//                    }
//
//                    current_nodes[fileIndex] = nextBlock.first;
//                    remaining_lengths[fileIndex] = nextBlock.second - 1;
//                    currentTag = current_nodes[fileIndex];
//                } else {
//                    std::cerr << "Emptyyy block encountered in file: " << files[fileIndex] << std::endl;
//                }
//            }
//        }
//
//        {
//            // Output current processing details
//            std::lock_guard <std::mutex> lock(output_mutex);
////            std::cerr << "Thread " << thread_id << " File " << fileIndex
////                      << " Block " << currentTag
////                      << " Remaining Length " << int(remaining_lengths[fileIndex])
////                      << std::endl;
//        }
//
//        notifyNext(thread_id); // Notify the next thread

        return currentTag;
    }

private:
    void initializeFiles() {
        std::cerr << "Initializing files" << std::endl;
        for (size_t i = 0; i < files.size(); i++) {
            sdsl::int_vector_buffer<8> in(files[i], std::ios::in);
//            std::ifstream in(files[i]);
            if (!in.is_open()) {
                throw std::runtime_error("Cannot open file: " + files[i]);
            }
            file_buffers[i] = std::move(in);

            // print the size of the files
//            std::cerr << "File size: " << file_buffers[i].size() << std::endl;

            // want to read batch_size of the tags from each file and store them in the tag_batches

            for (size_t j = 0; j < batch_size; j++) {
                if (file_positions[i] >= file_buffers[i].size()) {
                    std::cerr << "End of file reached for: " << files[i] << std::endl;
                    file_end_reached[i] = true;
                    break; // Stop reading if end of file is reached
                }

                auto tag_block = panindexer::TagArray::decode_run(
                        gbwt::ByteCode::read(file_buffers[i], file_positions[i])
                );
                tag_batches[i][j] = tag_block;
                current_length_tag_batch[i]++;

            }

        }
    }

    // This function checks for all the tag files it has that if the current number of batch files are less than batch_size/3 it will refill the tags
    void refill_tags() {
//        std::cerr << "REFILL" << std::endl;


        for (size_t i = 0; i < files.size(); i++){
            if (current_length_tag_batch[i] < batch_size / 3 && !file_end_reached[i]) {
                size_t remaining_tags = current_length_tag_batch[i];
                std::vector<std::pair<pos_t, uint8_t>> new_tags;
                size_t new_tags_size = 0;
                while (new_tags_size + remaining_tags < batch_size) {
                    if (file_positions[i] >= file_buffers[i].size()) {
                        // End of file reached, stop reading
                        file_end_reached[i] = true;
                        std::cerr << "REFILL End of file reached for: " << files[i] << std::endl;
                        break;
                    }

                    auto tag_block = panindexer::TagArray::decode_run(
                            gbwt::ByteCode::read(file_buffers[i], file_positions[i])
                    );

//                    new_tags.push_back(tag_block);
                    tag_batches[i][(current_length_tag_batch[i] + current_index_tag_batch[i]) % batch_size] = tag_block;
                    current_length_tag_batch[i]++;
                    new_tags_size++;
                }

//                if (new_tags_size > 0) {
//                    tag_batches[i].insert(tag_batches[i].end(), new_tags.begin(), new_tags.end());
//                }
            }
        }
//        std::cerr << "Refilling tags" << std::endl;
//        for (size_t i = 0; i < files.size(); i++) {
//            // Check if the remaining tags in the batch are less than batch_size / 3
//            if (tag_batches[i].size() < batch_size / 3 && !file_end_reached[i]) {
//                size_t remaining_tags = tag_batches[i].size();
//                // TODO: if there were memory issues, removing this might help
//                std::vector<std::pair<pos_t, uint8_t>> new_tags;
//
//                // Read more tags until we fill the batch or reach the end of the file
//                while (new_tags.size() + remaining_tags < batch_size) {
//                    if (file_positions[i] >= file_buffers[i].size()) {
//                        // End of file reached, stop reading
//                        file_end_reached[i] = true;
//                        std::cerr << "End of file reached for: " << files[i] << std::endl;
//                        break;
//                    }
//
//                    auto tag_block = panindexer::TagArray::decode_run(
//                            gbwt::ByteCode::read(file_buffers[i], file_positions[i])
//                    );
//                    new_tags.push_back(tag_block);
//                }
////
////                std::cerr << "file pos in refill " << file_positions[i] << std::endl;
////                std::cerr << "file size " << file_buffers[i].buffersize() << std::endl;
//
//
//                // if new tags is not empty Add new tags to the existing batch
//                if (!new_tags.empty()) {
//                    tag_batches[i].insert(tag_batches[i].end(), new_tags.begin(), new_tags.end());
//                }
//            }
//        }
    }




    void waitForTurn(size_t thread_id) {
        std::unique_lock <std::mutex> lock(thread_mutex);
        thread_cv.wait(lock, [this, thread_id] {
            return thread_id == current_thread_id;
        });
    }

    void notifyNext(size_t thread_id) {
        std::lock_guard <std::mutex> lock(thread_mutex);

        // Advance to the next thread
        current_thread_id++;
        if (current_thread_id == n_threads) {
            current_thread_id = 0; // Reset to thread 0
        }

        thread_cv.notify_all(); // Notify waiting threads
    }

    // Data Members
    std::vector <std::string> files;          // List of file paths
    std::vector <bool> file_end_reached;      // Flag to indicate end of file

    size_t n_threads;                        // Number of threads
    size_t current_thread_id;                // Current thread ID to execute
    std::vector <std::unique_ptr<std::mutex>> mutexes; // Mutex for each file
    std::mutex thread_mutex;                 // Mutex for thread coordination
    std::condition_variable thread_cv;       // Condition variable for thread coordination
    std::mutex output_mutex;                 // Mutex for synchronized output

    // File Handling
    std::vector <gbwt::size_type> file_positions;      // Current file positions
    std::vector <sdsl::int_vector_buffer<8>> file_buffers; // File buffers for reading
    size_t batch_size;                      // Number of tags to read in a batch
    std::vector <size_t> current_index_tag_batch; // Current index of tag batch for each file
    std::vector <std::vector<std::pair < pos_t, uint8_t>>> tag_batches; // buffer of tags of each tag_file
    std::vector <size_t> current_length_tag_batch; // Remaining length of tag batch for each file

};


// This function extract the tags starting from the starting_run first position to the starting_run + runs_per_thread last position
void extract_tags_batch(const FastLocate &r_index, FileReader &reader, size_t thread_id, std::vector<int> comp_to_file,
                        unordered_map <size_t, size_t> seq_id_to_comp_id, std::vector<std::pair<pos_t, uint8_t>> &buffer,
                        size_t starting_run, size_t runs_per_thread) {

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
//    std::chrono::duration<double> duration1 = time2 - time1;
//    std::cerr << "Indexing unique kmers took " << duration1.count() << " seconds" << std::endl;
#endif

//    std::cerr << "Extracting tags starting run " << starting_run << std::endl;


    std::vector<size_t> request(reader.get_file_num(), 0);
    std::vector<int> index_to_file;
    // get the sample at the start of the starting_run
    if (starting_run >= r_index.tot_runs()){
        std::cerr << "Starting run is greater than total runs" << std::endl;
        return;
    }
    size_t index = r_index.getSample(starting_run);

//    std::cerr << "Index: " << index << std::endl;
    size_t end_index;
    size_t last_run_index = -1;



    // TODO: check if this is correct
    if (starting_run + runs_per_thread >= r_index.tot_runs()){
        std::cerr << "End of runs reached" << std::endl;
        // what is the last run first position
        last_run_index = r_index.getSample(r_index.tot_runs() - 1);
        std::cerr << "Last run index: " << last_run_index << std::endl;
        end_index = -1;
    } else {
        end_index = r_index.getSample(starting_run + runs_per_thread);
    }


    int kk = 0;
//    std::cerr << "the end index is " << end_index << std::endl;
//    std::cerr << "Starting creating the request list " << std::endl;
    while (index != end_index && index != last_run_index){
//        kk++;
//        std::cerr << "kk " << kk << std::endl;

        auto seq_id = r_index.seqId(index);
        // want to get the file number that is associated with the seq id
        auto current_file = comp_to_file[seq_id_to_comp_id[seq_id]];
        index_to_file.push_back(current_file);
        request[current_file] += 1;

        index = r_index.locateNext(index);
//        std::cerr << index << std::endl;
    }


    if (index == last_run_index){
        // close all open files in reader
//        reader.closeAllFiles();


        std::cerr << "hit last run in block " << r_index.blocks.size() << std::endl;
        // what is the last run size
        auto last_block = r_index.blocks[r_index.blocks.size() - 1];
        if (last_block.runs.empty()){
            std::cerr << "last block is empty " << std::endl;
        } else {
            int last_run_size = r_index.blocks.back().last_run_size();

            // have to just traverse the last run
//            std::cerr << "LAST RUN" << std::endl;
            for (int i = 0; i < last_run_size; i++){
//                std::cerr << "index " << index;
                auto seq_id = r_index.seqId(index);
                // want to get the file number that is associated with the seq id
                auto current_file = comp_to_file[seq_id_to_comp_id[seq_id]];
                index_to_file.push_back(current_file);
                request[current_file] += 1;

                index = r_index.locateNext(index);
            }

        }


    }

//#if TIME
//    auto time2 = chrono::high_resolution_clock::now();
//    std::chrono::duration<double> duration1 = time2 - time1;
//    std::cerr << "The total times getting sequence ids take in extract tag " << duration1.count() << " seconds" << "for thread " << thread_id << std::endl;
//#endif

    // print the items in request
//    for (size_t i = 0; i < request.size(); i++){
//        std::cerr << "Request " << request[i] << std::endl;
//    }

//    std::cerr << "Starting extracting the tags" << std::endl;
    // we have the tags we want to extract from each file. now we extract the tags for each index
    std::vector<std::vector<std::pair<pos_t, uint8_t>>> run_tags = reader.extract_requested_tags(thread_id, request);
//
//#if TIME
//    auto time3 = chrono::high_resolution_clock::now();
//    std::chrono::duration<double> duration2 = time3 - time2;
//    std::cerr << "The total times reader.extract_tags take " << duration2.count() << " seconds" << "for thread " << thread_id  << std::endl;
//#endif
//    std::cerr << "Extracted the tags" << std::endl;
    pos_t current_tag;
    size_t current_index = 0;
    std::vector<std::pair<pos_t, uint8_t>> temp_buffer;

    // traverse the index_to_file and get the tags from the run_tags
    for (size_t i = 0; i < index_to_file.size(); i++){
        current_tag = run_tags[index_to_file[i]][0].first;
        // print the current tag
//        std::cerr << "Current tag: " << current_tag << std::endl;
        run_tags[index_to_file[i]][0].second--;
        if (run_tags[index_to_file[i]][0].second == 0){
            run_tags[index_to_file[i]].erase(run_tags[index_to_file[i]].begin());
        }
        if (i == 0) {
            temp_buffer.push_back(make_pair(current_tag, 1));
            current_index++;
        } else {
            if (current_tag == temp_buffer[current_index - 1].first) {
                temp_buffer[current_index - 1].second++;
            } else {
                temp_buffer.push_back(make_pair(current_tag, 1));
                current_index++;
            }
        }
    }

//#if TIME
//    auto time4 = chrono::high_resolution_clock::now();
//    std::chrono::duration<double> duration3 = time4 - time3;
//    std::cerr << "The total times converting back tags take " << duration3.count() << " seconds" << "for thread " << thread_id  << std::endl;
//#endif
    buffer = temp_buffer;

    // print the temp_buffer
//    for (size_t i = 0; i < temp_buffer.size(); i++){
//        std::cerr << "Temp buffer: " << temp_buffer[i].first << " " << int(temp_buffer[i].second) << std::endl;
//    }



}

        // get the number of threads we set using omp_set_num_threads(threads);

void extract_tag(const FastLocate &r_index, FileReader &reader, size_t thread_id, int &current, bool first,
                 std::vector<int> comp_to_file, unordered_map <size_t, size_t> seq_id_to_comp_id, pos_t &buffer,
                 int thread_num) {
    // current is the previous position which we have to run locate_next_nth on, this is when current is none zero
    // current equal to zero means this is the first time we are calling this thread_id
    if (first) {

        current = r_index.locateFirst();
        for (size_t i = 0; i < r_index.tot_strings() + thread_id; i++) {
            current = r_index.locateNext(current);
        }

    } else {
        // get the number of threads we set using omp_set_num_threads(threads);
//        int threads = omp_get_num_threads();
//        std::cerr << thread_num << std::endl;
        current = r_index.locate_next_nth(current, thread_num);
//        current = r_index.locate()
//        if (current == 0) std::cerr << "current" << current << std::endl;
    }


    auto seq_id = r_index.seqId(current);

    // want to get the file number that is associated with the seq id
    auto current_file = comp_to_file[seq_id_to_comp_id[seq_id]];


    pos_t block = reader.getNextBlock(thread_id, current_file);
    buffer = block;
    // print the block
//    std::cerr << "Block " << block << "from thread " << thread_id << std::endl;
}


std::vector <std::string> get_files_in_dir(const std::string &directoryPath) {
    std::vector <std::string> files;

    // Check if the path is a valid directory
    if (!fs::is_directory(directoryPath)) {
        std::cerr << "Path is not a directory: " << directoryPath << std::endl;
        return files;
    }

    // Iterate over entries in the directory
    for (const auto &entry: fs::directory_iterator(directoryPath)) {
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
    if (!sdsl::load_from_file(r_index, r_index_file)) {
        std::cerr << "Cannot load the r-index from " << r_index_file << std::endl;
        std::exit(EXIT_FAILURE);
    }


    std::cerr << "Getting the lists of tag files" << std::endl;
    // get the list of files in the directory
    std::vector <std::string> files = get_files_in_dir(tag_array_index_dir);

    int number_of_file = files.size();

    std::cerr << "The list of files are: " << std::endl;
    for (auto &file: files) {
        std::cerr << file << std::endl;
    }


    std::cerr << "Finding the node to component mapping" << std::endl;
    // finding the components
    auto node_to_comp_map = node_to_component(gbz);




//    std::cerr << "Node to comp map: " << std::endl;
//    // print all node_to_comp_map
//    for (auto it = node_to_comp_map.begin(); it != node_to_comp_map.end(); ++it){
//        std::cerr << it->first << " " << it->second << std::endl;
//    }


    std::vector<int> file_to_comp(number_of_file);
    std::vector<int> comp_to_file(number_of_file);

//    std::vector <pos_t> current_node_of_files(number_of_file);
//    std::vector <uint8_t> current_remaining_length_of_files(number_of_file);
//    std::vector <size_t> current_file_positions(number_of_file, 0);

    std::cerr << "Initializing the reader" << std::endl;
    FileReader reader(files, threads, 1000000);
    std::cerr << "Creating the mapping from comp to tag files" << std::endl;
    // for each tag block files, we read the first block and read the first node
    for (auto i = 0; i < number_of_file; i++) {
//        std::ifstream in(files[i]);
//        std::pair <pos_t, uint8_t> first_block = panindexer::TagArray::load_block_at(in, current_file_positions[i]);
//        std::cerr << "The first block of file " << files[i] << " is " << first_block.first << " " << int(first_block.second) << std::endl;
//        current_node_of_files[i] = first_block.first;
//        current_remaining_length_of_files[i] = first_block.second;

//        std::cerr << "The first block of file " << files[i] << " is " << first_block.first << " " << int(first_block.second) << std::endl;

        // get the component of the node of the first block
        size_t comp = node_to_comp_map[id(reader.get_first_tag(i))];

//        std::cerr << "The component of the first block of file is " << comp << std::endl;
        file_to_comp[i] = comp;
        comp_to_file[comp] = i;

    }


    std::cerr << "The mapping from comp to tag files is done" << std::endl;


    // print all file_to_comp
//    for (size_t i = 0; i < number_of_file; i++){
//        std::cerr << "comp " << i << " file " << comp_to_file[i] << std::endl;
//    }

    unordered_map <size_t, size_t> seq_id_to_comp_id;

    // get the first node of each path and get the component id of the node
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < r_index.tot_strings(); i++) {
        auto seq_graph_nodes = gbz.index.extract(i * 2);

        // extracting the first node of the sequence component id
        seq_id_to_comp_id[i] = node_to_comp_map[gbwt::Node::id(seq_graph_nodes[0])];
//        std::cerr << "seq id " << i << " comp id " << seq_id_to_comp_id[i] << std::endl;
    }
    std::cerr << "The mapping from seq id to comp id is done" << std::endl;

    // reading the r-index file and traversing the r-index and getting the sequence number of each BWT position

    auto current = r_index.locateFirst();
    // have to bypass the first #seqs because there are no tags for them
//    for (size_t i = 0; i < r_index.tot_strings(); i++){
////        cerr << "The first node is " << current << endl;
//        current = r_index.locateNext(current);
//    }

//    cerr << "The first nodes " << current << endl;
    pos_t current_tag;


    // note that the first #num_seq tags are correspond to the ENDMARKERs
    auto total_tags_count = r_index.get_sequence_size() - r_index.tot_strings();
    std::cerr << "Total tags count " << total_tags_count << std::endl;
//    vector<pair<pos_t, uint8_t>> tag_runs;
    size_t current_run = 0;


    cerr << "Creating the whole genome tag array indexing" << endl;
    TagArray tag_array;

#if TIME
    auto time2 = chrono::high_resolution_clock::now();
#endif
//    tag_array.serialize_run_by_run(tag_runs, "whole_genome_tag_array.tags");


    std::cerr << "Thread lists " << std::endl;
    std::vector <std::thread> threads_list;
    threads_list.reserve(threads);
    std::vector<std::vector<std::pair<pos_t, uint8_t>>> thread_buffers(threads);
    std::vector<int> current_position(threads, 0);
    int batch_size = 2000;
    vector <pos_t> tags_batch(batch_size);
//    vector <pair<pos_t, uint8_t>> tag_runs2;
//    vector<pos_t> tags2(total_tags_count);
    current_run = 0;
    vector <pair<pos_t, uint8_t>> tag_runs_check;
    uint8_t one = 1;

//    extract_tags_batch(r_index, reader, 0, comp_to_file, seq_id_to_comp_id, thread_buffers[0], threads, 0, 5);
//    extract_tags_batch(r_index, reader, 1, comp_to_file, seq_id_to_comp_id, thread_buffers[0], threads, 5, 5);

    const std::string filename = "whole_genome_tag_array_parallel_sdsl.tags";

    // Check if the file exists and delete it
    if (std::filesystem::exists(filename)) {
        if (std::remove(filename.c_str()) != 0) {
            std::cerr << "Error: Unable to delete the existing file.\n";
            return 1; // Exit with error
        } else {
            std::cout << "Existing file deleted successfully.\n";
        }
    }



    // Open the file for writing
    std::ofstream out(filename, std::ios::binary | std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open file for writing.\n";
        return 1; // Exit with error
    }

    size_t starting_run = 0;
    size_t run_per_thread = 100;


    std::cerr << "We will handle " << r_index.tot_runs() << " runs" << std::endl;
    cerr << "Merging tags and creating the whole genome tag array indexing" << endl;
    for (size_t to_read = 0; to_read < threads && to_read < r_index.tot_runs(); to_read++) {
        threads_list.emplace_back(extract_tags_batch, std::ref(r_index), std::ref(reader), to_read,
                                  comp_to_file, seq_id_to_comp_id,
                                  std::ref(thread_buffers[to_read]), starting_run, run_per_thread);
        starting_run += run_per_thread;
    }




    size_t number_of_jobs = (r_index.tot_runs() + run_per_thread - 1) / run_per_thread;
    std::pair<pos_t, uint8_t> previous_last_run;

    for (size_t to_write = 0; to_write < number_of_jobs; to_write++) {
        if (to_write % 10000 == 0) {
            std::cerr << "Writing job " << to_write << std::endl;
        }
        size_t thread_id = to_write % threads;
        threads_list[thread_id].join();
        std::vector<std::pair<pos_t, uint8_t>> current_tags = thread_buffers[thread_id];

        if (to_write == 0) {
            // write everything but the last pair
            previous_last_run = current_tags.back();
            current_tags.pop_back();
        } else {
            if (current_tags[0].first == previous_last_run.first) {
                current_tags[0].second += previous_last_run.second;
            } else {
                std::vector<std::pair<pos_t, uint8_t>> temp = {previous_last_run};
                tag_array.serialize_run_by_run(out, temp);
            }
            if (to_write < number_of_jobs - 1) {
                previous_last_run = current_tags.back();
                current_tags.pop_back();
            }
        }
        tag_array.serialize_run_by_run(out, current_tags);




        if (to_write + threads < number_of_jobs) {
            threads_list[thread_id] = std::thread(extract_tags_batch, std::ref(r_index), std::ref(reader), thread_id,
                                                  comp_to_file, seq_id_to_comp_id,
                                                  std::ref(thread_buffers[thread_id]), starting_run, run_per_thread);
            starting_run += run_per_thread;
        }


    }
    for (auto &thread: threads_list) {
        if (thread.joinable()) {
            thread.join();
        }
    }

//
//    for (size_t to_read = 0; to_read < threads && to_read < total_tags_count; to_read++) {
//        threads_list.emplace_back(extract_tag, std::ref(r_index), std::ref(reader), to_read,
//                                  std::ref(current_position[to_read]), true, comp_to_file, seq_id_to_comp_id,
//                                  std::ref(thread_buffers[to_read]), threads);
//    }
//
//
//    for (size_t to_write = 0; to_write < total_tags_count; to_write++) {
//        if (to_write % 1000 == 0) {
//            std::cerr << "Writing tag " << to_write << std::endl;
//        }
//        size_t thread_id = to_write % threads;
//        threads_list[thread_id].join();
//        pos_t current_tag = thread_buffers[thread_id];
//        tags_batch[to_write % batch_size] = current_tag;
//
//
//        if (current_run == 0) {
//            tag_runs2.push_back(make_pair(current_tag, one));
//            current_run++;
//        } else {
//            if (current_tag == tag_runs2[current_run - 1].first) {
//                tag_runs2[current_run - 1].second++;
//            } else {
//                tag_runs2.push_back(make_pair(current_tag, one));
//                current_run++;
//            }
//        }
//
//        // whenever the batch is full, we serialize the batch
//
//        if ((tag_runs2.size() == batch_size || to_write == total_tags_count - 1) && tag_runs2.size() != 0) {
//            // pop the last element
//            auto temp = tag_runs2.back();
//            if (to_write < total_tags_count - 1) {
//                tag_runs2.pop_back();
//            }
//            tag_array.serialize_run_by_run(out, tag_runs2);
//            tag_runs2.clear();
//            current_run = 1;
//            tag_runs2.push_back(temp);
//        }
//
//
//        if (to_write + threads < total_tags_count) {
//            threads_list[thread_id] = std::thread(extract_tag, std::ref(r_index), std::ref(reader), to_write % threads,
//                                                  std::ref(current_position[thread_id]), false, comp_to_file,
//                                                  seq_id_to_comp_id, std::ref(thread_buffers[thread_id]), threads);
//        }
//
//
//    }
//    for (auto &thread: threads_list) {
//        if (thread.joinable()) {
//            thread.join();
//        }
//    }

    out.close();


#if TIME
    auto time3 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = time3 - time2;
    std::cerr << "Converting tags using multiple threads took " << duration2.count() << " seconds" << std::endl;
#endif




//    TagArray tt;
//    tt.deserialize_run_by_run("whole_genome_tag_array_parallel_sdsl.tags");
//    auto tag_runs2 = tt.get_tag_runs();
//    // compare tag runs from test and tt
//
//    for (auto i = 0; i < tag_runs2.size() - 1; i++){
//        if (tag_runs2[i].first == tag_runs2[i+1].first){
//            std::cerr << "Wrong in index " << i << " " << tag_runs2[i].first << " " << tag_runs2[i+1].first << std::endl;
//        }
//    }


    return 0;


}


