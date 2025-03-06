//
// Created by seeskand on 3/4/25.
//


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

    void print_all_tags(){
        for (size_t i = 0; i < files.size(); i++){
            size_t tag_tot = 0;
            std::cerr << "File: " << files[i] << std::endl;
            for (size_t j = 0; j < current_length_tag_batch[i]; j++){
                std::cerr << "Tag: " << tag_batches[i][j].first << " Length: " << int(tag_batches[i][j].second) << std::endl;
                tag_tot += tag_batches[i][j].second;
            }
            std::cerr << "=============Total tags: " << tag_tot << std::endl;
        }
    }

    int get_file_num(){
        return files.size();
    }


    pos_t get_first_tag(int fileIndex) {
//        std::cerr << "Starting tag runs of file " << files[fileIndex] << std::endl;
//        for (size_t i = 0; i < 3; i++) {
//            std::cerr << tag_batches[fileIndex][i].first << " " << int(tag_batches[fileIndex][i].second) << std::endl;
//        }
        return tag_batches[fileIndex][0].first;
    }

    void print_first_n_item(int n, int fileIndex){
        for (int i = 0; i < n; i++){
            std::cerr << "Tag: " << tag_batches[fileIndex][i].first << " Length: " << int(tag_batches[fileIndex][i].second) << std::endl;
        }
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

    pos_t get_next_tag(int fileIndex) {
        if (current_length_tag_batch[fileIndex] == 0) {
            std::cerr << "GET NEXT TAG: No more tags is possible to read from file: " << files[fileIndex] << std::endl;
        }

        pos_t res = tag_batches[fileIndex][current_index_tag_batch[fileIndex]].first;
        tag_batches[fileIndex][current_index_tag_batch[fileIndex]].second--;

        if (tag_batches[fileIndex][current_index_tag_batch[fileIndex]].second == 0){
            current_index_tag_batch[fileIndex] = (current_index_tag_batch[fileIndex] + 1) % batch_size;
            current_length_tag_batch[fileIndex]--;
        }

        return res;
    }





    std::vector<std::vector<std::pair<pos_t, uint8_t>>> extract_requested_tags(
            size_t thread_id, const std::vector<size_t>& requests) {

        waitForTurn(thread_id);

        std::vector<std::vector<std::pair<pos_t, uint8_t>>> extracted_tags(files.size());

        for (size_t i = 0; i < files.size(); i++) {
            size_t current_extracted = 0;


            while (current_extracted < requests[i]){
                if (current_length_tag_batch[i] == 0) {
                    if (file_end_reached[i]) {
                        std::cerr << "No more tags is possible to read from file: " << files[i] << std::endl;
                        std::cerr << "Needed " << requests[i] << " but only extracted " << current_extracted << std::endl;
                        break;
                    } else {
                        refill_tags();
                    }
                }
                if (current_extracted + tag_batches[i][current_index_tag_batch[i]].second <= requests[i]){
                    // add the whole tag run to the extracted tags and delete it from the tag_batches
                    extracted_tags[i].push_back(tag_batches[i][current_index_tag_batch[i]]);
                    current_extracted += tag_batches[i][current_index_tag_batch[i]].second;
                    current_index_tag_batch[i] = (current_index_tag_batch[i] + 1) % batch_size;
                    current_length_tag_batch[i]--;

                } else {
                    // add the first part of the tag run to the extracted tags and update the tag run in the tag_batches
                    extracted_tags[i].push_back(std::make_pair(tag_batches[i][current_index_tag_batch[i]].first, requests[i] - current_extracted));
                    tag_batches[i][current_index_tag_batch[i]].second -= (requests[i] - current_extracted);
                    current_extracted = requests[i];
                }

            }

        }


        refill_tags();

        notifyNext(thread_id);
        return extracted_tags;
    }



private:
    void initializeFiles() {
        std::cerr << "Initializing files" << std::endl;
        for (size_t i = 0; i < files.size(); i++) {
            sdsl::int_vector_buffer<8> in(files[i], std::ios::in);
            if (!in.is_open()) {
                throw std::runtime_error("Cannot open file: " + files[i]);
            }
            file_buffers[i] = std::move(in);


            // want to read batch_size of the tags from each file and store them in the tag_batches

            for (size_t j = 0; j < batch_size; j++) {
                if (file_positions[i] >= file_buffers[i].size()) {
                    std::cerr << "End of file reached for: " << files[i] << std::endl;
                    file_end_reached[i] = true;
                    break;
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
        for (size_t i = 0; i < files.size(); i++) {
            if (current_length_tag_batch[i] < batch_size / 3 && !file_end_reached[i]) {
                size_t remaining_tags = current_length_tag_batch[i];
                std::vector <std::pair<pos_t, uint8_t>> new_tags;
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

                    tag_batches[i][(current_length_tag_batch[i] + current_index_tag_batch[i]) % batch_size] = tag_block;
                    current_length_tag_batch[i]++;
                    new_tags_size++;
                }
            }
        }
    }




    void waitForTurn(size_t thread_id) {
        std::unique_lock <std::mutex> lock(thread_mutex);
        thread_cv.wait(lock, [this, thread_id] {
            return thread_id == current_thread_id;
        });
    }

    void notifyNext(size_t thread_id) {
        std::lock_guard <std::mutex> lock(thread_mutex);

        current_thread_id++;
        if (current_thread_id == n_threads) {
            current_thread_id = 0;
        }

        thread_cv.notify_all();
    }

    std::vector <std::string> files;          // List of file paths
    std::vector <bool> file_end_reached;      // Flag to indicate end of file

    size_t n_threads;                        // Number of threads
    size_t current_thread_id;                // Current thread ID to execute
    std::vector <std::unique_ptr<std::mutex>> mutexes; // Mutex for each file
    std::mutex thread_mutex;                 // Mutex for thread coordination
    std::condition_variable thread_cv;       // Condition variable for thread coordination

    std::vector <gbwt::size_type> file_positions;      // Current file positions
    std::vector <sdsl::int_vector_buffer<8>> file_buffers; // File buffers for reading
    size_t batch_size;                      // Number of tags to read in a batch
    std::vector <size_t> current_index_tag_batch; // Current index of tag batch for each file
    std::vector <std::vector<std::pair < pos_t, uint8_t>>> tag_batches; // buffer of tags of each tag_file
    std::vector <size_t> current_length_tag_batch; // Remaining length of tag batch for each file

};








std::vector <std::string> get_files_in_dir(const std::string &directoryPath) {
    std::vector <std::string> files;
    if (!fs::is_directory(directoryPath)) {
        std::cerr << "Path is not a directory: " << directoryPath << std::endl;
        return files;
    }

    for (const auto &entry: fs::directory_iterator(directoryPath)) {
        if (fs::is_regular_file(entry.status())) {
            files.push_back(entry.path().string());
        }
    }

    return files;
}

int main(int argc, char **argv) {

    std::string gbz_graph = std::string(argv[1]);
    std::string r_index_file = std::string(argv[2]);
    std::string tag_array_index_dir = std::string(argv[3]);
    int threads = 8;

    GBZ gbz;
    cerr << "Loading the graph file" << endl;
    sdsl::simple_sds::load_from(gbz, gbz_graph);

    std::cerr << "Getting the lists of tag files" << std::endl;
    // get the list of files in the directory
    std::vector <std::string> files = get_files_in_dir(tag_array_index_dir);

    int number_of_file = files.size();

    std::cerr << "The list of files are: " << std::endl;
    for (auto &file: files) {
        std::cerr << file << std::endl;
    }


    cerr << "Reading the whole genome r-index file" << endl;
    FastLocate r_index;
    if (!sdsl::load_from_file(r_index, r_index_file)) {
        std::cerr << "Cannot load the r-index from " << r_index_file << std::endl;
        std::exit(EXIT_FAILURE);
    }


    vector<size_t> tag_count_per_file;
    tag_count_per_file.resize(number_of_file, 0);


    for (size_t i = 0; i < files.size(); i++) {
        sdsl::int_vector_buffer<8> in(files[i], std::ios::in);
        gbwt::size_type file_pos = 0;
        if (!in.is_open()) {
            throw std::runtime_error("Cannot open file: " + files[i]);
        }
        while (file_pos < in.size()) {
            auto tag_block = panindexer::TagArray::decode_run(
                    gbwt::ByteCode::read(in, file_pos)
            );

            tag_count_per_file[i] += tag_block.second;
        }


    }



    for (size_t i = 0; i < number_of_file; i++) {
        std::cerr << "File: " << files[i] << " Tags: " << tag_count_per_file[i] << std::endl;
    }



    std::cerr << "Finding the node to component mapping" << std::endl;
    // finding the components
    std::unordered_map<nid_t, size_t> node_to_comp_map = node_to_component(gbz);

    std::vector<int> file_to_comp(number_of_file);
    std::vector<int> comp_to_file(number_of_file);

    std::cerr << "Initializing the reader" << std::endl;
    FileReader reader(files, threads, 1000000);
    std::cerr << "Creating the mapping from comp to tag files" << std::endl;
    // for each tag block files, we read the first block and read the first node
    for (auto i = 0; i < number_of_file; i++) {
        // get the component of the node of the first block
        size_t comp = node_to_comp_map[id(reader.get_first_tag(i))];
        std::cerr << "The component of the first block of file " << files[i] << " is " << comp << " first tag " << reader.get_first_tag(i)  << " node id is " << id(reader.get_first_tag(i)) << std::endl;
        file_to_comp[i] = comp;
        comp_to_file[comp] = i;

    }

    std::cerr << "The mapping from comp to tag files is done" << std::endl;


    auto total_strings = r_index.tot_strings();
    std::vector<size_t> seq_id_to_comp_id;
    seq_id_to_comp_id.resize(total_strings);



    // TODO: make this multithreaded
    // get the first node of each path and get the component id of the node
#pragma omp parallel for
    for (size_t i = 0; i < total_strings; i++) {
        auto seq_graph_nodes = gbz.index.extract(i * 2);
        if (!seq_graph_nodes.empty()) {
            size_t node_id = gbwt::Node::id(seq_graph_nodes[0]);
            seq_id_to_comp_id[i] = node_to_comp_map.at(node_id);
        }
    }
    std::cerr << "The mapping from seq id to comp id is done" << std::endl;


    // note that the first #num_seq tags are correspond to the ENDMARKERs
    auto total_tags_count = r_index.get_sequence_size() - total_strings;
    std::cerr << "Total tags count " << total_tags_count << std::endl;



//    auto total_strings = r_index.tot_strings();
    auto first = r_index.locateFirst();

    for (int i = 0; i < total_strings - 1 ; i++){
        first = r_index.locateNext(first);
    }

    auto total_bwt_size = r_index.get_sequence_size();

    vector <size_t> rindex_per_file;
    rindex_per_file.resize(number_of_file, 0);

    std::cerr << "Total string " << total_strings << " Total bwt size " << total_bwt_size << std::endl;

    for (size_t i = total_strings; i < total_bwt_size; i++){
        auto seq_id = r_index.seqId(first);
//        std::cerr << "Seq id " << seq_id << std::endl;
        // want to get the file number that is associated with the seq id
        auto current_file = comp_to_file[seq_id_to_comp_id[seq_id]];
        rindex_per_file[current_file] += 1;
        first = r_index.locateNext(first);
    }


    for (size_t i = 0; i < number_of_file; i++) {
        std::cerr << "File: " << files[i] << " R-index: " << rindex_per_file[i] << " Tag count from file " << tag_count_per_file[i] << std::endl;
    }


}