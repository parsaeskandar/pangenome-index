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



// if not define Time define it
#ifndef TIME
#define TIME 1
#endif



//using namespace gbwtgraph;
using namespace panindexer;
using namespace std;
using namespace gbwtgraph;

namespace fs = std::filesystem;

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
    for (auto it = node_to_comp_map.begin(); it != node_to_comp_map.end(); ++it){
        std::cerr << it->first << " " << it->second << std::endl;
    }


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
    for (size_t i = 0; i < number_of_file; i++){
        std::cerr << "comp " << i << " file " << comp_to_file[i] << std::endl;
    }

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

    vector<int> each_file(2, 0);

    // TODO: note that the first #num_seq tags are correspond to the ENDMARKERs
    vector<pos_t> tags(r_index.get_sequence_size() - r_index.tot_strings());
    vector<pair<pos_t, uint8_t>> tag_runs;
    size_t current_run = 0;
    for (size_t i = r_index.tot_strings() ; i < r_index.get_sequence_size(); i++){
        // seq_id is the same as GBWT path identifier
        auto seq_id = r_index.seqId(current);

        // want to get the file number that is associated with the seq id
        auto current_file = comp_to_file[seq_id_to_comp_id[seq_id]];
        each_file[current_file]++;

        // reading the next tag from the current file
        if (current_remaining_length_of_files[current_file] == 0){
            std::cerr << "Error: reaching the end of file too soon " << std::endl;
        } else {
            // Just getting the node and decrementing the remaining length
            current_tag = current_node_of_files[current_file];
            current_remaining_length_of_files[current_file]--;

            if (current_remaining_length_of_files[current_file] == 0){
                // read the next block
                std::ifstream in(files[current_file]);
                std::pair<pos_t, uint8_t> next_block = panindexer::TagArray::load_block_at(in, current_file_positions[current_file]);

                current_node_of_files[current_file] = next_block.first;
                current_remaining_length_of_files[current_file] = next_block.second;

                if (next_block.second == 0){

                    std::cerr << "Error: the next block is empty " << seq_id << std::endl;
                    current_remaining_length_of_files[current_file] = 0;
                    current_node_of_files[current_file] = next_block.first;
                }

            }


            if (current_run == 0){
                tag_runs.push_back(make_pair(current_tag, 1));
                current_run++;
            } else {
                if (current_tag == tag_runs[current_run - 1].first){
                    tag_runs[current_run - 1].second++;
                } else {
                    tag_runs.push_back(make_pair(current_tag, 1));
                    current_run++;
                }
            }
            tags[ i - r_index.tot_strings()] = current_tag;
        }

        if (i != r_index.get_sequence_size() - 1){
            current = r_index.locateNext(current);
        }
    }


    // print the tags
    cerr << "tags size " << tags.size() << " r index size " << r_index.get_sequence_size() << endl;

    // print the tag_runs
//    for (size_t i = 0; i < tag_runs.size(); i++){
//        cerr << "tag " << tag_runs[i].first << " length " << int(tag_runs[i].second) << endl;
//    }

    TagArray tag_array;
    tag_array.serialize_run_by_run(tag_runs, "whole_genome_tag_array.tags");

























}
