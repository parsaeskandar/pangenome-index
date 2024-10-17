//
// Created by seeskand on 9/13/24.
//
//
// Created by Parsa Eskandar on 7/9/24.
//

#include <iostream>
#include "pangenome_index/tag_arrays.hpp"

namespace panindexer {


    TagArray::TagArray() {
        // Initialize the run_starts bit vector with the given size and set all bits to 0

    }


    void TagArray::load_bptree(BplusTree <Run> &bptree, size_t bwt_size) {
        // have to iterate over the bptree and add the runs to the encoded_runs vector and the bit-vectors
        auto number_of_runs = bptree.get_bpt_run_count();
        encoded_runs_starts = sdsl::bit_vector(7 * number_of_runs, 0);
        sdsl::sd_vector_builder builder(bwt_size + 1, number_of_runs);
        size_t encoded_runs_current_size = 0;


        for (auto it = bptree.begin(); it != bptree.end(); ++it) {
            Run current_item = *it;
            if (current_item.graph_position.value != 0) {
                // we have an actual run with graph position
                auto next_it = it;
                ++next_it;
                if (next_it != bptree.end()) {
                    Run next_item = *next_it;
                    // adding the pair (graph_value, length) as ByteCode to the encoded_runs vector
                    pos_t current_pos = current_item.graph_position.decode();
                    gbwt::size_type encoded =
                            (gbwtgraph::offset(current_pos)) | (gbwtgraph::is_rev(current_pos) << 10) |
                            ((next_item.start_position - current_item.start_position) << 11) |
                            (gbwtgraph::id(current_pos) << 19);
                    encoded_runs_starts[encoded_runs.size()] = 1;

//                    cout << encoded_runs.size() << endl;
//                    gbwt::size_type s = encoded_runs.size();
                    gbwt::ByteCode::write(encoded_runs, encoded);
                    // setting the bit in the run_starts bit vector
                    encoded_runs_current_size = encoded_runs.size();


//                    cout << encoded_runs.size() << endl;

                    // setting the bit in the bwt_intervals bit vector
                    builder.set(current_item.start_position);
//                    cout << "bwt " << current_item.start_position << endl;


                }
            }
        }
//        encoded_runs_starts[encoded_runs.size()] = 1;

        // resizing the
        encoded_runs_starts.resize(encoded_runs_current_size + 1);

        bwt_intervals = sdsl::sd_vector<>(builder);

//        cout << "encoded_runs size: " << encoded_runs.size() << endl;
//        cout << "encoded_runs_starts size: " << encoded_runs_starts.size() << endl;
//        cout << "bwt_intervals size: " << bwt_intervals.size() << endl;



    }


    void TagArray::serialize(std::ostream &out) {
        // Serialize encoded_runs_starts
        sdsl::serialize(encoded_runs_starts, out);

        // Serialize bwt_intervals
        sdsl::serialize(bwt_intervals, out);

        // Serialize encoded_runs
        size_t size = encoded_runs.size();
        out.write(reinterpret_cast<const char *>(&size), sizeof(size));
        out.write(reinterpret_cast<const char *>(encoded_runs.data()), size * sizeof(gbwt::byte_type));

        cerr << "Finished serializing" << endl;
    }


    void TagArray::load(std::istream &in) {
        try {
            // Deserialize encoded_runs_starts
            sdsl::load(encoded_runs_starts, in);

            // print the first 10 bits
//            for (size_t i = 0; i < encoded_runs_starts.size(); i++) {
//                cerr << encoded_runs_starts[i] << " ";
//            }
//            cerr << endl;

            if (in.fail()) {
                throw std::runtime_error("Failed to load encoded_runs_starts");
            }

            // Deserialize bwt_intervals
            sdsl::load(bwt_intervals, in);
            if (in.fail()) {
                throw std::runtime_error("Failed to load bwt_intervals");
            }

//            for (size_t i = 0; i < bwt_intervals.size(); i++) {
//                cerr << bwt_intervals[i] << " ";
//            }
//            cerr << endl;

            // Deserialize encoded_runs
            size_t size;
            in.read(reinterpret_cast<char *>(&size), sizeof(size));
            if (in.fail()) {
                throw std::runtime_error("Failed to read size of encoded_runs");
            }

            // print the size
            cerr << "size: " << size << endl;


            encoded_runs.resize(size);
            if (size > 0) {
                in.read(reinterpret_cast<char *>(encoded_runs.data()), size * sizeof(gbwt::byte_type));
                if (in.fail()) {
                    throw std::runtime_error("Failed to read data of encoded_runs");
                }
            }


            // initialize the select and rank data structures
            encoded_runs_starts_select = sdsl::bit_vector::select_1_type(&encoded_runs_starts);
            bwt_intervals_rank = sdsl::sd_vector<>::rank_1_type(&bwt_intervals);


            cerr << "Finished deserializing" << endl;

        } catch (const std::bad_alloc &e) {
            cerr << "Memory allocation failed: " << e.what() << endl;
            throw;
        } catch (const std::exception &e) {
            cerr << "Deserialization error: " << e.what() << endl;
            throw;
        }
    }


    void TagArray::query(size_t start, size_t end) {

        // first have to find ranks of start and end in the bwt_intervals which are number of 1's less than start and end

//        sdsl::sd_vector<>::rank_1_type bwt_intervals_rank(&bwt_intervals);
        size_t first_bit_index = bwt_intervals_rank(start);
        if (start > 0 && bwt_intervals[start] == 1) {
            // if the start is 1, then the tags start from start position and we don't need the last 1 before start
            first_bit_index++;
        }

        size_t end_bit_index = bwt_intervals_rank(end + 1);

        size_t number_of_runs = end_bit_index - first_bit_index + 1;

        // where the first run starts
        std::uint64_t bit_location = encoded_runs_starts_select(first_bit_index);
        std::unordered_set <std::uint64_t> unique_positions;

        while (number_of_runs > 0) {
            gbwt::size_type decc;
            // the read function changes the bit_location to the next bit location
            decc = gbwt::ByteCode::read(encoded_runs, bit_location);

            size_t decoded_offset = decc & ((1LL << 10) - 1);
            bool decoded_flag = (decc >> 10) & 0x1;
            uint8_t decoded_length = (decc >> 11) & 0xFF;
            int64_t decoded_node_id = (decc >> 19);

            number_of_runs--;
            unique_positions.insert(
                    gbwtgraph::Position::encode(pos_t(decoded_node_id, decoded_offset, decoded_flag)).value);
        }

//        std::cout << (end_bit_index - first_bit_index + 1) << '\t' << unique_positions.size() << '\t' << (double)(end_bit_index - first_bit_index + 1)/((double)unique_positions.size()) << std::endl;
//        std::cerr << "Number of unique positions in the interval: " << unique_positions.size() << std::endl;


        // now have to find the locations of the bits first_bit_index + 1 to end_bit_index in the bwt_intervals


//        for(size_t i = first_bit_index; i <= end_bit_index; i++){
////            cout << "i: " << i << endl;
//
//            std::uint64_t bit_location = encoded_runs_starts_select(i);
////            cerr << "bit location: " << bit_location << endl;
//
//            // for each bit location, decoding the encoded runs that start from that location
//            gbwt::size_type decc;
//            gbwt::size_type ind = 0;
//
//            decc = gbwt::ByteCode::read(encoded_runs, bit_location);
//
//            size_t decoded_offset = decc & ((1LL << 10) - 1);
//            bool decoded_flag = (decc >> 10) & 0x1;
//            uint8_t decoded_length = (decc >> 11) & 0xFF;
//            int64_t decoded_node_id = (decc >> 19);
//
//
//                // Print the decoded values
////            std::cerr << "Decoded offset: " << decoded_offset << std::endl;
////            std::cerr << "Decoded flag: " << decoded_flag << std::endl;
////            std::cerr << "Decoded length: " << static_cast<int>(decoded_length) << std::endl;
////            std::cerr << "Decoded node ID: " << decoded_node_id << std::endl;
//
//            unique_positions.insert(gbwtgraph::Position::encode(vg::pos_t(decoded_node_id, decoded_offset, decoded_flag)).value);
//
//        }

//        std::cout << (end_bit_index - first_bit_index + 1) << '\t' << unique_positions.size() << '\t' << (double)(end_bit_index - first_bit_index + 1)/((double)unique_positions.size()) << std::endl;
//        std::cerr << "Number of unique positions in the interval: " << unique_positions.size() << std::endl;






    }

}



