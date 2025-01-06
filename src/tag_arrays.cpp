//
// Created by seeskand on 9/13/24.
//
//
// Created by Parsa Eskandar on 7/9/24.
//


#include "pangenome_index/tag_arrays.hpp"

namespace panindexer {

    using nid_t = handlegraph::nid_t;
    using offset_t = handlegraph::offset_t;


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

    void TagArray::serialize_run_by_run(std::ofstream& out, const std::vector<std::pair<pos_t, uint8_t>>& tag_runs) {
        for (const auto& [value, run_length] : tag_runs) {
            out.write(reinterpret_cast<const char*>(&value), sizeof(pos_t));
            out.write(reinterpret_cast<const char*>(&run_length), sizeof(uint8_t));
        }
    }

    void TagArray::deserialize_run_by_run(const std::string& file_name) {
        std::ifstream in(file_name, std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "Error: Cannot open file for reading.\n";
        }

        pos_t value;
        uint8_t run_length;

        while (in.read(reinterpret_cast<char*>(&value), sizeof(pos_t))) {
            in.read(reinterpret_cast<char*>(&run_length), sizeof(uint8_t));
            this->tag_runs.emplace_back(value, run_length);
        }

        in.close();
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
//            cerr << "size: " << size << endl;


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

            // print
            cerr << "Decoded offset: " << decoded_offset << endl;
            cerr << "Decoded flag: " << decoded_flag << endl;
            cerr << "Decoded length: " << static_cast<int>(decoded_length) << endl;
            cerr << "Decoded node ID: " << decoded_node_id << endl;

            number_of_runs--;
            unique_positions.insert(
                    gbwtgraph::Position::encode(pos_t(decoded_node_id, decoded_offset, decoded_flag)).value);
        }

    }

    // this function just store the encoded_runs vector in the output stream
    void TagArray::store_blocks(std::ostream &out) {
        size_t size = encoded_runs.size();
        out.write(reinterpret_cast<const char *>(encoded_runs.data()), size * sizeof(gbwt::byte_type));



        cerr << "Finished storing blocks" << endl;

    }

    void TagArray::store_blocks_sdsl(std::string filename) {

        sdsl::int_vector_buffer<8> out(filename, std::ios::out);
        for (gbwt::size_type value : encoded_runs) {
            gbwt::ByteCode::write(out, value);
        }
        out.close();


        cerr << "Finished storing blocks" << endl;
    }

//    void TagArray::store_blocks_sdsl_run_by_run(std::string filename, const std::vector<std::pair<pos_t, uint8_t>>& tag_runs) {
//
//
//        // TODO: Implement this function
//        sdsl::int_vector_buffer<8> out(filename, std::ios::out);
//        for (gbwt::size_type value : encoded_runs) {
//            gbwt::ByteCode::write(out, value);
//        }
//        out.close();
//
//    }


    std::pair<pos_t, uint8_t> TagArray::decode_run(gbwt::size_type decc) {
        size_t decoded_offset = decc & ((1LL << 10) - 1);
        bool decoded_flag = (decc >> 10) & 0x1;
        uint8_t decoded_length = (decc >> 11) & 0xFF;
        int64_t decoded_node_id = (decc >> 19);

        pos_t graph_pos = make_pos_t(decoded_node_id, decoded_flag, decoded_offset);
        std::pair<pos_t, uint8_t> result = std::make_pair(graph_pos, decoded_length);

        return result;
    }



    


    // Function to load a run starting from a specified position in the file and return the start of the next run
    std::pair<pos_t, uint8_t> TagArray::load_block_at(std::istream &in, size_t &next_block_start) {
        // Seek to the start position of the current run
        in.seekg(next_block_start, std::ios::beg);
        if (in.fail()) {
            cerr << "Failed to seek to run start position" << endl;
            throw std::runtime_error("Failed to seek to run start position");
        }

        std::vector<gbwt::byte_type> current_run;
        gbwt::byte_type byte;

        size_t run_start_position = next_block_start;
        size_t run_end_position = run_start_position;

        // Read bytes until we reach the end of the current run
        while (in.read(reinterpret_cast<char*>(&byte), sizeof(gbwt::byte_type))) {
            current_run.push_back(byte);
            ++run_end_position;
            if (!(byte & 0x80)) { // Check if the last bit is 1, marking the end of the run
                break;
            }
        }

        // Update next_block_start for the next call
        next_block_start = run_end_position;

        if (run_start_position == run_end_position) {
            std::cerr << "Reached the end of file " << std::endl;
            // returning an empty pair
            return std::make_pair(pos_t{0, 0, 0}, 0);
        }

//        cerr << "Loaded run from position " << run_start_position << " to " << run_end_position << endl;

        // Decode the current run data using gbwt::ByteCode::read
        gbwt::size_type decc;
        std::uint64_t bit_location = 0; // Starting bit location in current_run

        decc = gbwt::ByteCode::read(current_run, bit_location); // Decode at bit location

        // Extract values from the decoded byte
        size_t decoded_offset = decc & ((1LL << 10) - 1);
        bool decoded_flag = (decc >> 10) & 0x1;
        uint8_t decoded_length = (decc >> 11) & 0xFF;
        int64_t decoded_node_id = (decc >> 19);

        pos_t graph_pos = make_pos_t(decoded_node_id, decoded_flag, decoded_offset);
        std::pair<pos_t, uint8_t> result = std::make_pair(graph_pos, decoded_length);


//         Output the decoded information for debugging
//        cerr << "Decoded offset: " << decoded_offset << endl;
//        cerr << "Decoded flag: " << decoded_flag << endl;
//        cerr << "Decoded length: " << static_cast<int>(decoded_length) << endl;
//        cerr << "Decoded node ID: " << decoded_node_id << endl;


        if (result.second == 0) {
            cerr << "Decoded length is 0" << endl;
        }
        return result;
    }



}



