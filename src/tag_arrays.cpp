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


//    TagArray::TagArray(int length_bits) : length_bits(length_bits) {
//        encoded_runs = std::vector<gbwt::byte_type>();
//        encoded_runs_starts = sdsl::bit_vector();
//        bwt_intervals = sdsl::sd_vector<>();
//    }


    TagArray::TagArray() {
    }

    gbwt::size_type TagArray::encode_run_length(size_t offset, bool is_rev, uint16_t length, int64_t node_id){
        gbwt::size_type encoded = (offset & 0x3FF) |  // 10 bits for offset
                           ((is_rev & 0x1) << 10) | // 1 bit for reverse flag
                           ((length & ((1LL << length_bits) - 1)) << 11) | // Dynamically allocated bits for length
                           (static_cast<uint64_t>(node_id) << (11 + length_bits)); // Remaining bits for node ID

        return encoded;

    }



    std::pair<pos_t, uint16_t> TagArray::decode_run(gbwt::size_type encoded) {
        size_t decoded_offset = encoded & 0x3FF; // Extract 10-bit offset
        bool decoded_flag = (encoded >> 10) & 0x1; // Extract 1-bit reverse flag
        uint16_t decoded_length = (encoded >> 11) & ((1LL << length_bits) - 1); // Extract length using dynamic bits
        int64_t decoded_node_id = (encoded >> (11 + length_bits)); // Extract node ID


        pos_t graph_pos = make_pos_t(decoded_node_id, decoded_flag, decoded_offset);
        std::pair<pos_t, uint16_t> result = std::make_pair(graph_pos, decoded_length);

        return result;
    }


    void TagArray::load_bptree_lite(BplusTree <Run> &bptree){
        // have to iterate over the bptree and add the runs to the encoded_runs vector
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
                    int total_length = next_item.start_position - current_item.start_position;
                    int max_tag_len = 1 << length_bits;
                    while (total_length >= max_tag_len){
                        gbwt::size_type encoded = encode_run_length(offset(current_pos), is_rev(current_pos), max_tag_len - 1, id(current_pos));
                        gbwt::ByteCode::write(encoded_runs, encoded);
                        total_length -= (max_tag_len - 1);
                    }

                    if (total_length > 0){
                        gbwt::size_type encoded = encode_run_length(offset(current_pos), is_rev(current_pos), total_length, id(current_pos));
                        gbwt::ByteCode::write(encoded_runs, encoded);
                    }


                }
            }
        }
    }

    void TagArray::serialize_run_by_run_batch(sdsl::int_vector_buffer<8>& out, const std::vector<std::pair<gbwtgraph::Position, uint16_t>>& tag_runs){

        for (auto it = tag_runs.begin(); it != tag_runs.end(); ++it) {
            std::vector<gbwt::byte_type> encoded_runs_temp;

            // adding the pair (graph_value, length) as ByteCode to the encoded_runs vector
            pos_t current_pos = it->first.decode();
            uint16_t total_length = it->second;
            int max_tag_len = 1 << length_bits;
            while (total_length >= max_tag_len){
                gbwt::size_type encoded = encode_run_length(offset(current_pos), is_rev(current_pos), max_tag_len - 1, id(current_pos));
                gbwt::ByteCode::write(out, encoded);
//                        out.push_back(encoded);
                total_length -= (max_tag_len - 1);
            }

            if (total_length > 0){
                gbwt::size_type encoded = encode_run_length(offset(current_pos), is_rev(current_pos), total_length, id(current_pos));
//                        out.push_back(encoded);
                gbwt::ByteCode::write(out, encoded);
            }

        }
    }

    void TagArray::serialize_bptree_lite(std::string filename, BplusTree <Run> &bptree){

        sdsl::int_vector_buffer<8> out(filename, std::ios::out | std::ios::trunc);

        // have to iterate over the bptree and add the runs to the encoded_runs vector
        for (auto it = bptree.begin(); it != bptree.end(); ++it) {
            std::vector<gbwt::byte_type> encoded_runs_temp;
            Run current_item = *it;
            if (current_item.graph_position.value != 0) {
                // we have an actual run with graph position
                auto next_it = it;
                ++next_it;
                if (next_it != bptree.end()) {
                    Run next_item = *next_it;
                    // adding the pair (graph_value, length) as ByteCode to the encoded_runs vector
                    pos_t current_pos = current_item.graph_position.decode();
                    int total_length = next_item.start_position - current_item.start_position;
                    int max_tag_len = 1 << length_bits;
                    while (total_length >= max_tag_len){
                        gbwt::size_type encoded = encode_run_length(offset(current_pos), is_rev(current_pos), max_tag_len - 1, id(current_pos));
                        gbwt::ByteCode::write(out, encoded);
//                        out.push_back(encoded);
                        total_length -= (max_tag_len - 1);
                    }

                    if (total_length > 0){
                        gbwt::size_type encoded = encode_run_length(offset(current_pos), is_rev(current_pos), total_length, id(current_pos));
//                        out.push_back(encoded);
                        gbwt::ByteCode::write(out, encoded);
                    }

                    // serialize the encoded_runs_temp
//                    main_out.write(reinterpret_cast<const char *>(encoded_runs.data()), size * sizeof(gbwt::byte_type));


                }
            }
        }

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
//                    gbwt::size_type encoded =
////                            (gbwtgraph::offset(current_pos)) | (gbwtgraph::is_rev(current_pos) << 10) |
////                            ((next_item.start_position - current_item.start_position) << 11) |
////                            (gbwtgraph::id(current_pos) << 19);
//                    gbwt::size_type encoded =
//                            (gbwtgraph::offset(current_pos)) | (gbwtgraph::is_rev(current_pos) << 10) |
//                            ((next_item.start_position - current_item.start_position) << 11) |  // Now uses 10 bits
//                            (gbwtgraph::id(current_pos) << 21);  // Adjust the shift to 21 instead of 20
                    gbwt::size_type encoded = encode_run_length(offset(current_pos), is_rev(current_pos), next_item.start_position - current_item.start_position, id(current_pos));


                    encoded_runs_starts[encoded_runs.size()] = 1;

//                    cout << encoded_runs.size() << endl;
//                    gbwt::size_type s = encoded_runs.size();
                    gbwt::ByteCode::write(encoded_runs, encoded);
                    // setting the bit in the run_starts bit vector
                    encoded_runs_current_size = encoded_runs.size();


//                    cout << encoded_runs.size() << endl;


//                    cout << "bwt " << current_item.start_position << endl;
                }
                // setting the bit in the bwt_intervals bit vector
                builder.set(current_item.start_position);
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

    void TagArray::serialize_run_by_run(std::ofstream& out, const std::vector<std::pair<pos_t, uint16_t>>& tag_runs) {
        for (const auto& [value, run_length] : tag_runs) {
            out.write(reinterpret_cast<const char*>(&value), sizeof(pos_t));
            out.write(reinterpret_cast<const char*>(&run_length), sizeof(uint16_t));
        }
    }


    void TagArray::deserialize_run_by_run(const std::string& file_name) {
        std::ifstream in(file_name, std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "Error: Cannot open file for reading.\n";
        }

        pos_t value;
        uint16_t run_length;
//        int i = 0;

        while (in.read(reinterpret_cast<char*>(&value), sizeof(pos_t))) {
            in.read(reinterpret_cast<char*>(&run_length), sizeof(uint16_t));
            // print the value and run_length
//            std::cerr << i << " " << value << " " << int(run_length) << std::endl;
//            i++;
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

//            size_t decoded_offset = decc & ((1LL << 10) - 1);
//            bool decoded_flag = (decc >> 10) & 0x1;
//            uint16_t decoded_length = (decc >> 11) & 0xFF;
//            int64_t decoded_node_id = (decc >> 19);

            auto decoded = decode_run(decc);

//            size_t decoded_offset = decc & ((1LL << 10) - 1);
//            bool decoded_flag = (decc >> 10) & 0x1;
//            uint16_t decoded_length = (decc >> 11) & 0x3FF;  // Extract 10 bits (0x3FF = 1023)
//            int64_t decoded_node_id = (decc >> 21);  // Adjust shift to match encoding
//
//
//
//            // print
//            cerr << "Decoded offset: " << decoded_offset << endl;
//            cerr << "Decoded flag: " << decoded_flag << endl;
//            cerr << "Decoded length: " << static_cast<int>(decoded_length) << endl;
//            cerr << "Decoded node ID: " << decoded_node_id << endl;

            number_of_runs--;
            unique_positions.insert(
                    gbwtgraph::Position::encode(decoded.first).value);
        }

    }

    // this function just store the encoded_runs vector in the output stream
    void TagArray::store_blocks(std::ostream &out) {
        size_t size = encoded_runs.size();
        out.write(reinterpret_cast<const char *>(encoded_runs.data()), size * sizeof(gbwt::byte_type));



        cerr << "Finished storing blocks" << endl;

    }

    void TagArray::store_blocks_sdsl(std::string filename) {


//        std::cerr << "encoded runs size " << encoded_runs.size() << std::endl;
        sdsl::int_vector_buffer<8> out(filename, std::ios::out | std::ios::trunc);
        for (gbwt::size_type value : encoded_runs) {
            out.push_back(value);
//            std::cerr << "value: " << value << std::endl;
//            gbwt::ByteCode::write(out, value);


        }
//        out.close();


        cerr << "Finished storing blocks" << endl;
    }





    // Function to load a run starting from a specified position in the file and return the start of the next run
    std::pair<pos_t, uint16_t> TagArray::load_block_at(std::istream &in, size_t &next_block_start) {
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
//        size_t decoded_offset = decc & ((1LL << 10) - 1);
//        bool decoded_flag = (decc >> 10) & 0x1;
//        uint16_t decoded_length = (decc >> 11) & 0x3FF;  // Extract 10 bits (0x3FF = 1023)
//        int64_t decoded_node_id = (decc >> 21);  // Adjust shift to match encoding

        auto result = decode_run(decc);


//        pos_t graph_pos = make_pos_t(decoded_node_id, decoded_flag, decoded_offset);
//        std::pair<pos_t, uint16_t> result = std::make_pair(graph_pos, decoded_length);




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

    // ============================================= compressed implementation =============================================

    void TagArray::merge_compressed_files(const std::string filename, const std::string encoded_starts_file, const std::string bwt_intervals_file){


        // merging the encoded_starts file with the main index file
        std::ifstream in_encoded_starts(encoded_starts_file, std::ios::binary);
        size_t start_pos_read;
        sdsl::sd_vector_builder builder(this->start_pos + 1, this->encoded_start_ones);

        while (in_encoded_starts.read(reinterpret_cast<char*>(&start_pos_read), sizeof(start_pos_read))){
            //        std::cerr << "Read start pos " << start_pos_read << std::endl;
            builder.set(start_pos_read);
        }
        in_encoded_starts.close();


        this->encoded_runs_starts_sd = sdsl::sd_vector<>(builder);
        // print some stats from the sd vector
        std::cerr << "The size of the encoded_runs_starts vector is: " << this->encoded_runs_starts_sd.size() << std::endl;
        std::cerr << "Number of 1s in the encoded_runs_starts vector is: " << this->encoded_start_ones << std::endl;

        std::ofstream main_out(filename, std::ios::binary | std::ios::app);
        if (!main_out.is_open()) {
            std::cerr << "Error: Cannot open file for writing.\n";
        }
        sdsl::serialize(this->encoded_runs_starts_sd, main_out);
        // delete the encoded_starts file
        std::remove(encoded_starts_file.c_str());


        // merging the bwt_intervals file with the main index file
        std::ifstream in_bwt_intervals(bwt_intervals_file, std::ios::binary);
        size_t bwt_pos_read;

        std::cerr << "Building the bwt_intervals vector with size " << this->cumulative_run_bwt_position + 1 << " and the number of 1s is " << this->remaining_run_to_write_start << std::endl;

        sdsl::sd_vector_builder builder_bwt(this->cumulative_run_bwt_position + 1, this->remaining_run_to_write_start);

        while (in_bwt_intervals.read(reinterpret_cast<char*>(&bwt_pos_read), sizeof(bwt_pos_read))){
            builder_bwt.set(bwt_pos_read);
        }
        in_bwt_intervals.close();

        this->bwt_intervals = sdsl::sd_vector<>(builder_bwt);

        std::cerr << "The size of the bwt_intervals vector is: " << this->bwt_intervals.size() << std::endl;
        std::cerr << "Number of 1s in the bwt_intervals vector is: " << this->remaining_run_to_write_start << std::endl;

        sdsl::serialize(this->bwt_intervals, main_out);
        std::remove(bwt_intervals_file.c_str());

        main_out.close();


        std::cerr << "The number of encoded runs vector size is " << this->cumulative_starts << std::endl;

        std::fstream out(filename, std::ios::in | std::ios::out | std::ios::binary);
        // write the size of the encoded_runs vector at the beginning of the file
        out.flush();
        out.seekp(0, std::ios::beg);
        out.write(reinterpret_cast<const char*>(&this->cumulative_starts), sizeof(size_t));

        out.close();
    }

    void TagArray::compressed_serialize(std::ostream &main_out, std::ostream &encoded_starts_file, std::ostream &bwt_intervals_file, std::vector<std::pair<pos_t, uint16_t>> &tag_runs){
        std::vector<gbwt::byte_type> encoded_runs;
        if (tag_runs.size() > 0) {
            for (auto [value, run_length] : tag_runs) {

                int max_tag_len = 1 << length_bits;
                // have to write run_length/(max_tag_len - 1) of the run (value, max_tag_len)

                while (run_length >= max_tag_len){

                    bwt_intervals_file.write(reinterpret_cast<const char*>(&this->cumulative_run_bwt_position), sizeof(this->cumulative_run_bwt_position));
                    this->cumulative_run_bwt_position += (max_tag_len - 1);


                    // keeping the data for creating the encoded starts
                    if (this->remaining_run_to_write_start % this->encoded_start_every_k_run == 0){
                        this->start_pos = this->cumulative_starts + encoded_runs.size();
                        // write the encoded start in file
                        encoded_starts_file.write(reinterpret_cast<const char*>(&this->start_pos), sizeof(this->start_pos));
                        this->encoded_start_ones++;
                    }

//                gbwt::size_type encoded =
//                        (gbwtgraph::offset(value)) | (gbwtgraph::is_rev(value) << 10) |
//                        (run_length << 11) |  // Now uses 10 bits
//                        (gbwtgraph::id(value) << 21);  // Adjust the shift to 21 instead of 20


                    gbwt::size_type encoded = encode_run_length(offset(value), is_rev(value), (max_tag_len - 1), id(value));


                    gbwt::ByteCode::write(encoded_runs, encoded);

                    this->remaining_run_to_write_start++;

                    run_length -= (max_tag_len - 1);

                }

                if (run_length > 0){
                    bwt_intervals_file.write(reinterpret_cast<const char*>(&this->cumulative_run_bwt_position), sizeof(this->cumulative_run_bwt_position));
                    this->cumulative_run_bwt_position += run_length;


                    // keeping the data for creating the encoded starts
                    if (this->remaining_run_to_write_start % this->encoded_start_every_k_run == 0){
                        this->start_pos = this->cumulative_starts + encoded_runs.size();
                        // write the encoded start in file
                        encoded_starts_file.write(reinterpret_cast<const char*>(&this->start_pos), sizeof(this->start_pos));
                        this->encoded_start_ones++;
                    }

//                gbwt::size_type encoded =
//                        (gbwtgraph::offset(value)) | (gbwtgraph::is_rev(value) << 10) |
//                        (run_length << 11) |  // Now uses 10 bits
//                        (gbwtgraph::id(value) << 21);  // Adjust the shift to 21 instead of 20


                    gbwt::size_type encoded = encode_run_length(offset(value), is_rev(value), run_length, id(value));


                    gbwt::ByteCode::write(encoded_runs, encoded);

                    this->remaining_run_to_write_start++;

                }

            }

            size_t size = encoded_runs.size();
//            out.write(reinterpret_cast<const char *>(&size), sizeof(size));
            main_out.write(reinterpret_cast<const char *>(encoded_runs.data()), size * sizeof(gbwt::byte_type));
            this->cumulative_starts += size;

        } else {
            std::cerr << "No tags to compress and write" << std::endl;
        }

    }


    void TagArray::load_compressed_tags(std::istream &in) {
        // Read the size of the encoded_runs vector
        size_t encoded_runs_size;
        in.read(reinterpret_cast<char*>(&encoded_runs_size), sizeof(size_t));

        std::cerr << "Loading encoded runs of size: " << encoded_runs_size << std::endl;

        // Read the encoded runs vector
        this->encoded_runs.resize(encoded_runs_size);
        in.read(reinterpret_cast<char*>(this->encoded_runs.data()), encoded_runs_size * sizeof(gbwt::byte_type));

        // Read the encoded_runs_starts_sd vector
        this->encoded_runs_starts_sd.load(in);

        // Read the bwt_intervals vector
        this->bwt_intervals.load(in);

        std::cerr << "Loaded encoded runs_starts_sd of size: " << this->encoded_runs_starts_sd.size() << std::endl;
        std::cerr << "Loaded bwt_intervals of size: " << this->bwt_intervals.size() << std::endl;

        this->encoded_runs_sd_starts_select = sdsl::sd_vector<>::select_1_type(&this->encoded_runs_starts_sd);
        this->bwt_intervals_rank = sdsl::sd_vector<>::rank_1_type(&this->bwt_intervals);

        std::cerr << "Select and rank data structures initialized" << std::endl;
    }



    void TagArray::query_compressed(size_t start, size_t end, size_t &number_of_runs) {


        // first have to find ranks of start and end in the bwt_intervals which are number of 1's less than start and end

//        sdsl::sd_vector<>::rank_1_type bwt_intervals_rank(&bwt_intervals);
//        std::cerr << "Finding the tags between " << start << " and " << end << std::endl;
        size_t first_bit_index = this->bwt_intervals_rank(start + 1);
//        if (start > 0 && bwt_intervals[start] == 1) {
//            // if the start is 1, then the tags start from start position and we don't need the last 1 before start
//            first_bit_index++;
//        }
//        std::cerr << "First bit index: " << first_bit_index << std::endl;

        // TODO: can handle this without the rank data structure
        size_t end_bit_index = this->bwt_intervals_rank(end + 1);

//        std::cerr << "End bit index: " << end_bit_index << std::endl;



        auto run_nums = end_bit_index - first_bit_index + 1;
        number_of_runs = run_nums;


        std::vector<std::uint64_t> unique_positions;
        unique_positions.reserve(run_nums);

        size_t current_tag_run_index = first_bit_index - (first_bit_index % this->encoded_start_every_k_run);
        size_t move_tags = first_bit_index % this->encoded_start_every_k_run;

        // where the first run starts
        // this is the actual bit location of latest % encoded_start_every_k_run == 0 in the encoded_runs vector
        std::uint64_t bit_location = this->encoded_runs_sd_starts_select(current_tag_run_index / this->encoded_start_every_k_run + 1);
        gbwt::size_type decc;


        while (move_tags > 1){
            // the read function changes the bit_location to the next bit location
            decc = gbwt::ByteCode::read(this->encoded_runs, bit_location);
            move_tags--;
        }


        while (run_nums > 0) {
            // the read function changes the bit_location to the next bit location
            decc = gbwt::ByteCode::read(this->encoded_runs, bit_location);

//            size_t decoded_offset = decc & ((1LL << 10) - 1);
//            bool decoded_flag = (decc >> 10) & 0x1;
//            uint16_t decoded_length = (decc >> 11) & 0x3FF;  // Extract 10 bits (0x3FF = 1023)
//            int64_t decoded_node_id = (decc >> 21);  // Adjust shift to match encoding

            auto decoded = decode_run(decc);
            // print
//            cerr << "Decoded offset: " << decoded_offset << endl;
//            cerr << "Decoded flag: " << decoded_flag << endl;
//            cerr << "Decoded length: " << static_cast<int>(decoded_length) << endl;
//            cerr << "Decoded node ID: " << decoded_node_id << endl;

            run_nums--;
            unique_positions.push_back(gbwtgraph::Position::encode(decoded.first).value);
//            unique_positions.insert(
//                    gbwtgraph::Position::encode(pos_t(decoded_node_id, decoded_offset, decoded_flag)).value);
        }
        std::sort(unique_positions.begin(), unique_positions.end());
        unique_positions.erase(std::unique(unique_positions.begin(), unique_positions.end()), unique_positions.end());

        std::cout << "Number of unique positions: " << unique_positions.size() << std::endl;
        // print the unique positions
        for (auto pos : unique_positions) {
            std::cout << pos << std::endl;
        }

    }



}



