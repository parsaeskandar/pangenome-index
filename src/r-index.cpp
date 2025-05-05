//
// Created by seeskand on 9/21/24.
// This implementation is based on the fast_locate implementation in GBWT from https://github.com/jltsiren/gbwt.git
//

#include "../include/pangenome_index/r-index.hpp"


//TODO: MAKE THE copy, swap, ... function with the new variables

namespace panindexer {
    /*
  Copyright (c) 2020, 2021, 2022 Jouni Siren

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

//------------------------------------------------------------------------------

// Numerical class constants.

    constexpr std::uint32_t
    FastLocate::Header::TAG;
    constexpr std::uint32_t
    FastLocate::Header::VERSION;

    constexpr size_type
    FastLocate::NO_POSITION;


//------------------------------------------------------------------------------

// Other class variables.

    const std::string FastLocate::EXTENSION = ".ri";

//------------------------------------------------------------------------------

    FastLocate::Header::Header() :
            tag(TAG), version(VERSION),
            max_length(1),
            flags(0) {
    }

    size_type
    FastLocate::Header::serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += sdsl::write_member(this->tag, out, child, "tag");
        written_bytes += sdsl::write_member(this->version, out, child, "version");
        written_bytes += sdsl::write_member(this->max_length, out, child, "max_length");
        written_bytes += sdsl::write_member(this->flags, out, child, "flags");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void
    FastLocate::Header::load(std::istream &in) {
        sdsl::read_member(this->tag, in);
        sdsl::read_member(this->version, in);
        sdsl::read_member(this->max_length, in);
        sdsl::read_member(this->flags, in);
    }

    void
    FastLocate::Header::check() const {
        if (this->tag != TAG) {
            throw sdsl::simple_sds::InvalidData("FastLocate: Invalid tag");
        }

        if (this->version != VERSION) {
            std::string msg =
                    "FastLocate: Expected v" + std::to_string(VERSION) + ", got v" + std::to_string(this->version);
            throw sdsl::simple_sds::InvalidData(msg);
        }

        std::uint64_t mask = 0;
        switch (this->version) {
            case VERSION:
                mask = 0;
                break;
        }
        if ((this->flags & mask) != this->flags) {
            throw sdsl::simple_sds::InvalidData("FastLocate: Invalid flags");
        }
    }

//------------------------------------------------------------------------------

    FastLocate::FastLocate() :
            buff_reader(nullptr) {
        this->initialize_complement_table();
    }

    FastLocate::FastLocate(const FastLocate &source) {
        this->copy(source);
    }

    FastLocate::FastLocate(FastLocate &&source) {
        *this = std::move(source);
    }

    FastLocate::~FastLocate() {
    }

    void
    FastLocate::swap(FastLocate &another) {
        if (this != &another) {
            std::swap(this->buff_reader, another.buff_reader);
            std::swap(this->header, another.header);
            this->samples.swap(another.samples);
            this->last.swap(another.last);
            this->last_to_run.swap(another.last_to_run);
        }
    }

    FastLocate &
    FastLocate::operator=(const FastLocate &source) {
        if (this != &source) { this->copy(source); }
        return *this;
    }

    FastLocate &
    FastLocate::operator=(FastLocate &&source) {
        if (this != &source) {
            this->buff_reader = source.buff_reader;
            this->header = std::move(source.header);
            this->samples = std::move(source.samples);
            this->last = std::move(source.last);
            this->last_to_run = std::move(source.last_to_run);
        }
        return *this;
    }

    size_type
    FastLocate::serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += this->header.serialize(out, child, "header");
        written_bytes += this->samples.serialize(out, child, "samples");
        written_bytes += this->last.serialize(out, child, "last");
        written_bytes += this->last_to_run.serialize(out, child, "last_to_run");


        written_bytes += this->sym_map.serialize(out, child, "sym_map");
        written_bytes += this->C.serialize(out, child, "C");

        written_bytes += this->blocks_start_pos.serialize(out, child, "blocks_start_pos");
        written_bytes += sdsl::write_member(sequence_size, out, child, "sequence_size");


        // write the number of blocks
        size_t blocks_size = this->blocks.size();
        written_bytes += sdsl::write_member(blocks_size, out, child, "blocks_size");

        for (const auto &block : blocks) {
            written_bytes += block.serialize(out, child, "block");
        }

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void
    FastLocate::load(std::istream &in) {
        this->header.load(in);
        this->header.check();
        this->header.setVersion(); // Update to the current version.


        this->samples.load(in);
        this->last.load(in);
        this->last_to_run.load(in);

        this->sym_map.load(in);
        this->C.load(in);

        sdsl::load(this->blocks_start_pos, in);

        sdsl::read_member(this->sequence_size, in);

        size_t blocks_size;
        sdsl::read_member(blocks_size, in);
        this->blocks.resize(blocks_size);

        for (size_t i = 0; i < blocks_size; ++i) {
            this->blocks[i].load(in);
        }

    }

    void
    FastLocate::copy(const FastLocate &source) {
        this->header = source.header;
        this->samples = source.samples;
        this->last = source.last;
        this->last_to_run = source.last_to_run;
    }

//------------------------------------------------------------------------------

    // This function return the run_id and the run_range of the run that the idx belongs to
    void
    FastLocate::bwt_index_run_id(unsigned long idx, gbwt::range_type &run, size_type &run_id) {
        size_type offset = 0;
        size_type runs_seen = 0;


        for (size_t i = 0; i < this->buff_reader->size(); i++) {
            size_t sym, freq;
            this->buff_reader->read_run(i, sym, freq);
            offset += freq;
            runs_seen += (sym == NENDMARKER ? freq : 1); // each endmarker is a separate run
            if (offset > idx) {
                if (sym == NENDMARKER) {
                    run.first = run.second = idx;
                    run_id = runs_seen - (offset - idx);
                } else {
                    run.first = offset - freq;
                    run.second = offset - 1;
                    run_id = runs_seen - 1;
                }
                break;
            }
        }
    }


    size_t FastLocate::bwt_char_at(size_t idx) {
        auto iter = this->blocks_start_pos.predecessor(idx);
        auto res = this->blocks[iter->first].bwt_char_at(idx - iter->second);
        return res;
    }

    // This function provide backward navigation in the BWT
    std::pair <size_t, size_t> FastLocate::psi(size_t idx) {
        size_t symbol = this->bwt_char_at(idx);
        return {symbol, this->C[this->sym_map[symbol]] + this->rankAt(idx, symbol)};
    }

    std::pair <size_t, size_t> FastLocate::psi_and_run_id(size_t idx, size_t &run_id, size_t &current_position) {
        run_id = 0;
        current_position = 0;
        size_t symbol = this->bwt_char_at(idx);
        size_t next_pos = this->C[this->sym_map[symbol]] + this->rankAt(idx, symbol, run_id, current_position);
        return {symbol, next_pos};
    }


    // current_position Tracks how far we've gone through the BWT
    // The number does not include the index symbol
    size_t FastLocate::rankAt(size_t pos, size_t symbol, size_t &run_id, size_t &current_position) const {

        auto iter = this->blocks_start_pos.predecessor(pos);
        size_t run_num = 0;
        size_t cur_pos = 0;
        auto cum_rank = this->blocks[iter->first].rankAt(pos - iter->second, symbol, run_num, cur_pos);
        run_id = iter->first * this->block_size + run_num;
        current_position = iter->second + cur_pos;
        auto cumulative_rank = cum_rank + this->blocks[iter->first].get_cum_ranks()[this->sym_map[symbol]];
        return cumulative_rank;
    }


    std::vector<size_t> FastLocate::rank_at_cached(size_t pos) const {
        auto iter = this->blocks_start_pos.predecessor(pos);
        std::vector<size_t> rank_vec(this->C.size(), 0);
        for (size_t i = 0; i < this->C.size(); i++) {
            size_t run_num = 0;
            size_t cur_pos = 0;
            auto cum_rank = this->blocks[iter->first].rankAt(pos - iter->second, nuc[i], run_num, cur_pos);
            rank_vec[i] = cum_rank + this->blocks[iter->first].get_cum_ranks()[this->sym_map[nuc[i]]];
        }
        return rank_vec;
    }





    // TODO: optimize for the ranges that are in the same or adjacent blocks
    range_type FastLocate::LF(range_type range, size_t sym, bool &starts_with_to, size_t &first_run) const {


        if(!this->sym_map[sym]) return {1,0};

        if (range.first > range.second) {
            return {1, 0};
        }
        starts_with_to = false;
        first_run = std::numeric_limits<size_type>::max();

        // Initialize run tracking variables
        size_t run_id = 0;                 // Track the current run id
        size_t current_position = 0;       // Track the current BWT position

        // Rank for the start of the range
        size_t start_offset = range.first;
        bool first_run_found = false;      // Track whether we've found the first run

        range.first = this->rankAt(range.first, sym, run_id, current_position);
//        std::cerr << "LF function range first " << range.first << std::endl;

        // Check if the first run contains the symbol and overlaps with range.first
//        if (current_position > start_offset) {
//            first_run_found = true;
//            if (sym == this->buff_reader->read_sym(run_id)) {  // Compare the symbol of the current run
//                starts_with_to = true;  // The range starts with the symbol
//                first_run = run_id; // Set the first run id
//            }
//        }

        // Rank for the end of the range
//        while (run_id < this->buff_reader->size() && current_position <= range.second) {
//            size_t symbol_in_run = 0;
//            size_t freq_in_run = 0;
//            this->buff_reader->read_run(run_id, symbol_in_run, freq_in_run);  // Read the run symbol and its frequency
//
//            ++run_id;  // Move to the next run
//            current_position += freq_in_run;  // Increment the BWT position
//
//            if (symbol_in_run == sym && first_run == std::numeric_limits<size_type>::max()) {
//                if (!first_run_found) {
//                    starts_with_to = true;
//                }
//                first_run = run_id - 1;  // Set the first run if not found
//            }
//            first_run_found = true;
//        }

        // Rank the second part of the range
        run_id = 0;
        current_position = 0;


        size_t sym_inside = this->rankAt(range.second + 1, sym, run_id, current_position) - range.first;

        if (sym_inside == 0) {
            return {1, 0};
        }
        range.first += this->C[this->sym_map[sym]];
        range.second = range.first + sym_inside - 1;


        return range;  // Return the updated range
    }


    // This function returns the exact number of runs, considering that each NENDMARKER is a separate run
    size_type FastLocate::total_runs() {
        size_t runs = 0;
        for (size_t i = 0; i < this->buff_reader->size(); i++) {
            size_t sym, freq;
            this->buff_reader->read_run(i, sym, freq);
            runs += (sym == NENDMARKER ? freq : 1);
        }
        return runs;
    }

    FastLocate::FastLocate(std::string source) :
            grlbwt_file(source)
//                buff_reader(&source)
    {
        this->initialize_complement_table();
        double start = readTimer();


        this->buff_reader = new bwt_buff_reader(source);
        this->calculate_C();

        if (this->buff_reader->size() == 0) {
            if (Verbosity::level >= Verbosity::FULL) {
                std::cerr << "FastLocate::FastLocate(): The input grlBWT is empty" << std::endl;
            }
            return;
        }

        // Determine the number of logical runs before each record.
        size_type total_runs = this->total_runs();
        auto n_seq = this->tot_strings();

        std::cerr << total_runs / this->block_size << std::endl;
        this->blocks.resize((total_runs / this->block_size) + 1);


        size_t run_iterator = 0;

        size_t run_nums = 0;
        size_t current_block_id = 0;
//        size_t start_buff = 0;
        std::vector <size_t> cumulative_freq;
        cumulative_freq.resize(this->C.size(), 0);
        std::cerr << cumulative_freq.size() << std::endl;

        size_t start_offset = 0;

        std::vector <std::pair<size_t, size_t>> run_buff;

        sdsl::sd_vector_builder block_sd_builder(this->bwt_size(), (total_runs / this->block_size) +
                                                                   (total_runs % this->block_size != 0));
        std::cerr << (total_runs / this->block_size) + (total_runs % this->block_size != 0) << std::endl;

        if (Verbosity::level >= Verbosity::FULL) {
            std::cerr << "FastLocate::FastLocate(): calculating the BWT information" << std::endl;
        }
        while (run_iterator < this->buff_reader->size()) {

            if (run_nums == 0) {
                // adding the cumulative freq of the previous runs to the current block

                this->blocks[current_block_id].set_character_cum_ranks(cumulative_freq);
                block_sd_builder.set_unsafe(start_offset);
            }

            size_t sym, freq;
            this->buff_reader->read_run(run_iterator, sym, freq);

//            std::cerr << "The run is: " << sym << " " << freq << std::endl;


            // handling each endmarker as a separate run
            if (sym == NENDMARKER) {
                // in this case we just add the endmarker run to the block
                if (run_nums + freq < this->block_size) {
                    run_nums += freq;
                    cumulative_freq[this->sym_map[sym]] += freq;
                    start_offset += freq;
                    for (size_t i = 0; i < freq; i++) {
                        run_buff.push_back({sym, 1});
                    }
                } else {
                    // in this case we got to add the endmarkers to the next block too
                    // adding the endmarkers to the current block
//                    std::cerr << "current freq of endmarker " << freq << std::endl;
//                    std::cerr << this->block_size - run_nums << " extra runs in the ENDMARKER case" << std::endl;
                    for (size_t i = 0; i < (this->block_size - run_nums); i++) {
//                        std::cerr << "adding the endmarker to the current block" << std::endl;
                        cumulative_freq[this->sym_map[sym]] += 1;
                        start_offset += 1;
                        run_buff.push_back({sym, 1});
                    }

                    auto remaining_freq = freq - (this->block_size - run_nums);

//                    std::cerr << "remaining freq " << remaining_freq << std::endl;

                    run_nums += this->block_size - run_nums;
                    assert(run_nums == this->block_size);
                    // now the current block is full
                    this->blocks[current_block_id].set_runs(run_buff);


                    auto block_fill = remaining_freq / this->block_size;

//                    std::cerr << "block fill " << block_fill << std::endl;
                    for (size_t i = 0; i < block_fill + 1; i++) {
                        run_buff.clear();
                        current_block_id++;
                        run_nums = 0;
                        if (i == block_fill && remaining_freq % this->block_size == 0) {
                            break;
                        }
                        block_sd_builder.set_unsafe(start_offset);
                        this->blocks[current_block_id].set_character_cum_ranks(cumulative_freq);

                        if (i == block_fill) {
                            for (size_t j = 0; j < remaining_freq % this->block_size; j++) {
                                cumulative_freq[this->sym_map[sym]] += 1;
                                start_offset++;
                                run_buff.push_back({sym, 1});
                                run_nums += 1;
                            }
                            assert(run_nums == remaining_freq % this->block_size);
                            assert(run_nums < this->block_size);


                        } else {

                            for (size_t j = 0; j < this->block_size; j++) {
                                cumulative_freq[this->sym_map[sym]] += 1;
                                start_offset++;
                                run_buff.push_back({sym, 1});
                                run_nums += 1;
                            }
                            this->blocks[current_block_id].set_runs(run_buff);
                        }

                    }



//                    current_block_id++; // moving to the next block
//                    run_buff.clear();
//
//
//                    run_nums = 0;
//
//                    block_sd_builder.set_unsafe(start_offset);
//                    this->blocks[current_block_id].set_character_cum_ranks(cumulative_freq);
//
//                    // now have to handle the remaining runs of the endmarker
//
//                    for (size_t i = 0; i < remaining_freq; i++){
//                        cumulative_freq[this->sym_map[sym]] += 1;
//                        start_offset++;
//                        run_buff.push_back({sym, 1});
//                        run_nums += 1;
//                    }

                }

            } else { // when the run is not and endmarker
                cumulative_freq[this->sym_map[sym]] += freq;
                run_buff.push_back({sym, freq});
                run_nums++;
                start_offset += freq;


            }




            // the block is full
            if (run_nums == this->block_size) {
                this->blocks[current_block_id].set_runs(run_buff);

                run_buff.clear();
                current_block_id++;
                run_nums = 0;
            }

            run_iterator++;

        }

        // have to check if the last block is added
        if (run_nums != 0) {
            this->blocks[current_block_id].set_runs(run_buff);
        }

        if (Verbosity::level >= Verbosity::FULL) {
            std::cerr << "FastLocate::FastLocate(): " << this->blocks.size() << " blocks of size " << this->block_size
                      << " calculated" << std::endl;
        }

        this->blocks_start_pos = sdsl::sd_vector<>(block_sd_builder);


        if (Verbosity::level >= Verbosity::FULL) {
            std::cerr << "FastLocate::FastLocate(): " << total_runs << " logical runs in the grlBWT" << std::endl;
        }

        // Global sample buffers.
        struct sample_record {
            size_type seq_id, seq_offset, run_id;

            // Sort by text position.
            bool operator<(const sample_record &another) const {
                return (this->seq_id < another.seq_id ||
                        (this->seq_id == another.seq_id && this->seq_offset < another.seq_offset));
            }
        };
        std::vector <sample_record> head_samples, tail_samples;
        head_samples.reserve(total_runs);
        tail_samples.reserve(total_runs);

        // Run identifier for each offset in the endmarker. We cannot get this
        // information efficiently with random access.
        if (Verbosity::level >= Verbosity::FULL) {
            std::cerr << "FastLocate::FastLocate(): Processing the endmarker record" << std::endl;
        }


        std::vector <size_type> endmarker_runs(n_seq, 0);
        {
            size_type run_id = 0;
            auto prev = this->psi(0);

            for (size_type i = 1; i < n_seq; i++) {
                auto curr = this->psi(i);

                if (curr.first == NENDMARKER || curr.first != prev.first) {
                    run_id++;
                    prev = curr;
                }
                endmarker_runs[i] = run_id;
            }

        }

        std::cerr << "endmarker runs size " << endmarker_runs.size() << std::endl;
//        // printing the endmarker runs
//        for (size_t i = 0; i < endmarker_runs.size(); i++) {
//            std::cerr << endmarker_runs[i] << " ";
//        }


        // Extract the samples from each sequence.
        double extract_start = readTimer();
        if (Verbosity::level >= Verbosity::FULL) {
            std::cerr << "FastLocate::FastLocate(): Extracting head/tail samples" << std::endl;
        }

        // running with how many threads
        std::cerr << omp_get_max_threads() << std::endl;
#pragma omp parallel for schedule(dynamic, 1)
        for (size_type i = 0; i < n_seq; i++) {

//            std::cerr << "Starting a sub job for seq" << i << std::endl;


            std::vector <sample_record> head_buffer, tail_buffer;
            size_type seq_offset = 0, run_id = endmarker_runs[i];
            if (i == 0 || run_id != endmarker_runs[i - 1]) {
                head_buffer.push_back({i, seq_offset, run_id});
            }
            if (i + 1 >= n_seq || run_id != endmarker_runs[i + 1]) {
                tail_buffer.push_back({i, seq_offset, run_id});
            }
            auto curr = this->psi(i);
            seq_offset++;
            range_type run(0, 0);

            while (curr.first != NENDMARKER) {


//                auto next1 = this->psi(curr.second); // TODO: make a function that does these two lines at once
//                range_type run1(0, 0);
//                size_type run_id1 = 0;
                size_t temp_run_id = 0;
                size_t tem_cur = 0;
//                this->bwt_index_run_id(curr.second, run1, run_id1);
                auto next = this->psi_and_run_id(curr.second, temp_run_id, tem_cur);
                std::pair tt = this->blocks[(temp_run_id / this->block_size)].get_run(temp_run_id % this->block_size);
                run.first = tem_cur;
                run.second = tem_cur + tt.second - 1;

                run_id = temp_run_id;

//                assert(1 == 2);
//                assert(static_cast<size_t>(run_id1) == temp_run_id);
////                std::cerr << static_cast<size_t>(run1.first) << " " << run.first << std::endl;
//                assert(static_cast<size_t>(run1.first) == run.first);
//                std::cerr << static_cast<size_t>(run1.second) << " " << run.second << std::endl;
//                assert(static_cast<size_t>(run1.second) == run.second);

                // have to find which run the next.first belongs to
                if (curr.second == run.first) {
                    head_buffer.push_back({i, seq_offset, run_id});
                }
                if (curr.second == run.second) {
                    tail_buffer.push_back({i, seq_offset, run_id});
                }
                curr = next;
                seq_offset++;

                if (seq_offset % 100000000 == 0) {
                    std::cerr << "seq " << i << " offset " << seq_offset << std::endl;
                }
            }


            for (sample_record &record: head_buffer) { record.seq_offset = seq_offset - 1 - record.seq_offset; }
            for (sample_record &record: tail_buffer) { record.seq_offset = seq_offset - 1 - record.seq_offset; }

            //
//            std::cerr << "completed a sub job for seq" << i  << " with len " << seq_offset << std::endl;

#pragma omp critical
            {
                this->header.max_length = std::max(this->header.max_length, seq_offset);
                head_samples.insert(head_samples.end(), head_buffer.begin(), head_buffer.end());
                tail_samples.insert(tail_samples.end(), tail_buffer.begin(), tail_buffer.end());
            }
        }

        sdsl::util::clear(endmarker_runs);
        if (Verbosity::level >= Verbosity::BASIC) {
            double seconds = readTimer() - extract_start;
            std::cerr << "FastLocate::FastLocate(): Extracted " << head_samples.size() << " / " << tail_samples.size()
                      << " head/tail samples in " << seconds << " seconds" << std::endl;
        }
//
        // Store the head samples.
        if (Verbosity::level >= Verbosity::FULL) {
            std::cerr << "FastLocate::FastLocate(): Storing the head samples" << std::endl;
        }
        parallelQuickSort(head_samples.begin(), head_samples.end(), [](const sample_record &a, const sample_record &b) {
            return (a.run_id < b.run_id);
        });
        this->samples.width(sdsl::bits::length(this->pack(n_seq - 1, this->header.max_length - 1)));
        this->samples.resize(total_runs);
        for (size_type i = 0; i < total_runs; i++) {
            this->samples[i] = this->pack(head_samples[i].seq_id, head_samples[i].seq_offset);
        }
        sdsl::util::clear(head_samples);

        // Store the tail samples.
        if (Verbosity::level >= Verbosity::FULL) {
            std::cerr << "FastLocate::FastLocate(): Storing the tail samples" << std::endl;
        }
        parallelQuickSort(tail_samples.begin(), tail_samples.end());
        sdsl::sd_vector_builder builder(n_seq * this->header.max_length, total_runs);
        this->last_to_run.width(sdsl::bits::length(total_runs - 1));
        this->last_to_run.resize(total_runs);
        for (size_type i = 0; i < total_runs; i++) {
            builder.set_unsafe(this->pack(tail_samples[i].seq_id, tail_samples[i].seq_offset));
            this->last_to_run[i] = tail_samples[i].run_id;
        }
        sdsl::util::clear(tail_samples);
        this->last = sdsl::sd_vector<>(builder);
//
        if (Verbosity::level >= Verbosity::BASIC) {
            double seconds = readTimer() - start;
            std::cerr << "FastLocate::FastLocate(): Processed " << n_seq << " sequences of total length "
                      << this->get_sequence_size()
                      << " in " << seconds << " seconds" << std::endl;
        }

    }


//------------------------------------------------------------------------------
//
//    range_type FastLocate::find(size_t sym, size_type &first) const {
//        if (this->sym_map[sym] == 0) { return SearchState(); }
//
//        CompressedRecord record = this->index->record(node);
//        if (!(record.empty())) {
//            first = this->getSample(node, 0);
//        }
//        return SearchState(node, 0, record.size() - 1);
//    }

    size_t FastLocate::tot_runs() const {
        return this->block_size * (this->blocks.size() - 1) + this->blocks.back().get_run_nums();
    }

    range_type FastLocate::extend(range_type state, size_t sym, size_type &first) const {

        if (state.second < state.first || (this->sym_map[sym] == 0)) {
            range_type empty(0, -1);
            return empty;
        }

//        CompressedRecord record = this->index->record(state.node);
        bool starts_with_node = false;
        size_t run_id = std::numeric_limits<size_type>::max();
        state = this->LF(state, sym, starts_with_node, run_id);
//    state.range = record.LF(state.range, node, starts_with_node, run_id);
        if (state.first >= state.second) {
            // The position at the start of the resulting range is the successor of the
            // first occurrence of the query node in the query range. We decrement the
            // offset instead of incrementing it, because sequence offset is the distance
            // to the end of the sequence.
            if (starts_with_node) { first++; }
            else {
                first = this->getSample(run_id) - 1;
            }
        }

        return state;
    }



    std::vector <size_type>
    FastLocate::locate(range_type state, size_type first) const {
        std::vector <size_type> result;
        if (state.second < state.first) { return result; }
        result.reserve(state.second - state.first + 1);

        // Find the nearest run start and use the sample there as the first hit,
        // if the caller did not provide it.
        size_type offset_of_first = state.first;
        if (first == NO_POSITION) // TODO: check this
        {
            auto iter = this->blocks_start_pos.predecessor(state.first);
            size_type run_id = 0;
            size_t offset = 0;

            size_t cur_pos = 0;

            auto run_num = this->blocks[iter->first].run_id_at(state.first - iter->second, cur_pos);
            run_id = iter->first * this->block_size + run_num;
            offset_of_first = iter->second + cur_pos;
            first = this->getSample(run_id);

        }



        // Iterate until the start of the range and locate the first occurrence.
        while (offset_of_first < state.first) {
            first = this->locateNext(first);
            offset_of_first++;
        }
        result.push_back(this->seqId(first));

        // Locate the remaining occurrences.
        for (size_type i = state.first + 1; i <= state.second; i++) {
            first = this->locateNext(first);
            result.push_back(this->seqId(first));
        }


        // TODO: change this sort to parallel ones like in gbwt
        std::sort(result.begin(), result.end());
        result.resize(std::unique(result.begin(), result.end()) - result.begin());

        return result;
    }

    std::vector <size_type>
    FastLocate::decompressSA() const {
        std::vector <size_type> result;

        result.reserve(this->get_sequence_size());
        result.push_back(this->locateFirst());
        for (size_type i = 1; i < this->get_sequence_size(); i++) {
            result.push_back(this->locateNext(result.back()));
        }

        return result;
    }

    std::vector <size_type>
    FastLocate::decompressDA() const {
        std::vector <size_type> result = this->decompressSA();
        for (size_type i = 0; i < result.size(); i++) {
            result[i] = this->seqId(result[i]);
        }
        return result;
    }

////------------------------------------------------------------------------------
//
    size_type FastLocate::locateNext(size_type prev) const {
    auto iter = this->last.predecessor(prev);
    return this->samples[this->last_to_run[iter->first] + 1] + (prev - iter->second);
}

//
////------------------------------------------------------------------------------
//
void
printStatistics(const FastLocate &index, const std::string &name) {
    printHeader("Runs");
    std::cout << index.size() << std::endl;
    printHeader("Size");
    std::cout << inMegabytes(sdsl::size_in_bytes(index)) << " MB" << std::endl;
    std::cout << std::endl;
}

std::string
indexType(const FastLocate &) {
    return "R-index";
}

//------------------------------------------------------------------------------
// Functionalities for FMD-index to be able to handle bidirectional BWT

// Backward extension from Algorithm 2 from Li's paper: doi:10.1093/bioinformatics/bts280
FastLocate::bi_interval FastLocate::backward_extend(const bi_interval& bint, size_t a) {
    size_t k = bint.forward;
    size_t k_prime = bint.reverse;
    int64_t s = bint.size;
    int64_t b = 0;

    auto rank_cache_ks = rank_at_cached(k + s);
    auto rank_cache_k = rank_at_cached(k);

    while (nuc[b] < this->complement(a)) {

        k_prime += (rank_cache_ks[this->sym_map[this->complement(nuc[b])]] - rank_cache_k[this->sym_map[this->complement(nuc[b])]]);

//        k_prime += (this->rankAt(k + s , this->complement(nuc[b])) - this->rankAt(k, this->complement(nuc[b])));
        b++;
    }

//    auto rank_ks = this->rankAt(k + s, a);
//    auto rank_k = this->rankAt(k, a);

    auto rank_ks = rank_cache_ks[this->sym_map[a]];
    auto rank_k = rank_cache_k[this->sym_map[a]];

    if (rank_k >= rank_ks) {
        return bi_interval(0, 0, 0);
    }

    s = rank_ks - rank_k;

    k = rank_k + this->C[this->sym_map[a]];
    return bi_interval(k, k_prime, s);
}

//FastLocate::bi_interval FastLocate::backward_extend(const bi_interval& bint, size_t a) {
//    size_t k = bint.forward;
//    size_t k_prime = bint.reverse;
//    int64_t s = bint.size;
//    int64_t b = 0;
//
////    auto rank_cache_ks = rank_at_cached(k + s);
////    auto rank_cache_k = rank_at_cached(k);
//
//    while (nuc[b] < this->complement(a)) {
//
////        k_prime += (rank_cache_ks[b] - rank_cache_k[b]);
//
//        k_prime += (this->rankAt(k + s , this->complement(nuc[b])) - this->rankAt(k, this->complement(nuc[b])));
//        b++;
//    }
//
//    auto rank_ks = this->rankAt(k + s, a);
//    auto rank_k = this->rankAt(k, a);
//
////    auto rank_ks = rank_cache_ks[this->sym_map[a]];
////    auto rank_k = rank_cache_k[this->sym_map[a]];
//
//    if (rank_k >= rank_ks) {
//        return bi_interval(0, 0, 0);
//    }
//
//    s = rank_ks - rank_k;
//
//    k = rank_k + this->C[this->sym_map[a]];
//    return bi_interval(k, k_prime, s);
//}

//FastLocate::bi_interval FastLocate::backward_extend(const bi_interval& bint, size_t a) {
//    size_t k = bint.forward;
//    size_t l = bint.reverse;
//    int64_t s = bint.size;
//    std::vector<char> nuc = {NENDMARKER, 'A', 'C', 'G', 'N', 'T'};
//
//
//    std::vector<size_t> k_vec;
//    k_vec.resize(6);
//    std::vector<size_t> s_vec;
//    s_vec.resize(6);
//    std::vector<size_t> l_vec;
//    l_vec.resize(6);
//
//
//
////    std::cerr << this->C[this->sym_map[nuc[5]]] << std::endl;
//// TODO: this doesn't go higher than 6?
//    for (int b = 0; b < this->C.size() ; b++){
//        size_t occ_b_k = this->rankAt(k, nuc[b]);
//        k_vec[b] = this->C[this->sym_map[nuc[b]]] + occ_b_k;
//        size_t occ_b_ks = this->rankAt(k + s, nuc[b]);
//        s_vec[b] = occ_b_ks - occ_b_k;
//
//    }
//    l_vec[0] = l;
//    l_vec[4] = l_vec[0] + s_vec[0];
//
//    for (int b = 3; b > 0; b--){
//        l_vec[b] = l_vec[b+1] + s_vec[b+1];
//    }
//    l_vec[5] = l_vec[1] + s_vec[1];
//    return bi_interval(k_vec[sym_map[a]], l_vec[sym_map[a]], s_vec[sym_map[a]]);
//}
//
//
//
FastLocate::bi_interval FastLocate::forward_extend(const bi_interval& bint, size_t symbol) {
    bi_interval tmp = bi_interval(bint.reverse, bint.forward, bint.size);
    // print tmp
    // std::cerr << "tmp: " << tmp.forward << " " << tmp.reverse << " " << tmp.size << std::endl;
    size_t comp = this->complement(symbol);
    // std::cerr << symbol << " " << comp << std::endl;
    tmp = this->backward_extend(tmp, comp);
    // std::cerr << "tmp: " << tmp.forward << " " << tmp.reverse << " " << tmp.size << std::endl;
    return bi_interval(tmp.reverse, tmp.forward, tmp.size);
}


void FastLocate::initialize_complement_table() {
    for (size_t i = 0; i < 256; ++i) {
        complement_table[i] = i;
    }
    complement_table['A'] = 'T';
    complement_table['C'] = 'G';
    complement_table['G'] = 'C';
    complement_table['T'] = 'A';
    complement_table['a'] = 't';
    complement_table['c'] = 'g';
    complement_table['g'] = 'c';
    complement_table['t'] = 'a';
    complement_table['N'] = 'N';
    complement_table['n'] = 'n';
    complement_table['$'] = '$';
    complement_table['\n'] = '\n';
    complement_table[NENDMARKER] = NENDMARKER;
}













}