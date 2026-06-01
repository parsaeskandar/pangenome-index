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
#include <algorithm>
#include <chrono>
#include <limits>
#include <memory>


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





    std::vector<std::vector<std::pair<pos_t, uint16_t>>> extract_requested_tags(
            size_t thread_id, const std::vector<size_t>& requests) {

        waitForTurn(thread_id);

        std::vector<std::vector<std::pair<pos_t, uint16_t>>> extracted_tags(files.size());

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
                std::vector <std::pair<pos_t, uint16_t>> new_tags;
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
    std::vector <std::vector<std::pair < pos_t, uint16_t>>> tag_batches; // buffer of tags of each tag_file
    std::vector <size_t> current_length_tag_batch; // Remaining length of tag batch for each file

};


// =============================================================================
// In-memory tag store (used by --in-memory mode).
//
// Decodes every tag file into RAM at startup (in parallel across files), then
// exposes:
//   * the same minimal interface FileReader exposes for main()
//     (get_file_num/get_first_tag/get_next_tag), so the shared setup code can
//     drive either reader via lambdas;
//   * random-access reads of the underlying run vector for workers, which
//     bypass the global round-robin lock that bottlenecked the streaming path.
// =============================================================================
class InMemoryTagStore {
public:
    InMemoryTagStore(const std::vector<std::string>& files, size_t n_threads_for_load = 0)
        : files_(files) {
        if (files.empty()) {
            throw std::invalid_argument("File list cannot be empty.");
        }
        runs_.resize(files.size());
        cursor_run_.assign(files.size(), 0);
        cursor_intra_.assign(files.size(), 0);

        std::cerr << "[in-memory] Loading " << files.size()
                  << " tag files fully into RAM..." << std::endl;
        auto t0 = std::chrono::high_resolution_clock::now();

        size_t nt = (n_threads_for_load == 0)
                        ? static_cast<size_t>(omp_get_max_threads())
                        : n_threads_for_load;
        nt = std::max<size_t>(1, std::min(nt, files.size()));

        std::atomic<bool> failed{false};
        std::string fail_msg;

        #pragma omp parallel for num_threads(nt) schedule(dynamic)
        for (size_t i = 0; i < files.size(); ++i) {
            if (failed.load()) continue;
            try {
                sdsl::int_vector_buffer<8> in(files_[i], std::ios::in);
                if (!in.is_open()) {
                    failed.store(true);
                    #pragma omp critical
                    { fail_msg = "Cannot open file: " + files_[i]; }
                    continue;
                }
                std::vector<std::pair<pos_t, uint16_t>>& runs = runs_[i];
                gbwt::size_type pos = 0;
                const gbwt::size_type sz = in.size();
                // Heuristic reserve: typical ByteCode entry is ~2-5 bytes/run.
                runs.reserve(std::max<size_t>(1024, sz / 3));
                while (pos < sz) {
                    auto blk = panindexer::TagArray::decode_run(
                            gbwt::ByteCode::read(in, pos));
                    runs.push_back(blk);
                }
                runs.shrink_to_fit();
            } catch (std::exception& e) {
                failed.store(true);
                #pragma omp critical
                { fail_msg = std::string("Exception while loading ") + files_[i] +
                             ": " + e.what(); }
            }
        }
        if (failed.load()) throw std::runtime_error(fail_msg);

        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = t1 - t0;

        size_t total_runs = 0;
        size_t total_tags = 0;
        for (size_t i = 0; i < files.size(); ++i) {
            size_t file_tags = 0;
            for (const auto& p : runs_[i]) file_tags += p.second;
            total_runs += runs_[i].size();
            total_tags += file_tags;
            std::cerr << "  " << files_[i] << ": " << runs_[i].size()
                      << " runs / " << file_tags << " tags" << std::endl;
        }
        std::cerr << "[in-memory] Loaded " << total_runs << " runs / "
                  << total_tags << " tags in " << dt.count() << " s" << std::endl;
    }

    // -------- FileReader-compatible interface (used by shared setup) --------
    int get_file_num() const { return static_cast<int>(files_.size()); }

    pos_t get_first_tag(int f) const {
        return runs_[f][0].first;
    }

    // Streaming-style cursor advance: returns the tag at the current cursor and
    // advances by one. Used by main() during the partial-first-run prefix.
    pos_t get_next_tag(int f) {
        if (cursor_run_[f] >= runs_[f].size()) {
            std::cerr << "InMemoryTagStore::get_next_tag: out of tags in "
                      << files_[f] << std::endl;
            return pos_t{0, 0, 0};
        }
        pos_t res = runs_[f][cursor_run_[f]].first;
        cursor_intra_[f]++;
        if (cursor_intra_[f] >= runs_[f][cursor_run_[f]].second) {
            cursor_run_[f]++;
            cursor_intra_[f] = 0;
        }
        return res;
    }

    // -------- Random-access interface (used by workers) --------
    const std::vector<std::pair<pos_t, uint16_t>>& runs(int f) const {
        return runs_[f];
    }

    // Returns the (run_index, intra_run_offset) cursor as left by the
    // streaming get_next_tag prefix consumption. Workers use these as the
    // starting point for computing per-batch start positions.
    size_t cursor_run(int f)   const { return cursor_run_[f]; }
    size_t cursor_intra(int f) const { return cursor_intra_[f]; }

private:
    std::vector<std::string> files_;
    // runs_[f] is the full decoded run list of file f.
    std::vector<std::vector<std::pair<pos_t, uint16_t>>> runs_;
    // Streaming cursor for get_next_tag (one per file).
    std::vector<size_t> cursor_run_;
    std::vector<size_t> cursor_intra_;
};


// This function extract the tags starting from the starting_run first position to the starting_run + runs_per_thread last position
void extract_tags_batch(const FastLocate &r_index, FileReader &reader, size_t thread_id, std::vector<int> comp_to_file,
                        std::vector<size_t> seq_id_to_comp_id, std::vector<std::pair<pos_t, uint16_t>> &buffer,
                        size_t starting_run, size_t runs_per_thread) {
    buffer.clear();

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
//    std::chrono::duration<double> duration1 = time2 - time1;
//    std::cerr << "Indexing unique kmers took " << duration1.count() << " seconds" << std::endl;
#endif


    std::vector<size_t> request(reader.get_file_num(), 0);
    std::vector<int> index_to_file;
    // get the sample at the start of the starting_run
    if (starting_run >= r_index.tot_runs()){
        std::cerr << "Starting run is greater than total runs" << std::endl;
        return;
    }
    size_t index = r_index.getSample(starting_run);
    size_t end_index;
    size_t last_run_index = -1;

    if (starting_run + runs_per_thread >= r_index.tot_runs()){
        std::cerr << "End of runs reached in thread " << thread_id << std::endl;
        // what is the last run first position
        last_run_index = r_index.getSample(r_index.tot_runs() - 1);
        std::cerr << "Last run index: " << last_run_index << std::endl;
        end_index = -1;
    } else {
        end_index = r_index.getSample(starting_run + runs_per_thread);
    }


    int kk = 0;

    while (index != end_index && index != last_run_index){
        auto seq_id = r_index.seqId(index);
        // want to get the file number that is associated with the seq id
        auto current_file = comp_to_file.at(seq_id_to_comp_id.at(seq_id));
        index_to_file.push_back(current_file);
        request[current_file] += 1;

        index = r_index.locateNext(index);
    }


    if (index == last_run_index){

        std::cerr << "hit last run" << std::endl;
        int last_run_size = r_index.last_run_size_global();


            // have to just traverse the last run
            std::cerr << "LAST RUN size " << last_run_size << std::endl;
            for (int i = 0; i < last_run_size; i++){
                std::cerr << "last run _ index " << index << std::endl;
                auto seq_id = r_index.seqId(index);
                // want to get the file number that is associated with the seq id
                auto current_file = comp_to_file.at(seq_id_to_comp_id.at(seq_id));
                index_to_file.push_back(current_file);
                request[current_file] += 1;

                if (i == last_run_size - 1){
                    break;
                }
                index = r_index.locateNext(index);
            }

        }


    // (continue inside function)

    // we have the tags we want to extract from each file. now we extract the tags for each index
    std::vector<std::vector<std::pair<pos_t, uint16_t>>> run_tags = reader.extract_requested_tags(thread_id, request);


    pos_t current_tag;
    size_t current_index = 0;
    std::vector<std::pair<pos_t, uint16_t>> temp_buffer;

    vector<size_t> current_index_run_tags_file;
    current_index_run_tags_file.resize(run_tags.size(), 0);


    // traverse the index_to_file and get the tags from the run_tags
    for (size_t i = 0; i < index_to_file.size(); i++){
        if (run_tags[index_to_file[i]].empty()){
            std::cerr << "ERROR: Run tags is empty for file " << index_to_file[i] << "skipping for now " << std::endl;
            continue;
        }
        if (run_tags[index_to_file[i]][current_index_run_tags_file[index_to_file[i]]].second == 0){
            if (current_index_run_tags_file[index_to_file[i]] + 1 >= run_tags[index_to_file[i]].size()){
                std::cerr << "ERROR: Must have reached the end of run tag and if there are more tags requested there is a problem" << std::endl;
                std::cerr << "Index " << i << " index to file  " << index_to_file[i] << " and size " << index_to_file.size() << std::endl;
                continue;
            }
        }
        current_tag = run_tags[index_to_file[i]][current_index_run_tags_file[index_to_file[i]]].first;
        run_tags[index_to_file[i]][current_index_run_tags_file[index_to_file[i]]].second--;
        if (run_tags[index_to_file[i]][current_index_run_tags_file[index_to_file[i]]].second == 0 && current_index_run_tags_file[index_to_file[i]] + 1 < run_tags[index_to_file[i]].size()){
            current_index_run_tags_file[index_to_file[i]]++;
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


    buffer = std::move(temp_buffer);


}


// =============================================================================
// In-memory merge pipeline.
//
// Replaces the streaming worker pool (extract_tags_batch + FileReader's
// strict round-robin) with a 3-stage parallel pipeline that obeys the only
// real ordering constraints: per-file FIFO consumption across batches, and
// batch-ID-ordered output. Everything else runs free.
//
//   Phase 1 (parallel, no locks): each batch walks its BWT range via
//     locateNext, producing request[b][f] and index_to_file[b].
//   Coordinator (sequential per file, parallel across files): walks each
//     file's run list once, computing the (run_idx, intra_offset) starting
//     position for every batch from the per-batch request counts.
//   Phase 2+3 (parallel, no locks): each batch initializes its per-file
//     cursors from the coordinator output, walks its index_to_file array,
//     reads tags straight out of RAM, and builds a local run-length buffer.
//   Writer (sequential, batch-ID order): mirrors the existing writer logic,
//     including the previous_last_run merge/flush bookkeeping.
//
// Caller contract: previous_last_run has been seeded by the partial-first-run
// prefix in main(), tag_array is already streaming-encoding via the same
// sidecar files used by the streaming path, starting_run points at the first
// run beyond the prefix, and the sidecar/output streams are open.
// =============================================================================
void merge_in_memory_pipeline(
        const FastLocate& r_index,
        InMemoryTagStore& store,
        const std::vector<int>& comp_to_file,
        const std::vector<size_t>& seq_id_to_comp_id,
        panindexer::TagArray& tag_array,
        std::ofstream& out_encoded_starts,
        std::ofstream& out_bwt_intervals,
        std::pair<pos_t, uint16_t>& previous_last_run,
        size_t& tag_run_count,
        size_t starting_run,
        size_t run_per_thread,
        int threads,
        size_t chunk_size)
{
    const size_t F = static_cast<size_t>(store.get_file_num());
    const size_t total_runs_r = r_index.tot_runs();
    if (starting_run >= total_runs_r) {
        std::cerr << "[in-memory] starting_run >= tot_runs; nothing to do" << std::endl;
        return;
    }
    const size_t num_batches =
        (total_runs_r - starting_run + run_per_thread - 1) / run_per_thread;

    // Chunk the pipeline so transient Phase-1 / coordinator / result-buffer
    // memory is bounded by a window of batches, not the entire input. The tag
    // store itself (the actual tags) stays resident the whole time -- this
    // bound only affects the per-batch metadata. See the OOM analysis in the
    // commit message / discussion: with millions of batches, holding every
    // batch's index_to_file[] simultaneously is what blows up.
    if (chunk_size == 0) chunk_size = 16384;
    chunk_size = std::min(chunk_size, num_batches);
    const size_t num_chunks = (num_batches + chunk_size - 1) / chunk_size;

    std::cerr << "[in-memory] " << num_batches << " batches across "
              << threads << " threads (run_per_thread=" << run_per_thread
              << ", F=" << F << ", chunk_size=" << chunk_size
              << ", num_chunks=" << num_chunks << ")" << std::endl;

    // Use uint8_t for file ids if we have <=255 files; uint16_t otherwise.
    // This halves the memory footprint of index_to_file across all batches.
    const bool small_file_ids = (F <= 255);

    // Per-file cursors carried ACROSS chunks. Start from wherever main()'s
    // partial-first-run prefix left the InMemoryTagStore cursor.
    std::vector<size_t> file_run_cursor(F);
    std::vector<size_t> file_intra_cursor(F);
    for (size_t f = 0; f < F; ++f) {
        file_run_cursor[f]   = store.cursor_run(static_cast<int>(f));
        file_intra_cursor[f] = store.cursor_intra(static_cast<int>(f));
    }

    // Precompute the sample at the start of the last run; used by the tail
    // batch (only) to detect when it has reached the absolute end of the BWT.
    const size_t last_run_first_sample = r_index.getSample(total_runs_r - 1);

    // Cumulative timing across chunks.
    double tot_p1 = 0.0, tot_co = 0.0, tot_p2 = 0.0, tot_w = 0.0;

    for (size_t chunk_idx = 0; chunk_idx < num_chunks; ++chunk_idx) {
        const size_t batch_lo = chunk_idx * chunk_size;
        const size_t batch_hi = std::min(num_batches, batch_lo + chunk_size);
        const size_t M = batch_hi - batch_lo;
        if (M == 0) continue;

        if (chunk_idx % 10 == 0 || chunk_idx + 1 == num_chunks) {
            std::cerr << "[in-memory] Chunk " << (chunk_idx + 1) << "/" << num_chunks
                      << ": batches [" << batch_lo << ", " << batch_hi << ")"
                      << std::endl;
        }

        // Per-chunk transient state. Allocated and freed inside the loop so
        // peak memory is bounded by M, not num_batches.
        std::vector<std::vector<uint8_t>>  itf8;
        std::vector<std::vector<uint16_t>> itf16;
        if (small_file_ids) itf8.resize(M);
        else                itf16.resize(M);
        std::vector<std::vector<size_t>> request(M, std::vector<size_t>(F, 0));

        // ------------------------------ Phase 1 ------------------------------
        // Parallel locateNext walks for the M batches in this chunk.
        auto t_p1 = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (size_t bi = 0; bi < M; ++bi) {
            const size_t b = batch_lo + bi;
            const size_t s = starting_run + b * run_per_thread;
            const size_t e = std::min(s + run_per_thread, total_runs_r);

            size_t end_sentinel;
            size_t last_run_sentinel;
            if (e >= total_runs_r) {
                end_sentinel      = static_cast<size_t>(-1);
                last_run_sentinel = last_run_first_sample;
            } else {
                end_sentinel      = r_index.getSample(e);
                last_run_sentinel = static_cast<size_t>(-1);
            }

            size_t idx = r_index.getSample(s);
            std::vector<size_t>& req = request[bi];

            const size_t reserve_hint = (e - s) * 16;
            if (small_file_ids) itf8[bi].reserve(reserve_hint);
            else                itf16[bi].reserve(reserve_hint);

            while (idx != end_sentinel && idx != last_run_sentinel) {
                const size_t sid = r_index.seqId(idx);
                const int    fid = comp_to_file[seq_id_to_comp_id[sid]];
                if (small_file_ids) itf8[bi].push_back(static_cast<uint8_t>(fid));
                else                itf16[bi].push_back(static_cast<uint16_t>(fid));
                req[fid] += 1;
                idx = r_index.locateNext(idx);
            }

            if (idx == last_run_sentinel) {
                const int last_run_size = static_cast<int>(r_index.last_run_size_global());
                for (int i = 0; i < last_run_size; ++i) {
                    const size_t sid = r_index.seqId(idx);
                    const int    fid = comp_to_file[seq_id_to_comp_id[sid]];
                    if (small_file_ids) itf8[bi].push_back(static_cast<uint8_t>(fid));
                    else                itf16[bi].push_back(static_cast<uint16_t>(fid));
                    req[fid] += 1;
                    if (i == last_run_size - 1) break;
                    idx = r_index.locateNext(idx);
                }
            }
        }

        tot_p1 += std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_p1).count();

        // ----------------------------- Coordinator ----------------------------
        // For each (batch, file) in this chunk, record the per-file
        // (run_idx, intra) starting position. Sequential per file, parallel
        // across files. The cursor state at the end of this chunk is carried
        // over to the next chunk via file_{run,intra}_cursor.
        auto t_co = std::chrono::high_resolution_clock::now();

        std::vector<std::vector<size_t>>   batch_start_run(M);
        std::vector<std::vector<uint32_t>> batch_start_intra(M);
        for (size_t bi = 0; bi < M; ++bi) {
            batch_start_run[bi].assign(F, 0);
            batch_start_intra[bi].assign(F, 0);
        }
        // End-of-chunk cursors per file; written by the parallel-for, then
        // copied into file_{run,intra}_cursor after the parallel section.
        std::vector<size_t> new_run_cursor(F);
        std::vector<size_t> new_intra_cursor(F);

        #pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (size_t f = 0; f < F; ++f) {
            const auto& runs_f = store.runs(static_cast<int>(f));
            size_t run_cursor = file_run_cursor[f];
            size_t intra      = file_intra_cursor[f];

            for (size_t bi = 0; bi < M; ++bi) {
                batch_start_run[bi][f]   = run_cursor;
                batch_start_intra[bi][f] = static_cast<uint32_t>(intra);

                size_t to_advance = request[bi][f];
                while (to_advance > 0) {
                    if (run_cursor >= runs_f.size()) {
                        #pragma omp critical
                        {
                            std::cerr << "[in-memory] WARNING: advancing past end of "
                                      << "file " << f << " at chunk " << chunk_idx
                                      << " batch " << (batch_lo + bi) << " ("
                                      << to_advance << " tags requested)"
                                      << std::endl;
                        }
                        break;
                    }
                    const size_t left = static_cast<size_t>(runs_f[run_cursor].second) - intra;
                    if (to_advance >= left) {
                        to_advance -= left;
                        run_cursor++;
                        intra = 0;
                    } else {
                        intra += to_advance;
                        to_advance = 0;
                    }
                }
            }
            new_run_cursor[f]   = run_cursor;
            new_intra_cursor[f] = intra;
        }

        // request[] is no longer needed within this chunk.
        std::vector<std::vector<size_t>>().swap(request);

        tot_co += std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_co).count();

        // ---------------------------- Phase 2 + 3 ----------------------------
        // Each batch builds its own run-length-encoded buffer from RAM. No
        // locks, no waiting -- the per-file cursors are independent across
        // batches because the coordinator already partitioned each file's tag
        // stream by batch.
        auto t_p2 = std::chrono::high_resolution_clock::now();

        std::vector<std::vector<std::pair<pos_t, uint16_t>>> results(M);

        #pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (size_t bi = 0; bi < M; ++bi) {
            const size_t total_pos = small_file_ids ? itf8[bi].size() : itf16[bi].size();
            if (total_pos == 0) continue;

            std::vector<size_t> run_idx(F);
            std::vector<size_t> intra(F);
            for (size_t f = 0; f < F; ++f) {
                run_idx[f] = batch_start_run[bi][f];
                intra[f]   = batch_start_intra[bi][f];
            }

            auto& out_buf = results[bi];
            out_buf.reserve(total_pos / 8 + 16);

            pos_t cur{};
            uint16_t cur_len = 0;
            bool first = true;

            auto consume_one = [&](size_t f) -> pos_t {
                const auto& runs_f = store.runs(static_cast<int>(f));
                const pos_t t = runs_f[run_idx[f]].first;
                intra[f]++;
                if (intra[f] >= runs_f[run_idx[f]].second) {
                    run_idx[f]++;
                    intra[f] = 0;
                }
                return t;
            };

            auto step = [&](size_t f) {
                const pos_t t = consume_one(f);
                if (first) {
                    cur = t;
                    cur_len = 1;
                    first = false;
                } else if (t == cur) {
                    if (cur_len == std::numeric_limits<uint16_t>::max()) {
                        out_buf.emplace_back(cur, cur_len);
                        cur_len = 0;
                    }
                    cur_len++;
                } else {
                    out_buf.emplace_back(cur, cur_len);
                    cur = t;
                    cur_len = 1;
                }
            };

            if (small_file_ids) {
                const auto& itf = itf8[bi];
                for (size_t i = 0; i < itf.size(); ++i) step(static_cast<size_t>(itf[i]));
            } else {
                const auto& itf = itf16[bi];
                for (size_t i = 0; i < itf.size(); ++i) step(static_cast<size_t>(itf[i]));
            }
            if (!first) out_buf.emplace_back(cur, cur_len);

            // Free this batch's metadata immediately.
            if (small_file_ids) std::vector<uint8_t>().swap(itf8[bi]);
            else                std::vector<uint16_t>().swap(itf16[bi]);
            std::vector<size_t>().swap(batch_start_run[bi]);
            std::vector<uint32_t>().swap(batch_start_intra[bi]);
        }

        tot_p2 += std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_p2).count();

        // ------------------------------- Writer ------------------------------
        // Sequential, in batch-ID order within this chunk. previous_last_run
        // is the caller's reference; it carries from chunk to chunk and back
        // to main() so the boundary merge logic works across the whole run.
        auto t_w = std::chrono::high_resolution_clock::now();

        for (size_t bi = 0; bi < M; ++bi) {
            const size_t b = batch_lo + bi;
            auto& current_tags = results[bi];
            if (current_tags.empty()) {
                std::cerr << "No tags extracted for batch " << b << std::endl;
                continue;
            }

            if (current_tags.front().first == previous_last_run.first) {
                current_tags.front().second += previous_last_run.second;
            } else {
                tag_array.append_compact_run_streamed(
                        previous_last_run.first, previous_last_run.second,
                        out_encoded_starts, out_bwt_intervals);
            }

            // Pop & carry the tail run, unless this is the absolute last
            // batch of the entire run (mirrors the streaming-path writer).
            if (b + 1 < num_batches) {
                previous_last_run = current_tags.back();
                current_tags.pop_back();
            }

            tag_run_count += current_tags.size();
            for (auto& p : current_tags) {
                tag_array.append_compact_run_streamed(
                        p.first, p.second, out_encoded_starts, out_bwt_intervals);
            }
            std::vector<std::pair<pos_t, uint16_t>>().swap(current_tags);
        }

        tot_w += std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_w).count();

        // Carry per-file cursors forward to the next chunk.
        file_run_cursor   = std::move(new_run_cursor);
        file_intra_cursor = std::move(new_intra_cursor);
    }

    std::cerr << "[in-memory] Totals: Phase 1=" << tot_p1
              << "s, Coordinator=" << tot_co
              << "s, Phase 2+3=" << tot_p2
              << "s, Writer=" << tot_w << "s" << std::endl;
}


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
    // Parse positional args + the --in-memory flag (can appear anywhere).
    // --in-memory swaps the streaming FileReader for an InMemoryTagStore and
    // runs the parallel in-memory pipeline. Faster but uses much more RAM
    // because every tag file is fully decoded into RAM up front.
    bool in_memory = false;
    // 0 means "let the pipeline pick a default" (currently 16384 batches/chunk).
    // Tune this when in-memory mode OOMs; lower = less transient RAM, slightly
    // more chunk-overhead per pass.
    size_t in_memory_chunk_size = 0;
    std::vector<std::string> pos_args;
    pos_args.reserve(argc);
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--in-memory" || a == "--in-ram") {
            in_memory = true;
        } else if (a == "--chunk" || a == "--chunk-size") {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << a << " requires a value" << std::endl;
                exit(1);
            }
            try {
                in_memory_chunk_size = std::stoull(argv[++i]);
            } catch (...) {
                std::cerr << "Error: invalid chunk size '" << argv[i] << "'" << std::endl;
                exit(1);
            }
        } else {
            pos_args.push_back(std::move(a));
        }
    }
    if (pos_args.size() != 3) {
        std::cerr << "usage: merge_tags [--in-memory] [--chunk N] "
                     "<gbz_graph> <r_index> <tag_array_dir>" << std::endl;
        exit(0);
    }

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
#endif


    int threads = omp_get_max_threads();
    omp_set_num_threads(threads);

    std::string gbz_graph = pos_args[0];
    std::string r_index_file = pos_args[1];
    std::string tag_array_index_dir = pos_args[2];

    std::cerr << "Mode: " << (in_memory ? "in-memory (fast, RAM-heavy)"
                                        : "streaming (default)") << std::endl;


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


    cerr << "Reading the whole genome r-index file (encoded)" << endl;
    FastLocate r_index;
    {
        std::ifstream rin(r_index_file, std::ios::binary);
        if (!rin) {
            std::cerr << "Cannot open r-index: " << r_index_file << std::endl;
            std::exit(EXIT_FAILURE);
        }
        r_index.load_encoded(rin);
    }


    std::cerr << "Finding the node to component mapping" << std::endl;
    // finding the components
    std::unordered_map<nid_t, size_t> node_to_comp_map = node_to_component(gbz);

    std::vector<int> file_to_comp(number_of_file);
    std::vector<int> comp_to_file(number_of_file);

    std::cerr << "Initializing the reader" << std::endl;
    // Only one of these is non-null depending on `in_memory`. We dispatch the
    // small number of reader calls in shared setup through lambdas so the rest
    // of main() doesn't have to care which one is active.
    std::unique_ptr<FileReader>         stream_reader;
    std::unique_ptr<InMemoryTagStore>   mem_reader;
    if (in_memory) {
        mem_reader = std::make_unique<InMemoryTagStore>(files);
    } else {
        stream_reader = std::make_unique<FileReader>(files, threads, 50000000);
    }
    auto reader_get_first_tag = [&](int f) -> pos_t {
        return in_memory ? mem_reader->get_first_tag(f)
                         : stream_reader->get_first_tag(f);
    };
    auto reader_get_next_tag = [&](int f) -> pos_t {
        return in_memory ? mem_reader->get_next_tag(f)
                         : stream_reader->get_next_tag(f);
    };

    std::cerr << "Creating the mapping from comp to tag files" << std::endl;
    // for each tag block files, we read the first block and read the first node
    for (auto i = 0; i < number_of_file; i++) {
        // get the component of the node of the first block
        pos_t first_tag = reader_get_first_tag(i);
        size_t comp = node_to_comp_map[id(first_tag)];
        std::cerr << "The component of the first block of file" << files[i] << " is " << comp << " first tag " << first_tag  << " node id is " << id(first_tag) << std::endl;
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
        auto seq_graph_nodes = gbz.index.extract(i);
        if (!seq_graph_nodes.empty()) {
            size_t node_id = gbwt::Node::id(seq_graph_nodes[0]);
            seq_id_to_comp_id[i] = node_to_comp_map.at(node_id);
        }
    }
    std::cerr << "The mapping from seq id to comp id is done" << std::endl;


    // note that the first #num_seq tags are correspond to the ENDMARKERs
    auto total_tags_count = r_index.get_sequence_size() - total_strings;
    std::cerr << "Total tags count " << total_tags_count << std::endl;


    cerr << "Creating the whole genome tag array indexing" << endl;
    TagArray tag_array;

#if TIME
    auto time2 = chrono::high_resolution_clock::now();
#endif


    std::cerr << "Thread lists " << std::endl;
    std::vector <std::thread> threads_list;
    threads_list.reserve(threads);
    std::vector<std::vector<std::pair<pos_t, uint16_t>>> thread_buffers(threads);


    const std::string filename = "whole_genome_tag_array_compressed.tags";
    const std::string encoded_runs_path = filename + ".encoded_runs.tmp";
    const std::string encoded_starts_file = "encoded_starts.bin";
    const std::string bwt_intervals_file = "bwt_intervals.bin";

    // Check if the file exists and delete it
    if (std::filesystem::exists(filename)) {
        if (std::remove(filename.c_str()) != 0) {
            std::cerr << "Error: Unable to delete the existing file.\n";
            return 1; // Exit with error
        } else {
            std::cerr << "Existing file deleted successfully.\n";
        }
    }

    if (std::filesystem::exists(encoded_starts_file)) {
        if (std::remove(encoded_starts_file.c_str()) != 0) {
            std::cerr << "Error: Unable to delete the existing file.\n";
            return 1; // Exit with error
        } else {
            std::cerr << "Existing file deleted successfully.\n";
        }
    }

    if (std::filesystem::exists(bwt_intervals_file)) {
        if (std::remove(bwt_intervals_file.c_str()) != 0) {
            std::cerr << "Error: Unable to delete the existing file.\n";
            return 1; // Exit with error
        } else {
            std::cerr << "Existing file deleted successfully.\n";
        }
    }



    // Open the file for writing
    std::ofstream out(filename, std::ios::binary | std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open file for writing.\n";
        return 1;
    }

    std::ofstream out_encoded_starts(encoded_starts_file, std::ios::binary | std::ios::app);
    if (!out_encoded_starts.is_open()) {
        std::cerr << "Error: Cannot open file for writing.\n";
        return 1;
    }

    // encoded_runs will be written as sdsl int_vector; no size header needed.


    std::ofstream out_bwt_intervals(bwt_intervals_file, std::ios::binary | std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open file for writing.\n";
        return 1;
    }


    // ############################################## variables to change
    size_t run_per_thread = 500;
    int encoded_start_every_k_run = 10;




    // ############################################## variables
    int remaining_run_to_write_start = 0;
    size_t cumulative_starts = 0;
    size_t encoded_start_ones = 0;
    size_t start_pos = 0;


//    TagArray tag_array;



    size_t tag_count = 0;
    size_t tag_run_count = 0;
    std::pair<pos_t, uint16_t> previous_last_run;
    std::vector<std::pair<pos_t, uint16_t>> temp_tag_runs;

    // have to calculate the number of ENDMARKERS at the beginning of the r-index and put special value 0 for them
    // Store in multiple runs so each run length fits in uint16_t (max 65535)
    auto num_endmarkers = total_strings;
    constexpr uint16_t max_run_len = 65535;
    pos_t endmarker_pos = pos_t{0, 0, 0};
    size_t remaining = num_endmarkers;
    while (remaining > 0) {
        uint16_t chunk = static_cast<uint16_t>(std::min(remaining, static_cast<size_t>(max_run_len)));
        temp_tag_runs.push_back(std::make_pair(endmarker_pos, chunk));
        remaining -= chunk;
    }
    tag_count += num_endmarkers;
    tag_run_count += temp_tag_runs.size();




    // Determine maximum node id to compute width for encoded runs
    nid_t max_node_id = 0;
    auto weak_components = gbwtgraph::weakly_connected_components(gbz.graph);
    for (const auto& comp : weak_components) {
        for (nid_t nid : comp) { if (nid > max_node_id) max_node_id = nid; }
    }
    size_t node_bits = sdsl::bits::hi(max_node_id) + 1;
    size_t width_bits = 10 + 1 + node_bits;
    tag_array.begin_encoded_runs_sdsl(encoded_runs_path, width_bits);

    // Now have to find the run and index of the first index that is not an ENDMARKER
    // having to find the run and index of the BWT position num_endmarkers
    // iter.first is the block number num_endmarkers is in
    // iter.second is the offset of the beginning of the block
    size_t run_id = 0; size_t offset_of_first = 0;
    r_index.run_id_and_offset_at(num_endmarkers, run_id, offset_of_first);

    // we want a job handling the remaining of the run_id run
    auto first = r_index.getSample(run_id);

    std::cerr << "start of the run first " << first << std::endl;

    // Iterate until the start of the range and locate the first occurrence.
    while (offset_of_first < num_endmarkers) {
        std::cerr << "Not in here" << std::endl;
        first = r_index.locateNext(first);
        offset_of_first++;
    }


    std::cerr << "run id " << run_id << " offset of first " << offset_of_first << std::endl;
    auto end = r_index.getSample(run_id + 1);
    // we have to handle the indexes that are between first and end


    while (first != end){

        auto seq_id = r_index.seqId(first);
//        std::cerr << "Seq id " << seq_id << std::endl;
        // want to get the file number that is associated with the seq id
        auto current_file = comp_to_file[seq_id_to_comp_id[seq_id]];

        auto temp_tag = reader_get_next_tag(current_file);

        if (temp_tag_runs.back().first == temp_tag){
            temp_tag_runs.back().second += 1;
            if (first == end){
                tag_count += temp_tag_runs.back().second;
            }
        } else {
            tag_count += temp_tag_runs.back().second;
            temp_tag_runs.push_back(std::make_pair(temp_tag, 1));
        }
        first = r_index.locateNext(first);
    }

    std::cerr << "Writing " << temp_tag_runs.size() << " tags before running actual jobs" << std::endl;
    // total size of the tags before running actual jobs
    std::cerr << "Total size of the tags before running actual jobs is " << tag_count << std::endl;

    tag_run_count += temp_tag_runs.size();


    previous_last_run = temp_tag_runs.back(); // have to check with next one and merge them if needed
    temp_tag_runs.pop_back();
    // TODO: handle the case in compressed version - also the last element might be needed for merging later


    for (auto &p : temp_tag_runs) {
        tag_array.append_compact_run_streamed(p.first, p.second, out_encoded_starts, out_bwt_intervals);
    }
//    std::vector<gbwt::byte_type> temp_encoded_runs;
//    if (temp_tag_runs.size() > 0) {
//        for (const auto& [value, run_length] : temp_tag_runs) {
//            if (remaining_run_to_write_start % encoded_start_every_k_run == 0){
//                start_pos = cumulative_starts + temp_encoded_runs.size();
//                // write the encoded start in file
//                out_encoded_starts.write(reinterpret_cast<const char*>(&start_pos), sizeof(start_pos));
//                encoded_start_ones++;
//            }
//
//            gbwt::size_type encoded1 =
//                    (gbwtgraph::offset(value)) | (gbwtgraph::is_rev(value) << 10) |
//                    (run_length << 11) |
//                    (gbwtgraph::id(value) << 19);
//
//            gbwt::ByteCode::write(temp_encoded_runs, encoded1);
//
//
//
//            remaining_run_to_write_start++;
//        }
//
//        size_t size = temp_encoded_runs.size();
////            out.write(reinterpret_cast<const char *>(&size), sizeof(size));
//        out.write(reinterpret_cast<const char *>(temp_encoded_runs.data()), size * sizeof(gbwt::byte_type));
//        cumulative_starts += size;
//
//    }
//    tag_array.serialize_run_by_run(out, temp_tag_runs);


    size_t starting_run = run_id + 1;








    size_t number_of_jobs = (r_index.tot_runs() - starting_run + run_per_thread - 1) / run_per_thread;


    std::cerr << "Total bwt indexes to find tags for is " << r_index.get_sequence_size() << std::endl;
    std::cerr << "We will handle " << r_index.tot_runs() << " runs" << std::endl;
    cerr << "Merging tags and creating the whole genome tag array indexing" << endl;

    if (in_memory) {
        // Fast RAM-heavy path: parallel locateNext walks + per-file
        // coordinator prefix sum + parallel RLE buffer builds + serial writer.
        merge_in_memory_pipeline(r_index, *mem_reader, comp_to_file, seq_id_to_comp_id,
                                 tag_array, out_encoded_starts, out_bwt_intervals,
                                 previous_last_run, tag_run_count,
                                 starting_run, run_per_thread, threads,
                                 in_memory_chunk_size);
    } else {
    for (size_t to_read = 0; to_read < threads && starting_run < r_index.tot_runs(); to_read++) {
        threads_list.emplace_back(extract_tags_batch, std::ref(r_index), std::ref(*stream_reader), to_read,
                                  comp_to_file, seq_id_to_comp_id,
                                  std::ref(thread_buffers[to_read]), starting_run, run_per_thread);
        starting_run += run_per_thread;
    }







    for (size_t to_write = 0; to_write < number_of_jobs; to_write++) {
        if (to_write % 10000 == 0) {
            std::cerr << "Writing job " << to_write << std::endl;
        }
        size_t thread_id = to_write % threads;
        threads_list[thread_id].join();
        std::vector<std::pair<pos_t, uint16_t>> current_tags;
        current_tags.swap(thread_buffers[thread_id]);

        if (to_write + threads < number_of_jobs) {
            threads_list[thread_id] = std::thread(extract_tags_batch, std::ref(r_index), std::ref(*stream_reader), thread_id,
                                                  comp_to_file, seq_id_to_comp_id,
                                                  std::ref(thread_buffers[thread_id]), starting_run, run_per_thread);
            starting_run += run_per_thread;
        }

        if (current_tags.size() == 0) {
            std::cerr << "No tags extracted for thread " << thread_id << std::endl;
            continue;
        }


        std::vector<gbwt::byte_type> encoded_runs;

        if (current_tags[0].first == previous_last_run.first) {
            current_tags[0].second += previous_last_run.second;
        } else {
            std::vector<std::pair<pos_t, uint16_t>> temp = {previous_last_run};
            for (auto &p : temp) {
                tag_array.append_compact_run_streamed(p.first, p.second, out_encoded_starts, out_bwt_intervals);
            }

//            std::cerr << "Handling the previous last run" << std::endl;

//            if (remaining_run_to_write_start % encoded_start_every_k_run == 0){
//                start_pos = cumulative_starts + encoded_runs.size();
//                // write the encoded start in file
//                out_encoded_starts.write(reinterpret_cast<const char*>(&start_pos), sizeof(start_pos));
//                encoded_start_ones++;
//            }
//
//            gbwt::size_type temp_encoded =
//                    (gbwtgraph::offset(previous_last_run.first)) | (gbwtgraph::is_rev(previous_last_run.first) << 10) |
//                    (previous_last_run.second << 11) |
//                    (gbwtgraph::id(previous_last_run.first) << 19);
//
//            gbwt::ByteCode::write(encoded_runs, temp_encoded);
//
//            remaining_run_to_write_start++;

//                std::vector<std::pair<pos_t, uint16_t>> temp = {previous_last_run};
//
//                tag_array.serialize_run_by_run(out, temp);

        }
        if (to_write < number_of_jobs - 1) {
            previous_last_run = current_tags.back();
            current_tags.pop_back();
        }


        tag_run_count += current_tags.size();
        for (auto &p : current_tags) {
            tag_array.append_compact_run_streamed(p.first, p.second, out_encoded_starts, out_bwt_intervals);
        }


    }
    } // end of `if (in_memory) ... else { ...streaming worker pool... }`
    // When there are no jobs, the main thread's last run was popped into previous_last_run and never written.
    // (Applies to both modes; in in-memory mode threads_list is empty so the join loop is a no-op.)
    if (number_of_jobs == 0) {
        tag_array.append_compact_run_streamed(previous_last_run.first, previous_last_run.second, out_encoded_starts, out_bwt_intervals);
    }
    for (auto &thread: threads_list) {
        if (thread.joinable()) {
            thread.join();
        }
    }


    // Finish encoded runs iv and serialize it to the main index file
    tag_array.end_encoded_runs_sdsl();
    {
        std::ifstream iv_in(encoded_runs_path, std::ios::binary);
        tag_array.load_encoded_runs_sdsl(iv_in);
        tag_array.serialize_encoded_runs_sdsl(out);
        iv_in.close();
        std::remove(encoded_runs_path.c_str());
    }

    out_encoded_starts.close();
    out_bwt_intervals.close();
    out.close();

    std::cerr << "Total tags count run " << tag_run_count << std::endl;
    std::cerr << "Total length of encoded runs " << start_pos + 1 << std::endl;
    std::cerr << "Saving the start of runs every " << encoded_start_every_k_run << " which lead to " << encoded_start_ones << " 1s in the start sd_vector" << std::endl;

    tag_array.merge_compressed_files_sdsl(filename, encoded_starts_file, bwt_intervals_file);
    std::cerr << "Index files merged and ready to use!" << std::endl;



//    out.close();
#if TIME
    auto time3 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = time3 - time2;
    std::cerr << "Converting tags using multiple threads took " << duration2.count() << " seconds" << std::endl;
#endif


    return 0;


}


