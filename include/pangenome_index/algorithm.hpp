//
// Created by seeskand on 9/16/24.
//
//
#ifndef PANGENOME_INDEX_ALGORITHM_HPP
#define PANGENOME_INDEX_ALGORITHM_HPP


#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <omp.h>
#include <hash_map.hpp>
#include "bplus_tree.hpp"
#include "gbwtgraph/gbz.h"
#include <hash_map.hpp>
#include <gbwt/internal.h>
#include "r-index.hpp"
#include "gbwtgraph/algorithms.h"
#include "utils.hpp"




using namespace std;
using namespace gbwtgraph;


namespace panindexer {



//    template<typename T>
//    class ThreadSafeQueue {
//    public:
//        // Push item to the queue
//        void push(const T &item);
//
//        // Try to pop item from the queue, return true if successful
//        bool try_pop(T &item);
//
//    private:
//        std::queue <T> queue_;
//        std::mutex mutex_;
//        std::condition_variable cond_var_;
//    };
//
//    template<typename T>
//    void ThreadSafeQueue<T>::push(const T &item) {
//        std::unique_lock <std::mutex> lock(mutex_);
//        queue_.push(item);
//        lock.unlock();
//        cond_var_.notify_one();
//    }
//
//    template<typename T>
//    bool ThreadSafeQueue<T>::try_pop(T &item) {
//        std::unique_lock <std::mutex> lock(mutex_);
//        if (queue_.empty()) {
//            return false;
//        }
//        item = queue_.front();
//        queue_.pop();
//        return true;
//    }

template<typename T>
class ThreadSafeQueue {
public:
    // Push item to the queue
    void push(const T &item) {
        std::unique_lock<std::mutex> lock(mutex_);
        queue_.push(item);
        cond_var_.notify_one();  // Notify while holding the lock
    }

    // Try to pop an item (non-blocking)
    bool try_pop(T &item) {
        std::unique_lock<std::mutex> lock(mutex_);
        if (queue_.empty()) {
            return false;
        }
        item = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    // Blocking pop (waits until an item is available)
    void wait_and_pop(T &item) {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_var_.wait(lock, [this]() { return !queue_.empty(); }); // Wait until not empty
        item = std::move(queue_.front());
        queue_.pop();
    }

    // Check if the queue is empty (thread-safe)
    bool empty() {
        std::unique_lock<std::mutex> lock(mutex_);
        return queue_.empty();
    }

    size_t size() {
        std::unique_lock<std::mutex> lock(mutex_);
        return queue_.size();
    }

private:
    std::queue<T> queue_;
    std::mutex mutex_;
    std::condition_variable cond_var_;
};

// This function input is the OCC vector of the end of the sequences and it returns the sorted end_of_seq vector
// which is the sorted vector of pairs (i, SA[i]) for the end of each sequence which is (ISA[j], j)
vector <pair<uint64_t, uint64_t>> sort_end_of_seq(vector <pair<uint64_t, uint64_t>> &OCC) {

    // Sort the end_of_seq vector by the second element of each pair
    sort(OCC.begin(), OCC.end(),
         [](const pair <uint64_t, uint64_t> &a, const pair <uint64_t, uint64_t> &b) {
             return a.second < b.second;
         });

    return OCC;
}




void kmers_to_bplustree_worker(FastLocate &idx, std::vector<std::pair<Run, size_t>> &queue,
                                   hash_map <gbwtgraph::Key64::value_type, gbwtgraph::Position> &index,
                                   size_t k, gbwt::range_type interval, const string &current_kmer) {


    if (current_kmer.length() == k && interval.first <= interval.second) {
    if (interval.second < interval.first) std::cerr << interval.first << " " << interval.second << " " << current_kmer << std::endl;

        // creating the kmer with the key type
        gbwtgraph::Key64 kmer_key = gbwtgraph::Key64::encode(current_kmer);
        // check if the kmer_key is in the index and if it is add the run to the queue
        auto it = index.find(kmer_key.get_key());
        if (it != index.end()) {

            Run run = {interval.first, it->second};

            queue.push_back({run, interval.second - interval.first + 1});

//            queue.push( {run, interval.second - interval.first + 1});
        }
        return;
    }

    for ( char base : {'A', 'C', 'G', 'T'}) {
        if (interval.first <= interval.second) {
            kmers_to_bplustree_worker(idx, queue, index, k, idx.LF(interval, base), base + current_kmer);
        }
    }
}

void parallel_kmers_to_bplustree(FastLocate &idx, BplusTree <Run> &bptree,
                                 hash_map <gbwtgraph::Key64::value_type, gbwtgraph::Position> &index, size_t k,
                                 gbwt::range_type interval) {

    int threads = omp_get_max_threads();
    // Thread-safe queue to collect results
//    ThreadSafeQueue <std::pair<Run, size_t>> queue;
    std::vector<std::vector<std::pair<Run, size_t>>> batches(threads);

    // number of threads

    // splitting the starting range (0, idx.bwt_size() - 1) into parts and call the kmers_to_bplustree_worker function
    // for each part
    size_t part_size = idx.bwt_size();
    part_size = (part_size - 1) / threads;

    std::cerr << "running each part " << part_size << " using threads " << threads << std::endl;

    std::vector<string> starting_kmers = {"A", "C", "G", "T"};
#pragma omp parallel for
    for (int i = 0; i < threads; i++) {
        size_t start = i * part_size;
        size_t end = (i + 1) * part_size - 1;
        if (i == threads - 1) {
            end = idx.bwt_size() - 1;
        }

        batches[i].reserve(end - start);
        kmers_to_bplustree_worker(idx, batches[i], index, k, {start, end}, "");
//        kmers_to_bplustree_worker(idx, queue, index, k, interval, starting_kmers[i]);
    }





    // Single-threaded insertion into BPlusTree
//    std::pair <Run, size_t> result;
//    std::vector<std::pair<Run, size_t>> merged;
//    for (const auto &batch : batches) {
//        merged.insert(merged.end(), batch.begin(), batch.end());
//    }

    // Step 2: Sort based on Run.start_position
//    std::sort(merged.begin(), merged.end(), [](const auto &a, const auto &b) {
//        return a.first.start_position < b.first.start_position;
//    });

//    std::cerr << "finished calculating the intervals with size " << merged.size() << std::endl;


    for (const auto &batch : batches) {
        for (const auto &item : batch) {
            bptree.insert(item.first, item.second);
        }
    }
//    for (const auto &item: merged){
//        bptree.insert(item.first, item.second);
//    }
//    while (queue.try_pop(result)) {
//        bptree.insert(result.first, result.second);
//    }
}




vector <pair<Run, size_t>>
extend_kmers_bfs_parallel(GBWTGraph &graph, FastLocate &idx, BplusTree <Run> &bptree, int batch_size) {
    vector <pair<Run, size_t>> extension_candidates;

    int num_threads = omp_get_max_threads();
    vector < std::queue < pair < Run, size_t > >> bfs_queues(num_threads);
    vector < vector < pair < Run, size_t>>> batches(num_threads);

    // Mutex for synchronizing access to bptree
    std::mutex bptree_mutex;

    // Add initial runs to the queues
    int item_count = 0;
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        Run current_item = *it;
        if (current_item.graph_position.value != 0) {
            auto next_it = it;
            ++next_it;
            if (next_it != bptree.end()) {
                Run next_item = *next_it;
                bfs_queues[item_count % num_threads].push(make_pair(current_item, next_item.start_position));
                ++item_count;
            }
        }
    }

    // Parallel BFS to extend kmers
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        while (!bfs_queues[thread_num].empty()) {
            pair <Run, size_t> current_pair = bfs_queues[thread_num].front();
            bfs_queues[thread_num].pop();

            Run current_item = current_pair.first;
            size_t interval_end = current_pair.second;

            if (current_item.graph_position.value != 0) {
                auto current_starting_pos = current_item.start_position;
                auto current_graph_pos = current_item.graph_position;

                // Decode the Position to a pos_t
                pos_t current_pos = current_graph_pos.decode();

                // Get the traversal of the graph position
                handle_t current_handle = graph.get_handle(id(current_pos), false);

                // The BWT interval of the current kmer
                gbwt::range_type current_interval = {current_starting_pos, interval_end - 1};

                if (!is_rev(current_pos)) {
                    if (offset(current_pos) > 0) {
                        auto prev_base = graph.get_base(graph.get_handle(id(current_pos), is_rev(current_pos)),
                                                        offset(current_pos) - 1);
                        auto prev_graph_pos = pos_t{id(current_pos), is_rev(current_pos),
                                                    offset(current_pos) - 1};

                        auto new_range = idx.LF(current_interval, prev_base);

                        if (new_range.first <= new_range.second) {
                            Run temp_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
                            batches[thread_num].push_back(make_pair(temp_run, new_range.second - new_range.first + 1));

                            if (batches[thread_num].size() >= batch_size) {
                                std::lock_guard <std::mutex> lock(bptree_mutex);
                                for (const auto &item: batches[thread_num]) {
//                                    bptree.insert(item.first, item.second);
//                                    bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                                    if (bptree.insert_success(item.first, item.second)) {
                                        bfs_queues[thread_num].push(
                                                make_pair(item.first, item.second + item.first.start_position));
                                    }
                                }
                                batches[thread_num].clear();
                            }
                        }
                    } else {
                        int prev_bases_num = 0;
                        handle_t prev_node;
                        char prev_base;
                        pos_t prev_graph_pos;

                        graph.follow_edges(current_handle, true, [&](const handle_t &prev) {
                            prev_bases_num++;
                            if (prev_bases_num != 1) return false;
                            prev_base = graph.get_base(prev, graph.get_length(prev) - 1);
                            prev_node = prev;
                            return true;
                        });

                        if (prev_bases_num == 1) {
                            prev_graph_pos = pos_t{graph.get_id(prev_node), graph.get_is_reverse(prev_node),
                                                   graph.get_length(prev_node) - 1};
                            auto new_range = idx.LF(current_interval, prev_base);
                            if (new_range.first <= new_range.second) {
                                Run new_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
                                batches[thread_num].push_back({new_run, new_range.second - new_range.first + 1});

                                if (batches[thread_num].size() >= batch_size) {
                                    std::lock_guard <std::mutex> lock(bptree_mutex);
                                    for (const auto &item: batches[thread_num]) {
//                                        bptree.insert(item.first, item.second);
//                                        bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                                        if (bptree.insert_success(item.first, item.second)) {
                                            bfs_queues[thread_num].push(
                                                    make_pair(item.first, item.second + item.first.start_position));
                                        }
                                    }
                                    batches[thread_num].clear();
                                }
                            }
                        }
                    }
                }
            }
        }

        // Process remaining items in the local batch
        if (!batches[thread_num].empty()) {
            std::lock_guard <std::mutex> lock(bptree_mutex);
            for (const auto &item: batches[thread_num]) {
//                bptree.insert(item.first, item.second);
//                bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                if (bptree.insert_success(item.first, item.second)) {
                    bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                }
            }
            batches[thread_num].clear();
        }
    }

    return extension_candidates;
}


void traverse_sequences_parallel(GBZ &gbz, BplusTree <Run> &bptree, FastLocate &idx,
                                 vector <pair<uint64_t, uint64_t>> &end_of_seq) {
    auto number_of_sequences = end_of_seq.size();

    std::vector<Run> tmp1;
    omp_lock_t lock;
    omp_init_lock(&lock);

    cerr << "Filling the gaps on the bptree" << endl;

    // Using `num_threads` for efficient merging
    int num_threads = omp_get_max_threads();
    std::vector<std::vector<Run>> thread_tmp(num_threads); // Each thread gets its own vector

#pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::vector<Run>& local_tmp1 = thread_tmp[thread_id]; // Thread-local storage

#pragma omp for
        for (int seq_num = 0; seq_num < number_of_sequences; ++seq_num) {

            auto seq_graph_nodes = gbz.index.extract(seq_num);
            auto bwt_index = end_of_seq[seq_num].first;

            auto current_nodes_index = seq_graph_nodes.size() - 1;
            auto current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
            auto in_node_index = gbz.graph.get_length(current_node) - 1;

            while (true) {
                auto temp = idx.psi(bwt_index);
                bwt_index = temp.second;
                auto first = temp.first;

                if (first == NENDMARKER) {
                    cerr << "The end of the sequence at BWT index " << bwt_index << endl;
                    break;
                }

                if (in_node_index == -1) {
                    if (current_nodes_index == 0) {
                        break;
                    }
                    current_nodes_index--;
                    current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
                    in_node_index = gbz.graph.get_length(current_node) - 1;
                }

                in_node_index--;

                auto bptree_search = bptree.search(bwt_index);

                if (bptree_search.graph_position.value == 0) {
                    pos_t current_pos = pos_t{
                            gbz.graph.get_id(current_node),
                            gbz.graph.get_is_reverse(current_node),
                            in_node_index + 1
                    };

                    Run current_run = {bwt_index, gbwtgraph::Position::encode(current_pos)};
                    local_tmp1.push_back(current_run);

                } else {
                    pos_t current_pos = pos_t{
                            gbz.graph.get_id(current_node),
                            gbz.graph.get_is_reverse(current_node),
                            in_node_index + 1
                    };
                    assert(gbwtgraph::Position::encode(current_pos).value == bptree_search.graph_position.value);
                }
            }
        }
    }

    cerr << "The traverse completed" << endl;
    std::cerr << "Adding runs with size 1 to the BPlusTree" << std::endl;

    // Merge all thread-local results into `tmp1` (single-threaded step)
    for (auto& local_tmp : thread_tmp) {
        for (auto &i : local_tmp) {
            bptree.insert(i, 1);
        }
//        tmp1.insert(tmp1.end(), local_tmp.begin(), local_tmp.end());
    }


    omp_destroy_lock(&lock);
}



// This function returns a mapping from each node id to its component id
std::unordered_map<nid_t, size_t> node_to_component(GBZ &gbz){
    auto weakly_connected_components = gbwtgraph::weakly_connected_components(gbz.graph);
    std::unordered_map<nid_t, size_t> node_to_component;
    size_t node_count = gbz.graph.get_node_count();
//    std::cerr << "number of weak components " << weakly_connected_components.size() << std::endl;
//    std::cerr << "The number of nodes in the graph is " << node_count << std::endl;
    node_to_component.reserve(node_count);

    for (size_t i = 0; i < weakly_connected_components.size(); i++) {
//        cerr << "weak comp " << i << endl;
        const std::vector<nid_t>& component = weakly_connected_components[i];

        for(nid_t node_id : component) {
//            cerr << node_id << endl;
            node_to_component[node_id] = i; }
    }

    return node_to_component;

}




// This function find the long MEMs in the pattern using the algorithm2 from https://arxiv.org/abs/2403.02008
void find_long_mems(const std::string& pattern, size_t min_len, FastLocate& fmd_index) {
    size_t m = pattern.length();
    size_t i = 0;

    while (i <= m - min_len) {
        size_t j = i + min_len - 1;
        size_t len_fwd = 0;

        FastLocate::bi_interval bint = {0, 0, fmd_index.bwt_size()};

        for (int k = j; k >= (int)i; --k) {
            bint = fmd_index.backward_extend(bint, pattern[k]);
            if (bint.size == 0) {
                break;
            }
            ++len_fwd;
        }

        if (len_fwd < (j - i + 1)) {
            // Failed to find a MEM of length ≥ L starting at i, skip ahead
            size_t k = j - len_fwd + 1;
            i = k;
            continue;
        }

        std::string rev_pattern = pattern.substr(i);
        std::reverse(rev_pattern.begin(), rev_pattern.end());
        for (char& c : rev_pattern) {
            c = fmd_index.complement(c);
        }

        FastLocate::bi_interval rev_bint = {0, 0, fmd_index.bwt_size()};
        size_t len_rev = 0;
        for (char c : rev_pattern) {
            rev_bint = fmd_index.backward_extend(rev_bint, c);
            if (rev_bint.size == 0) {
                break;
            }
            ++len_rev;
        }

        size_t end = i + len_rev - 1;

        // Report MEM: pattern[i..end]
        std::string mem = pattern.substr(i, len_rev);
        std::cout << "MEM: [" << i << ", " << end << "] = " << mem << "\n";

        // Step 3: prepare i for next iteration
        if (end < m - 1) {
            FastLocate::bi_interval skip_bint = {0, 0, fmd_index.bwt_size() };
            size_t len_skip = 0;
            for (int k = end + 1; k >= 0; --k) {
                skip_bint = fmd_index.backward_extend(skip_bint, pattern[k]);
                if (skip_bint.size == 0) break;
                ++len_skip;
            }
            i = end - len_skip + 2;
        } else {
            break;
        }
    }
}




}

#endif //PANGENOME_INDEX_ALGORITHM_HPP
