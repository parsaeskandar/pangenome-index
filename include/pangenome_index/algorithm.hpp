//
// Created by seeskand on 9/16/24.
//

#ifndef PANGENOME_INDEX_ALGORITHM_HPP
#define PANGENOME_INDEX_ALGORITHM_HPP


#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <omp.h>
#include <hash_map.hpp>
#include "../r-index/internal/r_index.hpp"
#include "bplus_tree.hpp"
#include "gbwtgraph/gbz.h"
#include <hash_map.hpp>
#include <gbwt/internal.h>
#include "r-index.hpp"




using namespace std;
using namespace ri;
using namespace gbwtgraph;


namespace panindexer {


    template<typename T>
    class ThreadSafeQueue {
    public:
        // Push item to the queue
        void push(const T &item);

        // Try to pop item from the queue, return true if successful
        bool try_pop(T &item);

    private:
        std::queue <T> queue_;
        std::mutex mutex_;
        std::condition_variable cond_var_;
    };

    template<typename T>
    void ThreadSafeQueue<T>::push(const T &item) {
        std::unique_lock <std::mutex> lock(mutex_);
        queue_.push(item);
        lock.unlock();
        cond_var_.notify_one();
    }

    template<typename T>
    bool ThreadSafeQueue<T>::try_pop(T &item) {
        std::unique_lock <std::mutex> lock(mutex_);
        if (queue_.empty()) {
            return false;
        }
        item = queue_.front();
        queue_.pop();
        return true;
    }

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

//void kmers_to_bplustree_worker(FastLocate &idx, ThreadSafeQueue<std::pair < Run, size_t>> &queue,
//                                   hash_map <gbwtgraph::Key64::value_type, gbwtgraph::Position> &index,
//                                   size_t k, range_t interval, const string &current_kmer) {
//
////    std::cerr << interval.first << " " << interval.second << " " << current_kmer << std::endl;
//    if (current_kmer.length() == k && interval.first <= interval.second) {
//
//        // creating the kmer with the key type
//        gbwtgraph::Key64 kmer_key = gbwtgraph::Key64::encode(current_kmer);
//        // check if the kmer_key is in the index and if it is add the run to the queue
//        auto it = index.find(kmer_key.get_key());
//        if (it != index.end()) {
//            Run run = {interval.first, it->second};
//
//            queue.push( {run, interval.second - interval.first + 1});
//        }
//        return;
//    }
//
//    for ( char base : {'A', 'C', 'G', 'T'}) {
//        if (interval.first <= interval.second) {
//            kmers_to_bplustree_worker(idx, queue, index, k, idx.LF(interval, base), base + current_kmer);
//        }
//    }
//}
//
//void parallel_kmers_to_bplustree(FastLocate &idx, BplusTree <Run> &bptree,
//                                 hash_map <gbwtgraph::Key64::value_type, gbwtgraph::Position> &index, size_t k,
//                                 range_t interval) {
//    // Thread-safe queue to collect results
//    ThreadSafeQueue <std::pair<Run, size_t>> queue;
//
//    // number of threads
//    int threads = omp_get_max_threads();
//    // splitting the starting range (0, idx.bwt_size() - 1) into parts and call the kmers_to_bplustree_worker function
//    // for each part
//    size_t part_size = idx.bwt_size();
//    part_size = (part_size - 1) / threads;
//
//    std::cerr << "running each part " << part_size << " using threads " << threads << std::endl;
//
//    std::vector<string> starting_kmers = {"A", "C", "G", "T"};
//#pragma omp parallel for
//    for (int i = 0; i < threads; i++) {
////        size_t start = i * part_size;
////        size_t end = (i + 1) * part_size - 1;
////        if (i == threads - 1) {
////            end = idx.bwt_size() - 1;
////        }
////        kmers_to_bplustree_worker(idx, queue, index, k, {start, end}, "");
//        kmers_to_bplustree_worker(idx, queue, index, k, interval, starting_kmers[i]);
//    }
//
//    // Single-threaded insertion into BPlusTree
//    std::pair <Run, size_t> result;
//    while (queue.try_pop(result)) {
//        bptree.insert(result.first, result.second);
//    }
//}

    void kmers_to_bplustree_worker(FastLocate &idx, ThreadSafeQueue<std::tuple<range_t, std::string>> &task_queue,
    hash_map<gbwtgraph::Key64::value_type, gbwtgraph::Position> &index,
            ThreadSafeQueue<std::pair<Run, size_t>> &result_queue, size_t k) {
        std::tuple<range_t, std::string> task;

        // Process tasks from the task queue
        while (task_queue.try_pop(task)) {
            range_t interval = std::get<0>(task);
            std::string current_kmer = std::get<1>(task);

            if (current_kmer.length() == k && interval.first <= interval.second) {
                // k-mer is complete, now handle the result
                gbwtgraph::Key64 kmer_key = gbwtgraph::Key64::encode(current_kmer);
                auto it = index.find(kmer_key.get_key());
                if (it != index.end()) {
                Run run = {interval.first, it->second};
                result_queue.push({run, interval.second - interval.first + 1});
            }
            continue;
            }

        // Try to extend the k-mer if the interval is valid
            for (char base : {'A', 'C', 'G', 'T'}) {
                if (interval.first <= interval.second) {
                    range_t new_interval = idx.LF(interval, base);
                    if (new_interval.first <= new_interval.second) {
                        // Push the new task to the task queue
                        task_queue.push(std::make_tuple(new_interval, base + current_kmer));
                    }
                }
            }
    }
}

void parallel_kmers_to_bplustree(FastLocate &idx, BplusTree<Run> &bptree,
                                 hash_map<gbwtgraph::Key64::value_type, gbwtgraph::Position> &index,
                                 size_t k, range_t initial_interval) {
    // Thread-safe queues for tasks and results
    ThreadSafeQueue<std::tuple<range_t, std::string>> task_queue;
    ThreadSafeQueue<std::pair<Run, size_t>> result_queue;

    // Initial tasks for each base
    for (char base : {'A', 'C', 'G', 'T'}) {
        range_t interval = idx.LF(initial_interval, base);
        if (interval.first <= interval.second) {
            task_queue.push(std::make_tuple(interval, std::string(1, base)));
        }
    }

    // Parallel workers to process tasks
    int threads = omp_get_max_threads();
#pragma omp parallel num_threads(threads)
    {
        kmers_to_bplustree_worker(idx, task_queue, index, result_queue, k);
    }

    // Single-threaded insertion into BPlusTree
    std::pair<Run, size_t> result;
    while (result_queue.try_pop(result)) {
        bptree.insert(result.first, result.second);
    }
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
                range_t current_interval = {current_starting_pos, interval_end - 1};

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
    int traverse = 0;

    vector<int> tmp;

    vector <Run> tmp1;
    omp_lock_t lock;
    omp_init_lock(&lock);

    cerr << "Filling the gaps on the bptree" << endl;

#pragma omp parallel for
    for (int seq_num = 0; seq_num < number_of_sequences; ++seq_num) {
//        cerr << "running for sequence number " << seq_num << endl;
//        cerr << traverse << endl;

        auto seq_graph_nodes = gbz.index.extract(seq_num * 2);
        auto bwt_index = end_of_seq[seq_num].first;

        auto current_nodes_index = seq_graph_nodes.size() - 1;
        auto current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
        auto in_node_index = gbz.graph.get_length(current_node) - 1;

        vector <Run> local_tmp1;

        // traversing the RLBWT of a sequence
        while (true) {
#pragma omp atomic
            traverse++;


            // moving backwards
            bwt_index = idx.LF(bwt_index);
            auto first = idx.F_at(bwt_index);
            if (first == ENDMARKER) { // TODO: check this
//                cerr << "The end of the sequence at bwt index " << bwt_index << endl;
                break;
            }

            // make sure there is still room to traverse on the current node
            if (in_node_index == -1) {
                if (current_nodes_index == 0) {
                    break;
                }
                current_nodes_index--;
                current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
                in_node_index = gbz.graph.get_length(current_node) - 1;
            }

            // traverse the nodes on the graph to get the same base
            assert(first == gbz.graph.get_base(current_node, in_node_index));

            in_node_index--;

            // have to check if the current bwt position is in the bptree or not
            // calling search function if it return a gap run then this position is not currently in the tree and have to add it
            // if it returns a non-gap run then this position is already in the tree and we can continue traversing the bwt

            auto bptree_search = bptree.search(bwt_index);

            // the current position is not in the tree
            if (bptree_search.graph_position.value == 0) {
                // adding the current position to the tree
                pos_t current_pos = pos_t{gbz.graph.get_id(current_node), gbz.graph.get_is_reverse(current_node),
                                          in_node_index + 1};

                Run current_run = {bwt_index, gbwtgraph::Position::encode(current_pos)};
                local_tmp1.push_back(current_run);

            } else {
                // not adding the current position to the tree however checking if the tree position and the current graph
                // positions are the same
                pos_t current_pos = pos_t{gbz.graph.get_id(current_node), gbz.graph.get_is_reverse(current_node),
                                          in_node_index + 1};

                assert(gbwtgraph::Position::encode(current_pos).value == bptree_search.graph_position.value);
            }
        }

        // Lock to update shared data
        omp_set_lock(&lock);
        tmp1.insert(tmp1.end(), local_tmp1.begin(), local_tmp1.end());
        omp_unset_lock(&lock);
    }

    cerr << "The traverse completed" << endl;

    // Sort before inserting into the bptree
    sort(tmp1.begin(), tmp1.end());

    // insert all the items in the tmp1 to the bptree
    for (auto &i: tmp1) {
        bptree.insert(i, 1);
    }

    omp_destroy_lock(&lock);
}

}

#endif //PANGENOME_INDEX_ALGORITHM_HPP


//
// Created by seeskand on 9/16/24.
//
//
//#ifndef PANGENOME_INDEX_ALGORITHM_HPP
//#define PANGENOME_INDEX_ALGORITHM_HPP
//
//
//#include <vector>
//#include <queue>
//#include <mutex>
//#include <condition_variable>
//#include <omp.h>
//#include <hash_map.hpp>
//#include "../r-index/internal/r_index.hpp"
//#include "bplus_tree.hpp"
//#include "gbwtgraph/gbz.h"
//#include <hash_map.hpp>
//
//
//using namespace std;
//using namespace ri;
//using namespace gbwtgraph;
//
//
//namespace panindexer {
//
//
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
//
//// This function input is the OCC vector of the end of the sequences (#) and it returns the sorted end_of_seq vector
//// which is the sorted vector of pairs (i, SA[i]) for the end of each sequence which is (ISA[j], j)
//    vector <pair<uint64_t, uint64_t>> sort_end_of_seq(vector <pair<uint64_t, uint64_t>> &OCC) {
//
//        // Sort the end_of_seq vector by the second element of each pair
//        sort(OCC.begin(), OCC.end(),
//             [](const pair <uint64_t, uint64_t> &a, const pair <uint64_t, uint64_t> &b) {
//                 return a.second < b.second;
//             });
//
//        return OCC;
//    }
//
//    void kmers_to_bplustree_worker(r_index<> &idx, ThreadSafeQueue<std::pair < Run, size_t>> &queue,
//    hash_map <gbwtgraph::Key64::value_type, gbwtgraph::Position> &index,
//            size_t k, range_t interval, const string &current_kmer) {
//    if (current_kmer.length() == k && interval.first <= interval.second) {
//
//    // creating the kmer with the key type
//    gbwtgraph::Key64 kmer_key = gbwtgraph::Key64::encode(current_kmer);
//    // check if the kmer_key is in the index and if it is add the run to the queue
//    auto it = index.find(kmer_key.get_key());
//    if (it != index.end()) {
//    Run run = {interval.first, it->second};
//    queue.push( {run, interval.second - interval.first + 1});
//}
//return;
//}
//
//for ( char base : {'A', 'C', 'G', 'T'}) {
//if (interval.first <= interval.second) {
//kmers_to_bplustree_worker(idx, queue, index, k, idx.LF(interval, base), base + current_kmer);
//}
//}
//}
//
//void parallel_kmers_to_bplustree(r_index<> &idx, BplusTree <Run> &bptree,
//                                 hash_map <gbwtgraph::Key64::value_type, gbwtgraph::Position> &index, size_t k,
//                                 range_t interval) {
//    // Thread-safe queue to collect results
//    ThreadSafeQueue <std::pair<Run, size_t>> queue;
//
//    // number of threads
//    int threads = omp_get_max_threads();
//    // splitting the starting range (0, idx.bwt_size() - 1) into parts and call the kmers_to_bplustree_worker function
//    // for each part
//    size_t part_size = (idx.bwt_size() - 1) / threads;
//#pragma omp parallel for
//    for (int i = 0; i < threads; i++) {
//        size_t start = i * part_size;
//        size_t end = (i + 1) * part_size - 1;
//        if (i == threads - 1) {
//            end = idx.bwt_size() - 1;
//        }
//        kmers_to_bplustree_worker(idx, queue, index, k, {start, end}, "");
//    }
//
//    // Single-threaded insertion into BPlusTree
//    std::pair <Run, size_t> result;
//    while (queue.try_pop(result)) {
//        bptree.insert(result.first, result.second);
//    }
//}
//
//vector <pair<Run, size_t>>
//extend_kmers_bfs_parallel(GBWTGraph &graph, r_index<> &idx, BplusTree <Run> &bptree, int batch_size) {
//    vector <pair<Run, size_t>> extension_candidates;
//
//    int num_threads = omp_get_max_threads();
//    vector < std::queue < pair < Run, size_t > >> bfs_queues(num_threads);
//    vector < vector < pair < Run, size_t>>> batches(num_threads);
//
//    // Mutex for synchronizing access to bptree
//    std::mutex bptree_mutex;
//
//    // Add initial runs to the queues
//    int item_count = 0;
//    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
//        Run current_item = *it;
//        if (current_item.graph_position.value != 0) {
//            auto next_it = it;
//            ++next_it;
//            if (next_it != bptree.end()) {
//                Run next_item = *next_it;
//                bfs_queues[item_count % num_threads].push(make_pair(current_item, next_item.start_position));
//                ++item_count;
//            }
//        }
//    }
//
//    // Parallel BFS to extend kmers
//#pragma omp parallel
//    {
//        int thread_num = omp_get_thread_num();
//        while (!bfs_queues[thread_num].empty()) {
//            pair <Run, size_t> current_pair = bfs_queues[thread_num].front();
//            bfs_queues[thread_num].pop();
//
//            Run current_item = current_pair.first;
//            size_t interval_end = current_pair.second;
//
//            if (current_item.graph_position.value != 0) {
//                auto current_starting_pos = current_item.start_position;
//                auto current_graph_pos = current_item.graph_position;
//
//                // Decode the Position to a pos_t
//                pos_t current_pos = current_graph_pos.decode();
//
//                // Get the traversal of the graph position
//                handle_t current_handle = graph.get_handle(id(current_pos), false);
//
//                // The BWT interval of the current kmer
//                range_t current_interval = {current_starting_pos, interval_end - 1};
//
//                if (!is_rev(current_pos)) {
//                    if (offset(current_pos) > 0) {
//                        auto prev_base = graph.get_base(graph.get_handle(id(current_pos), is_rev(current_pos)),
//                                                        offset(current_pos) - 1);
//                        auto prev_graph_pos = pos_t{id(current_pos), is_rev(current_pos),
//                                                    offset(current_pos) - 1};
//
//                        auto new_range = idx.LF(current_interval, prev_base);
//
//                        if (new_range.first <= new_range.second) {
//                            Run temp_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
//                            batches[thread_num].push_back(make_pair(temp_run, new_range.second - new_range.first + 1));
//
//                            if (batches[thread_num].size() >= batch_size) {
//                                std::lock_guard <std::mutex> lock(bptree_mutex);
//                                for (const auto &item: batches[thread_num]) {
////                                    bptree.insert(item.first, item.second);
////                                    bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
//                                    if (bptree.insert_success(item.first, item.second)) {
//                                        bfs_queues[thread_num].push(
//                                                make_pair(item.first, item.second + item.first.start_position));
//                                    }
//                                }
//                                batches[thread_num].clear();
//                            }
//                        }
//                    } else {
//                        int prev_bases_num = 0;
//                        handle_t prev_node;
//                        char prev_base;
//                        pos_t prev_graph_pos;
//
//                        graph.follow_edges(current_handle, true, [&](const handle_t &prev) {
//                            prev_bases_num++;
//                            if (prev_bases_num != 1) return false;
//                            prev_base = graph.get_base(prev, graph.get_length(prev) - 1);
//                            prev_node = prev;
//                            return true;
//                        });
//
//                        if (prev_bases_num == 1) {
//                            prev_graph_pos = pos_t{graph.get_id(prev_node), graph.get_is_reverse(prev_node),
//                                                   graph.get_length(prev_node) - 1};
//                            auto new_range = idx.LF(current_interval, prev_base);
//                            if (new_range.first <= new_range.second) {
//                                Run new_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
//                                batches[thread_num].push_back({new_run, new_range.second - new_range.first + 1});
//
//                                if (batches[thread_num].size() >= batch_size) {
//                                    std::lock_guard <std::mutex> lock(bptree_mutex);
//                                    for (const auto &item: batches[thread_num]) {
////                                        bptree.insert(item.first, item.second);
////                                        bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
//                                        if (bptree.insert_success(item.first, item.second)) {
//                                            bfs_queues[thread_num].push(
//                                                    make_pair(item.first, item.second + item.first.start_position));
//                                        }
//                                    }
//                                    batches[thread_num].clear();
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//
//        // Process remaining items in the local batch
//        if (!batches[thread_num].empty()) {
//            std::lock_guard <std::mutex> lock(bptree_mutex);
//            for (const auto &item: batches[thread_num]) {
////                bptree.insert(item.first, item.second);
////                bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
//                if (bptree.insert_success(item.first, item.second)) {
//                    bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
//                }
//            }
//            batches[thread_num].clear();
//        }
//    }
//
//    return extension_candidates;
//}
//
//
//void traverse_sequences_parallel(GBZ &gbz, BplusTree <Run> &bptree, r_index<> &idx,
//                                 vector <pair<uint64_t, uint64_t>> &end_of_seq) {
//    auto number_of_sequences = end_of_seq.size();
//    int traverse = 0;
//
//    vector<int> tmp;
//
//    vector <Run> tmp1;
//    omp_lock_t lock;
//    omp_init_lock(&lock);
//
//    cerr << "Filling the gaps on the bptree" << endl;
//
//#pragma omp parallel for
//    for (int seq_num = 0; seq_num < number_of_sequences; ++seq_num) {
////        cerr << "running for sequence number " << seq_num << endl;
////        cerr << traverse << endl;
//
//        auto seq_graph_nodes = gbz.index.extract(seq_num * 2);
//        auto bwt_index = end_of_seq[seq_num].first;
//
//        auto current_nodes_index = seq_graph_nodes.size() - 1;
//        auto current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
//        auto in_node_index = gbz.graph.get_length(current_node) - 1;
//
//        vector <Run> local_tmp1;
//
//        // traversing the RLBWT of a sequence
//        while (true) {
//#pragma omp atomic
//            traverse++;
//
//            // moving backwards
//            bwt_index = idx.LF(bwt_index);
//            auto first = idx.F_at(bwt_index);
//            if (first == '$') {
////                cerr << "The end of the sequence at bwt index " << bwt_index << endl;
//                break;
//            }
//
//            // make sure there is still room to traverse on the current node
//            if (in_node_index == -1) {
//                if (current_nodes_index == 0) {
//                    break;
//                }
//                current_nodes_index--;
//                current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
//                in_node_index = gbz.graph.get_length(current_node) - 1;
//            }
//
//            // traverse the nodes on the graph to get the same base
//            assert(first == gbz.graph.get_base(current_node, in_node_index));
//
//            in_node_index--;
//
//            // have to check if the current bwt position is in the bptree or not
//            // calling search function if it return a gap run then this position is not currently in the tree and have to add it
//            // if it returns a non-gap run then this position is already in the tree and we can continue traversing the bwt
//
//            auto bptree_search = bptree.search(bwt_index);
//
//            // the current position is not in the tree
//            if (bptree_search.graph_position.value == 0) {
//                // adding the current position to the tree
//                pos_t current_pos = pos_t{gbz.graph.get_id(current_node), gbz.graph.get_is_reverse(current_node),
//                                          in_node_index + 1};
//
//                Run current_run = {bwt_index, gbwtgraph::Position::encode(current_pos)};
//                local_tmp1.push_back(current_run);
//
//            } else {
//                // not adding the current position to the tree however checking if the tree position and the current graph
//                // positions are the same
//                pos_t current_pos = pos_t{gbz.graph.get_id(current_node), gbz.graph.get_is_reverse(current_node),
//                                          in_node_index + 1};
//
//                assert(gbwtgraph::Position::encode(current_pos).value == bptree_search.graph_position.value);
//            }
//        }
//
//        // Lock to update shared data
//        omp_set_lock(&lock);
//        tmp1.insert(tmp1.end(), local_tmp1.begin(), local_tmp1.end());
//        omp_unset_lock(&lock);
//    }
//
//    cerr << "The traverse completed" << endl;
//
//    // Sort before inserting into the bptree
//    sort(tmp1.begin(), tmp1.end());
//
//    // insert all the items in the tmp1 to the bptree
//    for (auto &i: tmp1) {
//        bptree.insert(i, 1);
//    }
//
//    omp_destroy_lock(&lock);
//}
//
//}
//
//#endif //PANGENOME_INDEX_ALGORITHM_HPP