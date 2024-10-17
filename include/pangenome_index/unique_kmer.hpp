//
// Created by seeskand on 9/16/24.
//

#ifndef PANGENOME_INDEX_UNIQUE_KMER_HPP
#define PANGENOME_INDEX_UNIQUE_KMER_HPP

#include <gbwtgraph/index.h>
#include <omp.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/utils.h>
#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/minimizer.h>
#include <hash_map.hpp>


namespace panindexer {
    using namespace std;
    using handlegraph::pos_t;

// Function to extract unique kmers from the graph in parallel
    // This function is a version of canonical_kmers function from gbwtgraph library but in this version we separate the kmer
// and its reverse complement. So this function returns the sorted vector of kmers
    template<class KeyType>
    std::vector <Kmer<KeyType>>
    forward_strand_kmers(std::string::const_iterator begin, std::string::const_iterator end, size_t k) {
        std::vector <Kmer<KeyType>> result;
        if (k == 0 || k > KeyType::KMER_MAX_LENGTH) {
            std::cerr << "The maximum kmer size is " << KeyType::KMER_MAX_LENGTH << std::endl;
            return result;
        }

        size_t valid_chars = 0, offset = 0;
        KeyType forward_key;
        std::string::const_iterator iter = begin;
        while (iter != end) {
            forward_key.forward(k, *iter, valid_chars);
            if (valid_chars >= k) {
                result.push_back({forward_key, forward_key.hash(), offset_type(offset - (k - 1)), false});
            }
            ++iter;
            offset++;
        }

        std::sort(result.begin(), result.end());
        return result;
    }

    template<class KeyType>
    std::vector <Kmer<KeyType>>
    forward_strand_kmers(const std::string &seq, size_t k) {
        return forward_strand_kmers<KeyType>(seq.begin(), seq.end(), k);
    }


// This function returns the unique kmers in the graph and stores them in the index
    template<class KeyType>
    void
    unique_kmers_parallel(const GBWTGraph &graph,
                          hash_map <gbwtgraph::Key64::value_type, gbwtgraph::Position> &index,
                          size_t k) {
        typedef KeyType key_type;
        typedef Kmer <key_type> kmer_type;
        constexpr
        size_t KMER_CACHE_SIZE = 1024;  // Adjust cache size as needed

        int threads = omp_get_max_threads();
        std::vector < std::vector < std::pair < gbwtgraph::Key64::value_type, gbwtgraph::Position>>> cache(threads);
        // make a duplicates cache to handle the duplicates faster
        hash_set <gbwtgraph::Key64::value_type> duplicates;


        // Lambda to flush the thread-local cache to the shared index
        auto flush_cache = [&](int thread_number) {
            auto current_cache = cache[thread_number];

            // Sort the cache by key
//        std::sort(current_cache.begin(), current_cache.end(), [](const auto &a, const auto &b) {
//            return a.first < b.first;
//        });
//        // Remove duplicates
//        auto last = std::unique(current_cache.begin(), current_cache.end(), [](const auto &a, const auto &b) {
//            return a.first == b.first;
//        });
//        current_cache.erase(last, current_cache.end());
#pragma omp critical
            {
                for (auto entry: current_cache) {
                    if (duplicates.find(entry.first) != duplicates.end()) {
                        // Key is a known duplicate, skip it
                        continue;
                    }
                    auto it = index.find(entry.first);
                    if (it == index.end()) {
                        index[entry.first] = entry.second;
                    } else if (it->second != entry.second) {
                        // remove the kmer from the index
                        index.erase(it);
                        duplicates.insert(entry.first);
                    }
                }
            }

            cache[thread_number].clear();  // Clear the cache after flushing
        };

        // Main lambda function to process kmers
        auto hash_kmers = [&](const std::vector <handle_t> &traversal, const std::string &seq) {
            int thread_id = omp_get_thread_num();
            std::vector <kmer_type> kmers = forward_strand_kmers<key_type>(seq, k);
//        cerr << kmers.size() << endl;
            auto iter = traversal.begin();
            size_t node_start = 0;

            for (auto kmer: kmers) {
                if (kmer.empty()) continue;


                size_t node_length = graph.get_length(*iter);
                while (node_start + node_length <= kmer.offset) {
                    node_start += node_length;
                    ++iter;
                    node_length = graph.get_length(*iter);
                }

                pos_t pos{graph.get_id(*iter), graph.get_is_reverse(*iter), kmer.offset - node_start};
//            cerr << " The pos is: " << gbwtgraph::id(pos) << " " << gbwtgraph::offset(pos) << " " << gbwtgraph::is_rev(pos) << endl;

//                if (debug) {
//                    if (gbwtgraph::Key64::encode("GACAAATCTGGGTTCAAATCCTCACTTTG") == kmer.key) {
//                        cerr << "The key is: " << kmer.key << " " << kmer.offset << " id " << graph.get_id(*iter)
//                             << " rev "
//                             << graph.get_is_reverse(*iter) << " kmer offset " << kmer.offset << " node_start "
//                             << node_start << " node_length " << node_length << endl;
//                    }
//
//                }

                if (!gbwtgraph::Position::valid_offset(pos)) {
#pragma omp critical (cerr)
                    {
                        std::cerr << "index_haplotypes(): Node offset " << offset(pos) << " is too large"
                                  << std::endl;
                    }
                    std::exit(EXIT_FAILURE);
                }
                // Use thread-local cache to reduce contention on the shared index
                cache[thread_id].emplace_back(kmer.key.get_key(), gbwtgraph::Position::encode(pos));



                // Flush the cache if it reaches the size limit
                if (cache[thread_id].size() >= KMER_CACHE_SIZE) {
                    flush_cache(thread_id);
                }
            }

        };

        // Parallel execution
//    gbwtgraph::for_each_nonredundant_window(graph, k, hash_kmers, true);
        gbwtgraph::for_each_haplotype_window(graph, k, hash_kmers, true);
        for (int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }

    }


}

#endif //PANGENOME_INDEX_UNIQUE_KMER_HPP
