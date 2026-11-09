//
// Created by seeskand on 3/26/25.
//
#include <gtest/gtest.h>
#include <algorithm>
#include <random>
#include <vector>

#include "../include/pangenome_index/bplus_tree.hpp"


using namespace panindexer;
using namespace gbwtgraph;


namespace {
    BplusTree<panindexer::Run> bptree(5);

    TEST(BplusTreeTest, TestOne) {
        size_t test_nums = 1000000;
        pos_t current_pos = pos_t{1, false, 0};

        // one run of length 1 at every 10th position
        for (size_t i = 0; i < test_nums; i++) {
            panindexer::Run current_run = {i * 10, gbwtgraph::Position(current_pos)};
            bptree.insert(current_run, 1);
        }

        // then fill the positions in between in a random order, all with the same graph position
        std::vector<size_t> indices;
        for (size_t i = 0; i < test_nums * 10 - 1; ++i) {
            if (i % 10 != 0) {
                indices.push_back(i);
            }
        }
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(indices.begin(), indices.end(), g);

        for (size_t i : indices) {
            panindexer::Run current_run = {i, gbwtgraph::Position(current_pos)};
            bptree.insert(current_run, 1);
        }

        // everything should be merged into one run and the gap after it
        std::vector<panindexer::Run> runs;
        for (auto it = bptree.begin(); it != bptree.end(); ++it) {
            runs.push_back(*it);
        }
        std::vector<panindexer::Run> results;
        results.push_back({0, gbwtgraph::Position(current_pos)});
        results.push_back({test_nums * 10 - 1, gbwtgraph::Position::no_pos()});
        ASSERT_EQ(runs, results) << "Invalid BPlusTree results";
    }
}
