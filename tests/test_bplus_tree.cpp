//
// Created by seeskand on 3/26/25.
//
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>

#include "../include/pangenome_index/bplus_tree.hpp"


using namespace panindexer;
using namespace gbwtgraph;


namespace {
    BplusTree<panindexer::Run> bptree(5);

    TEST(BplusTreeTest, TestOne) {



        size_t test_nums = 1000000;
        pos_t current_pos = pos_t{1, false, 0};
        for (int i = 0; i < test_nums; i++) {
            size_t first = i * 10;


            panindexer::Run current_run = {first, gbwtgraph::Position::encode(current_pos)};

            bptree.insert(current_run, 1);
        }
        std::vector<int> indices;
        for (int i = 0; i < test_nums * 10 - 1; ++i) {
            if (i % 10 != 0) {
                indices.push_back(i);
            }
        }

        // Shuffle the list
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(indices.begin(), indices.end(), g);

        // Insert in randomized order
        for (int i : indices) {
            pos_t current_pos = pos_t{1, false, 0};
            panindexer::Run current_run = {static_cast<size_t>(i), gbwtgraph::Position::encode(current_pos)};
//            std::cerr << "Inserting " << i << std::endl;
            bptree.insert(current_run, 1);
//            bptree.print_whole_tree();
        }
        vector<panindexer::Run> runs;
        for (auto it = bptree.begin(); it != bptree.end(); ++it) {
            panindexer::Run current_item = *it;
            runs.push_back(current_item);
        }
        vector<panindexer::Run> results;
        panindexer::Run run1 = {0, gbwtgraph::Position::encode(current_pos)};
        panindexer::Run run2 = {test_nums * 10 - 1, gbwtgraph::Position::no_value()};
        results.push_back(run1);
        results.push_back(run2);
        ASSERT_EQ(runs, results) << "Invalid BPlusTree results";

    }
}
