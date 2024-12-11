//
// Created by Parsa Eskandar on 10/10/24.
//


#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "../include/pangenome_index/r-index.hpp"
#include "../deps/grlBWT/scripts/fm_index.h"


using namespace panindexer;

namespace {

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

    struct Rotation {
        std::string rotation;
        int seq_index;
    };

// Function to create BWT from a given string and track the sequence index
    std::pair <std::string, std::vector<unsigned long long>>
    createBWTWithSequenceInfo(const std::string &input, const std::vector<int> &seq_indices) {
        int n = input.size();
        std::vector <Rotation> rotations;

        // Generate all rotations of the input string and associate each rotation with a sequence index
        for (int i = 0; i < n; ++i) {
            rotations.push_back({input.substr(i) + input.substr(0, i), seq_indices[i]});
        }

        std::sort(rotations.begin(), rotations.end(), [](const Rotation &a, const Rotation &b) {
            return a.rotation < b.rotation;
        });

        // Create BWT by taking the last column of sorted rotations
        std::string bwt;
        std::vector<unsigned long long> result_indices;

        for (const auto &rotation: rotations) {
            bwt += rotation.rotation.back();             // Collect BWT
            result_indices.push_back(rotation.seq_index); // Collect corresponding sequence index
        }


        return {bwt, result_indices};
    }

// Function to read strings from a file, concatenate them using $_i, and create BWT with sequence index tracking
    std::vector<unsigned long long> readFileAndCreateBWTWithIndices(const std::string &filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open the file!" << std::endl;
            return {"", {}};
        }

        std::string line;
        std::string concatenatedString;
        std::vector<int> seq_indices;
        int sequence_number = 0;


        char added_char = '$';
        while (std::getline(file, line)) {


            // Mark the characters of the current string as belonging to this sequence
            for (size_t i = 0; i < line.size(); ++i) {
                seq_indices.push_back(sequence_number);
            }

            // Add $ between strings

            seq_indices.push_back(sequence_number);
            concatenatedString += line;
            concatenatedString += added_char;


            added_char += 1;
            ++sequence_number;
        }

        file.close();

        return createBWTWithSequenceInfo(concatenatedString, seq_indices).second;
    }



    TEST(RINDEX_Test, Locate_small_test) {
        std::string filename = "../test_data/small_test_nl.txt";
        auto sequence_indices = readFileAndCreateBWTWithIndices(filename);


        std::string rlbwt_file = "../test_data/small_test_nl.rl_bwt";
        FastLocate r_index(rlbwt_file);


        auto x = r_index.decompressDA();
        ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

    }

    TEST(RINDEX_Test, Locate_medium_test) {
        std::string filename = "../test_data/med_test.txt";
        auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

        std::string rlbwt_file = "../test_data/med_test.rl_bwt";
        FastLocate r_index(rlbwt_file);

        auto x = r_index.decompressDA();
        ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

    }

    TEST(RINDEX_Test, Locate_big_test) {
        std::string filename = "../test_data/x.newline_separated";
        auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

        std::string rlbwt_file = "../test_data/x.rl_bwt";
        FastLocate r_index(rlbwt_file);

        auto x = r_index.decompressDA();
        ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

    }

    TEST(RINDEX_Test, Locate_N_test) {
        std::string filename = "../test_data/N_test.txt";
        auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

        std::string rlbwt_file = "../test_data/N_test.rl_bwt";
        FastLocate r_index(rlbwt_file);

        auto x = r_index.decompressDA();
        ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

    }

    TEST(FMINDEX_Test, small_test){
        std::string filename = "../test_data/x.rl_bwt";
        fm_index index(filename);
        FastLocate r_index(filename);

        for (int j = 0; j < 3; j++){
            auto start_fm = j;
            auto start_r = j;
            auto a = index.lf(start_fm);
            auto b = r_index.psi(start_r);

            start_fm = a.second;
            auto fm_char = a.first;
            start_r = b.second;
            auto r_char = b.first;

            int i = 0;
            while (r_char != NENDMARKER && fm_char != NENDMARKER){
                std::cerr << i << " " << start_fm << " " << start_r << std::endl;
                ASSERT_EQ(fm_char, r_char) << "Invalid LF results from the FM-index and the r-index";
                ASSERT_EQ(start_fm, start_r) << "Invalid LF results from the FM-index and the r-index";
                a = index.lf(start_fm);
                b = r_index.psi(start_r);
                start_fm = a.second;
                fm_char = a.first;
                start_r = b.second;
                r_char = b.first;
                i++;
            }



        }












    }



}