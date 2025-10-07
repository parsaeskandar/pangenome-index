//
// Created by Parsa Eskandar on 10/10/24.
//


#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>

#include "../include/pangenome_index/r-index.hpp"
#include "../deps/grlBWT/scripts/fm_index.h"
#include "../include/pangenome_index/algorithm.hpp"


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

TEST(RINDEX_Test, Locate_small_test_encoded) {
    std::string filename = "../test_data/small_test_nl.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/small_test_nl.rl_bwt";
    // Build legacy, serialize encoded, load encoded
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_small_test.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    auto x = r_index.decompressDA_encoded();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the encoded r-index";
}

TEST(RINDEX_Test, Locate_medium_test) {
    std::string filename = "../test_data/med_test.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/med_test.rl_bwt";
    FastLocate r_index(rlbwt_file);

    auto x = r_index.decompressDA();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

}

TEST(RINDEX_Test, Locate_medium_test_encoded) {
    std::string filename = "../test_data/med_test.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/med_test.rl_bwt";
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_med_test.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    auto x = r_index.decompressDA_encoded();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the encoded r-index";
}

TEST(RINDEX_Test, Locate_big_test) {
    std::string filename = "../test_data/x.newline_separated";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/x.rl_bwt";
    FastLocate r_index(rlbwt_file);

    auto x = r_index.decompressDA();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

}

TEST(RINDEX_Test, Locate_big_test_encoded) {
    std::string filename = "../test_data/x.newline_separated";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/x.rl_bwt";
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_big_test.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    auto x = r_index.decompressDA_encoded();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the encoded r-index";
}

TEST(RINDEX_Test, Locate_N_test) {
    std::string filename = "../test_data/N_test.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/N_test.rl_bwt";
    FastLocate r_index(rlbwt_file);

    auto x = r_index.decompressDA();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

}

TEST(RINDEX_Test, Locate_N_test_encoded) {
    std::string filename = "../test_data/N_test.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/N_test.rl_bwt";
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_N_test.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    auto x = r_index.decompressDA_encoded();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the encoded r-index";
}

// TEST(FMINDEX_Test, small_test){
//     std::string filename = "../big_test/merged_info.rl_bwt";
//     fm_index index(filename);
//     FastLocate r_index(filename);

//     for (int j = 0; j < 3; j++){
//         auto start_fm = j;
//         auto start_r = j;
//         auto a = index.lf(start_fm);
//         auto b = r_index.psi(start_r);

//         start_fm = a.second;
//         auto fm_char = a.first;
//         start_r = b.second;
//         auto r_char = b.first;

//         int i = 0;
//         while (r_char != NENDMARKER && fm_char != NENDMARKER){
// //                std::cerr << i << " " << start_fm << " " << start_r << std::endl;
//             ASSERT_EQ(fm_char, r_char) << "Invalid LF results from the FM-index and the r-index";
//             ASSERT_EQ(start_fm, start_r) << "Invalid LF results from the FM-index and the r-index";
//             a = index.lf(start_fm);
//             b = r_index.psi(start_r);
//             start_fm = a.second;
//             fm_char = a.first;
//             start_r = b.second;
//             r_char = b.first;
//             i++;
//         }



//     }



// }





TEST(FMDINDEX_Test, BackwardExtensionMatchesLF) {
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";
    FastLocate r_index(rlbwt_file);
//    r_index.initialize_complement_table();

    std::string kmer = "ATCAAAGAAAAAAGCCCAACATATCCATTACCATTACTAGTTACACATAGCATCAGGAACCAGAGAGTTGGA";
    std::string revcomp = kmer;
//    std::reverse(revcomp.begin(), revcomp.end());
//    for (char& c : revcomp) {
//        c = r_index.complement(c);
//    }

//    std::cerr << "Testing kmer: " << kmer << " and its reverse complement: " << revcomp << std::endl;

    // Initial full-range interval
    panindexer::FastLocate::bi_interval bint = {0, 0, r_index.bwt_size()};

    for (int i = kmer.size() - 1; i >= 0; --i) {
        char fwd_char = kmer[i];
        char rev_char = revcomp[i];

        // Compute manually via LF
        range_type fwd_expected = r_index.LF({bint.forward, bint.forward + bint.size - 1}, fwd_char);
        range_type rev_expected = r_index.LF({bint.reverse, bint.reverse + bint.size}, rev_char);

        size_t expected_size = 0;
        if (fwd_expected.first <= fwd_expected.second) {
            expected_size = fwd_expected.second - fwd_expected.first + 1;
        }

        panindexer::FastLocate::bi_interval extended = r_index.backward_extend(bint, fwd_char);


//        std::cerr << "char fwd: " << fwd_char << " char rev: " << rev_char << ", Interval size: " << extended.size << " Interval reverse " <<  extended.reverse << " Interval forward " << extended.forward << "\n";
        // Debug
//        std::cerr << "Char: " << fwd_char
//        << " | Expected size: " << expected_size
//        << " | Actual size: " << extended.size << std::endl;

        if (expected_size > 0){
            ASSERT_EQ(extended.size, expected_size) << "Mismatch in size for char " << fwd_char;
            ASSERT_EQ(extended.forward, fwd_expected.first);
        }

//        ASSERT_EQ(extended.reverse, rev_expected.first);

        // Move to next iteration
        bint = extended;
    }
}

TEST(FMDINDEX_Test, BackwardExtensionMatchesLF_Encoded) {
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";
    // Build and load encoded
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_bwd_ext.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    std::string kmer = "ATCAAAGAAAAAAGCCCAACATATCCATTACCATTACTAGTTACACATAGCATCAGGAACCAGAGAGTTGGA";
    panindexer::FastLocate::bi_interval bint = {0, 0, r_index.bwt_size()};

    for (int i = kmer.size() - 1; i >= 0; --i) {
        char fwd_char = kmer[i];
        range_type fwd_expected = r_index.LF_encoded({bint.forward, bint.forward + bint.size - 1}, fwd_char);
        size_t expected_size = 0;
        if (fwd_expected.first <= fwd_expected.second) {
            expected_size = fwd_expected.second - fwd_expected.first + 1;
        }
        panindexer::FastLocate::bi_interval extended = r_index.backward_extend_encoded(bint, fwd_char);
        if (expected_size > 0){
            ASSERT_EQ(extended.size, expected_size) << "Mismatch in size for char " << fwd_char;
            ASSERT_EQ(extended.forward, fwd_expected.first);
        }
        bint = extended;
    }
}


TEST(FMDINDEX_Test, CompareSampledKmersWithReverseComplementsBIGTEST) {
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";
    std::string text_file = "../test_data/big_test/merged_info";

    FastLocate r_index(rlbwt_file);
    r_index.initialize_complement_table();

    std::ifstream in(text_file);
    ASSERT_TRUE(in.is_open()) << "Could not open input text file";

    std::string full_text;
    std::string line;
    while (std::getline(in, line)) {
    if (!line.empty()) full_text += line;
    }

    ASSERT_GE(full_text.size(), 20) << "Input text is too short to sample from";

    const size_t k = 12;
    const size_t num_samples = 100;
    std::mt19937 rng(42); // fixed seed for reproducibility
    std::uniform_int_distribution<size_t> dist(0, full_text.size() - k);

    for (size_t i = 0; i < num_samples; ++i) {
    size_t pos = dist(rng);
    std::string kmer = full_text.substr(pos, k);

    // Skip invalid kmers with non-ACGTN symbols
    if (kmer.find_first_not_of("ACGTN") != std::string::npos) {
    i--;
    continue;
    }

    std::string revcomp = kmer;
    std::reverse(revcomp.begin(), revcomp.end());
    for (char& c : revcomp) {
    c = r_index.complement(c);
    }

    // Backward extend original
    panindexer::FastLocate::bi_interval int_kmer = {0, 0, r_index.bwt_size()};
    for (int j = kmer.size() - 1; j >= 0; --j) {
    int_kmer = r_index.backward_extend(int_kmer, kmer[j]);
    }

    // Backward extend reverse complement
    panindexer::FastLocate::bi_interval int_rc = {0, 0, r_index.bwt_size()};
    for (int j = revcomp.size() - 1; j >= 0; --j) {
    int_rc = r_index.backward_extend(int_rc, revcomp[j]);
    }

//    std::cerr << "k-mer       : " << kmer << ", Interval size: " << int_kmer.size << " Interval reverse " <<  int_kmer.reverse << " Interval forward " << int_kmer.forward << "\n";
//    std::cerr << "RevComp     : " << revcomp << ", Interval size: " << int_rc.size << " Interval reverse " <<  int_rc.reverse << " Interval forward " << int_rc.forward << "\n";

    ASSERT_EQ(int_kmer.forward, int_rc.reverse) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
    ASSERT_EQ(int_kmer.reverse, int_rc.forward) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
    ASSERT_EQ(int_kmer.size, int_rc.size) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";

    }
}

TEST(FMDINDEX_Test, CompareSampledKmersWithReverseComplementsBIGTEST_Encoded) {
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_big_rc.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    r_index.initialize_complement_table();

    std::ifstream in_txt("../test_data/big_test/merged_info");
    ASSERT_TRUE(in_txt.is_open()) << "Could not open input text file";

    std::string full_text, line;
    while (std::getline(in_txt, line)) { if (!line.empty()) full_text += line; }

    ASSERT_GE(full_text.size(), 20) << "Input text is too short to sample from";

    const size_t k = 12;
    const size_t num_samples = 100;
    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> dist(0, full_text.size() - k);

    for (size_t i = 0; i < num_samples; ++i) {
        size_t pos = dist(rng);
        std::string kmer = full_text.substr(pos, k);
        if (kmer.find_first_not_of("ACGTN") != std::string::npos) { i--; continue; }
        std::string revcomp = kmer;
        std::reverse(revcomp.begin(), revcomp.end());
        for (char& c : revcomp) { c = r_index.complement(c); }

        panindexer::FastLocate::bi_interval int_kmer = {0, 0, r_index.bwt_size()};
        for (int j = kmer.size() - 1; j >= 0; --j) { int_kmer = r_index.backward_extend_encoded(int_kmer, kmer[j]); }

        panindexer::FastLocate::bi_interval int_rc = {0, 0, r_index.bwt_size()};
        for (int j = revcomp.size() - 1; j >= 0; --j) { int_rc = r_index.backward_extend_encoded(int_rc, revcomp[j]); }

        ASSERT_EQ(int_kmer.forward, int_rc.reverse) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
        ASSERT_EQ(int_kmer.reverse, int_rc.forward) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
        ASSERT_EQ(int_kmer.size, int_rc.size) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
    }
}



}