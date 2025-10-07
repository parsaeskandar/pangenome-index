#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

#include "pangenome_index/tag_arrays.hpp"

using namespace std;
using namespace panindexer;
using handlegraph::pos_t;

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <input_algorithm_file> <output_compressed_file> [encoded_starts_tmp] [bwt_intervals_tmp]\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    const string input_file = argv[1];
    const string output_file = argv[2];

    string encoded_starts_tmp;
    string bwt_intervals_tmp;

    if (argc >= 4) {
        encoded_starts_tmp = argv[3];
    } else {
        encoded_starts_tmp = output_file + ".encoded_starts.tmp";
    }

    if (argc >= 5) {
        bwt_intervals_tmp = argv[4];
    } else {
        bwt_intervals_tmp = output_file + ".bwt_intervals.tmp";
    }

    // Read the entire input file (algorithm.hpp format: raw gbwt::ByteCode runs) into memory
    vector<gbwt::byte_type> encoded_bytes;
    try {
        ifstream in(input_file, ios::binary);
        if (!in) {
            cerr << "Error: cannot open input file: " << input_file << "\n";
            return 1;
        }
        in.seekg(0, ios::end);
        streampos sz = in.tellg();
        if (sz < 0) {
            cerr << "Error: failed to stat input file size\n";
            return 1;
        }
        encoded_bytes.resize(static_cast<size_t>(sz));
        in.seekg(0, ios::beg);
        if (!encoded_bytes.empty()) {
            in.read(reinterpret_cast<char*>(encoded_bytes.data()), encoded_bytes.size());
            if (!in) {
                cerr << "Error: failed to read input file: " << input_file << "\n";
                return 1;
            }
        }
        in.close();
    } catch (const std::exception& e) {
        cerr << "Exception while reading input: " << e.what() << "\n";
        return 1;
    }

    // Open outputs
    ofstream main_out(output_file, ios::binary | ios::trunc);
    if (!main_out) {
        cerr << "Error: cannot open output file: " << output_file << "\n";
        return 1;
    }

    // Reserve space for encoded_runs size header; will be filled in merge step
    size_t zero_size = 0;
    main_out.write(reinterpret_cast<const char*>(&zero_size), sizeof(zero_size));

    ofstream encoded_starts_out(encoded_starts_tmp, ios::binary | ios::trunc);
    if (!encoded_starts_out) {
        cerr << "Error: cannot open tmp file: " << encoded_starts_tmp << "\n";
        return 1;
    }

    ofstream bwt_intervals_out(bwt_intervals_tmp, ios::binary | ios::trunc);
    if (!bwt_intervals_out) {
        cerr << "Error: cannot open tmp file: " << bwt_intervals_tmp << "\n";
        return 1;
    }

    TagArray tag_array;

    // Decode runs and stream them into compressed serializer in batches
    const size_t batch_runs = 8192;
    vector<pair<pos_t, uint16_t>> batch;
    batch.reserve(batch_runs);

    uint64_t loc = 0;
    while (loc < encoded_bytes.size()) {
        gbwt::size_type decc = gbwt::ByteCode::read(encoded_bytes, loc);
        auto decoded = TagArray::decode_run(decc);
        batch.emplace_back(decoded.first, decoded.second);
        if (batch.size() >= batch_runs) {
            tag_array.compressed_serialize_compact(main_out, encoded_starts_out, bwt_intervals_out, batch);
            batch.clear();
        }
    }
    if (!batch.empty()) {
        tag_array.compressed_serialize_compact(main_out, encoded_starts_out, bwt_intervals_out, batch);
        batch.clear();
    }

    main_out.flush();
    main_out.close();
    encoded_starts_out.flush();
    encoded_starts_out.close();
    bwt_intervals_out.flush();
    bwt_intervals_out.close();

    // Merge sidecar tmp files into the main compressed index and write header
    try {
        tag_array.merge_compressed_files(output_file, encoded_starts_tmp, bwt_intervals_tmp);
    } catch (const std::exception& e) {
        cerr << "Error while merging compressed files: " << e.what() << "\n";
        return 1;
    }

    return 0;
}


