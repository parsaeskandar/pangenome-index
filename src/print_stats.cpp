#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>

#include "pangenome_index/r-index.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include <sdsl/util.hpp>

using namespace std;
using namespace panindexer;

static void usage(const char* prog){
    cerr << "Usage: " << prog << " <r_index.ri> <compressed_tags.tags>\n";
}

static void print_human_size(const string& name, size_t bytes, size_t denom_bits = 0){
    double mb = static_cast<double>(bytes) / (1024.0 * 1024.0);
    cout << name << ": " << bytes << " bytes (" << mb << " MB)";
    if (denom_bits > 0) {
        double bpr = (bytes * 8.0) / static_cast<double>(denom_bits);
        cout << ", " << bpr << " bits/run";
    }
    cout << "\n";
}

int main(int argc, char** argv){
    if (argc != 3){
        usage(argv[0]);
        return 1;
    }

    const string ri_path = argv[1];
    const string tags_path = argv[2];

    // Load r-index
    FastLocate r_index;
    if (!sdsl::load_from_file(r_index, ri_path)){
        cerr << "Error: cannot load r-index from " << ri_path << "\n";
        return 1;
    }

    // Load compressed tags
    TagArray tag_array;
    ifstream tin(tags_path, ios::binary);
    if (!tin){
        cerr << "Error: cannot open tags file: " << tags_path << "\n";
        return 1;
    }
    tag_array.load_compressed_tags(tin);
    tin.close();

    // High-level summary
    size_t text_len = r_index.get_sequence_size();
    size_t r_runs = r_index.size();
    size_t t_runs = tag_array.number_of_runs_compressed();

    cout << "=== High-level ===\n";
    cout << "Total sequence length (BWT size): " << text_len << "\n";
    cout << "BWT runs (r-index): " << r_runs << "\n";
    cout << "Tag array runs: " << t_runs << "\n\n";

    // r-index components
    cout << "=== R-index components ===\n";
    size_t bytes_header = sizeof(r_index.header.tag) + sizeof(r_index.header.version)
                        + sizeof(r_index.header.max_length) + sizeof(r_index.header.flags);
    size_t bytes_samples = sdsl::size_in_bytes(r_index.samples);
    size_t bytes_last = sdsl::size_in_bytes(r_index.last);
    size_t bytes_last_to_run = sdsl::size_in_bytes(r_index.last_to_run);
    size_t bytes_sym_map = sdsl::size_in_bytes(r_index.sym_map);
    size_t bytes_C = sdsl::size_in_bytes(r_index.C);
    size_t bytes_blocks_start_pos = sdsl::size_in_bytes(r_index.blocks_start_pos);

    size_t bytes_blocks_char = 0;
    size_t bytes_blocks_runs = 0;
    for (const auto& blk : r_index.blocks){
        bytes_blocks_char += sdsl::size_in_bytes(blk.get_cum_ranks());
        bytes_blocks_runs += blk.get_run_nums() * sizeof(std::pair<size_t,size_t>);
    }
    size_t bytes_blocks_total = bytes_blocks_char + bytes_blocks_runs;

    size_t bytes_misc = sizeof(r_index.sequence_size) + sizeof(r_index.block_size) + sizeof(r_index.complement_table);

    size_t bytes_rindex_total = bytes_header + bytes_samples + bytes_last + bytes_last_to_run + bytes_sym_map + bytes_C
                              + bytes_blocks_start_pos + bytes_blocks_total + bytes_misc;

    print_human_size("header", bytes_header, r_runs);
    print_human_size("samples", bytes_samples, r_runs);
    print_human_size("last (sd_vector)", bytes_last, r_runs);
    print_human_size("last_to_run", bytes_last_to_run, r_runs);
    print_human_size("sym_map", bytes_sym_map, r_runs);
    print_human_size("C", bytes_C, r_runs);
    print_human_size("blocks_start_pos (sd_vector)", bytes_blocks_start_pos, r_runs);
    print_human_size("blocks.character_cum_ranks", bytes_blocks_char, r_runs);
    print_human_size("blocks.runs (pairs)", bytes_blocks_runs, r_runs);
    print_human_size("misc (seq_size, block_size, complement)", bytes_misc, r_runs);
    print_human_size("TOTAL r-index (approx)", bytes_rindex_total, r_runs);
    cout << "\n";

    // Tag arrays components (compressed)
    cout << "=== Tag arrays (compressed) components ===\n";
    size_t bytes_encoded = tag_array.bytes_encoded_runs();
    size_t bytes_starts_sd = tag_array.bytes_encoded_runs_starts_sd();
    size_t bytes_bwt_int = tag_array.bytes_bwt_intervals();
    size_t bytes_tags_total = tag_array.bytes_total_compressed();

    print_human_size("encoded_runs (ByteCode)", bytes_encoded, t_runs);
    print_human_size("encoded_runs_starts (sd_vector)", bytes_starts_sd, t_runs);
    print_human_size("bwt_intervals (sd_vector)", bytes_bwt_int, t_runs);
    print_human_size("TOTAL tag arrays (compressed)", bytes_tags_total, t_runs);

    return 0;
}


