#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include <sdsl/wavelet_trees.hpp>
#include <gbwtgraph/gbz.h>
#include <gbwt/metadata.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <sstream>
#include <limits>
#include <map>
#include <chrono>
#include <iomanip>
#include <atomic>

using namespace std;
using namespace panindexer;
using namespace std::chrono;
using namespace gbwtgraph;

static std::atomic<bool> debug{false};

struct TagResult {
    vector<size_t> query_offsets; // offsets within the query interval (relative to l)
    unordered_map<size_t, vector<size_t>> sequence_offsets; // seq_id -> list of offsets in that sequence
};

// Parse interval string like "1000..2000"
static pair<size_t, size_t> parse_interval(const string& interval_str) {
    size_t dot_pos = interval_str.find("..");
    if (dot_pos == string::npos) {
        throw runtime_error("Invalid interval format: expected 'start..end'");
    }
    size_t start = stoull(interval_str.substr(0, dot_pos));
    size_t end = stoull(interval_str.substr(dot_pos + 2));
    if (start > end) {
        throw runtime_error("Invalid interval: start > end");
    }
    return {start, end};
}

// Find sequence ID by direct path name (for generic paths like "x", "y", "chr1")
// Returns numeric_limits<size_t>::max() if not found
static size_t find_sequence_id_by_path_name(const gbwt::GBWT& gbwt_index,
                                             const string& path_name_str,
                                             bool debug_output = false) {
    if (!gbwt_index.hasMetadata()) {
        cerr << "Error: GBWT index does not have metadata" << endl;
        return numeric_limits<size_t>::max();
    }
    
    const gbwt::Metadata& metadata = gbwt_index.metadata;
    
    if (debug_output) {
        cerr << "Searching for path by name: '" << path_name_str << "'" << endl;
        cerr << "GBWT has " << metadata.paths() << " paths" << endl;
    }
    
    // For generic paths, the path name is typically stored as the contig name
    // First try to find a contig that matches the path name
    size_t contig_id = metadata.contig(path_name_str);
    
    if (contig_id < metadata.contigs()) {
        // Found a matching contig, now find paths with this contig
        if (debug_output) {
            cerr << "Found contig_id=" << contig_id << " for path name '" << path_name_str << "'" << endl;
        }
        
        vector<size_t> matching_paths;
        for (size_t path_id = 0; path_id < metadata.paths(); ++path_id) {
            gbwt::PathName pn = metadata.path(path_id);
            if (pn.contig == contig_id) {
                matching_paths.push_back(path_id);
                if (debug_output) {
                    cerr << "Found matching path: path_id=" << path_id 
                         << " (sample=" << pn.sample
                         << ", contig=" << pn.contig
                         << ", phase=" << pn.phase
                         << ", count=" << pn.count << ")" << endl;
                }
            }
        }
        
        if (!matching_paths.empty()) {
            size_t gbwt_path_id = matching_paths[0];
            if (matching_paths.size() > 1) {
                cerr << "Warning: Multiple paths (" << matching_paths.size() 
                     << ") found for contig '" << path_name_str << "'. Using first match (path_id=" 
                     << gbwt_path_id << ")" << endl;
            }
            // RLBWT seq 2i = forward strand of GBWT path i (see coordinate_translation PathSpec::resolve).
            size_t rlbwt_seq_id = 2 * gbwt_path_id;
            if (debug_output) {
                cerr << "Resolved path name '" << path_name_str << "': GBWT path_id=" << gbwt_path_id
                     << " -> RLBWT seq_id=" << rlbwt_seq_id << " (forward)" << endl;
            }
            return rlbwt_seq_id;
        }
    }
    
    // If contig lookup failed, list available paths
    cerr << "Error: Path '" << path_name_str << "' not found" << endl;
    cerr << "Available contigs:" << endl;
    for (size_t i = 0; i < metadata.contigs(); ++i) {
        cerr << "  " << metadata.contig(i) << endl;
    }
    
    return numeric_limits<size_t>::max();
}

// Find sequence ID from path metadata (sample, contig, haplotype)
// Returns numeric_limits<size_t>::max() if not found
static size_t find_sequence_id_from_metadata(const gbwt::GBWT& gbwt_index,
                                              const string& sample_name,
                                              const string& contig_name,
                                              size_t haplotype,
                                              bool debug_output = false) {
    if (!gbwt_index.hasMetadata()) {
        cerr << "Error: GBWT index does not have metadata" << endl;
        return numeric_limits<size_t>::max();
    }
    
    const gbwt::Metadata& metadata = gbwt_index.metadata;
    
    if (debug_output) {
        cerr << "Searching for path: sample='" << sample_name 
             << "', contig='" << contig_name 
             << "', haplotype=" << haplotype << endl;
        cerr << "GBWT has " << metadata.paths() << " paths, "
             << metadata.samples() << " samples, "
             << metadata.contigs() << " contigs" << endl;
    }
    
    // Find sample ID from sample name
    size_t sample_id = metadata.sample(sample_name);
    if (sample_id >= metadata.samples()) {
        cerr << "Error: Sample '" << sample_name << "' not found in metadata" << endl;
        cerr << "Available samples:" << endl;
        for (size_t i = 0; i < metadata.samples(); ++i) {
            cerr << "  " << i << ": " << metadata.sample(i) << endl;
        }
        return numeric_limits<size_t>::max();
    }
    
    // Find contig ID from contig name
    size_t contig_id = metadata.contig(contig_name);
    if (contig_id >= metadata.contigs()) {
        cerr << "Error: Contig '" << contig_name << "' not found in metadata" << endl;
        cerr << "Available contigs:" << endl;
        for (size_t i = 0; i < metadata.contigs(); ++i) {
            cerr << "  " << i << ": " << metadata.contig(i) << endl;
        }
        return numeric_limits<size_t>::max();
    }
    
    if (debug_output) {
        cerr << "Found sample_id=" << sample_id << " for '" << sample_name << "'" << endl;
        cerr << "Found contig_id=" << contig_id << " for '" << contig_name << "'" << endl;
    }
    
    // Search through all paths to find matching one
    vector<size_t> matching_paths;
    for (size_t path_id = 0; path_id < metadata.paths(); ++path_id) {
        gbwt::PathName path_name = metadata.path(path_id);
        
        if (path_name.sample == sample_id && 
            path_name.contig == contig_id &&
            path_name.phase == haplotype) {
            matching_paths.push_back(path_id);
            if (debug_output) {
                cerr << "Found matching path: path_id=" << path_id 
                     << " (sample=" << path_name.sample
                     << ", contig=" << path_name.contig
                     << ", phase=" << path_name.phase
                     << ", count=" << path_name.count << ")" << endl;
            }
        }
    }
    
    if (matching_paths.empty()) {
        cerr << "Error: No path found matching sample='" << sample_name 
             << "', contig='" << contig_name 
             << "', haplotype=" << haplotype << endl;
        return numeric_limits<size_t>::max();
    }
    
    size_t gbwt_path_id = matching_paths[0];
    if (matching_paths.size() > 1) {
        cerr << "Warning: Multiple paths (" << matching_paths.size() 
             << ") found matching the criteria. Using first match (path_id=" 
             << gbwt_path_id << ")" << endl;
    }
    
    // GBWT metadata path_id i maps to RLBWT forward sequence 2*i (reverse is 2*i+1).
    size_t seq_id = 2 * gbwt_path_id;
    
    if (debug_output) {
        cerr << "Resolved: GBWT path_id=" << gbwt_path_id
             << " -> RLBWT seq_id=" << seq_id << " (forward)" << endl;
    }
    
    return seq_id;
}

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <r_index.ri> <sampled.tags> [options]" << endl;
    cerr << endl;
    cerr << "Path specification (choose one):" << endl;
    cerr << "  --seq-id ID                  Direct RLBWT sequence ID (GBWT path i forward = 2*i)" << endl;
    cerr << "  --gbz FILE --path-name NAME  Specify path by name (for generic paths like 'x', 'chr1')" << endl;
    cerr << "  --gbz FILE --sample NAME --contig NAME [--haplotype N]" << endl;
    cerr << "                               Specify path by metadata (for haplotype paths)" << endl;
    cerr << endl;
    cerr << "Required:" << endl;
    cerr << "  --interval START..END        Query interval on the path" << endl;
    cerr << endl;
    cerr << "Examples:" << endl;
    cerr << "  " << prog << " index.ri tags.stags --seq-id 0 --interval 1000..2000" << endl;
    cerr << "  " << prog << " index.ri tags.stags --gbz graph.gbz --path-name x --interval 1000..2000" << endl;
    cerr << "  " << prog << " index.ri tags.stags --gbz graph.gbz --sample GRCh38 --contig chr1 --interval 1000..2000" << endl;
    cerr << "  " << prog << " index.ri tags.stags --gbz graph.gbz --sample HG002 --contig chr1 --haplotype 1 --interval 1000..2000" << endl;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }
    
    string r_index_file = argv[1];
    string sampled_tags_file = argv[2];
    
    // Parse optional arguments
    size_t seq_id = numeric_limits<size_t>::max();
    pair<size_t, size_t> interval = {0, 0};
    bool has_interval = false;
    
    // Path metadata options (alternative to --seq-id)
    string gbz_file;
    string path_name;     // For generic paths like "x", "y", "chr1"
    string sample_name;   // For haplotype paths
    string contig_name;   // For haplotype paths
    size_t haplotype = 0;
    bool has_path_name = false;
    bool has_sample = false;
    bool has_contig = false;
    bool has_haplotype = false;
    
    for (int i = 3; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--interval" && i + 1 < argc) {
            interval = parse_interval(argv[++i]);
            has_interval = true;
        } else if (arg == "--seq-id" && i + 1 < argc) {
            seq_id = stoull(argv[++i]);
        } else if (arg == "--gbz" && i + 1 < argc) {
            gbz_file = argv[++i];
        } else if (arg == "--path-name" && i + 1 < argc) {
            path_name = argv[++i];
            has_path_name = true;
        } else if (arg == "--sample" && i + 1 < argc) {
            sample_name = argv[++i];
            has_sample = true;
        } else if (arg == "--contig" && i + 1 < argc) {
            contig_name = argv[++i];
            has_contig = true;
        } else if (arg == "--haplotype" && i + 1 < argc) {
            haplotype = stoull(argv[++i]);
            has_haplotype = true;
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return 0;
        } else {
            cerr << "Unknown argument: " << arg << endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    // Validate arguments: need either --seq-id OR (--gbz + --path-name) OR (--gbz + --sample + --contig)
    bool use_path_name = has_path_name;
    bool use_metadata = has_sample || has_contig;
    bool use_seq_id = (seq_id != numeric_limits<size_t>::max());
    
    if ((use_seq_id && use_path_name) || (use_seq_id && use_metadata) || (use_path_name && use_metadata)) {
        cerr << "Error: Cannot mix --seq-id, --path-name, and --sample/--contig options" << endl;
        usage(argv[0]);
        return 1;
    }
    
    if (use_path_name) {
        if (gbz_file.empty()) {
            cerr << "Error: --gbz is required when using --path-name" << endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    if (use_metadata) {
        if (gbz_file.empty()) {
            cerr << "Error: --gbz is required when using --sample/--contig" << endl;
            usage(argv[0]);
            return 1;
        }
        if (!has_sample) {
            cerr << "Error: --sample is required when using --contig" << endl;
            usage(argv[0]);
            return 1;
        }
        if (!has_contig) {
            cerr << "Error: --contig is required when using --sample" << endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    if (!use_metadata && !use_seq_id && !use_path_name) {
        cerr << "Error: Must specify either --seq-id, --path-name, or (--sample + --contig)" << endl;
        usage(argv[0]);
        return 1;
    }
    
    if (!has_interval) {
        cerr << "Error: --interval is required" << endl;
        usage(argv[0]);
        return 1;
    }
    
    // If using path-name or metadata, load GBZ and resolve sequence ID
    if (use_path_name || use_metadata) {
        if (debug) cerr << "Loading GBZ file for path lookup: " << gbz_file << "..." << endl;
        
        GBZ gbz;
        try {
            sdsl::simple_sds::load_from(gbz, gbz_file);
        } catch (const exception& e) {
            cerr << "Error loading GBZ file: " << e.what() << endl;
            return 1;
        }
        
        if (debug) cerr << "GBZ loaded. Resolving path..." << endl;
        
        if (use_path_name) {
            // Direct path name lookup (for generic paths like "x", "y", "chr1")
            seq_id = find_sequence_id_by_path_name(gbz.index, path_name, debug);
            
            if (seq_id == numeric_limits<size_t>::max()) {
                cerr << "Error: Could not find path '" << path_name << "'" << endl;
                return 1;
            }
            
            if (debug) {
                cerr << "Resolved path name '" << path_name 
                     << "' to sequence ID: " << seq_id << endl;
            }
        } else {
            // Structured metadata lookup (for haplotype paths)
            seq_id = find_sequence_id_from_metadata(gbz.index, sample_name, contig_name, haplotype, debug);
            
            if (seq_id == numeric_limits<size_t>::max()) {
                cerr << "Error: Could not resolve path from metadata" << endl;
                return 1;
            }
            
            if (debug) {
                cerr << "Resolved path (sample='" << sample_name 
                     << "', contig='" << contig_name 
                     << "', haplotype=" << haplotype 
                     << ") to sequence ID: " << seq_id << endl;
            }
        }
    }
    
    // Timing variables
    auto start_total = high_resolution_clock::now();
    auto start_load_rindex = high_resolution_clock::now();
    
    // Load r-index
    FastLocate r_index;
    {
        if (debug) cerr << "Loading r-index: " << r_index_file << "..." << endl;
        ifstream rin(r_index_file, ios::binary);
        if (!rin) { 
            cerr << "Cannot open r-index: " << r_index_file << endl; 
            return 1; 
        }
        r_index.load_encoded(rin);
    }
    auto end_load_rindex = high_resolution_clock::now();
    auto duration_load_rindex = duration_cast<milliseconds>(end_load_rindex - start_load_rindex);
    if (debug) cerr << "R-index loaded successfully. BWT size: " << r_index.bwt_size() << endl;
    
    auto start_load_tags = high_resolution_clock::now();
    // Load sampled tag array
    SampledTagArray sampled;
    {
        if (debug) cerr << "Loading sampled tag array: " << sampled_tags_file << "..." << endl;
        ifstream sin(sampled_tags_file, ios::binary);
        if (!sin) { 
            cerr << "Cannot open sampled.tags: " << sampled_tags_file << endl; 
            return 1; 
        }
        sampled.load(sin);
    }
    auto end_load_tags = high_resolution_clock::now();
    auto duration_load_tags = duration_cast<milliseconds>(end_load_tags - start_load_tags);
    if (debug) cerr << "Sampled tag array loaded successfully. Total runs: " << sampled.total_runs() << endl;
    
    // Convert path interval to sequence interval
    size_t seq_start = interval.first;
    size_t seq_end = interval.second;
    
    // Step 1: Convert sequence interval to packed text positions in last space
    // last marks packed text positions (not BWT positions)
    if (debug) cerr << "Converting sequence interval to packed text positions..." << endl;
    size_t text_pos_i = r_index.pack(seq_id, seq_start);  // Start of interval in packed text space
    size_t text_pos_j = r_index.pack(seq_id, seq_end);    // End of interval in packed text space
    
    if (debug) cerr << "Sequence interval [" << seq_start << ", " << seq_end << "] maps to text positions [" 
         << text_pos_i << ", " << text_pos_j << "] in last space" << endl;
    
    // Step 2: Find successor position on last at or after j
    // Use last_successor(j) to find the first marked position x at or after j
    r_index.ensure_last_rank();
    r_index.ensure_last_select();
    
    auto successor_result = r_index.last_successor(text_pos_j);
    size_t text_pos_x = successor_result.first;   // The marked text position at or after j
    size_t rank_x = successor_result.second;      // 0-based index into last_to_run

    size_t tot_runs = r_index.tot_runs();
    size_t run_id = 0;
    if (rank_x < r_index.last_to_run.size()) {
        run_id = r_index.last_to_run[rank_x];
        // last_to_run is int_vector; reading into size_t can sign-extend if stored value had high bit set
        if (r_index.last_to_run.width() < 64 && (run_id >> r_index.last_to_run.width()) != 0) {
            run_id &= (1ULL << r_index.last_to_run.width()) - 1;
        }
        if (run_id >= tot_runs) {
            cerr << "Error: last_to_run[" << rank_x << "]=" << r_index.last_to_run[rank_x]
                 << " is invalid (tot_runs=" << tot_runs << "). Index may be corrupt or built with a bug." << endl;
            return 1;
        }
    }

    // Find BWT position ISA[x] at the end of this run
    size_t bwt_pos_x = r_index.bwt_end_position_of_run(run_id);

    if (debug) cerr << "Found marked text position x=" << text_pos_x << " (rank=" << rank_x
         << ", run_id=" << run_id << ", BWT_pos=" << bwt_pos_x << ")" << endl;
    
    // Check if text_pos_x is after the sequence end, and if so, start from BWT[seq_id]
    size_t bwt_seq_end = seq_id;  // BWT position of $ at end of sequence seq_id
    size_t text_pos_seq_end;
    {
        size_t rindex_run_id = 0;
        size_t run_start_pos = 0;
        r_index.run_id_and_offset_at(bwt_seq_end, rindex_run_id, run_start_pos);
        text_pos_seq_end = r_index.getSample(rindex_run_id);
        // Navigate from run start to bwt_seq_end
        for (size_t p = run_start_pos; p < bwt_seq_end; ++p) {
            text_pos_seq_end = r_index.locateNext(text_pos_seq_end);
        }
    }
    
    if (debug) cerr << "Sequence end check: text_pos_x=" << text_pos_x 
         << ", text_pos_seq_end=" << text_pos_seq_end << " (BWT[" << bwt_seq_end << "])" << endl;
    
    // Check if text_pos_x (BWT rank from last array) is greater than sequence end
    if (text_pos_x > text_pos_seq_end) {
        cerr << "Error: BWT rank calculated from last array (text_pos_x=" << text_pos_x 
             << ") is greater than sequence end (text_pos_seq_end=" << text_pos_seq_end << ")" << endl;
        return 1;
    }
    
    // Step 3: Iterate backwards: (x-1, LF(ISA[x])) until x reaches i
    // Collect tags from TAG[ISA[x]] for each x in [i, j]
    auto start_find_tags = high_resolution_clock::now();
    if (debug) cerr << "Iterating backwards through text positions [" << text_pos_i << ", " << text_pos_x << "]..." << endl;
    sampled.ensure_run_rank();
    sampled.ensure_run_select();
    
    unordered_map<uint64_t, TagResult> results;
    size_t current_text_pos = text_pos_x;
    size_t current_bwt_pos = bwt_pos_x;
    size_t positions_checked = 0;
    size_t positions_with_tags = 0;
    
    // Iterate backwards from x down to i
    // Only process positions in [i, j] (x might be > j since it's the successor)
    while (current_text_pos >= text_pos_i) {
        // Only process positions in our interval [i, j]
        if (current_text_pos <= text_pos_j) {
            positions_checked++;
            // Get tag for current BWT position
            size_t sampled_run_id = sampled.run_id_at(current_bwt_pos);
            uint64_t tag_val = sampled.run_value(sampled_run_id);


            if (debug) cerr << "  Text pos " << current_text_pos << " -> BWT pos " << current_bwt_pos 
                     << " -> tag=" << tag_val << endl;

            // Do nothing if tag_val == 0
            // If we found a tag, unpack current_text_pos to get (seq_id, offset)
            if (tag_val != 0) {
                positions_with_tags++;
                // Unpack current_text_pos directly (it's already a packed position)
                size_t seq_id_found = r_index.seqId(current_text_pos);
                size_t offset_from_start = r_index.seqOffset(current_text_pos);
                
                // Calculate offset relative to query start (only for query sequence)
                if (seq_id_found == seq_id && offset_from_start >= seq_start && offset_from_start <= seq_end) {
                    size_t offset_in_query = offset_from_start - seq_start;
                    results[tag_val].query_offsets.push_back(offset_in_query);
                }
                
                if (debug) cerr << "  Text pos " << current_text_pos << " (seq_id=" << seq_id_found 
                     << ", offset=" << offset_from_start << ") -> BWT pos " << current_bwt_pos 
                     << " -> tag=" << tag_val << endl;
            }
        }
        
        // Move to previous text position: x--
        if (current_text_pos > text_pos_i) {
            current_text_pos--;
            
            // Apply LF-mapping: ISA[x-1] = LF(ISA[x])
            current_bwt_pos = r_index.LF(current_bwt_pos);
        } else {
            break; // Reached i
        }
    }
    auto end_find_tags = high_resolution_clock::now();
    auto duration_find_tags = duration_cast<milliseconds>(end_find_tags - start_find_tags);
    if (debug) cerr << "Found " << results.size() << " unique tags from backward iteration" << endl;
    
    if (results.empty()) {
        auto end_total = high_resolution_clock::now();
        auto duration_total = duration_cast<milliseconds>(end_total - start_total);
        cout << "seq_id=" << seq_id << "\tinterval=" << seq_start << ".." << seq_end << "\tno_tags" << endl;
        cerr << "\n=== Query Statistics ===" << endl;
        cerr << "  Sequence ID: " << seq_id << endl;
        if (use_metadata) {
            cerr << "  Sample: " << sample_name << endl;
            cerr << "  Contig: " << contig_name << endl;
            cerr << "  Haplotype: " << haplotype << endl;
        }
        cerr << "  Interval: " << seq_start << ".." << seq_end << " (length: " << (seq_end - seq_start + 1) << ")" << endl;
        cerr << "  Positions checked: " << positions_checked << endl;
        cerr << "  Positions with tags: " << positions_with_tags << endl;
        cerr << "  Unique tags found: 0" << endl;
        cerr << "\n=== Timing Statistics ===" << endl;
        cerr << "  Load r-index: " << duration_load_rindex.count() << " ms" << endl;
        cerr << "  Load sampled tags: " << duration_load_tags.count() << " ms" << endl;
        cerr << "  Find tags in interval: " << duration_find_tags.count() << " ms" << endl;
        cerr << "  Total time: " << duration_total.count() << " ms" << endl;
        return 0;
    }
    
    // Step 3: For each tag, find all occurrences using wt_gmr select queries
    auto start_find_sequences = high_resolution_clock::now();
    if (debug) cerr << "Enumerating all occurrences for " << results.size() << " tags..." << endl;
    const auto& wm = sampled.values();

    
    size_t total_tag_runs = 0;
    size_t total_bwt_positions_processed = 0;
    size_t total_unique_sequences = 0;
    size_t total_occurrences = 0;
    
    for (auto& kv : results) {
        uint64_t tag_code = kv.first;
        TagResult& res = kv.second;
        
        // Deduplicate query offsets
        sort(res.query_offsets.begin(), res.query_offsets.end());
        res.query_offsets.erase(unique(res.query_offsets.begin(), res.query_offsets.end()), res.query_offsets.end());
        
        // Get total occurrences of this tag using wt_gmr rank
        size_t total_occ = wm.rank(wm.size(), tag_code);
        if (total_occ == 0) {
            if (debug) cerr << "Warning: Tag " << tag_code << " has 0 occurrences in wavelet matrix" << endl;
            continue;
        }
        
        total_tag_runs += total_occ;
        if (debug) cerr << "  Tag " << tag_code << ": Found " << total_occ << " tag runs containing this tag" << endl;
        
        // Step 4: For each run containing this tag, locate all (seq_id, offset) pairs
        // Use select queries to find all runs containing this tag
        // We want ALL sequences that have this tag, not just those in the query interval
        bool first_run_is_gap = sampled.is_first_run_gap();
        for (size_t j = 1; j <= total_occ; ++j) {
            size_t tag_run_id = wm.select(j, tag_code);  // This is a TAG run_id (0-indexed into wt_gmr, non-gap only)
            // Convert TAG run_id to BWT run_id (which includes gaps)
            // If first_run_is_gap=false: BWT run_id = 2 * TAG run_id (even = normal tags)
            // If first_run_is_gap=true: BWT run_id = 2 * TAG run_id + 1 (odd = normal tags)
            size_t bwt_run_id = 2 * tag_run_id + 1 - static_cast<size_t>(first_run_is_gap);
            auto run_span = sampled.run_span(bwt_run_id);
            
            if (debug) cerr << "    Tag Run " << j << "/" << total_occ << ": tag_run_id=" << tag_run_id 
                 << ", bwt_run_id=" << bwt_run_id
                 << ", BWT interval=[" << run_span.first << ", " << run_span.second 
                 << "] (size=" << (run_span.second - run_span.first + 1) << ")" << endl;
            
            // Process the entire run to get all sequences with this tag
            size_t locate_start = run_span.first;
            size_t locate_end = run_span.second;
            
            // IMPORTANT: Handle BWT size mismatch between SampledTagArray and r-index
            // SampledTagArray includes endmarkers, but r-index.bwt_size() doesn't
            // Skip entire run if it's completely out of bounds
            if (locate_start >= r_index.bwt_size()) {
                if (debug) cerr << "      Skipping entire run: locate_start=" << locate_start 
                     << " >= r_index.bwt_size()=" << r_index.bwt_size() << endl;
                continue;
            }
            
            // Clamp locate_end to valid BWT range
            if (locate_end >= r_index.bwt_size()) {
                if (debug) cerr << "      Clamping locate_end from " << locate_end 
                     << " to " << (r_index.bwt_size() - 1) << endl;
                locate_end = r_index.bwt_size() - 1;
            }
            
            size_t positions_processed = 0;
            size_t unique_sequences_found = 0;

            size_t packed_pos;
            {
                size_t rindex_run_id = 0;
                size_t run_start_pos = 0;
                r_index.run_id_and_offset_at(locate_start, rindex_run_id, run_start_pos);
                packed_pos = r_index.getSample(rindex_run_id);
                
                // Navigate from run start to locate_start
                for (size_t p = run_start_pos; p < locate_start; ++p) {
                    packed_pos = r_index.locateNext(packed_pos);
                }
            }
            
            // Process all positions sequentially
            for (size_t pos = locate_start; pos <= locate_end; ++pos) {
                if (pos > locate_start) {
                    packed_pos = r_index.locateNext(packed_pos);
                }
                
                // Unpack to get (seq_id, offset)
                auto pr = r_index.unpack(packed_pos);
                size_t seq_id_found = pr.first;
                size_t offset_from_start = pr.second; // offset is already distance from START
                
                // Each BWT position maps to a unique (seq_id, offset) pair
                res.sequence_offsets[seq_id_found].push_back(offset_from_start);
                unique_sequences_found++;
                
                // Print first few conversions for debugging
                if (debug && positions_processed < 5) {
                    cerr << "        BWT[" << pos << "] -> (seq_id=" << seq_id_found 
                         << ", offset_from_start=" << offset_from_start << ")" << endl;
                }
                positions_processed++;
            }
            
            total_bwt_positions_processed += positions_processed;
            total_unique_sequences += unique_sequences_found;
            total_occurrences += positions_processed;
            
            if (debug) cerr << "      Processed " << positions_processed << " BWT positions, found " 
                 << unique_sequences_found << " unique (seq_id, offset) pairs" << endl;
        }
    }
    auto end_find_sequences = high_resolution_clock::now();
    auto duration_find_sequences = duration_cast<milliseconds>(end_find_sequences - start_find_sequences);
    
    // Step 5: Output results
    // For each TAG, show the different sequence IDs that pass through those TAGS
    if (debug) cerr << "Writing results..." << endl;
    
    cout << "Query Sequence:" << endl;
    cout << "  Sequence ID: " << seq_id << endl;
    if (use_metadata) {
        cout << "  Sample: " << sample_name << endl;
        cout << "  Contig: " << contig_name << endl;
        cout << "  Haplotype: " << haplotype << endl;
    }
    cout << "  Interval: " << seq_start << ".." << seq_end << endl;
    cout << "  Number of tags found: " << results.size() << endl;
    cout << endl;
    
    // Sort tags by tag_code before printing
    vector<pair<uint64_t, TagResult>> sorted_results(results.begin(), results.end());
    sort(sorted_results.begin(), sorted_results.end(), 
         [](const pair<uint64_t, TagResult>& a, const pair<uint64_t, TagResult>& b) {
             return a.first < b.first;
         });
    
    // Print results for each tag
    for (const auto& kv : sorted_results) {
        uint64_t tag_code = kv.first;
        const TagResult& res = kv.second;
        
        // Decode tag code to get node_id and is_rev
        uint64_t decoded = tag_code - 1;
        int64_t node_id = (decoded >> 1) + 1;
        bool is_rev = (decoded & 1) != 0;
        
        cout << "Tag Code: " << tag_code << endl;
        cout << "  Tag Details:" << endl;
        cout << "    Node ID: " << node_id << endl;
        cout << "    Is Reverse: " << (is_rev ? "true" : "false") << endl;
        cout << "    Query offsets (relative to query start): ";
        for (size_t i = 0; i < res.query_offsets.size(); ++i) {
            if (i > 0) cout << ", ";
            cout << res.query_offsets[i];
        }
        cout << endl;
        
        // Show sequence IDs that pass through this tag
        cout << "  Sequence IDs passing through this tag: " << res.sequence_offsets.size() << endl;
        for (const auto& seq_kv : res.sequence_offsets) {
            size_t seq_id_found = seq_kv.first;
            const vector<size_t>& offsets = seq_kv.second;
            
            cout << "    Sequence ID: " << seq_id_found << " (occurrences: " << offsets.size() << ")" << endl;
            cout << "      Offsets: ";
            for (size_t i = 0; i < offsets.size() && i < 10; ++i) {
                if (i > 0) cout << ", ";
                cout << offsets[i];
            }
            if (offsets.size() > 10) {
                cout << ", ... (and " << (offsets.size() - 10) << " more)";
            }
            cout << endl;
        }
        cout << endl;
    }
    
    auto end_total = high_resolution_clock::now();
    auto duration_total = duration_cast<milliseconds>(end_total - start_total);
    
    // Print statistics
    cerr << "\n=== Query Statistics ===" << endl;
    cerr << "  Sequence ID: " << seq_id << endl;
    if (use_metadata) {
        cerr << "  Sample: " << sample_name << endl;
        cerr << "  Contig: " << contig_name << endl;
        cerr << "  Haplotype: " << haplotype << endl;
    }
    cerr << "  Interval: " << seq_start << ".." << seq_end << " (length: " << (seq_end - seq_start + 1) << ")" << endl;
    cerr << "  Positions checked: " << positions_checked << endl;
    cerr << "  Positions with tags: " << positions_with_tags << endl;
    cerr << "  Unique tags found: " << results.size() << endl;
    cerr << "  Total tag runs: " << total_tag_runs << endl;
    cerr << "  Total BWT positions processed: " << total_bwt_positions_processed << endl;
    cerr << "  Total unique sequences: " << total_unique_sequences << endl;
    cerr << "  Total occurrences: " << total_occurrences << endl;
    if (results.size() > 0) {
        cerr << "  Average occurrences per tag: " << fixed << setprecision(2) 
             << (double)total_occurrences / results.size() << endl;
        cerr << "  Average sequences per tag: " << fixed << setprecision(2) 
             << (double)total_unique_sequences / results.size() << endl;
    }
    
    cerr << "\n=== Timing Statistics ===" << endl;
    cerr << "  Load r-index: " << duration_load_rindex.count() << " ms" << endl;
    cerr << "  Load sampled tags: " << duration_load_tags.count() << " ms" << endl;
    cerr << "  Find tags in interval: " << duration_find_tags.count() << " ms" << endl;
    cerr << "  Find all sequences for tags: " << duration_find_sequences.count() << " ms" << endl;
    cerr << "  Total time: " << duration_total.count() << " ms" << endl;
    
    if (debug) cerr << "Query completed successfully" << endl;
    
    return 0;
}
