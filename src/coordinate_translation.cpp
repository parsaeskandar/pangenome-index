#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include "pangenome_index/translation_tables.hpp"
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/simple_sds.hpp>
#include <gbwt/gbwt.h>
#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>
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
#include <random>
#include <atomic>

using namespace std;
using namespace std::chrono;

// Use explicit namespace for panindexer types to avoid conflict with gbwt::FastLocate
using panindexer::FastLocate;
using panindexer::SampledTagArray;
using panindexer::TranslationTable1;
using panindexer::TranslationTable2;
using panindexer::PathInterval;
using panindexer::SubpathInfo;
using panindexer::TargetInterval;
using panindexer::IntervalMapping;

static std::atomic<bool> debug{false};

// Structure to store tag information with source offsets
struct TagInfo {
    uint64_t tag_code;                    // Encoded (node_id, is_rev)
    vector<size_t> source_offsets;        // Offsets on source haplotype where this tag appears
    vector<size_t> source_bwt_positions; // BWT positions where source visits this tag
    vector<size_t> source_packed_positions; // Packed text positions for source visits
};

// Structure to store a visit to a node by a haplotype
struct NodeVisit {
    size_t seq_id;
    size_t offset;              // Offset on the haplotype
    size_t bwt_pos;             // BWT position
    size_t packed_pos;          // Packed text position (SA value)
    uint64_t tag_code;          // Tag code for this node
};

// Diagnostics for find_sequences_for_tag's LF cost, accumulated across calls
// in a thread-local so callers (e.g. the anchor builder) can measure how much
// of a query is spent enumerating a node's pangenome-wide visits. Reset to {}
// before a measured region and read afterwards. Defined identically (file
// scope) in surject_anchor_builder.cpp / pangenome_server.cpp per this file's
// existing extern-struct convention (NodeVisit, TagInfo, …).
struct FindSeqStats {
    size_t calls = 0;             // find_sequences_for_tag invocations
    size_t runs = 0;              // total tag runs ("vectors") iterated
    size_t lf_steps = 0;          // total locateNext (LF) calls — navigation + run walk
    size_t visits = 0;            // total NodeVisit entries produced (node's pangenome usage)
    size_t last_run_nav_steps = 0;// LF steps to navigate from the sample to the LAST run's start
    size_t last_run_length = 0;   // number of positions in the LAST run iterated
};
thread_local FindSeqStats g_find_seq_stats;

// Structure for coordinate translation result
struct TranslationResult {
    size_t source_offset;       // Offset on source haplotype
    size_t target_offset;       // Offset on target haplotype
    size_t target_seq_id;       // Target sequence ID (for fragments)
    uint64_t tag_code;          // Node ID at this position
};

// Structure to store fragment information
struct FragmentInfo {
    size_t seq_id;              // Sequence ID of this fragment
    size_t start_offset;        // Start offset on this fragment
    size_t end_offset;          // End offset on this fragment
    vector<size_t> anchor_offsets; // Anchor offsets where this fragment aligns
};

// Optional timing breakdown for find_tags_in_interval (used with --benchmark)
struct FindTagsInIntervalTiming {
    size_t num_lf_before_phase1 = 0;     // number of LF steps in initial skip loop (text_pos_x - text_pos_j)
    double init_skip_lf_ms = 0;          // time spent in initial skip loop (LF until current_text_pos <= text_pos_j)
    size_t phase1_lf_count = 0;        // number of r_index.LF() calls in phase1
    double phase1_lf_ms = 0;            // time spent in r_index.LF() in phase1
    double phase1_tag_lookup_ms = 0;     // sampled.run_id_at() and sampled.run_value()
    double phase1_position_lookup_ms = 0; // r_index.seqId() and r_index.seqOffset()
    double phase1_tag_map_ms = 0;        // tag_map find/insert/push_back operations
    double phase1_fast_path_ms = 0;      // fast path: decode_tag, decompressSA, loops
    bool used_fast_path = false;         // true if fast path was enabled (gbwt + graph + target + bidirectional)
    bool found_last_common = false;     // true if we stopped at last common node and broke out of phase1 loop
    size_t last_common_source_base = 0;  // source sequence offset (base position) where last common node was found; 0 if not found
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
            if (matching_paths.size() > 1) {
                cerr << "Warning: Multiple paths (" << matching_paths.size() 
                     << ") found for contig '" << path_name_str << "'. Using first match (path_id=" 
                     << matching_paths[0] << ")" << endl;
            }
            
            if (debug_output) {
                cerr << "Resolved path name '" << path_name_str << "' to sequence ID: " << matching_paths[0] << endl;
            }
            return matching_paths[0];
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
    
    if (matching_paths.size() > 1) {
        cerr << "Warning: Multiple paths (" << matching_paths.size() 
             << ") found matching the criteria. Using first match (path_id=" 
             << matching_paths[0] << ")" << endl;
    }
    
    // The path_id in GBWT metadata corresponds to the sequence ID
    size_t seq_id = matching_paths[0];
    
    if (debug_output) {
        cerr << "Resolved to sequence ID: " << seq_id << endl;
    }
    
    return seq_id;
}

// Structure to hold path specification options for source or target
struct PathSpec {
    size_t seq_id = numeric_limits<size_t>::max();
    string path_name;
    string sample_name;
    string contig_name;
    size_t haplotype = 0;
    bool has_path_name = false;
    bool has_sample = false;
    bool has_contig = false;
    bool has_haplotype = false;
    bool reverse_strand = false;  // If true, use reverse complement strand (seq 2i+1 instead of 2i)
    
    bool use_path_name() const { return has_path_name; }
    bool use_metadata() const { return has_sample || has_contig; }
    bool use_seq_id() const { return seq_id != numeric_limits<size_t>::max(); }
    
    // Resolve the sequence ID using the provided GBWT index
    // For RLBWT, sequences are stored in pairs: seq 2i = forward, seq 2i+1 = reverse complement
    // When resolving from path name or metadata, we get the GBWT path_id and convert to RLBWT seq_id
    // Returns true on success, false on error
    bool resolve(const gbwt::GBWT& gbwt_index, const string& label, bool debug_output) {
        if (use_seq_id()) {
            // Already have seq_id (direct RLBWT seq_id), nothing to do
            // Note: when using --*-id directly, user provides the RLBWT seq_id
            return true;
        }
        
        size_t gbwt_path_id = numeric_limits<size_t>::max();
        
        if (use_path_name()) {
            gbwt_path_id = find_sequence_id_by_path_name(gbwt_index, path_name, debug_output);
            if (gbwt_path_id == numeric_limits<size_t>::max()) {
                cerr << "Error: Could not find " << label << " path '" << path_name << "'" << endl;
                return false;
            }
        } else if (use_metadata()) {
            gbwt_path_id = find_sequence_id_from_metadata(gbwt_index, sample_name, contig_name, haplotype, debug_output);
            if (gbwt_path_id == numeric_limits<size_t>::max()) {
                cerr << "Error: Could not resolve " << label << " path from metadata" << endl;
                return false;
            }
        } else {
            return false;
        }
        
        // Convert GBWT path_id to RLBWT seq_id
        // RLBWT stores both strands: seq 2i = forward, seq 2i+1 = reverse complement
        seq_id = 2 * gbwt_path_id + (reverse_strand ? 1 : 0);
        
        if (debug_output) {
            cerr << "Resolved " << label << ": GBWT path_id=" << gbwt_path_id 
                 << " -> RLBWT seq_id=" << seq_id 
                 << " (strand: " << (reverse_strand ? "reverse" : "forward") << ")" << endl;
            if (use_path_name()) {
                cerr << "  (from path name '" << path_name << "')" << endl;
            } else if (use_metadata()) {
                cerr << "  (from sample='" << sample_name 
                     << "', contig='" << contig_name 
                     << "', haplotype=" << haplotype << ")" << endl;
            }
        }
        
        return true;
    }
    
    // Validate that only one method is used
    bool validate(const string& label) const {
        int methods = (use_seq_id() ? 1 : 0) + (use_path_name() ? 1 : 0) + (use_metadata() ? 1 : 0);
        if (methods > 1) {
            cerr << "Error: Cannot mix --" << label << "-id, --" << label << "-path-name, and --" 
                 << label << "-sample/--" << label << "-contig options" << endl;
            return false;
        }
        if (use_metadata()) {
            if (!has_sample) {
                cerr << "Error: --" << label << "-sample is required when using --" << label << "-contig" << endl;
                return false;
            }
            if (!has_contig) {
                cerr << "Error: --" << label << "-contig is required when using --" << label << "-sample" << endl;
                return false;
            }
        }
        return true;
    }
    
    // Check if any specification method was provided
    bool is_specified() const {
        return use_seq_id() || use_path_name() || use_metadata();
    }
};

// Single-path --benchmark (without tables): fixed RLBWT pair — forward strand of GBWT path 0 vs path 1.
constexpr size_t kBenchmarkSourceSeqId = 0;
constexpr size_t kBenchmarkTargetSeqId = 2;

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <rlbwt_rindex.ri> <gbwt_rindex.ri> <sampled.tags> [options]" << endl;
    cerr << endl;
    cerr << "Required files:" << endl;
    cerr << "  rlbwt_rindex.ri         RLBWT R-index file" << endl;
    cerr << "  gbwt_rindex.ri          GBWT FastLocate R-index file (.ri)" << endl;
    cerr << "  sampled.tags            Sampled tag array file" << endl;
    cerr << endl;
    cerr << "Source specification (choose one):" << endl;
    cerr << "  --source-id ID                            Direct RLBWT sequence ID (already accounts for strand)" << endl;
    cerr << "  --source-path-name NAME [--source-reverse]" << endl;
    cerr << "                                            Specify by path name (e.g., 'x', 'chr1')" << endl;
    cerr << "  --source-sample NAME --source-contig NAME [--source-haplotype N] [--source-reverse]" << endl;
    cerr << "                                            Specify by metadata" << endl;
    cerr << endl;
    cerr << "Target specification (choose one):" << endl;
    cerr << "  --target-id ID                            Direct RLBWT sequence ID (already accounts for strand)" << endl;
    cerr << "  --target-path-name NAME [--target-reverse]" << endl;
    cerr << "                                            Specify by path name (e.g., 'x', 'chr1')" << endl;
    cerr << "  --target-sample NAME --target-contig NAME [--target-haplotype N] [--target-reverse]" << endl;
    cerr << "                                            Specify by metadata" << endl;
    cerr << endl;
    cerr << "Note on RLBWT sequence IDs:" << endl;
    cerr << "  RLBWT stores both strands: seq 2i = forward, seq 2i+1 = reverse complement" << endl;
    cerr << "  When using --*-path-name or --*-sample/--*-contig, the GBWT path_id is converted:" << endl;
    cerr << "    RLBWT seq_id = 2 * GBWT_path_id (forward) or 2 * GBWT_path_id + 1 (reverse)" << endl;
    cerr << "  Use --source-reverse or --target-reverse to select reverse complement strand" << endl;
    cerr << endl;
    cerr << "Required options:" << endl;
    cerr << "  --interval START..END   Interval on source haplotype (base offsets; not required if --benchmark)" << endl;
    cerr << "  --gbz FILE              GBZ file (.gbz); contains GBWT index + GBWTGraph (required)" << endl;
    cerr << endl;
    cerr << "Table-driven mode (use together; can be combined with --benchmark):" << endl;
    cerr << "  --table1 FILE           Translation Table 1 (path name + global interval -> path_id, local interval)" << endl;
    cerr << "  --table2 FILE           Translation Table 2 (src_path_id, tgt haplotype -> tgt_path_ids)" << endl;
    cerr << "  --source-haplotype-name NAME   Source haplotype name (e.g. GRCh38#0); interval is in its global coords" << endl;
    cerr << "  --target-haplotype-name NAME   Target haplotype name (e.g. HG002#1)" << endl;
    cerr << endl;
    cerr << "Other options:" << endl;
    cerr << "  --benchmark             Run benchmark: translate seq 0->2 at gene/megabase scales (excludes load time)" << endl;
    cerr << "  --benchmark-intervals L1,L2,...  Comma-separated interval lengths (default: 5000,10000,1000000)" << endl;
    cerr << "  --debug                 Enable debug output" << endl;
    cerr << "  --no-debug              Disable debug output" << endl;
    cerr << endl;
    cerr << "Examples:" << endl;
    cerr << "  # Using direct RLBWT sequence IDs" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --source-id 0 --target-id 2 --interval 1000..2000 --gbz graph.gbz" << endl;
    cerr << endl;
    cerr << "  # Using path names (forward strand)" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --source-path-name x --target-path-name y --interval 1000..2000 --gbz graph.gbz" << endl;
    cerr << endl;
    cerr << "  # Using path names with reverse complement strand for target" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --source-path-name x --target-path-name y --target-reverse --interval 1000..2000 --gbz graph.gbz" << endl;
    cerr << endl;
    cerr << "  # Using metadata" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --source-sample HG002 --source-contig chr1 --source-haplotype 1 \\" << endl;
    cerr << "       --target-sample GRCh38 --target-contig chr1 --interval 1000..2000 --gbz graph.gbz" << endl;
    cerr << endl;
    cerr << "  # Table-driven mode (haplotype names + translation tables)" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --table1 t1.bin --table2 t2.bin \\" << endl;
    cerr << "       --source-haplotype-name GRCh38#0 --target-haplotype-name HG002#1 --interval 1000..2000 --gbz graph.gbz" << endl;
    cerr << endl;
    cerr << "  # Benchmark with translation tables (random global intervals, no --interval needed)" << endl;
    cerr << "  " << prog << " idx.ri gbwt.ri tags.stags --benchmark --table1 t1.bin --table2 t2.bin \\" << endl;
    cerr << "       --source-haplotype-name GRCh38#0 --target-haplotype-name HG002#1 --gbz graph.gbz" << endl;
}

/// Build path_id -> (path_name, subpath_start) from Table 1 for converting path-local offsets to haplotype (global) coordinates.
static unordered_map<size_t, pair<string, size_t>> build_path_id_to_global(const TranslationTable1& t1) {
    unordered_map<size_t, pair<string, size_t>> out;
    vector<string> names = t1.names();
    for (const string& name : names) {
        vector<SubpathInfo> subpaths = t1.subpaths(name);
        for (const SubpathInfo& sp : subpaths) {
            out[sp.path_id] = {name, sp.subpath_start};
        }
    }
    return out;
}

/// Return path names in Table 1 whose base name starts with haplotype_prefix (e.g. "GRCh38#0" -> "GRCh38#0#chr1", ...).
static vector<string> path_names_for_haplotype(const TranslationTable1& t1, const string& haplotype_prefix) {
    vector<string> names = t1.names();
    vector<string> result;
    string prefix = haplotype_prefix;
    if (prefix.empty() || prefix.back() != '#') {
        prefix += '#';
    }
    for (const string& name : names) {
        if (name.size() >= prefix.size() && name.compare(0, prefix.size(), prefix) == 0) {
            result.push_back(name);
        }
    }
    return result;
}

// Forward declaration for helper used by find_tags_in_interval
pair<int64_t, bool> decode_tag(uint64_t tag_code);

// Recover packed text position (SA value) at a given BWT position using r-index.
// Used to verify that (text_pos, bwt_pos) pairs maintained during LF iteration are correct.
static size_t recover_text_pos_from_bwt(FastLocate& r_index, size_t bwt_pos) {
    size_t run_id = 0, bwt_run_start = 0;
    r_index.run_id_and_offset_at(bwt_pos, run_id, bwt_run_start);
    size_t packed = r_index.getSample(run_id);
    for (size_t k = bwt_run_start; k < bwt_pos; ++k) {
        packed = r_index.locateNext(packed);
    }
    return packed;
}

// Find all tags (nodes) visited by the source haplotype in the given interval.
// Optional fast path: when gbwt_index, gbwt_fast_locate, graph, and target_seq_id are provided
// and the GBWT is bidirectional, stops at the last common node with the target haplotype and
// completes the interval using GBWT inverseLF (backward walk) instead of r-index LF.
vector<TagInfo> find_tags_in_interval(FastLocate& r_index, SampledTagArray& sampled,
                                       size_t source_seq_id, size_t seq_start, size_t seq_end,
                                       const gbwt::GBWT* gbwt_index = nullptr,
                                       const gbwt::FastLocate* gbwt_fast_locate = nullptr,
                                       const gbwtgraph::GBWTGraph* graph = nullptr,
                                       size_t target_seq_id = numeric_limits<size_t>::max(),
                                       FindTagsInIntervalTiming* out_timing = nullptr) {
    auto t_find_tags_start = high_resolution_clock::now();
    vector<TagInfo> tags;
    unordered_map<uint64_t, TagInfo> tag_map;

    // Convert sequence interval to packed text positions
    size_t text_pos_i = r_index.pack(source_seq_id, seq_start);
    size_t text_pos_j = r_index.pack(source_seq_id, seq_end);
    
    if (debug) {
        cerr << "Finding tags in interval [" << seq_start << ", " << seq_end 
             << "] (text positions [" << text_pos_i << ", " << text_pos_j << "])" << endl;
    }
    
    auto successor_result = r_index.last_successor(text_pos_j);
    size_t text_pos_x = successor_result.first;
    size_t rank_x = successor_result.second;

    if (debug) {
        cerr << "rank_x=" << rank_x << ", last_to_run.size()=" << r_index.last_to_run.size() << endl;
    }

    size_t tot_runs = r_index.tot_runs();
    size_t run_id = 0;
    if (rank_x < r_index.last_to_run.size()) {
        run_id = r_index.last_to_run[rank_x];
        if (r_index.last_to_run.width() < 64 && (run_id >> r_index.last_to_run.width()) != 0) {
            run_id &= (1ULL << r_index.last_to_run.width()) - 1;
        }
        if (run_id >= tot_runs) {
            cerr << "Error: last_to_run[" << rank_x << "]=" << r_index.last_to_run[rank_x]
                 << " is invalid (tot_runs=" << tot_runs << "). Index may be corrupt or built with a bug." << endl;
            return tags;
        }
    }
    size_t bwt_pos_x = r_index.bwt_end_position_of_run(run_id);

    if (debug) {
        cerr << "Found marked text position x=" << text_pos_x << " (rank=" << rank_x
             << ", run_id=" << run_id << ", BWT_pos=" << bwt_pos_x << ")" << endl;
    }
    
    // Check if text_pos_x is after the sequence end
    size_t bwt_seq_end = source_seq_id;
    size_t text_pos_seq_end;
    {
        size_t rindex_run_id = 0;
        size_t run_start_pos = 0;
        r_index.run_id_and_offset_at(bwt_seq_end, rindex_run_id, run_start_pos);
        text_pos_seq_end = r_index.getSample(rindex_run_id);
        for (size_t p = run_start_pos; p < bwt_seq_end; ++p) {
            text_pos_seq_end = r_index.locateNext(text_pos_seq_end);
        }
    }
    
    if (debug) {
        cerr << "Sequence end check: text_pos_x=" << text_pos_x 
             << ", text_pos_seq_end=" << text_pos_seq_end << " (BWT[" << bwt_seq_end << "])" << endl;
    }
    
    if (text_pos_x > text_pos_seq_end) {
        cerr << "Error: BWT rank calculated from last array (text_pos_x=" << text_pos_x 
             << ") is greater than sequence end (text_pos_seq_end=" << text_pos_seq_end << ")" << endl;
        return tags;
    }
    
    size_t current_text_pos = text_pos_x;
    size_t current_bwt_pos = bwt_pos_x;

    // Number of LF steps required before phase1 (skip from text_pos_x down to text_pos_j)
    if (out_timing) {
        out_timing->num_lf_before_phase1 = (text_pos_x > text_pos_j) ? (text_pos_x - text_pos_j) : 0;
    }

    // Skip positions after the interval end
    auto t_init_skip_start = high_resolution_clock::now();
    while (current_text_pos > text_pos_j) {
        current_text_pos--;
        current_bwt_pos = r_index.LF(current_bwt_pos);
    }
    if (out_timing) {
        out_timing->init_skip_lf_ms = duration_cast<microseconds>(high_resolution_clock::now() - t_init_skip_start).count() / 1000.0;
    }

    if (debug) {
        cerr << "    is_first_run_gap=" << sampled.is_first_run_gap() << endl;
    }

    unordered_set<size_t> verify_text_positions;
    if (debug && text_pos_j >= text_pos_i) {
        const size_t num_verify = 10;
        for (size_t v = 0; v < num_verify; ++v) {
            size_t t = text_pos_i + (text_pos_j - text_pos_i) * v / (num_verify > 1 ? num_verify - 1 : 1);
            if (t >= text_pos_i && t <= text_pos_j) verify_text_positions.insert(t);
        }
        cerr << "    BWT<->Text verification: will check " << verify_text_positions.size() 
             << " positions (recover text_pos from bwt_pos via r-index)" << endl;
    }
    
    // Determine if we can use the fast path: GBWT inverseLF from last common node backward
    const bool use_fast_path = (gbwt_index != nullptr && gbwt_fast_locate != nullptr && graph != nullptr
                                && target_seq_id != numeric_limits<size_t>::max()
                                && gbwt_index->bidirectional());
    if (out_timing) {
        out_timing->used_fast_path = use_fast_path;
    }

    if (debug) {
        cerr << "use_fast_path=" << use_fast_path << endl;
    }
    
    // Last common node state (used only in fast path)
    bool found_last_common = false;
    gbwt::node_type last_common_node = gbwt::ENDMARKER;
    size_t last_common_source_base = 0;
    gbwt::edge_type source_gbwt_edge = { gbwt::ENDMARKER, 0 };  // sentinel for "not set"
    
    // ----- Phase 1: Traverse with r_index.last + LF backward, collect TAGs, stop at last common if fast path -----
    while (current_text_pos >= text_pos_i) {
        if (current_text_pos <= text_pos_j) {
            auto t_tag_lookup_start = high_resolution_clock::now();
            size_t sampled_run_id = sampled.run_id_at(current_bwt_pos);
            uint64_t tag_val = sampled.run_value(sampled_run_id);
            if (out_timing) {
                out_timing->phase1_tag_lookup_ms += duration_cast<microseconds>(high_resolution_clock::now() - t_tag_lookup_start).count() / 1000.0;
            }

            if (debug) {
                cerr << "  Text pos " << current_text_pos << " -> BWT pos " << current_bwt_pos 
                     << " -> tag=" << tag_val << " -> run_id=" << sampled_run_id << endl;
            }
            
            if (debug && verify_text_positions.count(current_text_pos)) {
                size_t recovered = recover_text_pos_from_bwt(r_index, current_bwt_pos);
                bool ok = (recovered == current_text_pos);
                cerr << "    [VERIFY] text_pos=" << current_text_pos << " bwt_pos=" << current_bwt_pos
                     << " -> recovered_text_pos=" << recovered << (ok ? " OK" : " MISMATCH") << endl;
            }
            
            if (tag_val != 0) {
                auto t_pos_lookup_start = high_resolution_clock::now();
                size_t seq_id_found = r_index.seqId(current_text_pos);
                size_t offset_from_start = r_index.seqOffset(current_text_pos);
                if (out_timing) {
                    out_timing->phase1_position_lookup_ms += duration_cast<microseconds>(high_resolution_clock::now() - t_pos_lookup_start).count() / 1000.0;
                }
                
                if (seq_id_found == source_seq_id && 
                    offset_from_start >= seq_start && offset_from_start <= seq_end) {
                    
                    auto t_tag_map_start = high_resolution_clock::now();
                    if (tag_map.find(tag_val) == tag_map.end()) {
                        tag_map[tag_val] = TagInfo{tag_val, {}, {}, {}};
                    }
                    
                    tag_map[tag_val].source_offsets.push_back(offset_from_start);
                    tag_map[tag_val].source_bwt_positions.push_back(current_bwt_pos);
                    tag_map[tag_val].source_packed_positions.push_back(current_text_pos);
                    if (out_timing) {
                        out_timing->phase1_tag_map_ms += duration_cast<microseconds>(high_resolution_clock::now() - t_tag_map_start).count() / 1000.0;
                    }
                }
                
                // Fast path: check if this TAG is the last common node (target haplotype also passes)
                if (use_fast_path) {
                    auto t_fast_path_start = high_resolution_clock::now();
                    auto [node_id, is_rev] = decode_tag(tag_val);
                    gbwt::node_type node = gbwt::Node::encode(node_id, is_rev);
                    std::vector<gbwt::size_type> sa_values = gbwt_fast_locate->decompressSA(node);
                    bool target_passes = false;
                    for (size_t i = 0; i < sa_values.size(); ++i) {
                        if (gbwt_fast_locate->seqId(sa_values[i]) == target_seq_id) {
                            target_passes = true;
                            break;
                        }
                    }
                    if (debug) {
                        size_t check_source_base = r_index.seqOffset(current_text_pos);
                        cerr << "  [fast path] check node_id=" << node_id
                             << " is_rev=" << is_rev
                             << " tag_code=" << tag_val
                             << " source_base=" << check_source_base
                             << " sa_values.size=" << sa_values.size()
                             << " target_seq_id=" << target_seq_id
                             << " target_passes=" << (target_passes ? "YES" : "no");
                        if (!target_passes) {
                            cerr << " (seq_ids visiting: ";
                            size_t max_print = std::min<size_t>(sa_values.size(), 8);
                            for (size_t i = 0; i < max_print; ++i) {
                                if (i > 0) cerr << ",";
                                cerr << gbwt_fast_locate->seqId(sa_values[i]);
                            }
                            if (sa_values.size() > max_print) cerr << ",...";
                            cerr << ")";
                        }
                        cerr << endl;
                    }
                    if (target_passes) {
                        // This is the last common node (first common we hit when walking backward from end)
                        found_last_common = true;
                        last_common_node = node;
                        last_common_source_base = r_index.seqOffset(current_text_pos);
                        
                        // Find GBWT handle for source haplotype on this node (same as check_common_node for last common):
                        // decompressSA(node) gives all visits; for last common we want source visit with SMALLEST
                        // GBWT seq_offset (latest in path). The index i in that array is the offset within node's record.
                        gbwt::size_type source_offset_in_node = gbwt::invalid_offset();
                        gbwt::size_type min_seq_offset = gbwt::invalid_offset();  // smallest = latest in path
                        for (size_t i = 0; i < sa_values.size(); ++i) {
                            if (gbwt_fast_locate->seqId(sa_values[i]) == source_seq_id) {
                                gbwt::size_type so = gbwt_fast_locate->seqOffset(sa_values[i]);
                                if (min_seq_offset == gbwt::invalid_offset() || so < min_seq_offset) {
                                    min_seq_offset = so;
                                    source_offset_in_node = i;
                                }
                            }
                        }
                        if (source_offset_in_node != gbwt::invalid_offset()) {
                            source_gbwt_edge = gbwt::edge_type(node, source_offset_in_node);
                        }
                        // Always break out of the while when we found the last common node (target_passes),
                        // so we stop the expensive LF walk. Phase 2 (inverseLF) runs only when source_gbwt_edge is set.
                        if (out_timing) {
                            out_timing->found_last_common = true;
                            out_timing->last_common_source_base = last_common_source_base;
                        }
                        if (debug) {
                            cerr << "  [fast path] Last common node: tag_code=" << tag_val
                                 << " (node_id=" << node_id << ", is_rev=" << is_rev << ")"
                                 << ", last_common_source_base=" << last_common_source_base
                                 << ", source_gbwt_edge=(" << source_gbwt_edge.first << "," << source_gbwt_edge.second << ")" << endl;
                        }
                        if (out_timing) {
                            out_timing->phase1_fast_path_ms += duration_cast<microseconds>(high_resolution_clock::now() - t_fast_path_start).count() / 1000.0;
                        }
                        break;
                    }
                    if (out_timing) {
                        out_timing->phase1_fast_path_ms += duration_cast<microseconds>(high_resolution_clock::now() - t_fast_path_start).count() / 1000.0;
                    }
                }
            }
        }
        
        if (found_last_common) break;
        
        if (current_text_pos > text_pos_i) {
            current_text_pos--;
            auto t_lf_start = high_resolution_clock::now();
            current_bwt_pos = r_index.LF(current_bwt_pos);
            if (out_timing) {
                out_timing->phase1_lf_count++;
                out_timing->phase1_lf_ms += duration_cast<microseconds>(high_resolution_clock::now() - t_lf_start).count() / 1000.0;
            }
        } else {
            break;
        }
    }

    // Fast path: if no common node was found inside the interval (e.g. interval is too small
    // to cover any sampled tag run), extend the BWT walk backward past text_pos_i down to the
    // haplotype start, looking for the first source tag whose node is also visited by the target.
    if (use_fast_path && !found_last_common) {
        size_t text_pos_haplotype_start = r_index.pack(source_seq_id, 0);
        if (debug) {
            cerr << "  [extend] No common node in interval; extending backward past interval start "
                 << "(text_pos_i=" << text_pos_i << ", haplotype_start=" << text_pos_haplotype_start << ")" << endl;
        }
        while (current_text_pos > text_pos_haplotype_start && !found_last_common) {
            auto t_lf_start = high_resolution_clock::now();
            current_text_pos--;
            current_bwt_pos = r_index.LF(current_bwt_pos);
            if (out_timing) {
                out_timing->phase1_lf_count++;
                out_timing->phase1_lf_ms += duration_cast<microseconds>(high_resolution_clock::now() - t_lf_start).count() / 1000.0;
            }

            size_t sampled_run_id = sampled.run_id_at(current_bwt_pos);
            uint64_t tag_val = sampled.run_value(sampled_run_id);

            if (debug) {
                cerr << "  [extend] Text pos " << current_text_pos << " -> BWT pos " << current_bwt_pos
                     << " -> tag=" << tag_val << " -> run_id=" << sampled_run_id << endl;
            }

            if (tag_val == 0) continue;

            size_t seq_id_found = r_index.seqId(current_text_pos);
            size_t offset_from_start = r_index.seqOffset(current_text_pos);
            if (seq_id_found != source_seq_id) continue;

            auto [ext_node_id, ext_is_rev] = decode_tag(tag_val);
            gbwt::node_type ext_node = gbwt::Node::encode(ext_node_id, ext_is_rev);
            std::vector<gbwt::size_type> sa_values = gbwt_fast_locate->decompressSA(ext_node);
            bool target_passes = false;
            for (size_t i = 0; i < sa_values.size(); ++i) {
                if (gbwt_fast_locate->seqId(sa_values[i]) == target_seq_id) {
                    target_passes = true;
                    break;
                }
            }

            if (debug) {
                cerr << "  [extend] check node_id=" << ext_node_id << " is_rev=" << ext_is_rev
                     << " tag_code=" << tag_val << " source_base=" << offset_from_start
                     << " sa_values.size=" << sa_values.size()
                     << " target_seq_id=" << target_seq_id
                     << " target_passes=" << (target_passes ? "YES" : "no");
                if (!target_passes) {
                    cerr << " (seq_ids visiting: ";
                    size_t max_print = std::min<size_t>(sa_values.size(), 8);
                    for (size_t i = 0; i < max_print; ++i) {
                        if (i > 0) cerr << ",";
                        cerr << gbwt_fast_locate->seqId(sa_values[i]);
                    }
                    if (sa_values.size() > max_print) cerr << ",...";
                    cerr << ")";
                }
                cerr << endl;
            }

            if (!target_passes) continue;

            // Found a common node before the interval start: record as last common anchor.
            found_last_common = true;
            last_common_node = ext_node;
            last_common_source_base = offset_from_start;

            if (tag_map.find(tag_val) == tag_map.end()) {
                tag_map[tag_val] = TagInfo{tag_val, {}, {}, {}};
            }
            tag_map[tag_val].source_offsets.push_back(offset_from_start);
            tag_map[tag_val].source_bwt_positions.push_back(current_bwt_pos);
            tag_map[tag_val].source_packed_positions.push_back(current_text_pos);

            // Source visit's offset_in_node = i with smallest seqOffset (latest in path) — same as Phase 1.
            gbwt::size_type source_offset_in_node = gbwt::invalid_offset();
            gbwt::size_type min_seq_offset = gbwt::invalid_offset();
            for (size_t i = 0; i < sa_values.size(); ++i) {
                if (gbwt_fast_locate->seqId(sa_values[i]) == source_seq_id) {
                    gbwt::size_type so = gbwt_fast_locate->seqOffset(sa_values[i]);
                    if (min_seq_offset == gbwt::invalid_offset() || so < min_seq_offset) {
                        min_seq_offset = so;
                        source_offset_in_node = i;
                    }
                }
            }
            if (source_offset_in_node != gbwt::invalid_offset()) {
                source_gbwt_edge = gbwt::edge_type(ext_node, source_offset_in_node);
            }
            if (out_timing) {
                out_timing->found_last_common = true;
                out_timing->last_common_source_base = last_common_source_base;
            }
            if (debug) {
                cerr << "  [extend] Last common anchor found before interval: tag_code=" << tag_val
                     << " (node_id=" << ext_node_id << ", is_rev=" << ext_is_rev << ")"
                     << ", source_base=" << offset_from_start
                     << " (interval starts at " << seq_start << ")" << endl;
            }
        }
        if (debug && !found_last_common) {
            cerr << "  [extend] No common source tag found before interval start (reached haplotype start)" << endl;
        }
    }

    // Fast path was enabled but no last common node found (even after extending): no mapping possible.
    if (use_fast_path && !found_last_common) {
        auto elapsed_ms = duration_cast<milliseconds>(high_resolution_clock::now() - t_find_tags_start).count();
        cerr << "No last common node found between source and target in interval [" << seq_start << ", " << seq_end
             << "] (extended search reached haplotype start); no mapping possible. (total time: " << elapsed_ms << " ms)" << endl;
        return {};
    }

    // ----- Phase 2 (fast path only): From last common node, walk backward with GBWT inverseLF until offset < seq_start -----
    // Skip Phase 2 if the anchor is already before seq_start (extended-search case): walking further
    // backward would only collect tags outside the interval; the single anchor is enough for downstream.
    if (use_fast_path && found_last_common && source_gbwt_edge.first != gbwt::ENDMARKER &&
        last_common_source_base >= seq_start) {
        size_t current_base = last_common_source_base;
        gbwt::edge_type current_edge = source_gbwt_edge;
        
        if (debug) {
            cerr << "  [fast path] Walking backward with inverseLF from base=" << current_base << endl;
        }
        
        while (true) {
            gbwt::edge_type prev_edge = gbwt_index->inverseLF(current_edge);
            if (prev_edge.first == gbwt::ENDMARKER) {
                if (debug) cerr << "  [fast path] inverseLF returned ENDMARKER, stopping" << endl;
                break;
            }
            
            gbwtgraph::handle_t prev_handle = graph->get_handle(gbwt::Node::id(prev_edge.first), gbwt::Node::is_reverse(prev_edge.first));
            size_t prev_node_len = graph->get_length(prev_handle);
            if (current_base < prev_node_len) {
                if (debug) cerr << "  [fast path] current_base=" << current_base << " < prev_node_len=" << prev_node_len << ", underflow, stopping" << endl;
                break;
            }
            size_t new_base = current_base - prev_node_len;
            
            int64_t prev_node_id = gbwt::Node::id(prev_edge.first);
            bool prev_is_rev = gbwt::Node::is_reverse(prev_edge.first);
            uint64_t prev_tag_code = SampledTagArray::encode_value(prev_node_id, prev_is_rev);
            
            if (tag_map.find(prev_tag_code) == tag_map.end()) {
                tag_map[prev_tag_code] = TagInfo{prev_tag_code, {}, {}, {}};
            }
            tag_map[prev_tag_code].source_offsets.push_back(new_base);
            tag_map[prev_tag_code].source_bwt_positions.push_back(0);
            tag_map[prev_tag_code].source_packed_positions.push_back(0);
            
            if (debug) {
                cerr << "  [fast path] prev node id=" << prev_node_id << " rev=" << prev_is_rev
                     << " tag_code=" << prev_tag_code << " new_base=" << new_base << endl;
            }
            
            if (new_base < seq_start) {
                if (debug) cerr << "  [fast path] new_base < seq_start, stopping backward walk" << endl;
                break;
            }
            
            current_base = new_base;
            current_edge = prev_edge;
        }
    } else if (!use_fast_path) {
        // ----- No fast path: continue backward past interval start for one more tag (legacy behavior) -----
        size_t text_pos_haplotype_start = r_index.pack(source_seq_id, 0);
        bool found_tag_before_interval = false;
        
        if (debug) {
            cerr << "Looking for first tag before interval start (text_pos_i=" << text_pos_i 
                 << ", haplotype_start=" << text_pos_haplotype_start << ")" << endl;
        }
        
        if (current_text_pos > text_pos_haplotype_start) {
            current_text_pos--;
            current_bwt_pos = r_index.LF(current_bwt_pos);
            
            while (current_text_pos >= text_pos_haplotype_start && !found_tag_before_interval) {
                size_t sampled_run_id = sampled.run_id_at(current_bwt_pos);
                uint64_t tag_val = sampled.run_value(sampled_run_id);
                
                if (debug) {
                    cerr << "  [before interval] Text pos " << current_text_pos << " -> BWT pos " 
                         << current_bwt_pos << " -> tag=" << tag_val << endl;
                }
                
                if (tag_val != 0) {
                    size_t seq_id_found = r_index.seqId(current_text_pos);
                    size_t offset_from_start = r_index.seqOffset(current_text_pos);
                    
                    if (seq_id_found == source_seq_id) {
                        if (tag_map.find(tag_val) == tag_map.end()) {
                            tag_map[tag_val] = TagInfo{tag_val, {}, {}, {}};
                        }
                        
                        tag_map[tag_val].source_offsets.push_back(offset_from_start);
                        tag_map[tag_val].source_bwt_positions.push_back(current_bwt_pos);
                        tag_map[tag_val].source_packed_positions.push_back(current_text_pos);
                        
                        found_tag_before_interval = true;
                        
                        if (debug) {
                            cerr << "  Found tag before interval: tag_code=" << tag_val 
                                 << " at offset=" << offset_from_start << endl;
                        }
                    }
                }
                
                if (current_text_pos > text_pos_haplotype_start && !found_tag_before_interval) {
                    current_text_pos--;
                    current_bwt_pos = r_index.LF(current_bwt_pos);
                } else {
                    break;
                }
            }
        }
        
        if (debug && !found_tag_before_interval) {
            cerr << "  No tag found before interval (reached haplotype start)" << endl;
        }
    }

    // Convert map to vector and sort by first source offset (earliest in interval)
    for (auto& kv : tag_map) {
        tags.push_back(kv.second);
    }
    
    sort(tags.begin(), tags.end(), [](const TagInfo& a, const TagInfo& b) {
        if (a.source_offsets.empty() || b.source_offsets.empty()) return false;
        return a.source_offsets[0] < b.source_offsets[0];
    });

    if (debug) {
        cerr << "Found " << tags.size() << " unique tags in source interval" << endl;
    }
    if (debug) {
        cerr << "Tags found in source interval:" << endl;
        for (const auto& tag : tags) {
            cerr << "Tag code: " << tag.tag_code << endl;
            cerr << "Source offsets: [";
            for (size_t i = 0; i < tag.source_offsets.size(); i++) {
                if (i > 0) cerr << ", ";
                cerr << tag.source_offsets[i];
            }
            cerr << "]" << endl;
        }
    }
    
    return tags;
}

// Find all sequences (and their offsets) that pass through a given tag
vector<NodeVisit> find_sequences_for_tag(FastLocate& r_index, SampledTagArray& sampled,
                                          uint64_t tag_code) {
    vector<NodeVisit> visits;

    g_find_seq_stats.calls++;   // diagnostics (see FindSeqStats)

    if (debug) {
        cerr << "[find_sequences_for_tag] Searching for tag_code=" << tag_code << endl;
    }

    const auto& wm = sampled.values();
    size_t total_occ = wm.rank(wm.size(), tag_code);
    
    if (debug) {
        cerr << "  wm.size()=" << wm.size() << ", total_occ=" << total_occ << endl;
    }
    
    if (total_occ == 0) {
        if (debug) {
            cerr << "  No occurrences found, returning empty" << endl;
        }
        return visits;
    }
    
    bool first_run_is_gap = sampled.is_first_run_gap();
    if (debug) {
        cerr << "  first_run_is_gap=" << first_run_is_gap << endl;
    }
    
    for (size_t j = 1; j <= total_occ; ++j) {
        if (debug) {
            cerr << "  Processing occurrence " << j << "/" << total_occ << endl;
        }
        
        size_t tag_run_id = wm.select(j, tag_code);
        if (debug) {
            cerr << "    tag_run_id=" << tag_run_id << endl;
        }
        
        size_t bwt_run_id = 2 * tag_run_id + static_cast<size_t>(first_run_is_gap);
        if (debug) {
            cerr << "    bwt_run_id=" << bwt_run_id << endl;
        }
        
        auto run_span = sampled.run_span(bwt_run_id);
        size_t locate_start = run_span.first;
        size_t locate_end = run_span.second;
        
        if (debug) {
            cerr << "    run_span=[" << locate_start << ", " << locate_end 
                 << "], size=" << (locate_end - locate_start + 1) << endl;
        }
        
        // IMPORTANT: Skip positions >= r_index.bwt_size() (those are endmarkers)
        // The sampled tag array includes endmarkers, but RLBWT r-index doesn't
        size_t r_index_size = r_index.bwt_size();
        if (locate_start >= r_index_size) {
            if (debug) {
                cerr << "    Skipping entire run: locate_start=" << locate_start 
                     << " >= r_index.bwt_size()=" << r_index_size << " (endmarker positions)" << endl;
            }
            continue;  // Skip this entire run
        }
        
        // Clamp locate_end to valid range
        if (locate_end >= r_index_size) {
            if (debug) {
                cerr << "    Clamping locate_end from " << locate_end << " to " << (r_index_size - 1) 
                     << " (r_index.bwt_size()=" << r_index_size << ")" << endl;
            }
            locate_end = r_index_size - 1;
        }
        
        // Get initial packed position
        size_t packed_pos;
        size_t stats_nav_steps = 0;   // LF steps to navigate to this run's start (diagnostics)
        {
            size_t rindex_run_id = 0;
            size_t run_start_pos = 0;
            
            if (debug) {
                cerr << "    Getting initial packed position at locate_start=" << locate_start << endl;
                cerr << "      r_index.bwt_size()=" << r_index.bwt_size() << endl;
                cerr << "      r_index.size() (num runs)=" << r_index.size() << endl;
            }
            
            r_index.run_id_and_offset_at(locate_start, rindex_run_id, run_start_pos);
            
            if (debug) {
                cerr << "      rindex_run_id=" << rindex_run_id << ", run_start_pos=" << run_start_pos << endl;
            }
            
            packed_pos = r_index.getSample(rindex_run_id);
            
            if (debug) {
                cerr << "      initial packed_pos=" << packed_pos << endl;
            }
            
            // Navigate from run start to locate_start
            size_t nav_steps = locate_start - run_start_pos;
            stats_nav_steps = nav_steps;   // diagnostics
            if (debug) {
                cerr << "      Need to navigate " << nav_steps << " steps from run_start_pos="
                     << run_start_pos << " to locate_start=" << locate_start << endl;
            }
            
            for (size_t p = run_start_pos; p < locate_start; ++p) {
                if (debug) {
                    cerr << "      Before locateNext: p=" << p << ", packed_pos=" << packed_pos << endl;
                }
                packed_pos = r_index.locateNext(packed_pos);
                if (debug) {
                    cerr << "      After locateNext: p=" << p << ", packed_pos=" << packed_pos << endl;
                }
            }
            
            if (debug) {
                cerr << "      After navigation: packed_pos=" << packed_pos << endl;
            }
        }
        
        // Process all positions in the run
        for (size_t pos = locate_start; pos <= locate_end; ++pos) {
            if (debug) {
                cerr << "    Processing pos=" << pos << " (offset=" << (pos - locate_start) 
                     << "/" << (locate_end - locate_start) << ")" << endl;
            }
            if (debug && ((pos - locate_start) % 100 == 0 || pos == locate_start || pos == locate_end)) {
                cerr << "    Processing pos=" << pos << " (offset=" << (pos - locate_start) 
                     << "/" << (locate_end - locate_start) << ")" << endl;
            }
            
            if (pos > locate_start) {
                if (debug) {
                    cerr << "      Before locateNext: pos=" << pos << ", packed_pos=" << packed_pos << endl;
                    cerr << "        Checking if packed_pos is valid (< bwt_size=" << r_index.bwt_size() << ")..." << endl;
                }
                
                // if (packed_pos >= r_index.bwt_size()) {
                //     cerr << "ERROR: packed_pos=" << packed_pos << " is >= bwt_size=" << r_index.bwt_size() << endl;
                //     cerr << "  This should not happen! packed_pos should be < bwt_size" << endl;
                //     return visits;  // Return what we have so far
                // }
                
                packed_pos = r_index.locateNext(packed_pos);
                
                if (debug) {
                    cerr << "      After locateNext: pos=" << pos << ", packed_pos=" << packed_pos << endl;
                }
                if (debug && ((pos - locate_start) % 100 == 0 || pos == locate_end)) {
                    cerr << "      After locateNext: packed_pos=" << packed_pos << endl;
                }
            }
            
            // Unpack to get (seq_id, offset)
            if (debug && ((pos - locate_start) % 100 == 0 || pos == locate_start || pos == locate_end)) {
                cerr << "      About to unpack packed_pos=" << packed_pos << endl;
            }
            auto pr = r_index.unpack(packed_pos);

            size_t seq_id_found = pr.first;
            size_t offset_from_start = pr.second;
            
            if (debug && ((pos - locate_start) % 100 == 0 || pos == locate_start || pos == locate_end)) {
                cerr << "      Unpacked: seq_id=" << seq_id_found << ", offset=" << offset_from_start << endl;
            }
            
            NodeVisit visit;
            visit.seq_id = seq_id_found;
            visit.offset = offset_from_start;
            visit.bwt_pos = pos;
            visit.packed_pos = packed_pos;
            visit.tag_code = tag_code;
            
            visits.push_back(visit);
            // std::cerr << "done finding sequences for tag" << endl;
        }

        // Diagnostics: LF cost of this run = navigation (sample → run start) +
        // run walk (one locateNext per position after the first). The "last
        // run" fields end up holding the values for the final run iterated.
        {
            const size_t run_len = locate_end - locate_start + 1;
            const size_t walk_steps = locate_end - locate_start;  // locateNext calls in the walk
            g_find_seq_stats.runs++;
            g_find_seq_stats.lf_steps += stats_nav_steps + walk_steps;
            g_find_seq_stats.visits += run_len;
            g_find_seq_stats.last_run_nav_steps = stats_nav_steps;
            g_find_seq_stats.last_run_length = run_len;
        }

        if (debug) {
            cerr << "    Finished occurrence " << j << ", total visits so far: " << visits.size() << endl;
        }
    }
    
    if (debug) {
        cerr << "  Total visits found: " << visits.size() << endl;
    }
    
    return visits;
}

// Find all common nodes between source and target haplotypes (handles fragments)
// Returns vector of anchor points: (source_offset, target_offset, tag_code, source_packed_pos, target_packed_pos, target_seq_id)
vector<tuple<size_t, size_t, uint64_t, size_t, size_t, size_t>> find_all_common_nodes(
    FastLocate& r_index, SampledTagArray& sampled,
    const vector<TagInfo>& source_tags, size_t target_seq_id) {
    
    vector<tuple<size_t, size_t, uint64_t, size_t, size_t, size_t>> anchors;
    
    for (const auto& tag_info : source_tags) {
        // Find all sequences that pass through this tag
        vector<NodeVisit> all_visits = find_sequences_for_tag(r_index, sampled, tag_info.tag_code);
        
        // Find visits by target haplotype (may include multiple fragments)
        vector<NodeVisit> target_visits;
        for (const auto& visit : all_visits) {
            if (visit.seq_id == target_seq_id) {
                target_visits.push_back(visit);
            }
        }
        
        if (!target_visits.empty()) {
            // Found common node(s)! Match visits by SA value ordering
            
            for (size_t src_idx = 0; src_idx < tag_info.source_offsets.size(); src_idx++) {
                size_t source_offset = tag_info.source_offsets[src_idx];
                size_t source_packed_pos = tag_info.source_packed_positions[src_idx];
                
                // Find target visit with closest SA value
                size_t best_target_idx = 0;
                size_t min_diff = numeric_limits<size_t>::max();
                
                for (size_t tgt_idx = 0; tgt_idx < target_visits.size(); tgt_idx++) {
                    size_t diff = (source_packed_pos > target_visits[tgt_idx].packed_pos) 
                        ? (source_packed_pos - target_visits[tgt_idx].packed_pos)
                        : (target_visits[tgt_idx].packed_pos - source_packed_pos);
                    
                    if (diff < min_diff) {
                        min_diff = diff;
                        best_target_idx = tgt_idx;
                    }
                }
                
                size_t target_offset = target_visits[best_target_idx].offset;
                size_t target_packed_pos = target_visits[best_target_idx].packed_pos;
                size_t target_fragment_seq_id = target_visits[best_target_idx].seq_id;
                
                if (debug) {
                    uint64_t decoded = tag_info.tag_code - 1;
                    int64_t node_id = (decoded >> 1) + 1;
                    bool is_rev = (decoded & 1) != 0;
                    cerr << "Found common node: node_id=" << node_id 
                         << ", is_rev=" << is_rev
                         << ", source_offset=" << source_offset
                         << ", target_seq_id=" << target_fragment_seq_id
                         << ", target_offset=" << target_offset << endl;
                }
                
                anchors.push_back(make_tuple(source_offset, target_offset, tag_info.tag_code, 
                                            source_packed_pos, target_packed_pos, target_fragment_seq_id));
            }
        }
    }
    
    return anchors;
}

// Decode tag code to get node ID and orientation
pair<int64_t, bool> decode_tag(uint64_t tag_code) {
    uint64_t decoded = tag_code - 1;
    int64_t node_id = (decoded >> 1) + 1;
    bool is_rev = (decoded & 1) != 0;
    return {node_id, is_rev};
}

// Convert node offset to base offset for a given GBWT sequence
// Returns the cumulative base offset at the start of the node at node_offset
// Requires GBWTGraph to get node lengths
size_t node_offset_to_base_offset(
    const gbwt::GBWT& gbwt_index,
    const gbwtgraph::GBWTGraph& graph,
    size_t seq_id,
    size_t node_offset) {
    
    // Extract the path for this sequence
    gbwt::vector_type path = gbwt_index.extract(gbwt::Path::encode(seq_id, false));
    
    if (path.empty()) {
        cerr << "Warning: Sequence " << seq_id << " not found in GBWT" << endl;
        return 0;
    }
    
    size_t cumulative_bases = 0;
    size_t current_node_offset = 0;
    
    for (gbwt::node_type node : path) {
        if (node == gbwt::ENDMARKER) break;
        if (current_node_offset >= node_offset) break;
        
        // Get node handle from graph
        gbwtgraph::handle_t handle = graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node));
        size_t node_length = graph.get_length(handle);
        
        cumulative_bases += node_length;
        current_node_offset++;
    }
    
    return cumulative_bases;
}

// Convert base offset to node offset for a given GBWT sequence
// Returns the node offset (number of nodes from start) that contains the base offset
// Requires GBWTGraph to get node lengths
size_t base_offset_to_node_offset(
    const gbwt::GBWT& gbwt_index,
    const gbwtgraph::GBWTGraph& graph,
    size_t seq_id,
    size_t base_offset) {
    
    if (debug) {
        cerr << "[base_offset_to_node_offset] seq_id=" << seq_id 
             << ", base_offset=" << base_offset << endl;
    }
    
    // Extract the path for this sequence
    gbwt::vector_type path = gbwt_index.extract(gbwt::Path::encode(seq_id, false));
    
    if (path.empty()) {
        cerr << "Warning: Sequence " << seq_id << " not found in GBWT" << endl;
        return 0;
    }
    
    if (debug) {
        cerr << "  Path length: " << path.size() << " nodes" << endl;
    }
    
    size_t cumulative_bases = 0;
    size_t node_offset = 0;
    
    for (gbwt::node_type node : path) {
        if (node == gbwt::ENDMARKER) break;
        
        // Get node handle from graph
        gbwtgraph::handle_t handle = graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node));
        size_t node_length = graph.get_length(handle);
        
        if (debug) {
            cerr << "  node_offset=" << node_offset << ", node_id=" << gbwt::Node::id(node)
                 << ", cumulative_bases=[" << cumulative_bases << ", " 
                 << (cumulative_bases + node_length) << ")" << endl;
        }
        
        // Check if base_offset is within this node
        // cumulative_bases is the start of current node, cumulative_bases + node_length is the end
        if (base_offset < cumulative_bases + node_length) {
            // Base offset is within this node
            if (debug) {
                cerr << "  Found! base_offset " << base_offset << " is in node_offset " << node_offset << endl;
            }
            return node_offset;
        }
        
        cumulative_bases += node_length;
        node_offset++;
    }
    
    // Base offset is beyond the path, return last node offset
    if (debug) {
        cerr << "  Warning: base_offset " << base_offset << " is beyond path end (cumulative_bases=" 
             << cumulative_bases << "), returning last node_offset=" << (node_offset > 0 ? node_offset - 1 : 0) << endl;
    }
    return node_offset > 0 ? node_offset - 1 : 0;
}

// Structure to hold first and last common nodes
struct CommonNodes {
    size_t first_source_offset;      // GBWT node offset
    size_t first_target_offset;      // GBWT node offset
    size_t first_source_base;        // RLBWT base offset
    size_t first_target_base;        // RLBWT base offset
    uint64_t first_tag_code;
    size_t last_source_offset;       // GBWT node offset
    size_t last_target_offset;      // GBWT node offset
    size_t last_source_base;         // RLBWT base offset
    size_t last_target_base;        // RLBWT base offset
    uint64_t last_tag_code;
    bool found;
};

// Structure to store a node visit along a path during traversal
struct PathNode {
    uint64_t tag_code;               // Tag code for this node
    gbwt::node_type node;            // GBWT node
    size_t node_offset;              // Node offset from start of path
    size_t base_offset;               // Base offset from start of path
    gbwt::edge_type gbwt_position;   // GBWT edge position for LF mapping
};

// Helper function to check if a tag is a common node and extract offsets
// use_largest_offset controls offset selection:
//   false: smallest RLBWT base + largest GBWT offset (both = earliest) -> for FIRST node
//   true:  largest RLBWT base + smallest GBWT offset (both = latest) -> for LAST node
bool check_common_node(
    const gbwt::FastLocate& gbwt_fast_locate,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const TagInfo& tag_info,
    size_t source_seq_id, size_t target_seq_id,
    size_t& source_offset, size_t& target_offset,
    size_t& source_base, size_t& target_base,
    bool use_largest_offset = false) {
    
    // Decode tag to get GBWT node
    auto [node_id, is_rev] = decode_tag(tag_info.tag_code);
    gbwt::node_type node = gbwt::Node::encode(node_id, is_rev);
    
    if (debug) {
        cerr << "  Checking tag_code=" << tag_info.tag_code 
             << " (node_id=" << node_id << ", is_rev=" << is_rev << ")" << endl;
    }
    
    // Use GBWT FastLocate to find all paths visiting this node
    std::vector<gbwt::size_type> sa_values = gbwt_fast_locate.decompressSA(node);
    
    if (debug) {
        cerr << "    Found " << sa_values.size() << " path occurrences on this node" << endl;
    }
    
    // Find source and target visits
    // Note: In GBWT, larger offset = earlier in path (opposite of RLBWT)
    vector<pair<gbwt::size_type, size_t>> source_visits;  // (sa_value, path_offset)
    vector<pair<gbwt::size_type, size_t>> target_visits;  // (sa_value, path_offset)
    
    for (size_t i = 0; i < sa_values.size(); i++) {
        gbwt::size_type seq_id = gbwt_fast_locate.seqId(sa_values[i]);
        gbwt::size_type seq_offset = gbwt_fast_locate.seqOffset(sa_values[i]);

        if (debug) {
            cerr << "    SA value: " << sa_values[i] << " (seq_id=" << seq_id << ", seq_offset=" << seq_offset << ")" << endl;
        }
        
        if (seq_id == source_seq_id) {
            source_visits.push_back({sa_values[i], seq_offset});
        }
        if (seq_id == target_seq_id) {
            target_visits.push_back({sa_values[i], seq_offset});
        }
    }
    
    if (debug) {
        cerr << "    Source visits: " << source_visits.size() << endl;
        cerr << "    Target visits: " << target_visits.size() << endl;
    }
    
    if (source_visits.empty() || target_visits.empty()) {
        return false;
    }
    
    // Found common node! Match visits based on offsets
    // Sort visits by offset (larger offset = earlier position in path)
    sort(source_visits.begin(), source_visits.end(),
         [](const pair<gbwt::size_type, size_t>& a, const pair<gbwt::size_type, size_t>& b) {
             return a.second > b.second;  // Larger offset first (earlier in path)
         });
    sort(target_visits.begin(), target_visits.end(),
         [](const pair<gbwt::size_type, size_t>& a, const pair<gbwt::size_type, size_t>& b) {
             return a.second > b.second;  // Larger offset first (earlier in path)
         });
    
    // Select GBWT offset based on whether we want first or last common node
    // IMPORTANT: In GBWT, LARGER offset = EARLIER in path (opposite of RLBWT!)
    if (use_largest_offset) {
        // use_largest_offset = true: LAST common node
        // Select SMALLEST GBWT offset (latest in path)
        source_offset = source_visits[source_visits.size() - 1].second;  // last index = smallest offset
        target_offset = target_visits[target_visits.size() - 1].second;
        if (debug) {
            cerr << "    Selected SMALLEST GBWT offsets (latest in path, for last common node): source=" 
                 << source_offset << ", target=" << target_offset << endl;
        }
    } else {
        // use_largest_offset = false: FIRST common node
        // Select LARGEST GBWT offset (earliest in path)
        source_offset = source_visits[0].second;  // index 0 = largest offset after sorting
        target_offset = target_visits[0].second;
        if (debug) {
            cerr << "    Selected LARGEST GBWT offsets (earliest in path, for first common node): source=" 
                 << source_offset << ", target=" << target_offset << endl;
        }
    }
    
    // Get base offsets from RLBWT r-index for this tag
    // Collect all source and target visits from RLBWT
    vector<NodeVisit> source_visits_rlbwt = find_sequences_for_tag(rlbwt_rindex, sampled, tag_info.tag_code);
    
    // Separate source and target visits
    vector<size_t> source_base_offsets;
    vector<size_t> target_base_offsets;
    
    for (const auto& visit : source_visits_rlbwt) {
        if (visit.seq_id == source_seq_id) {
            source_base_offsets.push_back(visit.offset);
        } else if (visit.seq_id == target_seq_id) {
            target_base_offsets.push_back(visit.offset);
        }
    }
    
    if (source_base_offsets.empty() || target_base_offsets.empty()) {
        cerr << "  No common nodes found in RLBWT But in GBWT we found both source and target visits" << endl;
        return false;  // Not common to both
    }
    
    // Sort offsets (ascending order)
    sort(source_base_offsets.begin(), source_base_offsets.end());
    sort(target_base_offsets.begin(), target_base_offsets.end());
    
    // For RLBWT: smaller offset = earlier in sequence, larger offset = later in sequence
    if (use_largest_offset) {
        // use_largest_offset = true: LAST common node
        // Select largest RLBWT base offset (latest in sequence)
        source_base = source_base_offsets.back();
        target_base = target_base_offsets.back();
        
        // Verify: the largest source offset from find_sequences_for_tag should match
        // the largest offset in tag_info.source_offsets (which came from find_tags_in_interval)
        if (!tag_info.source_offsets.empty()) {
            size_t expected_source_base = *max_element(tag_info.source_offsets.begin(), tag_info.source_offsets.end());
            assert(source_base == expected_source_base && 
                   "source_base from find_sequences_for_tag doesn't match expected from TagInfo");
        }
    } else {
        // use_largest_offset = false: FIRST common node
        // Select smallest RLBWT base offset (earliest in sequence)
        source_base = source_base_offsets[0];
        target_base = target_base_offsets[0];
        
        // Verify: the smallest source offset from find_sequences_for_tag should match
        // the smallest offset in tag_info.source_offsets (which came from find_tags_in_interval)
        if (!tag_info.source_offsets.empty()) {
            size_t expected_source_base = *min_element(tag_info.source_offsets.begin(), tag_info.source_offsets.end());
            assert(source_base == expected_source_base && 
                   "source_base from find_sequences_for_tag doesn't match expected from TagInfo");
        }
    }

    if (debug) {
        cerr << "    RLBWT base offsets:" << endl;
        if (use_largest_offset) {
            cerr << "      Selected LARGEST RLBWT offsets (latest in sequence, for last common node):" << endl;
        } else {
            cerr << "      Selected SMALLEST RLBWT offsets (earliest in sequence, for first common node):" << endl;
        }
        cerr << "      Source base: " << source_base << endl;
        cerr << "      Target base: " << target_base << endl;
    }

    
    return true;
}

// Find the first (smallest tag_code) and last (largest tag_code) common nodes
// between source and target haplotypes using GBWT FastLocate
// Optimized: finds first by iterating from beginning, finds last by iterating from end
// Also looks up base offsets from RLBWT r-index
// Returns CommonNodes structure with both first and last common nodes
CommonNodes find_first_and_last_common_nodes_gbwt(
    const gbwt::FastLocate& gbwt_fast_locate,
    FastLocate& rlbwt_rindex,
    SampledTagArray& sampled,
    const vector<TagInfo>& source_tags, 
    size_t source_seq_id, size_t target_seq_id) {
    
    CommonNodes result;
    result.found = false;
    
    // source_tags is already sorted by source_offsets[0] ascending (earliest in path first)
    // from find_tags_in_interval, so no need to sort again.
    
    if (debug) {
        cerr << "Finding first and last common nodes (optimized search):" << endl;
        cerr << "  Total tags: " << source_tags.size() << endl;
        cerr << "  Tags (already by source base offset):" << endl;
        for (size_t i = 0; i < source_tags.size(); i++) {
            auto [node_id, is_rev] = decode_tag(source_tags[i].tag_code);
            size_t offset = source_tags[i].source_offsets.empty() ? 0 : source_tags[i].source_offsets[0];
            cerr << "    [" << i << "] tag_code=" << source_tags[i].tag_code 
                 << " (node_id=" << node_id << "), base_offset=" << offset << endl;
        }
    }
    
    // Step 1: Find first common node by iterating from beginning (earliest in path)
    if (debug) {
        cerr << "  Searching for FIRST common node (earliest in path)..." << endl;
    }
    
    for (size_t i = 0; i < source_tags.size(); i++) {
        const auto& tag_info = source_tags[i];
        size_t source_offset, target_offset, source_base, target_base;
        
        // For first common node: smallest RLBWT offset + largest GBWT offset (both = earliest)
        if (check_common_node(gbwt_fast_locate, rlbwt_rindex, sampled, tag_info,
                              source_seq_id, target_seq_id,
                              source_offset, target_offset, source_base, target_base,
                              false)) {  // false = smallest RLBWT + largest GBWT (earliest)
            // Found first common node!
            result.first_source_offset = source_offset;
            result.first_target_offset = target_offset;
            result.first_source_base = source_base;
            result.first_target_base = target_base;
            result.first_tag_code = tag_info.tag_code;
            result.found = true;
            
            if (debug) {
                cerr << "  ✓ Found FIRST common node (earliest in path)!" << endl;
                cerr << "    (smallest RLBWT base, largest GBWT offset)" << endl;
                cerr << "    Source: GBWT_offset=" << source_offset 
                     << ", RLBWT_base=" << source_base << endl;
                cerr << "    Target: GBWT_offset=" << target_offset
                     << ", RLBWT_base=" << target_base << endl;
            }
            break;
        }
    }
    
    if (!result.found) {
        if (debug) {
            cerr << "  No common nodes found" << endl;
        }
        return result;
    }
    
    // Step 2: Find last common node by iterating from end (latest in path)
    if (debug) {
        cerr << "  Searching for LAST common node (latest in path)..." << endl;
    }
    
    for (size_t i = source_tags.size(); i > 0; i--) {
        const auto& tag_info = source_tags[i - 1];
        size_t source_offset, target_offset, source_base, target_base;
        
        // For last common node: largest RLBWT offset + smallest GBWT offset (both = latest)
        if (check_common_node(gbwt_fast_locate, rlbwt_rindex, sampled, tag_info,
                              source_seq_id, target_seq_id,
                              source_offset, target_offset, source_base, target_base,
                              true)) {  // true = largest RLBWT + smallest GBWT (latest)
            // Found last common node!
            result.last_source_offset = source_offset;
            result.last_target_offset = target_offset;
            result.last_source_base = source_base;
            result.last_target_base = target_base;
            result.last_tag_code = tag_info.tag_code;
            
            if (debug) {
                cerr << "  ✓ Found LAST common node (latest in path)!" << endl;
                cerr << "    (largest RLBWT base, smallest GBWT offset)" << endl;
                cerr << "    Source: GBWT_offset=" << source_offset 
                     << ", RLBWT_base=" << source_base << endl;
                cerr << "    Target: GBWT_offset=" << target_offset
                     << ", RLBWT_base=" << target_base << endl;
            }
            break;
        }
    }
    
    if (debug) {
        cerr << "\n  Summary:" << endl;
        cerr << "    First common node: tag_code=" << result.first_tag_code << endl;
        cerr << "      Source: GBWT_node_offset=" << result.first_source_offset 
             << ", RLBWT_base_offset=" << result.first_source_base << endl;
        cerr << "      Target: GBWT_node_offset=" << result.first_target_offset 
             << ", RLBWT_base_offset=" << result.first_target_base << endl;
        cerr << "    Last common node: tag_code=" << result.last_tag_code << endl;
        cerr << "      Source: GBWT_node_offset=" << result.last_source_offset 
             << ", RLBWT_base_offset=" << result.last_source_base << endl;
        cerr << "      Target: GBWT_node_offset=" << result.last_target_offset 
             << ", RLBWT_base_offset=" << result.last_target_base << endl;
    }
    
    return result;
}

// Trace paths using GBWT LF mapping (one node at a time)
// Uses GBWT's efficient forward traversal to map positions between source and target haplotypes
// FastLocate is used to find initial positions, then GBWT index LF() is used for traversal
// Each path (source and target) is traversed independently until it reaches its respective
// last common node at the correct base offset.
unordered_map<size_t, vector<size_t>> build_offset_mapping_with_gbwt_lf(
    const gbwt::GBWT& gbwt_index,
    const gbwt::FastLocate& gbwt_fast_locate,
    const gbwtgraph::GBWTGraph& graph,
    size_t source_seq_id, size_t source_start, size_t source_end,
    size_t target_seq_id, size_t anchor_source_offset, size_t anchor_target_offset,
    size_t anchor_source_base, size_t anchor_target_base,
    uint64_t anchor_tag_code, 
    size_t last_common_source_base, size_t last_common_target_base,
    uint64_t last_common_tag_code) {
    
    unordered_map<size_t, vector<size_t>> offset_map;
    
    if (debug) {
        cerr << "Building offset mapping using GBWT LF mapping:" << endl;
        cerr << "  Anchor source node offset: " << anchor_source_offset << endl;
        cerr << "  Anchor target node offset: " << anchor_target_offset << endl;
        cerr << "  Last common source base: " << last_common_source_base << endl;
        cerr << "  Last common target base: " << last_common_target_base << endl;
    }
    
    // Decode anchor node and encode as GBWT node
    auto [anchor_node_id, anchor_is_rev] = decode_tag(anchor_tag_code);
    gbwt::node_type anchor_node = gbwt::Node::encode(anchor_node_id, anchor_is_rev);
    
    if (debug) {
        cerr << "  Anchor node: id=" << anchor_node_id 
             << ", is_rev=" << anchor_is_rev 
             << ", encoded=" << anchor_node << endl;
    }
    
    // Decompress SA for just the anchor node (efficient - only this node's occurrences)
    std::vector<gbwt::size_type> sa_values = gbwt_fast_locate.decompressSA(anchor_node);
    
    if (debug) {
        cerr << "  SA values for anchor node: " << sa_values.size() << " occurrences" << endl;
    }
    
    // Find the offsets (indices) for source and target sequences in the SA result
    // The index in sa_values IS the offset needed for LF mapping
    gbwt::size_type source_offset_in_node = gbwt::invalid_offset();
    gbwt::size_type target_offset_in_node = gbwt::invalid_offset();
    
    for (size_t i = 0; i < sa_values.size(); i++) {
        gbwt::size_type seq_id = gbwt_fast_locate.seqId(sa_values[i]);
        gbwt::size_type seq_offset = gbwt_fast_locate.seqOffset(sa_values[i]);
        
        // Print all indexes where source sequence visits this node
        if (seq_id == source_seq_id) {
            if (debug) {
                cerr << "  Source seq " << source_seq_id << " visits node at index " << i 
                     << ", seq_offset=" << seq_offset << endl;
            }
            if (seq_offset == anchor_source_offset) {
                source_offset_in_node = i;
                if (debug) {
                    cerr << "    ^ This is the anchor (matches anchor_source_offset=" << anchor_source_offset << ")" << endl;
                }
            }
        }
        
        // Print all indexes where target sequence visits this node
        if (seq_id == target_seq_id) {
            if (debug) {
                cerr << "  Target seq " << target_seq_id << " visits node at index " << i 
                     << ", seq_offset=" << seq_offset << endl;
            }
            if (seq_offset == anchor_target_offset) {
                target_offset_in_node = i;
                if (debug) {
                    cerr << "    ^ This is the anchor (matches anchor_target_offset=" << anchor_target_offset << ")" << endl;
                }
            }
        }
    }
    
    // For now, we assume only one occurrence per sequence at this node
    if (source_offset_in_node == gbwt::invalid_offset()) {
        cerr << "Warning: Source sequence " << source_seq_id 
             << " not found at anchor node (offset " << anchor_source_offset << ")" << endl;
        return offset_map;
    }
    if (target_offset_in_node == gbwt::invalid_offset()) {
        cerr << "Warning: Target sequence " << target_seq_id 
             << " not found at anchor node (offset " << anchor_target_offset << ")" << endl;
        return offset_map;
    }
    
    // Create GBWT positions for LF mapping
    // edge_type = (node, offset_within_node's_record)
    gbwt::edge_type source_pos = {anchor_node, source_offset_in_node};
    gbwt::edge_type target_pos = {anchor_node, target_offset_in_node};
    
    if (debug) {
        cerr << "  Starting positions for LF traversal:" << endl;
        cerr << "    Source: (node=" << source_pos.first << ", offset=" << source_pos.second << ")" << endl;
        cerr << "    Target: (node=" << target_pos.first << ", offset=" << target_pos.second << ")" << endl;
    }
    
    // Track current offsets in each path (in terms of node count from start of path)
    size_t source_path_offset = anchor_source_offset;
    size_t target_path_offset = anchor_target_offset;
    
    // Track cumulative base offsets (starting from RLBWT base offsets, then adding node lengths)
    size_t source_base_offset = anchor_source_base;
    size_t target_base_offset = anchor_target_base;
    
    // Start with anchor point (using base offsets from RLBWT)
    if (anchor_source_base >= source_start && anchor_source_base <= source_end) {
        offset_map[anchor_source_base].push_back(anchor_target_base);
    }
    
    if (debug) {
        cerr << "  Starting from anchor (base offsets from RLBWT):" << endl;
        cerr << "    Anchor: source_node=" << anchor_source_offset << ", source_base=" << anchor_source_base << endl;
        cerr << "    Anchor: target_node=" << anchor_target_offset << ", target_base=" << anchor_target_base << endl;
        cerr << "    Input interval: [" << source_start << ", " << source_end << "] bases" << endl;
    }
    
    // ========== FORWARD TRAVERSAL using GBWT LF mapping ==========
    // GBWT's LF mapping moves forward along the path (one node at a time)
    // Unlike standard FM-index LF which moves backward, GBWT's LF moves forward because:
    //   1. GBWT stores paths as sequences of nodes (not bases)
    //   2. Each node's record stores outgoing edges (next nodes)
    //   3. LF(pos) returns the next node position along the path
    //
    // IMPORTANT: Paths can diverge and converge again
    // Example: Source: 1→2→3, Target: 1→4→5→6→2→3
    // We traverse each path INDEPENDENTLY, store all nodes, then map at common nodes
    
    // Store paths as vectors of PathNode entries
    vector<PathNode> source_path;
    vector<PathNode> target_path;
    
    // Add anchor node to both paths
    PathNode anchor_src_node;
    anchor_src_node.tag_code = anchor_tag_code;
    anchor_src_node.node = anchor_node;
    anchor_src_node.node_offset = anchor_source_offset;
    anchor_src_node.base_offset = anchor_source_base;
    anchor_src_node.gbwt_position = source_pos;
    source_path.push_back(anchor_src_node);
    
    PathNode anchor_tgt_node;
    anchor_tgt_node.tag_code = anchor_tag_code;
    anchor_tgt_node.node = anchor_node;
    anchor_tgt_node.node_offset = anchor_target_offset;
    anchor_tgt_node.base_offset = anchor_target_base;
    anchor_tgt_node.gbwt_position = target_pos;
    target_path.push_back(anchor_tgt_node);
    
    if (debug) {
        cerr << "  Starting independent path traversal..." << endl;
        cerr << "  Anchor node: tag_code=" << anchor_tag_code << endl;
        cerr << "    Decoded: node_id=" << anchor_node_id << ", is_rev=" << anchor_is_rev << endl;
        cerr << "    GBWT node (encoded): " << anchor_node << endl;
        gbwtgraph::handle_t anchor_handle = graph.get_handle(anchor_node_id, anchor_is_rev);
        cerr << "    Anchor node sequence: " << graph.get_sequence(anchor_handle) << endl;
        cerr << "    Anchor node length: " << graph.get_length(anchor_handle) << endl;
        cerr << "    Source: node_offset=" << anchor_source_offset << ", base_offset=" << anchor_source_base << endl;
        cerr << "    Target: node_offset=" << anchor_target_offset << ", base_offset=" << anchor_target_base << endl;
        cerr << "    Source edge_type: (node=" << source_pos.first << ", offset=" << source_pos.second << ")" << endl;
        cerr << "    Target edge_type: (node=" << target_pos.first << ", offset=" << target_pos.second << ")" << endl;
    }
    
    // ========== Traverse source path independently ==========
    gbwt::edge_type src_pos = source_pos;
    size_t src_node_offset = anchor_source_offset;
    size_t src_base_offset = anchor_source_base;
    
    if (debug) {
        cerr << "  Traversing source path..." << endl;
    }
    
    while (true) {
        gbwt::edge_type next_src_pos = gbwt_index.LF(src_pos);
        
        if (next_src_pos.first == gbwt::ENDMARKER) {
            if (debug) {
                cerr << "    Source path stopped: reached ENDMARKER after " << source_path.size() << " nodes" << endl;
            }
            break;
        }
        
        // Update base offset by adding length of node we just left
        gbwtgraph::handle_t curr_src_handle = graph.get_handle(
            gbwt::Node::id(src_pos.first), gbwt::Node::is_reverse(src_pos.first));
        size_t node_len = graph.get_length(curr_src_handle);
        src_base_offset += node_len;
        src_node_offset++;
        
        if (debug) {
            cerr << "    Left GBWT node: " << src_pos.first 
                 << " (id=" << gbwt::Node::id(src_pos.first) 
                 << ", rev=" << gbwt::Node::is_reverse(src_pos.first) << ")" << endl;
            cerr << "    Moved forward: added node_len=" << node_len 
                 << ", new src_base_offset=" << src_base_offset << endl;
            cerr << "    Left node sequence: " << graph.get_sequence(curr_src_handle) << endl;
        }
        
        // Get tag code for next node
        gbwt::node_type next_src_node = next_src_pos.first;
        int64_t next_src_node_id = gbwt::Node::id(next_src_node);
        bool next_src_is_rev = gbwt::Node::is_reverse(next_src_node);
        uint64_t src_decoded = ((static_cast<uint64_t>(next_src_node_id) - 1) << 1) | (next_src_is_rev ? 1 : 0);
        uint64_t next_src_tag_code = src_decoded + 1;
        
        // Store this node first (we need it even if it's the last one)
        PathNode src_node;
        src_node.tag_code = next_src_tag_code;
        src_node.node = next_src_node;
        src_node.node_offset = src_node_offset;
        src_node.base_offset = src_base_offset;
        src_node.gbwt_position = next_src_pos;
        source_path.push_back(src_node);
        
        if (debug) {
            gbwtgraph::handle_t next_src_handle = graph.get_handle(
                gbwt::Node::id(next_src_node), gbwt::Node::is_reverse(next_src_node));
            cerr << "    Next GBWT node: " << next_src_node 
                 << " (id=" << next_src_node_id 
                 << ", rev=" << next_src_is_rev << ")" << endl;
            cerr << "    Added source node: tag_code=" << next_src_tag_code 
                 << ", node_offset=" << src_node_offset 
                 << ", base_offset=" << src_base_offset << endl;
            cerr << "    Added node sequence: " << graph.get_sequence(next_src_handle) << endl;
        }
        
        // Check stop condition: stop when we've reached or passed the last common node's base offset
        // We store the node first, then check if we should stop
        if (src_base_offset >= last_common_source_base) {
            if (debug) {
                cerr << "    Source path stopped: reached last common base offset " 
                     << last_common_source_base << " (current: " << src_base_offset << ")"
                     << " after " << source_path.size() << " nodes" << endl;
            }
            break;
        }
        
        // Also stop if we've gone past the source interval
        if (src_base_offset > source_end) {
            if (debug) {
                cerr << "    Source path stopped: src_base_offset=" << src_base_offset 
                     << " > source_end=" << source_end 
                     << " after " << source_path.size() << " nodes" << endl;
            }
            break;
        }
        
        src_pos = next_src_pos;
    }
    
    // ========== Traverse target path independently ==========
    gbwt::edge_type tgt_pos = target_pos;
    size_t tgt_node_offset = anchor_target_offset;
    size_t tgt_base_offset = anchor_target_base;
    
    if (debug) {
        cerr << "  Traversing target path..." << endl;
    }
    
    while (true) {
        gbwt::edge_type next_tgt_pos = gbwt_index.LF(tgt_pos);
        
        if (next_tgt_pos.first == gbwt::ENDMARKER) {
            if (debug) {
                cerr << "    Target path stopped: reached ENDMARKER after " << target_path.size() << " nodes" << endl;
            }
            break;
        }
        
        // Update base offset by adding length of node we just left
        gbwtgraph::handle_t curr_tgt_handle = graph.get_handle(
            gbwt::Node::id(tgt_pos.first), gbwt::Node::is_reverse(tgt_pos.first));
        size_t node_len = graph.get_length(curr_tgt_handle);
        tgt_base_offset += node_len;
        tgt_node_offset++;
        
        if (debug) {
            cerr << "    Left GBWT node: " << tgt_pos.first 
                 << " (id=" << gbwt::Node::id(tgt_pos.first) 
                 << ", rev=" << gbwt::Node::is_reverse(tgt_pos.first) << ")" << endl;
            cerr << "    Moved forward: added node_len=" << node_len 
                 << ", new tgt_base_offset=" << tgt_base_offset << endl;
            cerr << "    Left node sequence: " << graph.get_sequence(curr_tgt_handle) << endl;
        }
        
        // Get tag code for next node
        gbwt::node_type next_tgt_node = next_tgt_pos.first;
        int64_t next_tgt_node_id = gbwt::Node::id(next_tgt_node);
        bool next_tgt_is_rev = gbwt::Node::is_reverse(next_tgt_node);
        uint64_t tgt_decoded = ((static_cast<uint64_t>(next_tgt_node_id) - 1) << 1) | (next_tgt_is_rev ? 1 : 0);
        uint64_t next_tgt_tag_code = tgt_decoded + 1;
        
        // Store this node first (we need it even if it's the last one)
        PathNode tgt_node;
        tgt_node.tag_code = next_tgt_tag_code;
        tgt_node.node = next_tgt_node;
        tgt_node.node_offset = tgt_node_offset;
        tgt_node.base_offset = tgt_base_offset;
        tgt_node.gbwt_position = next_tgt_pos;
        target_path.push_back(tgt_node);
        
        if (debug) {
            gbwtgraph::handle_t next_tgt_handle = graph.get_handle(
                gbwt::Node::id(next_tgt_node), gbwt::Node::is_reverse(next_tgt_node));
            cerr << "    Next GBWT node: " << next_tgt_node 
                 << " (id=" << next_tgt_node_id 
                 << ", rev=" << next_tgt_is_rev << ")" << endl;
            cerr << "    Added target node: tag_code=" << next_tgt_tag_code 
                 << ", node_offset=" << tgt_node_offset 
                 << ", base_offset=" << tgt_base_offset << endl;
            cerr << "    Added node sequence: " << graph.get_sequence(next_tgt_handle) << endl;
        }
        
        // Check stop condition: stop when we've reached or passed the last common node's base offset
        if (tgt_base_offset >= last_common_target_base) {
            if (debug) {
                cerr << "    Target path stopped: reached last common base offset " 
                     << last_common_target_base << " (current: " << tgt_base_offset << ")"
                     << " after " << target_path.size() << " nodes" << endl;
            }
            break;
        }
        
        tgt_pos = next_tgt_pos;
    }
    
    if (debug) {
        cerr << "  Path traversal completed:" << endl;
        cerr << "    Source path: " << source_path.size() << " nodes" << endl;
        cerr << "    Target path: " << target_path.size() << " nodes" << endl;
    }
    
    // ========== Find common nodes and map positions ==========
    // Create a map from tag_code to target path nodes for efficient lookup
    unordered_map<uint64_t, vector<size_t>> target_nodes_by_tag;
    for (size_t i = 0; i < target_path.size(); i++) {
        target_nodes_by_tag[target_path[i].tag_code].push_back(i);
    }
    
    if (debug) {
        cerr << "  Finding common nodes and mapping positions..." << endl;
    }
    
    // Iterate through source path. For each common node, emit per-base mappings for every base
    // of the node that lies within [source_start, source_end]. Source and target traverse the
    // same graph node, so within-node offsets are preserved 1:1.
    for (const auto& src_node : source_path) {
        size_t src_node_len = graph.get_length(
            graph.get_handle(gbwt::Node::id(src_node.node), gbwt::Node::is_reverse(src_node.node)));
        size_t src_node_end = src_node.base_offset + src_node_len;  // exclusive

        // Skip nodes whose coverage doesn't overlap [source_start, source_end].
        if (src_node_end <= source_start || src_node.base_offset > source_end) {
            continue;
        }

        auto it = target_nodes_by_tag.find(src_node.tag_code);
        if (it == target_nodes_by_tag.end()) continue;

        if (it->second.size() > 1 && debug) {
            cerr << "    Common node with MULTIPLE target occurrences: tag_code=" << src_node.tag_code
                 << " (" << it->second.size() << " occurrences)" << endl;
        }

        // Compute the slice of the source node that overlaps the query interval.
        size_t overlap_start = std::max(src_node.base_offset, source_start);
        size_t overlap_end   = std::min(src_node_end, source_end + 1);  // exclusive

        for (size_t tgt_idx : it->second) {
            const PathNode& tgt_node = target_path[tgt_idx];
            for (size_t src_off = overlap_start; src_off < overlap_end; ++src_off) {
                size_t within_node = src_off - src_node.base_offset;
                size_t tgt_off = tgt_node.base_offset + within_node;
                offset_map[src_off].push_back(tgt_off);
            }
            if (debug) {
                cerr << "    Common node: tag_code=" << src_node.tag_code << endl;
                cerr << "      Source: base_offset=" << src_node.base_offset
                     << ", node_offset=" << src_node.node_offset
                     << ", node_len=" << src_node_len << endl;
                cerr << "      Target: base_offset=" << tgt_node.base_offset
                     << ", node_offset=" << tgt_node.node_offset << endl;
                cerr << "      -> Mapped per-base over [" << overlap_start << "," << overlap_end << ") "
                     << "(" << (overlap_end - overlap_start) << " bases) to target ["
                     << (tgt_node.base_offset + (overlap_start - src_node.base_offset)) << ","
                     << (tgt_node.base_offset + (overlap_end - src_node.base_offset)) << ")" << endl;
            }
        }
    }
    
    if (debug) {
        cerr << "  Mapping complete: " << offset_map.size() << " positions mapped" << endl;
    }
    
    return offset_map;
}

// Legacy wrapper using panindexer::FastLocate (for compatibility)
// TODO: Remove this once all callers are updated to use GBWT directly
unordered_map<size_t, vector<size_t>> build_offset_mapping_with_gbwt_rindex(
    FastLocate& rlbwt_rindex, FastLocate& gbwt_rindex,
    size_t source_seq_id, size_t source_start, size_t source_end,
    size_t target_seq_id, size_t anchor_source_offset, size_t anchor_target_offset,
    uint64_t anchor_tag_code) {
    
    cerr << "Warning: Using legacy build_offset_mapping_with_gbwt_rindex. "
         << "Please update to use build_offset_mapping_with_gbwt_lf with GBWT objects." << endl;
    
    // Return just the anchor point mapping as fallback
    unordered_map<size_t, vector<size_t>> offset_map;
    offset_map[anchor_source_offset].push_back(anchor_target_offset);
    return offset_map;
}

// Trace coordinates forward from anchor point using GBWT LF mapping
vector<TranslationResult> trace_coordinates_gbwt(
    const gbwt::GBWT& gbwt_index,
    const gbwt::FastLocate& gbwt_fast_locate,
    const gbwtgraph::GBWTGraph& graph,
    size_t source_seq_id, size_t source_start, size_t source_end,
    size_t target_seq_id, size_t anchor_source_offset, size_t anchor_target_offset,
    size_t anchor_source_base, size_t anchor_target_base,
    uint64_t anchor_tag_code, 
    size_t last_common_source_base, size_t last_common_target_base,
    uint64_t last_common_tag_code) {
    
    vector<TranslationResult> translations;
    
    if (debug) {
        cerr << "Tracing coordinates from anchor using GBWT LF mapping:" << endl;
        cerr << "  Anchor source offset: " << anchor_source_offset << endl;
        cerr << "  Anchor target offset: " << anchor_target_offset << endl;
        cerr << "  Last common source base: " << last_common_source_base << endl;
        cerr << "  Last common target base: " << last_common_target_base << endl;
    }
    
    // Build offset mapping using GBWT LF-based tracing
    unordered_map<size_t, vector<size_t>> offset_map = build_offset_mapping_with_gbwt_lf(
        gbwt_index, gbwt_fast_locate, graph, source_seq_id, source_start, source_end,
        target_seq_id, anchor_source_offset, anchor_target_offset,
        anchor_source_base, anchor_target_base,
        anchor_tag_code, last_common_source_base, last_common_target_base, last_common_tag_code);
    
    if (debug) {
        cerr << "\n=== Final Offset Mapping ===" << endl;
        cerr << "Total source positions mapped: " << offset_map.size() << endl;
        for (const auto& [src, tgt_vec] : offset_map) {
            cerr << "  " << src << " -> [ ";
            for (size_t i = 0; i < tgt_vec.size(); i++) {
                if (i > 0) cerr << ", ";
                cerr << tgt_vec[i];
            }
            cerr << " ]" << endl;
        }
        cerr << "========================\n" << endl;
    }
    
    // Convert to translation results (create one result per source->target mapping)
    for (size_t src_off = source_start; src_off <= source_end; src_off++) {
        if (offset_map.find(src_off) != offset_map.end()) {
            // For each source offset, create a result for EACH target offset
            const vector<size_t>& target_offsets = offset_map[src_off];
            for (size_t tgt_off : target_offsets) {
                TranslationResult result;
                result.source_offset = src_off;
                result.target_seq_id = target_seq_id;
                result.target_offset = tgt_off;
                result.tag_code = anchor_tag_code;
                translations.push_back(result);
            }
        } else {
            // No mapping found for this source offset
            TranslationResult result;
            result.source_offset = src_off;
            result.target_seq_id = target_seq_id;
            result.target_offset = 0;
            result.tag_code = anchor_tag_code;
            translations.push_back(result);
            
            if (debug) {
                cerr << "Warning: Source offset " << src_off << " not mapped to target" << endl;
            }
        }
    }
    
    return translations;
}

// Trace coordinates forward and backward from anchor point using legacy r-index
vector<TranslationResult> trace_coordinates(
    FastLocate& rlbwt_rindex, FastLocate& gbwt_rindex,
    size_t source_seq_id, size_t source_start, size_t source_end,
    size_t target_seq_id, size_t anchor_source_offset, size_t anchor_target_offset,
    uint64_t anchor_tag_code) {
    
    vector<TranslationResult> translations;
    
    if (debug) {
        cerr << "Tracing coordinates from anchor using legacy r-index:" << endl;
        cerr << "  Anchor source offset: " << anchor_source_offset << endl;
        cerr << "  Anchor target offset: " << anchor_target_offset << endl;
    }
    
    // Build offset mapping using legacy r-index tracing
    unordered_map<size_t, vector<size_t>> offset_map = build_offset_mapping_with_gbwt_rindex(
        rlbwt_rindex, gbwt_rindex, source_seq_id, source_start, source_end,
        target_seq_id, anchor_source_offset, anchor_target_offset,
        anchor_tag_code);
    
    // Convert to translation results (create one result per source->target mapping)
    for (size_t src_off = source_start; src_off <= source_end; src_off++) {
        if (offset_map.find(src_off) != offset_map.end()) {
            // For each source offset, create a result for EACH target offset
            const vector<size_t>& target_offsets = offset_map[src_off];
            for (size_t tgt_off : target_offsets) {
                TranslationResult result;
                result.source_offset = src_off;
                result.target_seq_id = target_seq_id;
                result.target_offset = tgt_off;
                result.tag_code = anchor_tag_code;
                translations.push_back(result);
            }
        } else {
            // No mapping found for this source offset
            TranslationResult result;
            result.source_offset = src_off;
            result.target_seq_id = target_seq_id;
            result.target_offset = 0;
            result.tag_code = anchor_tag_code;
            translations.push_back(result);
        }
    }
    
    return translations;
}

#ifndef COORDINATE_TRANSLATION_NO_MAIN
int main(int argc, char** argv) {
    if (argc < 4) {
        usage(argv[0]);
        return 1;
    }
    
    string rlbwt_rindex_file = argv[1];
    string gbwt_rindex_file = argv[2];
    string sampled_tags_file = argv[3];
    
    // Parse arguments
    PathSpec source_spec;
    PathSpec target_spec;
    pair<size_t, size_t> interval = {0, 0};
    bool has_interval = false;
    string gbz_file = "";
    bool benchmark_mode = false;
    string benchmark_interval_lengths_str = "";  // e.g. "5000,10000,1000000"
    string table1_file = "";
    string table2_file = "";
    string source_haplotype_name = "";
    string target_haplotype_name = "";

    for (int i = 4; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--interval" && i + 1 < argc) {
            interval = parse_interval(argv[++i]);
            has_interval = true;
        } 
        // Source specification options
        else if (arg == "--source-id" && i + 1 < argc) {
            source_spec.seq_id = stoull(argv[++i]);
        } else if (arg == "--source-path-name" && i + 1 < argc) {
            source_spec.path_name = argv[++i];
            source_spec.has_path_name = true;
        } else if (arg == "--source-sample" && i + 1 < argc) {
            source_spec.sample_name = argv[++i];
            source_spec.has_sample = true;
        } else if (arg == "--source-contig" && i + 1 < argc) {
            source_spec.contig_name = argv[++i];
            source_spec.has_contig = true;
        } else if (arg == "--source-haplotype" && i + 1 < argc) {
            source_spec.haplotype = stoull(argv[++i]);
            source_spec.has_haplotype = true;
        }
        // Target specification options
        else if (arg == "--target-id" && i + 1 < argc) {
            target_spec.seq_id = stoull(argv[++i]);
        } else if (arg == "--target-path-name" && i + 1 < argc) {
            target_spec.path_name = argv[++i];
            target_spec.has_path_name = true;
        } else if (arg == "--target-sample" && i + 1 < argc) {
            target_spec.sample_name = argv[++i];
            target_spec.has_sample = true;
        } else if (arg == "--target-contig" && i + 1 < argc) {
            target_spec.contig_name = argv[++i];
            target_spec.has_contig = true;
        } else if (arg == "--target-haplotype" && i + 1 < argc) {
            target_spec.haplotype = stoull(argv[++i]);
            target_spec.has_haplotype = true;
        } else if (arg == "--target-reverse") {
            target_spec.reverse_strand = true;
        }
        // Source strand option
        else if (arg == "--source-reverse") {
            source_spec.reverse_strand = true;
        }
        // Other options
        else if (arg == "--gbz" && i + 1 < argc) {
            gbz_file = argv[++i];
        } else if (arg == "--benchmark") {
            benchmark_mode = true;
        } else if (arg == "--benchmark-intervals" && i + 1 < argc) {
            benchmark_interval_lengths_str = argv[++i];
        } else if (arg == "--table1" && i + 1 < argc) {
            table1_file = argv[++i];
        } else if (arg == "--table2" && i + 1 < argc) {
            table2_file = argv[++i];
        } else if (arg == "--source-haplotype-name" && i + 1 < argc) {
            source_haplotype_name = argv[++i];
        } else if (arg == "--target-haplotype-name" && i + 1 < argc) {
            target_haplotype_name = argv[++i];
        } else if (arg == "--debug") {
            debug = true;
        } else if (arg == "--no-debug") {
            debug = false;
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return 0;
        } else {
            cerr << "Unknown argument: " << arg << endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    bool table_mode = (!table1_file.empty() && !table2_file.empty() &&
                       !source_haplotype_name.empty() && !target_haplotype_name.empty());

    // Validate source/target specification (not required in table-driven or single-path benchmark mode)
    if (!table_mode && !benchmark_mode) {
        if (!source_spec.validate("source")) {
            usage(argv[0]);
            return 1;
        }
        if (!source_spec.is_specified()) {
            cerr << "Error: Source haplotype must be specified (--source-id, --source-path-name, or --source-sample/--source-contig)"
                 << " or use table mode (--table1 --table2 --source-haplotype-name --target-haplotype-name)" << endl;
            usage(argv[0]);
            return 1;
        }
        if (!target_spec.validate("target")) {
            usage(argv[0]);
            return 1;
        }
        if (!target_spec.is_specified()) {
            cerr << "Error: Target haplotype must be specified (--target-id, --target-path-name, or --target-sample/--target-contig)"
                 << " or use table mode (--table1 --table2 --source-haplotype-name --target-haplotype-name)" << endl;
            usage(argv[0]);
            return 1;
        }
    }
    
    if (!has_interval && !benchmark_mode) {
        cerr << "Error: --interval is required (unless --benchmark)" << endl;
        usage(argv[0]);
        return 1;
    }
    
    if (gbz_file.empty()) {
        cerr << "Error: --gbz is required (contains GBWT index + GBWTGraph)" << endl;
        usage(argv[0]);
        return 1;
    }
    
    size_t seq_start = interval.first;
    size_t seq_end = interval.second;
    
    auto start_total = high_resolution_clock::now();
    (void)start_total;  // used for total-including-load in debug; after load we use start_after_load
    
    // Load RLBWT r-index
    FastLocate rlbwt_rindex;
    {
        if (debug) cerr << "Loading RLBWT r-index: " << rlbwt_rindex_file << "..." << endl;
        ifstream rin(rlbwt_rindex_file, ios::binary);
        if (!rin) { 
            cerr << "Cannot open RLBWT r-index: " << rlbwt_rindex_file << endl; 
            return 1; 
        }
        rlbwt_rindex.load_encoded(rin);
        rlbwt_rindex.ensure_last_rank();
        rlbwt_rindex.ensure_last_select();
    }
    
    // Load GBZ (contains GBWT index + GBWTGraph) - required for LF mapping and node lengths
    if (debug) cerr << "Loading GBZ (GBWT + graph): " << gbz_file << "..." << endl;
    gbwtgraph::GBZ gbz;
    {
        try {
            sdsl::simple_sds::load_from(gbz, gbz_file);
            if (debug) cerr << "  GBZ loaded. GBWT sequences: " << gbz.index.sequences() << endl;
        } catch (const std::exception& e) {
            cerr << "Error loading GBZ file: " << e.what() << endl;
            return 1;
        }
    }
    const gbwt::GBWT* gbwt_index_ptr = &gbz.index;
    const gbwtgraph::GBWTGraph& graph = gbz.graph;
    
    // Load GBWT FastLocate (r-index) from .ri file
    if (debug) cerr << "Loading GBWT r-index (FastLocate): " << gbwt_rindex_file << "..." << endl;
    std::unique_ptr<gbwt::FastLocate> gbwt_rindex_ptr;
    {
        ifstream gin(gbwt_rindex_file, ios::binary);
        if (!gin) { 
            cerr << "Cannot open GBWT r-index: " << gbwt_rindex_file << endl; 
            return 1; 
        }
        
        // Check file size
        gin.seekg(0, ios::end);
        size_t file_size = gin.tellg();
        gin.seekg(0, ios::beg);
        if (debug) cerr << "  File size: " << file_size << " bytes" << endl;
        
        if (file_size == 0) {
            cerr << "Error: GBWT r-index file is empty" << endl;
            return 1;
        }
        
        try {
            gbwt_rindex_ptr = std::make_unique<gbwt::FastLocate>();
            gbwt_rindex_ptr->load(gin);
            if (debug) cerr << "  GBWT r-index loaded successfully" << endl;
        } catch (const std::exception& e) {
            cerr << "Error loading GBWT r-index: " << e.what() << endl;
            return 1;
        } catch (...) {
            cerr << "Unknown error loading GBWT r-index (possibly segfault)" << endl;
            return 1;
        }
    }
    
    // Associate FastLocate with the GBWT index (required for decompressSA and other operations)
    // FastLocate stores a pointer to the GBWT index, which must be set after loading from file
    if (gbwt_index_ptr && gbwt_rindex_ptr) {
        gbwt_rindex_ptr->setGBWT(gbz.index);
        if (debug) cerr << "  FastLocate associated with GBWT index" << endl;
    }
    
    // Resolve source and target path specifications to sequence IDs (skip in table-driven mode)
    if (!table_mode && !benchmark_mode &&
        (!source_spec.use_seq_id() || !target_spec.use_seq_id())) {
        if (!gbwt_index_ptr) {
            cerr << "Error: GBWT index is required for path name resolution" << endl;
            return 1;
        }
        
        if (!source_spec.resolve(*gbwt_index_ptr, "source", debug)) {
            return 1;
        }
        
        if (!target_spec.resolve(*gbwt_index_ptr, "target", debug)) {
            return 1;
        }
    }
    
    size_t source_seq_id;
    size_t target_seq_id;
    if (table_mode) {
        source_seq_id = 0;
        target_seq_id = 0;
    } else if (benchmark_mode) {
        source_seq_id = kBenchmarkSourceSeqId;
        target_seq_id = kBenchmarkTargetSeqId;
    } else {
        source_seq_id = source_spec.seq_id;
        target_seq_id = target_spec.seq_id;
    }
    
    if (debug) {
        cerr << "Resolved sequence IDs:" << endl;
        cerr << "  Source: " << source_seq_id;
        if (source_spec.use_path_name()) {
            cerr << " (from path name '" << source_spec.path_name << "')";
        } else if (source_spec.use_metadata()) {
            cerr << " (from sample='" << source_spec.sample_name 
                 << "', contig='" << source_spec.contig_name 
                 << "', haplotype=" << source_spec.haplotype << ")";
        }
        cerr << endl;
        
        cerr << "  Target: " << target_seq_id;
        if (target_spec.use_path_name()) {
            cerr << " (from path name '" << target_spec.path_name << "')";
        } else if (target_spec.use_metadata()) {
            cerr << " (from sample='" << target_spec.sample_name 
                 << "', contig='" << target_spec.contig_name 
                 << "', haplotype=" << target_spec.haplotype << ")";
        }
        cerr << endl;
    }
    
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
        sampled.ensure_run_rank();
        sampled.ensure_run_select();
    }
    
    if (debug) {
        cerr << "\n=== All files loaded successfully ===" << endl;
        cerr << "RLBWT R-index: loaded" << endl;
        if (gbwt_index_ptr && gbwt_rindex_ptr) {
            cerr << "GBWT index: loaded (" << gbwt_index_ptr->sequences() << " sequences)" << endl;
            cerr << "GBWT FastLocate: loaded" << endl;
        }
        cerr << "Sampled tag array: loaded" << endl;
        cerr << "=====================================\n" << endl;
    }
    
    auto start_after_load = high_resolution_clock::now();

    // Table-driven mode (single run): use Table 1 + Table 2, then translate and convert to haplotype coords
    if (table_mode && !benchmark_mode) {
        TranslationTable1 t1;
        TranslationTable2 t2;
        {
            ifstream t1in(table1_file, ios::binary);
            if (!t1in) { cerr << "Error: Cannot open Table 1 file: " << table1_file << endl; return 1; }
            t1.load(t1in);
            ifstream t2in(table2_file, ios::binary);
            if (!t2in) { cerr << "Error: Cannot open Table 2 file: " << table2_file << endl; return 1; }
            t2.load(t2in);
        }
        unordered_map<size_t, pair<string, size_t>> path_to_global = build_path_id_to_global(t1);
        vector<string> source_path_names = path_names_for_haplotype(t1, source_haplotype_name);
        // Optional: further restrict to a specific contig (e.g. chr10) if provided via --source-contig.
        if (source_spec.has_contig) {
            string suffix = "#" + source_spec.contig_name;
            vector<string> filtered;
            for (const string& name : source_path_names) {
                if (name.size() >= suffix.size() &&
                    name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
                    filtered.push_back(name);
                }
            }
            source_path_names.swap(filtered);
        }
        if (source_path_names.empty()) {
            cerr << "Error: No path names in Table 1 for source haplotype '" << source_haplotype_name << "'";
            if (source_spec.has_contig) {
                cerr << " (contig='" << source_spec.contig_name << "')";
            }
            cerr << endl;
            return 1;
        }
        // Restrict to the input interval only: Table 1 returns (path_id, local_interval) for the exact
        // overlap of [global_start, global_end] with each subpath; Table 2 and translation use that same interval.
        size_t global_start = interval.first;
        size_t global_end = interval.second;  // user interval is [start, end] inclusive
        vector<PathInterval> source_intervals;
        for (const string& name : source_path_names) {
            vector<PathInterval> pis = t1.lookup(name, global_start, global_end + 1);  // Table 1 uses [start, end) exclusive
            for (PathInterval& pi : pis) {
                source_intervals.push_back(pi);
            }
        }
        if (source_intervals.empty()) {
            cerr << "Error: No path intervals in Table 1 for source haplotype '" << source_haplotype_name
                 << "' and interval [" << global_start << ".." << global_end << "]" << endl;
            return 1;
        }
        struct HaplotypeTranslation {
            size_t source_haplotype_offset = 0;
            size_t target_haplotype_offset = 0;
            size_t target_path_id = 0;
        };
        vector<HaplotypeTranslation> all_results;
        for (const PathInterval& pi : source_intervals) {
            size_t src_path_id = pi.path_id;
            size_t local_start = pi.start;
            size_t local_end = pi.end;  // exclusive
            if (local_end <= local_start) continue;
            size_t seq_start_incl = local_start;
            size_t seq_end_incl = local_end - 1;
            vector<TargetInterval> tgt_intervals = t2.lookup(src_path_id, target_haplotype_name, local_start, local_end);
            cerr << "Table2 (single-run): src_path_id=" << src_path_id
                 << " local_interval=[" << local_start << "," << local_end << ")"
                 << " target_haplotype=\"" << target_haplotype_name << "\" -> "
                 << tgt_intervals.size() << " target path(s)";
            if (!tgt_intervals.empty()) {
                cerr << " {";
                size_t max_print = std::min<size_t>(tgt_intervals.size(), 8);
                for (size_t i = 0; i < max_print; i++) {
                    if (i > 0) cerr << ",";
                    cerr << tgt_intervals[i].tgt_path_id;
                }
                if (tgt_intervals.size() > max_print) cerr << ",...";
                cerr << "}";
            }
            cerr << endl;
            if (tgt_intervals.empty()) continue;
            size_t src_seq_id = 2 * src_path_id;  // forward strand
            auto it_src = path_to_global.find(src_path_id);
            if (it_src == path_to_global.end()) continue;
            size_t src_subpath_start = it_src->second.second;
            // Filter by contig (if any) and deduplicate by tgt_path_id (Table 2 returns one per segment).
            unordered_set<size_t> distinct_tgt_paths;
            for (const TargetInterval& ti : tgt_intervals) {
                size_t tgt_path_id = ti.tgt_path_id;
                if (target_spec.has_contig) {
                    auto it_tgt = path_to_global.find(tgt_path_id);
                    if (it_tgt == path_to_global.end()) continue;
                    const string& tgt_name = it_tgt->second.first;
                    string suffix = "#" + target_spec.contig_name;
                    if (tgt_name.size() < suffix.size() ||
                        tgt_name.compare(tgt_name.size() - suffix.size(), suffix.size(), suffix) != 0) {
                        continue;
                    }
                }
                distinct_tgt_paths.insert(tgt_path_id);
            }
            // Compute Table 2 extent per target (min..max of overlapping segments); use it to limit find_tags and trace.
            vector<IntervalMapping> segs_single = t2.segments(src_path_id, target_haplotype_name);
            for (size_t tgt_path_id : distinct_tgt_paths) {
                // Extent = [min_start, max_end) of Table 2 segments for this target overlapping the query.
                size_t extent_start = local_end, extent_end = local_start;
                for (const IntervalMapping& m : segs_single) {
                    if (m.tgt_path_id != tgt_path_id) continue;
                    size_t overlap_start = std::max(local_start, m.src_start);
                    size_t overlap_end = std::min(local_end, m.src_end);
                    if (overlap_start >= overlap_end) continue;
                    extent_start = std::min(extent_start, overlap_start);
                    extent_end = std::max(extent_end, overlap_end);
                }
                if (extent_start >= extent_end) continue;  // no overlapping segment
                size_t extent_start_incl = extent_start;
                size_t extent_end_incl = extent_end - 1;
                size_t tgt_seq_id = 2 * tgt_path_id;
                vector<TagInfo> tags = find_tags_in_interval(rlbwt_rindex, sampled, src_seq_id, extent_start_incl, extent_end_incl,
                                                             gbwt_index_ptr, gbwt_rindex_ptr.get(), &graph, tgt_seq_id);
                if (tags.empty()) continue;
                CommonNodes common = find_first_and_last_common_nodes_gbwt(*gbwt_rindex_ptr, rlbwt_rindex, sampled, tags, src_seq_id, tgt_seq_id);
                if (!common.found) continue;
                vector<TranslationResult> trans = trace_coordinates_gbwt(
                    *gbwt_index_ptr, *gbwt_rindex_ptr, graph, src_seq_id, extent_start_incl, extent_end_incl,
                    tgt_seq_id, common.first_source_offset, common.first_target_offset,
                    common.first_source_base, common.first_target_base, common.first_tag_code,
                    common.last_source_base, common.last_target_base, common.last_tag_code);
                auto it_tgt = path_to_global.find(tgt_path_id);
                if (it_tgt == path_to_global.end()) continue;
                size_t tgt_subpath_start = it_tgt->second.second;
                for (const TranslationResult& tr : trans) {
                    if (tr.target_offset == 0) continue;
                    HaplotypeTranslation ht;
                    ht.source_haplotype_offset = src_subpath_start + tr.source_offset;
                    ht.target_haplotype_offset = tgt_subpath_start + tr.target_offset;
                    ht.target_path_id = tgt_path_id;
                    all_results.push_back(ht);
                }
            }
        }
        cout << "Coordinate Translation Results (table-driven, haplotype coordinates):" << endl;
        cout << "Source haplotype: " << source_haplotype_name << endl;
        cout << "Target haplotype: " << target_haplotype_name << endl;
        cout << "Interval (global): " << global_start << ".." << global_end << endl;
        cout << endl;
        cout << "Source_haplotype_offset\tTarget_haplotype_offset\tTarget_path_id" << endl;
        std::sort(all_results.begin(), all_results.end(),
                  [](const HaplotypeTranslation& a, const HaplotypeTranslation& b) {
                      if (a.source_haplotype_offset != b.source_haplotype_offset)
                          return a.source_haplotype_offset < b.source_haplotype_offset;
                      if (a.target_haplotype_offset != b.target_haplotype_offset)
                          return a.target_haplotype_offset < b.target_haplotype_offset;
                      return a.target_path_id < b.target_path_id;
                  });
        for (const HaplotypeTranslation& ht : all_results) {
            cout << ht.source_haplotype_offset << "\t" << ht.target_haplotype_offset << "\t" << ht.target_path_id << endl;
        }
        return 0;
    }

    // Benchmark mode: run translation at gene/megabase scales (excluding load time). Supports table-driven mode.
    if (benchmark_mode) {
        vector<size_t> lengths;
        if (!benchmark_interval_lengths_str.empty()) {
            string s = benchmark_interval_lengths_str;
            for (size_t i = 0; i < s.size(); ) {
                size_t j = s.find(',', i);
                if (j == string::npos) j = s.size();
                string part = s.substr(i, j - i);
                if (!part.empty()) lengths.push_back(stoull(part));
                i = j + (j < s.size() ? 1 : 0);
            }
        }
        if (lengths.empty()) {
            lengths = { 5000, 10000, 1000000 };  // gene scale (5kb, 10kb), megabase (1Mb)
        }
        const int num_runs_per_interval = 20;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> dist(100000, 50000000);

        if (table_mode) {
            // Table-driven benchmark: load tables, resolve random global intervals via Table 1/2, time each translation
            TranslationTable1 t1;
            TranslationTable2 t2;
            {
                ifstream t1in(table1_file, ios::binary);
                if (!t1in) { cerr << "Error: Cannot open Table 1 file: " << table1_file << endl; return 1; }
                t1.load(t1in);
                ifstream t2in(table2_file, ios::binary);
                if (!t2in) { cerr << "Error: Cannot open Table 2 file: " << table2_file << endl; return 1; }
                t2.load(t2in);
            }
            unordered_map<size_t, pair<string, size_t>> path_to_global = build_path_id_to_global(t1);
            vector<string> source_path_names = path_names_for_haplotype(t1, source_haplotype_name);
            // Optional contig restriction (e.g. chr10) from --source-contig.
            if (source_spec.has_contig) {
                string suffix = "#" + source_spec.contig_name;
                vector<string> filtered;
                for (const string& name : source_path_names) {
                    if (name.size() >= suffix.size() &&
                        name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
                        filtered.push_back(name);
                    }
                }
                source_path_names.swap(filtered);
            }
            if (source_path_names.empty()) {
                cerr << "Error: No path names in Table 1 for source haplotype '" << source_haplotype_name << "'";
                if (source_spec.has_contig) {
                    cerr << " (contig='" << source_spec.contig_name << "')";
                }
                cerr << endl;
                return 1;
            }
            cout << "scale\tinterval_len\trun\tnum_source_intervals\tnum_source_path_ids\tnum_translations\tfind_tags_ms\tfind_common_ms\ttrace_gbwt_ms\ttotal_ms\tno_mapping_reason" << endl;
            for (size_t L : lengths) {
                for (int run_id = 0; run_id < num_runs_per_interval; run_id++) {
                    size_t global_start = dist(gen);
                    size_t global_end = global_start + L;
                    vector<PathInterval> source_intervals;
                    for (const string& name : source_path_names) {
                        vector<PathInterval> pis = t1.lookup(name, global_start, global_end + 1);
                        for (PathInterval& pi : pis) {
                            source_intervals.push_back(pi);
                        }
                    }
                    if (source_intervals.empty()) {
                        cerr << "Benchmark (table): SKIP L=" << L << " run=" << (run_id + 1)
                             << " global=[" << global_start << "," << global_end << "): no mapping — Table 1 returned 0 path intervals for source haplotype '"
                             << source_haplotype_name << "' (interval does not overlap any subpath). Using Table 1." << endl;
                        continue;
                    }
                    int64_t total_find_ms = 0, total_common_ms = 0, total_trace_ms = 0;
                    int num_translations = 0;
                    int num_src_with_t2_targets = 0;
                    int num_t2_targets_total = 0;
                    int num_t2_targets_after_contig = 0;
                    int num_skipped_no_tags = 0, num_skipped_no_common = 0;
                    // Count distinct source path_ids we will potentially consider
                    std::unordered_set<size_t> distinct_src_paths;
                    for (const PathInterval& pi : source_intervals) {
                        distinct_src_paths.insert(pi.path_id);
                    }
                    cerr << "Benchmark (table): L=" << L << " run=" << (run_id + 1)
                         << " global=[" << global_start << "," << global_end << ")"
                         << " -> " << source_intervals.size() << " source interval(s) on "
                         << distinct_src_paths.size() << " distinct source path_id(s) from Table 1." << endl;
                    for (const PathInterval& pi : source_intervals) {
                        size_t src_path_id = pi.path_id;
                        size_t local_start = pi.start;
                        size_t local_end = pi.end;
                        if (local_end <= local_start) continue;
                        size_t seq_start_incl = local_start;
                        size_t seq_end_incl = local_end - 1;
                        vector<TargetInterval> tgt_intervals = t2.lookup(src_path_id, target_haplotype_name, local_start, local_end);
                        cerr << "Table2 (benchmark): src_path_id=" << src_path_id
                             << " local_interval=[" << local_start << "," << local_end << ")"
                             << " target_haplotype=\"" << target_haplotype_name << "\" -> "
                             << tgt_intervals.size() << " target path(s) (before contig filter)" << endl;
                        if (tgt_intervals.empty()) continue;
                        num_src_with_t2_targets++;
                        num_t2_targets_total += static_cast<int>(tgt_intervals.size());
                        size_t src_seq_id = 2 * src_path_id;
                        // First, fully filter target paths by contig (if any), then deduplicate by tgt_path_id
                        // (Table 2 returns one entry per overlapping segment; many segments can map to the same path).
                        std::unordered_set<size_t> filtered_tgt_path_set;
                        for (const TargetInterval& ti : tgt_intervals) {
                            size_t tgt_path_id = ti.tgt_path_id;
                            if (target_spec.has_contig) {
                                auto it_tgt_name = path_to_global.find(tgt_path_id);
                                if (it_tgt_name == path_to_global.end()) continue;
                                const string& tgt_name = it_tgt_name->second.first;
                                string suffix = "#" + target_spec.contig_name;
                                if (tgt_name.size() < suffix.size() ||
                                    tgt_name.compare(tgt_name.size() - suffix.size(), suffix.size(), suffix) != 0) {
                                    continue;
                                }
                            }
                            filtered_tgt_path_set.insert(tgt_path_id);
                        }
                        std::vector<size_t> filtered_tgt_paths(filtered_tgt_path_set.begin(), filtered_tgt_path_set.end());
                        if (target_spec.has_contig) {
                            num_t2_targets_after_contig += static_cast<int>(filtered_tgt_paths.size());
                            cerr << "Table2 (benchmark): src_path_id=" << src_path_id
                                 << " local_interval=[" << local_start << "," << local_end << ")"
                                 << " target_haplotype=\"" << target_haplotype_name << "\" -> "
                                 << filtered_tgt_paths.size() << " distinct target path(s) (after contig='"
                                 << target_spec.contig_name << "' filter, deduplicated). Names: {";
                            bool first = true;
                            for (size_t tgt_path_id : filtered_tgt_paths) {
                                auto it_tgt_name = path_to_global.find(tgt_path_id);
                                if (it_tgt_name == path_to_global.end()) continue;
                                if (!first) cerr << ", ";
                                cerr << it_tgt_name->second.first << " (id=" << tgt_path_id << ")";
                                first = false;
                            }
                            cerr << "}" << endl;
                        }
                        // Table 2 segments: compute extent per target; run find_tags/trace on extent only.
                        vector<IntervalMapping> segs_bench = t2.segments(src_path_id, target_haplotype_name);
                        // Run translation once per distinct target path, using Table 2 extent only.
                        for (size_t tgt_path_id : filtered_tgt_paths) {
                            // Extent = [min_start, max_end) of Table 2 segments for this target overlapping the query.
                            size_t extent_start = local_end, extent_end = local_start;
                            for (const IntervalMapping& m : segs_bench) {
                                if (m.tgt_path_id != tgt_path_id) continue;
                                size_t overlap_start = std::max(local_start, m.src_start);
                                size_t overlap_end = std::min(local_end, m.src_end);
                                if (overlap_start >= overlap_end) continue;
                                extent_start = std::min(extent_start, overlap_start);
                                extent_end = std::max(extent_end, overlap_end);
                            }
                            if (extent_start >= extent_end) continue;
                            size_t extent_start_incl = extent_start;
                            size_t extent_end_incl = extent_end - 1;
                            size_t tgt_seq_id = 2 * tgt_path_id;
                            auto t0 = high_resolution_clock::now();
                            vector<TagInfo> tags = find_tags_in_interval(rlbwt_rindex, sampled, src_seq_id, extent_start_incl, extent_end_incl,
                                                                         gbwt_index_ptr, gbwt_rindex_ptr.get(), &graph, tgt_seq_id);
                            auto t1 = high_resolution_clock::now();
                            if (tags.empty()) {
                                num_skipped_no_tags++;
                                continue;
                            }
                            CommonNodes common = find_first_and_last_common_nodes_gbwt(*gbwt_rindex_ptr, rlbwt_rindex, sampled, tags, src_seq_id, tgt_seq_id);
                            auto t2 = high_resolution_clock::now();
                            if (!common.found) {
                                num_skipped_no_common++;
                                continue;
                            }
                            (void) trace_coordinates_gbwt(
                                *gbwt_index_ptr, *gbwt_rindex_ptr, graph, src_seq_id, extent_start_incl, extent_end_incl,
                                tgt_seq_id, common.first_source_offset, common.first_target_offset,
                                common.first_source_base, common.first_target_base, common.first_tag_code,
                                common.last_source_base, common.last_target_base, common.last_tag_code);
                            auto t3 = high_resolution_clock::now();
                            total_find_ms += duration_cast<milliseconds>(t1 - t0).count();
                            total_common_ms += duration_cast<milliseconds>(t2 - t1).count();
                            total_trace_ms += duration_cast<milliseconds>(t3 - t2).count();
                            num_translations++;
                        }
                    }
                    if (num_translations == 0) {
                        cerr << "Benchmark (table): SKIP L=" << L << " run=" << (run_id + 1)
                             << " global=[" << global_start << "," << global_end << "): no mapping — Table 1 gave "
                             << source_intervals.size() << " source interval(s); Table 2 gave " << num_t2_targets_total
                             << " target path(s) in " << num_src_with_t2_targets << " interval(s). Of those: "
                             << num_skipped_no_tags << " had no tags (find_tags), " << num_skipped_no_common
                             << " had no common nodes. ";
                        if (target_spec.has_contig) {
                            cerr << "After applying target contig filter ('" << target_spec.contig_name << "'), "
                                 << num_t2_targets_after_contig << " candidate target path(s) remained. ";
                        }
                        cerr << "Tables only indicate which (path,interval)×target to try; "
                             << "coordinate translation (find_tags/common_nodes) is the authority — no mapping possible for this run." << endl;
                        continue;
                    }
                    int64_t total_ms = total_find_ms + total_common_ms + total_trace_ms;
                    string scale = (L >= 1000000) ? "megabase" : "gene";
                    cout << scale << "\t" << L << "\t" << (run_id + 1) << "\t" << source_intervals.size()
                         << "\t" << distinct_src_paths.size() << "\t" << num_translations << "\t"
                         << total_find_ms << "\t" << total_common_ms << "\t" << total_trace_ms << "\t" << total_ms << "\t-"
                         << endl;
                }
            }
        } else {
            // Single-path benchmark (original behaviour); not using translation tables
            cout << "scale\tinterval_len\trun\tnum_tags\tfind_tags_ms\tused_fast_path\tfound_last_common\tlast_common_source_base\tnum_lf_before_phase1\tphase1_lf_count\tinit_skip_lf_ms\tfind_tags_lf_ms\tfind_tags_tag_lookup_ms\tfind_tags_position_lookup_ms\tfind_tags_tag_map_ms\tfind_tags_fast_path_ms\tfind_common_ms\ttrace_gbwt_ms\ttotal_ms\tno_mapping_reason" << endl;
            for (size_t L : lengths) {
                for (int run_id = 0; run_id < num_runs_per_interval; run_id++) {
                    size_t run_start = dist(gen);
                    size_t run_end = run_start + L;
                    FindTagsInIntervalTiming find_timing;
                    auto t0 = high_resolution_clock::now();
                    vector<TagInfo> run_tags = find_tags_in_interval(rlbwt_rindex, sampled, source_seq_id, run_start, run_end,
                                                                      gbwt_index_ptr, gbwt_rindex_ptr.get(), &graph, target_seq_id,
                                                                      &find_timing);
                    auto t1 = high_resolution_clock::now();
                    if (run_tags.empty()) {
                        cerr << "Benchmark (single-path): SKIP L=" << L << " run=" << (run_id + 1)
                             << " interval=[" << run_start << "," << run_end << "): no mapping — find_tags returned 0 tags. "
                             << "Not using translation tables. Interval is path-local on source path (source_seq_id="
                             << source_seq_id << ", target_seq_id=" << target_seq_id << "). "
                             << "If tables report a mapping for this range, they use global coords and different path(s); use --table1/--table2 benchmark to compare." << endl;
                        continue;
                    }
                    CommonNodes run_common;
                    if (gbwt_index_ptr && gbwt_rindex_ptr) {
                        run_common = find_first_and_last_common_nodes_gbwt(*gbwt_rindex_ptr, rlbwt_rindex, sampled, run_tags, source_seq_id, target_seq_id);
                    }
                    auto t2 = high_resolution_clock::now();
                    vector<TranslationResult> run_translations;
                    if (run_common.found && gbwt_index_ptr && gbwt_rindex_ptr) {
                        run_translations = trace_coordinates_gbwt(
                            *gbwt_index_ptr, *gbwt_rindex_ptr, graph, source_seq_id, run_start, run_end,
                            target_seq_id, run_common.first_source_offset, run_common.first_target_offset,
                            run_common.first_source_base, run_common.first_target_base, run_common.first_tag_code,
                            run_common.last_source_base, run_common.last_target_base, run_common.last_tag_code);
                    }
                    auto t3 = high_resolution_clock::now();
                    auto d_find = duration_cast<milliseconds>(t1 - t0).count();
                    auto d_common = duration_cast<milliseconds>(t2 - t1).count();
                    auto d_trace = duration_cast<milliseconds>(t3 - t2).count();
                    auto d_total = duration_cast<milliseconds>(t3 - t0).count();
                    string scale = (L >= 1000000) ? "megabase" : "gene";
                    cout << scale << "\t" << L << "\t" << (run_id + 1) << "\t" << run_tags.size() << "\t"
                         << d_find << "\t"
                         << (find_timing.used_fast_path ? 1 : 0) << "\t"
                         << (find_timing.found_last_common ? 1 : 0) << "\t"
                         << find_timing.last_common_source_base << "\t"
                         << find_timing.num_lf_before_phase1 << "\t" << find_timing.phase1_lf_count << "\t"
                         << fixed << setprecision(3)
                         << find_timing.init_skip_lf_ms << "\t" << find_timing.phase1_lf_ms << "\t"
                         << find_timing.phase1_tag_lookup_ms << "\t" << find_timing.phase1_position_lookup_ms << "\t"
                         << find_timing.phase1_tag_map_ms << "\t" << find_timing.phase1_fast_path_ms << "\t"
                         << setprecision(0)
                         << d_common << "\t" << d_trace << "\t" << d_total << "\t-"
                         << endl;
                }
            }
        }
        return 0;
    }
    
    // Step 1: Find all tags in source interval (using RLBWT r-index)
    auto start_find_tags = high_resolution_clock::now();
    vector<TagInfo> source_tags = find_tags_in_interval(rlbwt_rindex, sampled, source_seq_id, seq_start, seq_end,
                                                        gbwt_index_ptr, gbwt_rindex_ptr.get(), &graph, target_seq_id);
    auto end_find_tags = high_resolution_clock::now();
    auto duration_find_tags = duration_cast<milliseconds>(end_find_tags - start_find_tags);
    
    if (source_tags.empty()) {
        cerr << "No tags found in source interval" << endl;
        return 1;
    }
    
    // Verification: traverse main GBWT index on source sequence to find nodes in the interval
    if (debug && gbwt_index_ptr) {
        cerr << "\n========== Verification: nodes in source interval (from GBWT traversal) ==========" << endl;
        cerr << "Source path seq_id=" << source_seq_id << ", interval [" << seq_start << ", " << seq_end << ")" << endl;
        gbwt::vector_type path = gbwt_index_ptr->extract(gbwt::Path::encode(source_seq_id, false));
        if (path.empty()) {
            cerr << "  Warning: Source path not found in GBWT" << endl;
        } else {
            size_t cumulative_bases = 0;
            size_t node_offset = 0;
            for (gbwt::node_type node : path) {
                if (node == gbwt::ENDMARKER) break;
                int64_t node_id = gbwt::Node::id(node);
                bool node_rev = gbwt::Node::is_reverse(node);
                gbwtgraph::handle_t handle = graph.get_handle(node_id, node_rev);
                size_t node_len = graph.get_length(handle);
                size_t node_base_start = cumulative_bases;
                size_t node_base_end = cumulative_bases + node_len;
                // Node overlaps [seq_start, seq_end) iff [node_base_start, node_base_end) overlaps
                bool in_interval = (node_base_start < seq_end && node_base_end > seq_start);
                if (in_interval) {
                    cerr << "  node_offset=" << node_offset
                         << " node_id=" << node_id
                         << " rev=" << node_rev
                         << " base_range=[" << node_base_start << "," << node_base_end << ")"
                         << endl;
                }
                cumulative_bases += node_len;
                node_offset++;
            }
            cerr << "  (end of path traversal)" << endl;
        }
        cerr << "==================================================================================" << endl << endl;
    }
    
    // Step 2: Find first and last common nodes with target haplotype using GBWT FastLocate
    // Sorts tags by tag_code (smallest to largest) and uses decompressSA() to find matching paths
    auto start_find_common = high_resolution_clock::now();
    CommonNodes common_nodes;
    
    if (gbwt_index_ptr && gbwt_rindex_ptr) {
        common_nodes = find_first_and_last_common_nodes_gbwt(*gbwt_rindex_ptr, rlbwt_rindex, sampled, source_tags, source_seq_id, target_seq_id);
    } else {
        cerr << "Error: GBWT index and FastLocate required for finding common nodes" << endl;
        return 1;
    }
    
    auto end_find_common = high_resolution_clock::now();
    auto duration_find_common = duration_cast<milliseconds>(end_find_common - start_find_common);
    
    if (!common_nodes.found) {
        cerr << "No common nodes found between source and target haplotypes" << endl;
        return 1;
    }
    
    size_t anchor_source_offset = common_nodes.first_source_offset;
    size_t anchor_target_offset = common_nodes.first_target_offset;
    size_t anchor_source_base = common_nodes.first_source_base;
    size_t anchor_target_base = common_nodes.first_target_base;
    uint64_t anchor_tag_code = common_nodes.first_tag_code;
    size_t last_common_source_offset = common_nodes.last_source_offset;
    size_t last_common_target_offset = common_nodes.last_target_offset;
    size_t last_common_source_base = common_nodes.last_source_base;
    size_t last_common_target_base = common_nodes.last_target_base;
    uint64_t last_common_tag_code = common_nodes.last_tag_code;
    
    // ========== DEBUG: Verify anchor positions by traversing from start ==========
    // NOTE: GBWT seqOffset is counted from the END of the path (distance to end),
    // not from the start. We traverse from start and then convert.
    if (debug && gbwt_index_ptr) {
        cerr << "\n========== VERIFICATION: Traversing paths from start ==========" << endl;
        cerr << "NOTE: GBWT seqOffset = distance from END (larger = closer to start)" << endl;
        
        // Decode anchor tag to get node_id
        auto [anchor_node_id, anchor_is_rev] = decode_tag(anchor_tag_code);
        gbwt::node_type anchor_gbwt_node = gbwt::Node::encode(anchor_node_id, anchor_is_rev);
        
        auto [last_node_id, last_is_rev] = decode_tag(last_common_tag_code);
        gbwt::node_type last_gbwt_node = gbwt::Node::encode(last_node_id, last_is_rev);
        
        cerr << "First common node (anchor): tag_code=" << anchor_tag_code 
             << " (node_id=" << anchor_node_id << ", is_rev=" << anchor_is_rev << ")" << endl;
        cerr << "  Expected source: node_offset_from_end=" << anchor_source_offset << ", base_offset=" << anchor_source_base << endl;
        cerr << "  Expected target: node_offset_from_end=" << anchor_target_offset << ", base_offset=" << anchor_target_base << endl;
        
        cerr << "Last common node: tag_code=" << last_common_tag_code 
             << " (node_id=" << last_node_id << ", is_rev=" << last_is_rev << ")" << endl;
        cerr << "  Expected source: node_offset_from_end=" << last_common_source_offset << ", base_offset=" << last_common_source_base << endl;
        cerr << "  Expected target: node_offset_from_end=" << last_common_target_offset << ", base_offset=" << last_common_target_base << endl;
        
        // Traverse SOURCE path from start
        cerr << "\n--- Traversing SOURCE path (seq_id=" << source_seq_id << ") from start ---" << endl;
        gbwt::edge_type src_pos = gbwt_index_ptr->start(source_seq_id);
        size_t src_node_idx = 0;
        size_t src_base_offset_calc = 0;
        bool found_anchor_in_source = false;
        bool found_last_in_source = false;
        size_t anchor_src_node_idx_calc = 0;
        size_t anchor_src_base_calc = 0;
        size_t last_src_node_idx_calc = 0;
        size_t last_src_base_calc = 0;
        
        while (src_pos.first != gbwt::ENDMARKER) {
            gbwt::node_type curr_node = src_pos.first;
            int64_t curr_id = gbwt::Node::id(curr_node);
            bool curr_rev = gbwt::Node::is_reverse(curr_node);
            
            gbwtgraph::handle_t curr_handle = graph.get_handle(curr_id, curr_rev);
            size_t curr_len = graph.get_length(curr_handle);
            
            // Check if this is the anchor node
            if (curr_node == anchor_gbwt_node && !found_anchor_in_source) {
                found_anchor_in_source = true;
                anchor_src_node_idx_calc = src_node_idx;
                anchor_src_base_calc = src_base_offset_calc;
                cerr << "  FOUND ANCHOR at node_idx=" << src_node_idx 
                     << ", base_offset=" << src_base_offset_calc 
                     << " (node_id=" << curr_id << ", len=" << curr_len << ")" << endl;
            }
            
            // Check if this is the last common node
            if (curr_node == last_gbwt_node) {
                found_last_in_source = true;
                last_src_node_idx_calc = src_node_idx;
                last_src_base_calc = src_base_offset_calc;
                cerr << "  FOUND LAST at node_idx=" << src_node_idx 
                     << ", base_offset=" << src_base_offset_calc 
                     << " (node_id=" << curr_id << ", len=" << curr_len << ")" << endl;
            }
            
            // Move to next node
            src_base_offset_calc += curr_len;
            src_node_idx++;
            src_pos = gbwt_index_ptr->LF(src_pos);
            
            // Limit output for very long paths
            if (src_node_idx > 10000000) {
                cerr << "  (stopped after 10M nodes)" << endl;
                break;
            }
        }
        
        size_t src_total_nodes = src_node_idx;
        cerr << "  Source path total: " << src_total_nodes << " nodes, " << src_base_offset_calc << " bases" << endl;
        
        // GBWT seqOffset is from the END of the path, so convert:
        // offset_from_end = (total_nodes - 1) - offset_from_start
        if (found_anchor_in_source) {
            size_t anchor_src_offset_from_end = (src_total_nodes - 1) - anchor_src_node_idx_calc;
            cerr << "  ANCHOR verification:" << endl;
            cerr << "    Node offset from START (calculated): " << anchor_src_node_idx_calc << endl;
            cerr << "    Node offset from END (converted):    " << anchor_src_offset_from_end << endl;
            cerr << "    Expected (GBWT seqOffset from END):  " << anchor_source_offset 
                 << (anchor_source_offset == anchor_src_offset_from_end ? " ✓" : " ✗ MISMATCH!") << endl;
            cerr << "    Base offset: expected=" << anchor_source_base << ", calculated=" << anchor_src_base_calc 
                 << (anchor_source_base == anchor_src_base_calc ? " ✓" : " ✗ MISMATCH!") << endl;
        } else {
            cerr << "  WARNING: Anchor node NOT FOUND in source path!" << endl;
        }
        
        if (found_last_in_source) {
            size_t last_src_offset_from_end = (src_total_nodes - 1) - last_src_node_idx_calc;
            cerr << "  LAST verification:" << endl;
            cerr << "    Node offset from START (calculated): " << last_src_node_idx_calc << endl;
            cerr << "    Node offset from END (converted):    " << last_src_offset_from_end << endl;
            cerr << "    Expected (GBWT seqOffset from END):  " << last_common_source_offset 
                 << (last_common_source_offset == last_src_offset_from_end ? " ✓" : " ✗ MISMATCH!") << endl;
            cerr << "    Base offset: expected=" << last_common_source_base << ", calculated=" << last_src_base_calc 
                 << (last_common_source_base == last_src_base_calc ? " ✓" : " ✗ MISMATCH!") << endl;
        } else {
            cerr << "  WARNING: Last common node NOT FOUND in source path!" << endl;
        }
        
        // ========== Print first 8000 bases of SOURCE path with node details ==========
        cerr << "\n--- SOURCE path: First 8000 bases (node by node) ---" << endl;
        cerr << "Format: [node_idx] node_id (rev?) len=X | base_range: start..end | seq: SEQUENCE" << endl;
        
        gbwt::edge_type src_detail_pos = gbwt_index_ptr->start(source_seq_id);
        size_t src_detail_node_idx = 0;
        size_t src_detail_base_offset = 0;
        
        while (src_detail_pos.first != gbwt::ENDMARKER && src_detail_base_offset < 8000) {
            gbwt::node_type curr_node = src_detail_pos.first;
            int64_t curr_id = gbwt::Node::id(curr_node);
            bool curr_rev = gbwt::Node::is_reverse(curr_node);
            
            gbwtgraph::handle_t curr_handle = graph.get_handle(curr_id, curr_rev);
            size_t curr_len = graph.get_length(curr_handle);
            std::string curr_seq = graph.get_sequence(curr_handle);
            
            // Truncate sequence if too long for display
            std::string display_seq = curr_seq;
            if (display_seq.length() > 50) {
                display_seq = display_seq.substr(0, 25) + "..." + display_seq.substr(display_seq.length() - 22);
            }
            
            // Compute tag_code for this node (reverse of decode_tag)
            // decode_tag: decoded = tag_code - 1; node_id = (decoded >> 1) + 1; is_rev = (decoded & 1)
            // encode: decoded = ((node_id - 1) << 1) | (is_rev ? 1 : 0); tag_code = decoded + 1
            uint64_t tag_code = (((curr_id - 1) << 1) | (curr_rev ? 1 : 0)) + 1;
            
            cerr << "[" << src_detail_node_idx << "] node_id=" << curr_id 
                 << (curr_rev ? " (rev)" : "") 
                 << " len=" << curr_len 
                 << " | base_range: " << src_detail_base_offset << ".." << (src_detail_base_offset + curr_len - 1)
                 << " | tag_code=" << tag_code
                 << " | seq: " << display_seq << endl;
            
            // Move to next node
            src_detail_base_offset += curr_len;
            src_detail_node_idx++;
            src_detail_pos = gbwt_index_ptr->LF(src_detail_pos);
        }
        
        cerr << "--- End of first 8000 bases (stopped at base " << src_detail_base_offset << ") ---\n" << endl;
        
        // Traverse TARGET path from start
        cerr << "\n--- Traversing TARGET path (seq_id=" << target_seq_id << ") from start ---" << endl;
        gbwt::edge_type tgt_pos = gbwt_index_ptr->start(target_seq_id);
        size_t tgt_node_idx = 0;
        size_t tgt_base_offset_calc = 0;
        bool found_anchor_in_target = false;
        bool found_last_in_target = false;
        size_t anchor_tgt_node_idx_calc = 0;
        size_t anchor_tgt_base_calc = 0;
        size_t last_tgt_node_idx_calc = 0;
        size_t last_tgt_base_calc = 0;
        
        while (tgt_pos.first != gbwt::ENDMARKER) {
            gbwt::node_type curr_node = tgt_pos.first;
            int64_t curr_id = gbwt::Node::id(curr_node);
            bool curr_rev = gbwt::Node::is_reverse(curr_node);
            
            gbwtgraph::handle_t curr_handle = graph.get_handle(curr_id, curr_rev);
            size_t curr_len = graph.get_length(curr_handle);
            
            // Check if this is the anchor node
            if (curr_node == anchor_gbwt_node && !found_anchor_in_target) {
                found_anchor_in_target = true;
                anchor_tgt_node_idx_calc = tgt_node_idx;
                anchor_tgt_base_calc = tgt_base_offset_calc;
                cerr << "  FOUND ANCHOR at node_idx=" << tgt_node_idx 
                     << ", base_offset=" << tgt_base_offset_calc 
                     << " (node_id=" << curr_id << ", len=" << curr_len << ")" << endl;
            }
            
            // Check if this is the last common node
            if (curr_node == last_gbwt_node) {
                found_last_in_target = true;
                last_tgt_node_idx_calc = tgt_node_idx;
                last_tgt_base_calc = tgt_base_offset_calc;
                cerr << "  FOUND LAST at node_idx=" << tgt_node_idx 
                     << ", base_offset=" << tgt_base_offset_calc 
                     << " (node_id=" << curr_id << ", len=" << curr_len << ")" << endl;
            }
            
            // Move to next node
            tgt_base_offset_calc += curr_len;
            tgt_node_idx++;
            tgt_pos = gbwt_index_ptr->LF(tgt_pos);
            
            // Limit output for very long paths
            if (tgt_node_idx > 10000000) {
                cerr << "  (stopped after 10M nodes)" << endl;
                break;
            }
        }
        
        size_t tgt_total_nodes = tgt_node_idx;
        cerr << "  Target path total: " << tgt_total_nodes << " nodes, " << tgt_base_offset_calc << " bases" << endl;
        
        // GBWT seqOffset is from the END of the path, so convert:
        // offset_from_end = (total_nodes - 1) - offset_from_start
        if (found_anchor_in_target) {
            size_t anchor_tgt_offset_from_end = (tgt_total_nodes - 1) - anchor_tgt_node_idx_calc;
            cerr << "  ANCHOR verification:" << endl;
            cerr << "    Node offset from START (calculated): " << anchor_tgt_node_idx_calc << endl;
            cerr << "    Node offset from END (converted):    " << anchor_tgt_offset_from_end << endl;
            cerr << "    Expected (GBWT seqOffset from END):  " << anchor_target_offset 
                 << (anchor_target_offset == anchor_tgt_offset_from_end ? " ✓" : " ✗ MISMATCH!") << endl;
            cerr << "    Base offset: expected=" << anchor_target_base << ", calculated=" << anchor_tgt_base_calc 
                 << (anchor_target_base == anchor_tgt_base_calc ? " ✓" : " ✗ MISMATCH!") << endl;
        } else {
            cerr << "  WARNING: Anchor node NOT FOUND in target path!" << endl;
        }
        
        if (found_last_in_target) {
            size_t last_tgt_offset_from_end = (tgt_total_nodes - 1) - last_tgt_node_idx_calc;
            cerr << "  LAST verification:" << endl;
            cerr << "    Node offset from START (calculated): " << last_tgt_node_idx_calc << endl;
            cerr << "    Node offset from END (converted):    " << last_tgt_offset_from_end << endl;
            cerr << "    Expected (GBWT seqOffset from END):  " << last_common_target_offset 
                 << (last_common_target_offset == last_tgt_offset_from_end ? " ✓" : " ✗ MISMATCH!") << endl;
            cerr << "    Base offset: expected=" << last_common_target_base << ", calculated=" << last_tgt_base_calc 
                 << (last_common_target_base == last_tgt_base_calc ? " ✓" : " ✗ MISMATCH!") << endl;
        } else {
            cerr << "  WARNING: Last common node NOT FOUND in target path!" << endl;
        }
        
        // ========== Print first 1500 bases of TARGET path with node details ==========
        cerr << "\n--- TARGET path: First 1500 bases (node by node) ---" << endl;
        cerr << "Format: [node_idx] node_id (rev?) len=X | base_range: start..end | seq: SEQUENCE" << endl;
        
        gbwt::edge_type tgt_detail_pos = gbwt_index_ptr->start(target_seq_id);
        size_t detail_node_idx = 0;
        size_t detail_base_offset = 0;
        
        while (tgt_detail_pos.first != gbwt::ENDMARKER && detail_base_offset < 1500) {
            gbwt::node_type curr_node = tgt_detail_pos.first;
            int64_t curr_id = gbwt::Node::id(curr_node);
            bool curr_rev = gbwt::Node::is_reverse(curr_node);
            
            gbwtgraph::handle_t curr_handle = graph.get_handle(curr_id, curr_rev);
            size_t curr_len = graph.get_length(curr_handle);
            std::string curr_seq = graph.get_sequence(curr_handle);
            
            // Truncate sequence if too long for display
            std::string display_seq = curr_seq;
            if (display_seq.length() > 50) {
                display_seq = display_seq.substr(0, 25) + "..." + display_seq.substr(display_seq.length() - 22);
            }
            
            cerr << "[" << detail_node_idx << "] node_id=" << curr_id 
                 << (curr_rev ? " (rev)" : "") 
                 << " len=" << curr_len 
                 << " | base_range: " << detail_base_offset << ".." << (detail_base_offset + curr_len - 1)
                 << " | seq: " << display_seq << endl;
            
            // Move to next node
            detail_base_offset += curr_len;
            detail_node_idx++;
            tgt_detail_pos = gbwt_index_ptr->LF(tgt_detail_pos);
        }
        
        cerr << "--- End of first 1500 bases (stopped at base " << detail_base_offset << ") ---\n" << endl;
        
        // ========== CRITICAL: Show nodes ACTUALLY in the query interval ==========
        cerr << "\n--- NODES ACTUALLY IN SOURCE INTERVAL [" << seq_start << ", " << seq_end << "] ---" << endl;
        cerr << "These are the TAGs that SHOULD have been found by find_tags_in_interval:" << endl;
        
        gbwt::edge_type interval_pos = gbwt_index_ptr->start(source_seq_id);
        size_t interval_node_idx = 0;
        size_t interval_base_offset = 0;
        std::vector<uint64_t> actual_tags_in_interval;
        
        while (interval_pos.first != gbwt::ENDMARKER) {
            gbwt::node_type curr_node = interval_pos.first;
            int64_t curr_id = gbwt::Node::id(curr_node);
            bool curr_rev = gbwt::Node::is_reverse(curr_node);
            
            gbwtgraph::handle_t curr_handle = graph.get_handle(curr_id, curr_rev);
            size_t curr_len = graph.get_length(curr_handle);
            
            // Check if this node overlaps with the query interval
            size_t node_end = interval_base_offset + curr_len - 1;
            if (interval_base_offset <= seq_end && node_end >= seq_start) {
                // This node is in the interval!
                uint64_t tag_code = ((static_cast<uint64_t>(curr_id) - 1) << 1) | (curr_rev ? 1 : 0);
                tag_code += 1;  // +1 because 0 is reserved for gaps
                
                actual_tags_in_interval.push_back(tag_code);
                
                if (actual_tags_in_interval.size() <= 20) {  // Print first 20
                    cerr << "  [" << interval_node_idx << "] node_id=" << curr_id 
                         << (curr_rev ? " (rev)" : "")
                         << " | base_range: " << interval_base_offset << ".." << node_end
                         << " | tag_code=" << tag_code << endl;
                }
            }
            
            // Stop if we've passed the interval
            if (interval_base_offset > seq_end) {
                break;
            }
            
            interval_base_offset += curr_len;
            interval_node_idx++;
            interval_pos = gbwt_index_ptr->LF(interval_pos);
        }
        
        cerr << "Total nodes in interval: " << actual_tags_in_interval.size() << endl;
        
        // Compare with TAGs found by find_tags_in_interval
        cerr << "\n--- COMPARISON: TAGs found vs TAGs that should exist ---" << endl;
        cerr << "TAGs found by find_tags_in_interval: " << source_tags.size() << endl;
        for (size_t i = 0; i < std::min(source_tags.size(), (size_t)10); i++) {
            const auto& tag = source_tags[i];
            auto [tid, trev] = decode_tag(tag.tag_code);
            bool found_in_actual = std::find(actual_tags_in_interval.begin(), 
                                              actual_tags_in_interval.end(), 
                                              tag.tag_code) != actual_tags_in_interval.end();
            cerr << "  tag_code=" << tag.tag_code << " (node_id=" << tid << ")"
                 << (found_in_actual ? " ✓ EXISTS in interval" : " ✗ NOT in interval!") << endl;
        }
        
        cerr << "========== END VERIFICATION ==========\n" << endl;
    }
    // ========== END DEBUG ==========
    
    // Step 3: Trace coordinates using GBWT LF mapping
    auto start_trace = high_resolution_clock::now();
    vector<TranslationResult> translations;
    
    if (gbwt_index_ptr && gbwt_rindex_ptr) {
        // Use GBWT LF mapping for forward path traversal
        // FastLocate finds initial positions, GBWT index LF() does the traversal
        // Each path is traversed until it reaches its last common node at the correct base offset
        translations = trace_coordinates_gbwt(
            *gbwt_index_ptr, *gbwt_rindex_ptr, graph, source_seq_id, seq_start, seq_end,
            target_seq_id, anchor_source_offset, anchor_target_offset,
            anchor_source_base, anchor_target_base, anchor_tag_code,
            last_common_source_base, last_common_target_base, last_common_tag_code);
    } else {
        cerr << "Error: Both GBWT index and GBWT r-index are required" << endl;
        return 1;
    }
    
    auto end_trace = high_resolution_clock::now();
    auto duration_trace = duration_cast<milliseconds>(end_trace - start_trace);
    
    // Output results
    cout << "Coordinate Translation Results:" << endl;
    cout << "Source haplotype:" << endl;
    cout << "  Sequence ID: " << source_seq_id << endl;
    if (source_spec.use_path_name()) {
        cout << "  Path name: " << source_spec.path_name << endl;
        cout << "  Strand: " << (source_spec.reverse_strand ? "reverse" : "forward") << endl;
    } else if (source_spec.use_metadata()) {
        cout << "  Sample: " << source_spec.sample_name << endl;
        cout << "  Contig: " << source_spec.contig_name << endl;
        cout << "  Haplotype: " << source_spec.haplotype << endl;
        cout << "  Strand: " << (source_spec.reverse_strand ? "reverse" : "forward") << endl;
    }
    cout << "  Interval: " << seq_start << ".." << seq_end << endl;
    cout << "Target haplotype:" << endl;
    cout << "  Sequence ID: " << target_seq_id << endl;
    if (target_spec.use_path_name()) {
        cout << "  Path name: " << target_spec.path_name << endl;
        cout << "  Strand: " << (target_spec.reverse_strand ? "reverse" : "forward") << endl;
    } else if (target_spec.use_metadata()) {
        cout << "  Sample: " << target_spec.sample_name << endl;
        cout << "  Contig: " << target_spec.contig_name << endl;
        cout << "  Haplotype: " << target_spec.haplotype << endl;
        cout << "  Strand: " << (target_spec.reverse_strand ? "reverse" : "forward") << endl;
    }
    cout << "First common node (anchor): tag_code=" << anchor_tag_code << endl;
    cout << "  Source: node_offset=" << anchor_source_offset 
         << ", base_offset=" << anchor_source_base << endl;
    cout << "  Target: node_offset=" << anchor_target_offset 
         << ", base_offset=" << anchor_target_base << endl;
    cout << "Last common node (stop): tag_code=" << last_common_tag_code << endl;
    cout << "  Source: node_offset=" << last_common_source_offset 
         << ", base_offset=" << last_common_source_base << endl;
    cout << "  Target: node_offset=" << last_common_target_offset 
         << ", base_offset=" << last_common_target_base << endl;
    cout << endl;
    
    // Check for fragments
    // auto all_anchors = find_all_common_nodes(rlbwt_rindex, sampled, source_tags, target_seq_id);
    
    // Calculate the actual interval covered by the translation
    // (from anchor_source_base to last_common_source_base, which may be larger than the input interval)
    size_t actual_start = anchor_source_base;
    size_t actual_end = last_common_source_base;
    if (actual_start > actual_end) {
        std::swap(actual_start, actual_end);
    }
    
    cout << "Actual translated interval: " << actual_start << ".." << actual_end << endl;
    cout << endl;
    
    cout << "Translations:" << endl;
    cout << "Source_Offset\tTarget_Seq_ID\tTarget_Offset" << endl;
    for (const auto& trans : translations) {
        if (trans.target_offset != 0) {
            cout << trans.source_offset << "\t" << trans.target_seq_id << "\t" 
                 << trans.target_offset << endl;
        }
    }
    
    auto end_total = high_resolution_clock::now();
    auto duration_total = duration_cast<milliseconds>(end_total - start_total);
    auto duration_excl_load = duration_cast<milliseconds>(end_total - start_after_load);
    
    if (debug) {
        cerr << "\n=== Timing Statistics ===" << endl;
        cerr << "  Find tags in interval: " << duration_find_tags.count() << " ms" << endl;
        cerr << "  Find common node: " << duration_find_common.count() << " ms" << endl;
        cerr << "  Trace coordinates (GBWT): " << duration_trace.count() << " ms" << endl;
        cerr << "  Total (excluding load): " << duration_excl_load.count() << " ms" << endl;
        cerr << "  Total (including load): " << duration_total.count() << " ms" << endl;
    }

    return 0;
}
#endif // COORDINATE_TRANSLATION_NO_MAIN
