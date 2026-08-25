#ifndef PANGENOME_INDEX_TRANSLATION_TABLES_HPP
#define PANGENOME_INDEX_TRANSLATION_TABLES_HPP

/**
 * Translation Table 1: Named path + global interval → (path_id, local interval)
 *
 * Many reference paths in the graph are split into subpaths with different
 * start coordinates and lengths. For example, GRCh38#0#chr10 might appear as:
 *
 *   GRCh38#0#chr10[13060]      length 38515470  →  covers [13060, 38528540)
 *   GRCh38#0#chr10[38573338]   length 994246    →  covers [38573338, 39567584)
 *   GRCh38#0#chr10[42099786]   length 5680582   →  covers [42099786, 47780368)
 *   GRCh38#0#chr10[47881464]   length 85785427  →  covers [47881464, 133666891)
 *
 * "Global" coordinates for that name are these reference intervals. A query
 * like GRCh38#0#chr10[1234567-1300000] is in global coords; 1234567 lies
 * in the first subpath [13060, 38528540), so it maps to path_id of that
 * subpath with local interval [1234567 - 13060, 1300000 - 13060].
 *
 * Table 1 stores, for each named path (e.g. "GRCh38#0#chr10"), the list of
 * subpaths (path_id, subpath_start, length). Lookup(name, start, end) returns
 * all (path_id, local_start, local_end) pairs for subpaths overlapping [start,end].
 *
 * Translation Table 2: (source_path_id, target_haplotype) → sorted interval segment mappings
 *
 * For each (source_path_id, target_haplotype_name) pair that shares at least one
 * graph node, Table 2 stores a sorted list of IntervalMapping entries describing
 * which source intervals map to which target_path_id. Source intervals can overlap
 * and map to different target paths (e.g. src[0,100]->77, src[10,90]->78).
 *
 * Example:
 *   (source_path_id=12, "HG002#1") → [
 *     { src=[0, 38000),   tgt_path_id=77 },
 *     { src=[10, 90),     tgt_path_id=78 },  // overlaps
 *     { src=[45000, 110000), tgt_path_id=78 },
 *   ]
 *
 * Lookup returns tgt_path_id for each segment whose source interval overlaps
 * the query [local_start, local_end).
 */

#include <cstdint>
#include <limits>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iosfwd>

namespace panindexer {

/// One subpath: GBWT path id and its extent in global coordinates [subpath_start, subpath_start + length).
struct SubpathInfo {
    size_t path_id = 0;         ///< GBWT sequence/path id for this subpath
    size_t subpath_start = 0;  ///< Start position in global (reference) coordinates
    size_t length = 0;         ///< Length of this subpath in bases

    size_t end() const { return subpath_start + length; }
};

/// One (path_id, local interval) result from Table 1 lookup. Interval is [start, end) (end exclusive).
struct PathInterval {
    size_t path_id = 0;
    size_t start = 0;  ///< Start offset within the path (local coordinates)
    size_t end = 0;    ///< End offset within the path (exclusive)
};

/**
 * Translation Table 1: converts (named path, global interval) to (path_id, local interval) list.
 *
 * - Keys: named path (e.g. "GRCh38#0#chr10"). Each name has a list of SubpathInfo.
 * - Lookup: given (path_name, global_start, global_end), returns all PathInterval for
 *   subpaths overlapping that range, with coordinates in local path space.
 */
class TranslationTable1 {
public:
    TranslationTable1() = default;

    /// Add a subpath for a named path. Subpaths for the same name should be added in order of subpath_start.
    void add_subpath(const std::string& path_name, size_t path_id, size_t subpath_start, size_t length);

    /// Lookup: (path_name, global_start, global_end) → list of (path_id, local_start, local_end).
    /// Global interval is [global_start, global_end) (end exclusive).
    /// Returns empty if path_name is unknown or interval does not overlap any subpath.
    std::vector<PathInterval> lookup(const std::string& path_name,
                                    size_t global_start,
                                    size_t global_end) const;

    /// Lookup by name id (index into names()). Use when you have already resolved name → id.
    std::vector<PathInterval> lookup_by_name_id(size_t name_id,
                                               size_t global_start,
                                               size_t global_end) const;

    /// Number of named paths (path names) in the table.
    size_t num_names() const { return name_to_subpaths_.size(); }

    /// All path names; order matches name_id (index) used by lookup_by_name_id.
    std::vector<std::string> names() const;

    /// Get subpath list for a path name (for inspection/debug). Empty if unknown.
    std::vector<SubpathInfo> subpaths(const std::string& path_name) const;

    /// Serialize to a binary stream (format: version, then name→subpaths).
    void serialize(std::ostream& out) const;

    /// Load from a binary stream.
    void load(std::istream& in);

private:
    /// Path names in order of first add (used for lookup_by_name_id and stable serialize).
    std::vector<std::string> names_;
    /// path_name → list of subpaths
    std::map<std::string, std::vector<SubpathInfo>> name_to_subpaths_;
};

/**
 * One segment entry in Table 2.
 * States: source local interval [src_start, src_end) on source_path_id maps to
 *         target path tgt_path_id (no target coordinates stored).
 * Source intervals can overlap: e.g. src[0,100]->i and src[10,90]->j.
 * Segments are stored sorted by src_start within each (src_path_id, tgt_haplotype) key.
 */
struct IntervalMapping {
    size_t src_start   = 0;  ///< Start on source path (local, inclusive)
    size_t src_end     = 0;  ///< End on source path   (local, exclusive)
    size_t tgt_path_id = 0;  ///< GBWT path_id of the target subpath
    /// Target interval [tgt_start, tgt_end) on tgt_path_id, path-local.
    /// Populated only by "B2" tables (build_table2_b2); zero and meaningless
    /// when has_target_coords() is false. Always stored ascending, even when
    /// tgt_reverse is set — the flag records the direction of correspondence,
    /// not the storage order.
    size_t tgt_start   = 0;
    size_t tgt_end     = 0;
    /// True when the source traverses this block in the opposite orientation
    /// to the target, i.e. src_start corresponds to tgt_end rather than
    /// tgt_start. Callers that ignore this will mirror inverted blocks.
    bool   tgt_reverse = false;

    size_t src_len() const { return src_end - src_start; }
    size_t tgt_len() const { return tgt_end - tgt_start; }
};

/// Result of a Table 2 lookup. Target coordinates are present only for tables
/// built with target coordinates (see TranslationTable2::has_target_coords);
/// check has_coords before using tgt_start/tgt_end rather than testing them
/// against zero, which is a legitimate offset.
struct TargetInterval {
    size_t tgt_path_id = 0;
    size_t tgt_start   = 0;
    size_t tgt_end     = 0;
    bool   tgt_reverse = false;
    bool   has_coords  = false;
};

/**
 * Translation Table 2: (source_path_id, target_haplotype_name) → sorted IntervalMapping list.
 *
 * Sparse: only pairs that share at least one graph node are stored.
 * Within each key the list is sorted by src_start for binary-search lookup.
 */
class TranslationTable2 {
public:
    TranslationTable2() = default;

    /// Add one interval mapping for (src_path_id, tgt_haplotype).
    /// Mappings need not be added in sorted order; call finalize() when done.
    void add_mapping(size_t src_path_id,
                     const std::string& tgt_haplotype,
                     const IntervalMapping& mapping);

    /// Sort all segment lists by src_start. Must be called before lookup or serialize.
    void finalize();

    /**
     * Lookup: (src_path_id, tgt_haplotype, local_start, local_end) →
     *         list of TargetInterval (tgt_path_id only) for each overlapping segment.
     * Returns empty if no entry or no overlap.
     * finalize() must have been called first.
     */
    std::vector<TargetInterval> lookup(size_t src_path_id,
                                       const std::string& tgt_haplotype,
                                       size_t local_start,
                                       size_t local_end) const;

    /// Number of (src_path_id, tgt_haplotype) pairs stored.
    size_t num_entries() const { return entries_.size(); }

    /// True when the loaded table carries target intervals ("B2" form). False
    /// for tables whose entries only name a target path. Set by load(); for a
    /// table populated in memory via add_mapping() it reflects whether any
    /// mapping supplied a non-empty target interval.
    bool has_target_coords() const { return has_target_coords_; }

    /// Total number of IntervalMapping segments across all keys.
    size_t total_segments() const;

    /// All (src_path_id, tgt_haplotype) keys.
    ///
    /// WARNING: materializes one string per key. An all-pairs B2 table has tens
    /// of millions of keys, so this is a whole-table copy — fine for offline
    /// tools, never acceptable on a request path. Use target_haplotypes_for()
    /// to get one source path's targets.
    std::vector<std::pair<size_t, std::string>> keys() const;

    /// Target haplotype names stored for one source path.
    ///
    /// Keys are ordered by (src_path_id, tgt_haplotype), so one source path's
    /// keys are contiguous and this is a binary search plus a short walk:
    /// O(log n + k) in place of keys()'s O(n) copy of the entire table.
    std::vector<std::string> target_haplotypes_for(size_t src_path_id) const;

    /// Segment list for a specific key (for inspection/debug).
    std::vector<IntervalMapping> segments(size_t src_path_id,
                                          const std::string& tgt_haplotype) const;

    /// Serialize to binary stream.
    void serialize(std::ostream& out) const;

    /// Load from binary stream.
    void load(std::istream& in);

private:
    struct Key {
        size_t      src_path_id;
        std::string tgt_haplotype;
        bool operator<(const Key& o) const {
            if (src_path_id != o.src_path_id) return src_path_id < o.src_path_id;
            return tgt_haplotype < o.tgt_haplotype;
        }
    };
    std::map<Key, std::vector<IntervalMapping>> entries_;
    bool has_target_coords_ = false;
};

/**
 * Streaming writer for Table 2 files.
 *
 * TranslationTable2 keeps a std::map keyed by (path id, haplotype string). An
 * all-pairs B2 table has one key per (source path, target haplotype) — tens of
 * millions of keys holding one segment each — where that map costs several GB of
 * node, string and vector overhead on top of the payload. This writer takes
 * groups in sorted key order and emits the file directly, so a builder never
 * materializes the map.
 *
 * Usage: call begin_key() for each key in ascending (src_path_id, haplotype)
 * order, add_segment() for its segments in ascending src_start order, then
 * finish(). Keys must not repeat.
 */
class TranslationTable2Writer {
public:
    /// Writes to `out`; `haplotype_names` must list every name begin_key() will
    /// use (it is interned in the file, so names cost 4 bytes per key, not a
    /// full string). Nothing is written until finish().
    TranslationTable2Writer(std::ostream& out,
                            const std::vector<std::string>& haplotype_names);

    /// Start a new key. Returns false if the name is not in haplotype_names or
    /// the key is not greater than the previous one.
    bool begin_key(size_t src_path_id, const std::string& tgt_haplotype);

    /// Append a segment to the current key.
    void add_segment(const IntervalMapping& seg);

    /// Write header, name table and all buffered keys. Must be called once.
    void finish();

    size_t num_keys() const { return keys_.size(); }
    size_t num_segments() const { return segments_.size(); }

private:
    std::ostream& out_;
    std::vector<std::string> names_;
    std::map<std::string, uint32_t> name_ids_;
    /// (src_path_id, hap_id, first segment index, segment count)
    struct KeyRec { size_t src_path_id; uint32_t hap_id; size_t first; size_t count; };
    std::vector<KeyRec> keys_;
    std::vector<IntervalMapping> segments_;
    bool finished_ = false;
};

} // namespace panindexer

#endif // PANGENOME_INDEX_TRANSLATION_TABLES_HPP
