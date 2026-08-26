/**
 * Translation Table 1 implementation.
 * Converts (named path, global interval) → list of (path_id, local interval).
 */

#include "pangenome_index/translation_tables.hpp"
#include <algorithm>
#include <cassert>
#include <istream>
#include <ostream>
#include <stdexcept>

namespace panindexer {

void TranslationTable1::add_subpath(const std::string& path_name,
                                     size_t path_id,
                                     size_t subpath_start,
                                     size_t length) {
    if (name_to_subpaths_.find(path_name) == name_to_subpaths_.end()) {
        names_.push_back(path_name);
    }
    name_to_subpaths_[path_name].push_back(
        SubpathInfo{path_id, subpath_start, length});
}

std::vector<PathInterval> TranslationTable1::lookup(const std::string& path_name,
                                                    size_t global_start,
                                                    size_t global_end) const {
    auto it = name_to_subpaths_.find(path_name);
    if (it == name_to_subpaths_.end()) {
        return {};
    }
    const std::vector<SubpathInfo>& subpaths = it->second;
    std::vector<PathInterval> result;
    for (const SubpathInfo& sp : subpaths) {
        size_t subpath_end = sp.subpath_start + sp.length;
        // Overlap [global_start, global_end) with [sp.subpath_start, subpath_end)
        size_t overlap_start = std::max(global_start, sp.subpath_start);
        size_t overlap_end = std::min(global_end, subpath_end);
        if (overlap_start >= overlap_end) {
            continue;
        }
        PathInterval pi;
        pi.path_id = sp.path_id;
        pi.start = overlap_start - sp.subpath_start;
        pi.end = overlap_end - sp.subpath_start;
        result.push_back(pi);
    }
    return result;
}

std::vector<PathInterval> TranslationTable1::lookup_by_name_id(size_t name_id,
                                                               size_t global_start,
                                                               size_t global_end) const {
    if (name_id >= names_.size()) {
        return {};
    }
    return lookup(names_[name_id], global_start, global_end);
}

std::vector<std::string> TranslationTable1::names() const {
    return names_;
}

std::vector<SubpathInfo> TranslationTable1::subpaths(const std::string& path_name) const {
    auto it = name_to_subpaths_.find(path_name);
    if (it == name_to_subpaths_.end()) {
        return {};
    }
    return it->second;
}

namespace {

const uint32_t TABLE1_MAGIC = 0x54543100;  // "TT1\0"
const uint32_t TABLE1_VERSION = 1;

void write_uint32(std::ostream& out, uint32_t x) {
    out.put(static_cast<char>(x & 0xff));
    out.put(static_cast<char>((x >> 8) & 0xff));
    out.put(static_cast<char>((x >> 16) & 0xff));
    out.put(static_cast<char>((x >> 24) & 0xff));
}

void write_uint64(std::ostream& out, uint64_t x) {
    for (int i = 0; i < 8; ++i) {
        out.put(static_cast<char>(x & 0xff));
        x >>= 8;
    }
}

uint32_t read_uint32(std::istream& in) {
    uint32_t x = 0;
    for (int i = 0; i < 4; ++i) {
        int c = in.get();
        if (c == std::char_traits<char>::eof()) {
            throw std::runtime_error("TranslationTable1::load: unexpected EOF");
        }
        x |= static_cast<uint32_t>(static_cast<unsigned char>(c)) << (i * 8);
    }
    return x;
}

uint64_t read_uint64(std::istream& in) {
    uint64_t x = 0;
    for (int i = 0; i < 8; ++i) {
        int c = in.get();
        if (c == std::char_traits<char>::eof()) {
            throw std::runtime_error("TranslationTable1::load: unexpected EOF");
        }
        x |= static_cast<uint64_t>(static_cast<unsigned char>(c)) << (i * 8);
    }
    return x;
}

void write_string(std::ostream& out, const std::string& s) {
    write_uint64(out, s.size());
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
}

std::string read_string(std::istream& in) {
    uint64_t len = read_uint64(in);
    std::string s(len, '\0');
    in.read(&s[0], static_cast<std::streamsize>(len));
    if (!in) {
        throw std::runtime_error("TranslationTable1::load: failed to read string");
    }
    return s;
}

} // namespace

void TranslationTable1::serialize(std::ostream& out) const {
    write_uint32(out, TABLE1_MAGIC);
    write_uint32(out, TABLE1_VERSION);
    write_uint64(out, names_.size());
    for (const std::string& name : names_) {
        write_string(out, name);
        const std::vector<SubpathInfo>& subpaths = name_to_subpaths_.at(name);
        write_uint64(out, subpaths.size());
        for (const SubpathInfo& sp : subpaths) {
            write_uint64(out, sp.path_id);
            write_uint64(out, sp.subpath_start);
            write_uint64(out, sp.length);
        }
    }
}

void TranslationTable1::load(std::istream& in) {
    uint32_t magic = read_uint32(in);
    if (magic != TABLE1_MAGIC) {
        throw std::runtime_error("TranslationTable1::load: invalid magic (not a Table1 file?)");
    }
    uint32_t version = read_uint32(in);
    if (version != TABLE1_VERSION) {
        throw std::runtime_error("TranslationTable1::load: unsupported version");
    }
    name_to_subpaths_.clear();
    names_.clear();
    uint64_t num_names = read_uint64(in);
    for (uint64_t i = 0; i < num_names; ++i) {
        std::string name = read_string(in);
        uint64_t num_subpaths = read_uint64(in);
        std::vector<SubpathInfo> subpaths;
        subpaths.reserve(num_subpaths);
        for (uint64_t j = 0; j < num_subpaths; ++j) {
            SubpathInfo sp;
            sp.path_id = read_uint64(in);
            sp.subpath_start = read_uint64(in);
            sp.length = read_uint64(in);
            subpaths.push_back(sp);
        }
        names_.push_back(name);
        name_to_subpaths_[name] = std::move(subpaths);
    }
}

// ============================================================
// TranslationTable2 implementation
// ============================================================

void TranslationTable2::add_mapping(size_t src_path_id,
                                     const std::string& tgt_haplotype,
                                     const IntervalMapping& mapping) {
    Key k{src_path_id, tgt_haplotype};
    if (mapping.tgt_end > mapping.tgt_start) has_target_coords_ = true;
    entries_[k].push_back(mapping);
}

void TranslationTable2::finalize() {
    for (auto& kv : entries_) {
        auto& segs = kv.second;
        std::sort(segs.begin(), segs.end(),
                  [](const IntervalMapping& a, const IntervalMapping& b) {
                      return a.src_start < b.src_start;
                  });
        // Overlapping source intervals are allowed: src[a,b] and src[c,d] can both
        // map to different tgt_path_ids (e.g. src[0,100]->i, src[10,90]->j).
    }
}

std::vector<TargetInterval> TranslationTable2::lookup(size_t src_path_id,
                                                       const std::string& tgt_haplotype,
                                                       size_t local_start,
                                                       size_t local_end) const {
    Key k{src_path_id, tgt_haplotype};
    auto it = entries_.find(k);
    if (it == entries_.end()) {
        return {};
    }
    const auto& segs = it->second;
    if (segs.empty() || local_start >= local_end) {
        return {};
    }

    // Binary search for first segment whose src_end > local_start
    size_t lo = 0, hi = segs.size();
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (segs[mid].src_end <= local_start) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    std::vector<TargetInterval> result;
    for (size_t i = lo; i < segs.size(); ++i) {
        const IntervalMapping& seg = segs[i];
        if (seg.src_start >= local_end) break;  // past query range

        size_t clip_src_start = std::max(local_start, seg.src_start);
        size_t clip_src_end   = std::min(local_end,   seg.src_end);
        if (clip_src_start >= clip_src_end) continue;

        TargetInterval ti;
        ti.tgt_path_id = seg.tgt_path_id;
        if (has_target_coords_) {
            // Deliberately the WHOLE block, not a slice interpolated from the
            // clipped source range. A B2 block is coarse: its internal mapping
            // is unknown, so interpolating could exclude the true target
            // position whenever a sizeable indel sits inside the block — a
            // false negative, the one error class this table must never
            // introduce. Callers wanting a tighter bound can interpolate using
            // src_start/src_end, accepting that risk explicitly.
            ti.tgt_start   = seg.tgt_start;
            ti.tgt_end     = seg.tgt_end;
            ti.tgt_reverse = seg.tgt_reverse;
            ti.has_coords  = true;
        }
        result.push_back(ti);
    }
    return result;
}

size_t TranslationTable2::total_segments() const {
    size_t total = 0;
    for (const auto& kv : entries_) {
        total += kv.second.size();
    }
    return total;
}

std::vector<std::pair<size_t, std::string>> TranslationTable2::keys() const {
    std::vector<std::pair<size_t, std::string>> result;
    result.reserve(entries_.size());
    for (const auto& kv : entries_) {
        result.emplace_back(kv.first.src_path_id, kv.first.tgt_haplotype);
    }
    return result;
}

std::vector<std::string>
TranslationTable2::target_haplotypes_for(size_t src_path_id) const {
    std::vector<std::string> out;
    // Key ordering is (src_path_id, tgt_haplotype), so seeking to the first key
    // for this path and walking while the path id holds visits exactly the
    // entries that belong to it.
    Key lo{src_path_id, std::string()};
    for (auto it = entries_.lower_bound(lo);
         it != entries_.end() && it->first.src_path_id == src_path_id; ++it) {
        out.push_back(it->first.tgt_haplotype);
    }
    return out;
}

std::vector<IntervalMapping> TranslationTable2::segments(size_t src_path_id,
                                                          const std::string& tgt_haplotype) const {
    Key k{src_path_id, tgt_haplotype};
    auto it = entries_.find(k);
    if (it == entries_.end()) return {};
    return it->second;
}

namespace {

const uint32_t TABLE2_MAGIC   = 0x54543200;  // "TT2\0"
// v1: had tgt_start/tgt_end.  v2: dropped them (target path id only).
// v3: restores them ("B2" form) and interns haplotype names, which matters
//     because an all-pairs B2 table has one key per (source path, haplotype)
//     and a repeated name string per key would outweigh the payload.
// v4: interned names, NO target interval. This is the routing-only table the
// query path actually uses; v2 expressed the same content but repeated the
// haplotype name as a string in every key, which costs hundreds of MB across
// tens of millions of keys.
const uint32_t TABLE2_VERSION = 3;
const uint32_t TABLE2_VERSION_ROUTING = 4;
const uint32_t TABLE2_FLAG_TGT_REVERSE = 1u;

}  // namespace

void TranslationTable2::serialize(std::ostream& out) const {
    // Intern haplotype names so each key costs a 4-byte id, not a string.
    std::vector<std::string> names;
    std::map<std::string, uint32_t> ids;
    for (const auto& kv : entries_) {
        if (ids.emplace(kv.first.tgt_haplotype,
                        static_cast<uint32_t>(names.size())).second) {
            names.push_back(kv.first.tgt_haplotype);
        }
    }

    // Emit the LEGACY v2 layout when no mapping carries a target interval.
    // Writing v3 regardless would advertise target coordinates that are
    // uniformly zero, and load() would believe them — so a builder that never
    // sets them (build_translation_tables, build_table2_coarse) keeps producing
    // exactly the bytes it produced before, and only a table that really has
    // coordinates claims to.
    const bool v3 = has_target_coords_;
    write_uint32(out, TABLE2_MAGIC);
    write_uint32(out, v3 ? TABLE2_VERSION : 2u);
    if (v3) {
        write_uint64(out, static_cast<uint64_t>(names.size()));
        for (const std::string& n : names) write_string(out, n);
    }
    write_uint64(out, static_cast<uint64_t>(entries_.size()));
    for (const auto& kv : entries_) {
        write_uint64(out, static_cast<uint64_t>(kv.first.src_path_id));
        if (v3) write_uint32(out, ids[kv.first.tgt_haplotype]);
        else    write_string(out, kv.first.tgt_haplotype);
        const auto& segs = kv.second;
        write_uint64(out, static_cast<uint64_t>(segs.size()));
        for (const IntervalMapping& seg : segs) {
            write_uint64(out, static_cast<uint64_t>(seg.src_start));
            write_uint64(out, static_cast<uint64_t>(seg.src_end));
            write_uint64(out, static_cast<uint64_t>(seg.tgt_path_id));
            if (v3) {
                write_uint64(out, static_cast<uint64_t>(seg.tgt_start));
                write_uint64(out, static_cast<uint64_t>(seg.tgt_end));
                write_uint32(out, seg.tgt_reverse ? TABLE2_FLAG_TGT_REVERSE : 0u);
            }
        }
    }
}

void TranslationTable2::load(std::istream& in) {
    uint32_t magic = read_uint32(in);
    if (magic != TABLE2_MAGIC) {
        throw std::runtime_error("TranslationTable2::load: invalid magic (not a Table2 file?)");
    }
    uint32_t version = read_uint32(in);
    if (version != 1 && version != 2 && version != TABLE2_VERSION &&
        version != TABLE2_VERSION_ROUTING) {
        throw std::runtime_error("TranslationTable2::load: unsupported version");
    }
    // v1 had target intervals, v2 dropped them, v3 restored them, v4 drops them
    // again but keeps v3's interned names.
    has_target_coords_ = (version == 1 || version == TABLE2_VERSION);
    const bool interned = (version >= 3);

    entries_.clear();
    std::vector<std::string> names;
    if (interned) {
        uint64_t n_names = read_uint64(in);
        names.reserve(static_cast<size_t>(n_names));
        for (uint64_t i = 0; i < n_names; ++i) names.push_back(read_string(in));
    }

    uint64_t num_keys = read_uint64(in);
    for (uint64_t i = 0; i < num_keys; ++i) {
        Key k;
        k.src_path_id = static_cast<size_t>(read_uint64(in));
        if (interned) {
            uint32_t hid = read_uint32(in);
            if (hid >= names.size()) {
                throw std::runtime_error("TranslationTable2::load: haplotype id out of range");
            }
            k.tgt_haplotype = names[hid];
        } else {
            k.tgt_haplotype = read_string(in);
        }
        uint64_t num_segs = read_uint64(in);
        std::vector<IntervalMapping> segs;
        segs.reserve(static_cast<size_t>(num_segs));
        for (uint64_t j = 0; j < num_segs; ++j) {
            IntervalMapping seg;
            seg.src_start   = static_cast<size_t>(read_uint64(in));
            seg.src_end     = static_cast<size_t>(read_uint64(in));
            seg.tgt_path_id = static_cast<size_t>(read_uint64(in));
            if (has_target_coords_) {
                seg.tgt_start = static_cast<size_t>(read_uint64(in));
                seg.tgt_end   = static_cast<size_t>(read_uint64(in));
                if (version >= 3) {
                    seg.tgt_reverse = (read_uint32(in) & TABLE2_FLAG_TGT_REVERSE) != 0;
                }
            }
            segs.push_back(seg);
        }
        entries_[k] = std::move(segs);
    }
}


// --------------------------------------------------------------------------
// TranslationTable2Writer
// --------------------------------------------------------------------------

TranslationTable2Writer::TranslationTable2Writer(
        std::ostream& out, const std::vector<std::string>& haplotype_names,
        bool with_target_coords)
    : out_(out), names_(haplotype_names), with_coords_(with_target_coords) {
    for (size_t i = 0; i < names_.size(); ++i) {
        name_ids_[names_[i]] = static_cast<uint32_t>(i);
    }
}

bool TranslationTable2Writer::begin_key(size_t src_path_id,
                                        const std::string& tgt_haplotype) {
    auto it = name_ids_.find(tgt_haplotype);
    if (it == name_ids_.end()) return false;
    if (!keys_.empty()) {
        const KeyRec& prev = keys_.back();
        // Same ordering as TranslationTable2's Key: path id, then name. Compare
        // names by string, not by interned id — ids follow insertion order in
        // haplotype_names, which need not be alphabetical.
        const bool ordered =
            (src_path_id > prev.src_path_id) ||
            (src_path_id == prev.src_path_id &&
             tgt_haplotype > names_[prev.hap_id]);
        if (!ordered) return false;
    }
    keys_.push_back(KeyRec{src_path_id, it->second, segments_.size(), 0});
    return true;
}

void TranslationTable2Writer::add_segment(const IntervalMapping& seg) {
    if (keys_.empty()) return;
    segments_.push_back(seg);
    keys_.back().count++;
}

void TranslationTable2Writer::finish() {
    if (finished_) return;
    finished_ = true;
    write_uint32(out_, TABLE2_MAGIC);
    write_uint32(out_, with_coords_ ? TABLE2_VERSION : TABLE2_VERSION_ROUTING);
    write_uint64(out_, static_cast<uint64_t>(names_.size()));
    for (const std::string& n : names_) write_string(out_, n);
    write_uint64(out_, static_cast<uint64_t>(keys_.size()));
    for (const KeyRec& k : keys_) {
        write_uint64(out_, static_cast<uint64_t>(k.src_path_id));
        write_uint32(out_, k.hap_id);
        write_uint64(out_, static_cast<uint64_t>(k.count));
        for (size_t i = 0; i < k.count; ++i) {
            const IntervalMapping& seg = segments_[k.first + i];
            write_uint64(out_, static_cast<uint64_t>(seg.src_start));
            write_uint64(out_, static_cast<uint64_t>(seg.src_end));
            write_uint64(out_, static_cast<uint64_t>(seg.tgt_path_id));
            if (with_coords_) {
                write_uint64(out_, static_cast<uint64_t>(seg.tgt_start));
                write_uint64(out_, static_cast<uint64_t>(seg.tgt_end));
                write_uint32(out_, seg.tgt_reverse ? TABLE2_FLAG_TGT_REVERSE : 0u);
            }
        }
    }
}

} // namespace panindexer
