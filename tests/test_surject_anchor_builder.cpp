// test_surject_anchor_builder.cpp
//
// Tests for build_surject_anchors_for_path() using the pre-built test
// indexes under coord_translation_tests/test0_target_loop/.
//
// Test graph (graph.gfa):
//   Nodes:
//     1: ATCG  (len 4)
//     2: A     (len 1)
//     3: G     (len 1)
//     4: TTC   (len 3)
//     5: A     (len 1)
//     6: C     (len 1)
//     7: TAG   (len 3)
//
//   Paths:
//     x: 1+, 2+, 4+, 5+, 7+              (12 bp total)
//     y: 1+, 3+, 4+, 1+, 2+, 4+, 6+, 7+  (20 bp total, loops back through 1+, 4+)
//
// In our tests path x stands in for the giraffe-output graph alignment
// (source); path y is the target haplotype we surject onto. We construct
// synthetic SourceMappings that walk path x's nodes in order, then ask the
// builder to find anchors against y.

#include <gtest/gtest.h>

#include <gbwt/fast_locate.h>
#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>
#include <sdsl/simple_sds.hpp>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "pangenome_index/r-index.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include "pangenome_index/surject_anchor_builder.hpp"
#include "pangenome_index/translation_tables.hpp"

namespace {

// ---------------------------------------------------------------------------
// Test data and helper for loading the pre-built indexes
// ---------------------------------------------------------------------------

const std::string TEST_DIR =
    "../coord_translation_tests/test0_target_loop/";

struct LoadedIndexes {
    gbwtgraph::GBZ gbz;
    panindexer::FastLocate rlbwt;
    panindexer::SampledTagArray sampled;
    gbwt::FastLocate gbwt_fast_locate;
    panindexer::TranslationTable1 t1;
};

std::unique_ptr<LoadedIndexes> load_all() {
    auto idx = std::make_unique<LoadedIndexes>();

    // 1) GBZ
    sdsl::simple_sds::load_from(idx->gbz, TEST_DIR + "graph.gbz");

    // 2) RLBWT r-index
    {
        std::ifstream rin(TEST_DIR + "rlbwt_rindex.ri", std::ios::binary);
        if (!rin) throw std::runtime_error("cannot open rlbwt_rindex.ri");
        idx->rlbwt.load_encoded(rin);
        idx->rlbwt.ensure_last_rank();
        idx->rlbwt.ensure_last_select();
    }

    // 3) GBWT FastLocate
    {
        std::ifstream gin(TEST_DIR + "fastlocate.ri", std::ios::binary);
        if (!gin) throw std::runtime_error("cannot open fastlocate.ri");
        idx->gbwt_fast_locate.load(gin);
        idx->gbwt_fast_locate.setGBWT(idx->gbz.index);
    }

    // 4) Sampled tag array
    {
        std::ifstream sin(TEST_DIR + "sampled.tags", std::ios::binary);
        if (!sin) throw std::runtime_error("cannot open sampled.tags");
        idx->sampled.load(sin);
        idx->sampled.ensure_run_rank();
        idx->sampled.ensure_run_select();
    }

    // 5) Translation Table 1 (optional for these tests but loaded for completeness)
    {
        std::ifstream t1in(TEST_DIR + "output.t1", std::ios::binary);
        if (!t1in) throw std::runtime_error("cannot open output.t1");
        idx->t1.load(t1in);
    }

    return idx;
}

// Build a SourceMapping list from a list of (node_id, is_reverse, node_length)
// triples, accumulating read offsets.
std::vector<panindexer::SourceMapping> make_source_mappings(
    const std::vector<std::tuple<int64_t, bool, size_t>>& mappings)
{
    std::vector<panindexer::SourceMapping> out;
    out.reserve(mappings.size());
    size_t read_pos = 0;
    for (const auto& m : mappings) {
        panindexer::SourceMapping sm;
        sm.node_id = std::get<0>(m);
        sm.is_reverse = std::get<1>(m);
        sm.read_begin_offset = read_pos;
        sm.read_end_offset = read_pos + std::get<2>(m);
        sm.node_offset_in_node = 0;
        sm.mapping_from_length = std::get<2>(m);
        out.push_back(sm);
        read_pos += std::get<2>(m);
    }
    return out;
}

// Pretty-print one anchor.
void print_anchor(size_t i, const panindexer::PrecomputedAnchor& a) {
    std::cout
        << "  anchor[" << i << "]"
        << "  src_mappings=[" << a.source_mapping_begin << ", " << a.source_mapping_end << ")"
        << "  read=[" << a.read_begin_offset << ", " << a.read_end_offset << ")"
        << "  target_base=[" << a.path_offset_step_begin << ", " << a.path_offset_step_end << "]"
        << "  gbwt_edge_begin=(node=" << a.gbwt_edge_begin.first
            << ", offset=" << a.gbwt_edge_begin.second << ")"
        << "  gbwt_edge_end=(node=" << a.gbwt_edge_end.first
            << ", offset=" << a.gbwt_edge_end.second << ")"
        << std::endl;
}

void print_result_header(const std::string& label,
                         const panindexer::AnchorBuildResult& r) {
    using S = panindexer::AnchorBuildResult::Status;
    std::cout << "\n=== " << label << " ===\n";
    std::cout << "status            = ";
    switch (r.status) {
        case S::Ok:             std::cout << "Ok"; break;
        case S::EmptyAlignment: std::cout << "EmptyAlignment"; break;
        case S::UnknownPath:    std::cout << "UnknownPath"; break;
        case S::NoCommonNodes:  std::cout << "NoCommonNodes"; break;
    }
    std::cout << "\n"
              << "target_path_length = " << r.target_path_length << "\n"
              << "rev_strand         = " << (r.target_rev_strand ? "true" : "false") << "\n"
              << "num_anchors        = " << r.anchors.size() << std::endl;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST(SurjectAnchorBuilder, MetadataAndPathIds) {
    auto idx = load_all();
    // x = path 0, y = path 1 (GFA-declared order). Two paths total.
    EXPECT_EQ(idx->gbz.index.sequences() / 2, 2u);

    std::cout << "\n=== Graph metadata ===\n"
              << "num_paths     = " << (idx->gbz.index.sequences() / 2) << "\n"
              << "num_nodes     = " << idx->gbz.graph.get_node_count() << "\n";

    // Confirm path 0 (x) and path 1 (y) lengths by walking their GBWT paths.
    for (size_t pid : {0u, 1u}) {
        gbwt::vector_type nodes = idx->gbz.index.extract(gbwt::Path::encode(pid, false));
        size_t total = 0;
        std::vector<int64_t> node_ids;
        for (gbwt::node_type n : nodes) {
            if (n == gbwt::ENDMARKER) break;
            node_ids.push_back(gbwt::Node::id(n));
            total += idx->gbz.graph.get_length(
                idx->gbz.graph.get_handle(gbwt::Node::id(n), gbwt::Node::is_reverse(n)));
        }
        std::cout << "  path[" << pid << "]  length=" << total << "  nodes=[";
        for (size_t k = 0; k < node_ids.size(); ++k) {
            if (k) std::cout << ", ";
            std::cout << node_ids[k];
        }
        std::cout << "]" << std::endl;
    }
}

// Source = path x (nodes 1, 2, 4, 5, 7). Target = path y. Expect:
//   - Common nodes between x and y: {1, 2, 4, 7}  (node 5 is on x only).
//   - y visits node 1 twice (positions 0 and 8); 4 twice (positions 5 and 13).
//   - First anchor's expected target base = 0 (node 1 first visit).
//   - Last anchor's expected target base  = 17 (node 7 visit on y).
TEST(SurjectAnchorBuilder, PathXOntoPathY) {
    auto idx = load_all();

    auto src = make_source_mappings({
        {1, false, 4},  // node 1: ATCG -> read[0,4)
        {2, false, 1},  // node 2: A    -> read[4,5)
        {4, false, 3},  // node 4: TTC  -> read[5,8)
        {5, false, 1},  // node 5: A    -> read[8,9)    (off-target: not on y)
        {7, false, 3},  // node 7: TAG  -> read[9,12)
    });

    constexpr size_t TARGET_PATH_ID = 1;  // path y

    auto result = panindexer::build_surject_anchors_for_path(
        idx->gbz, idx->rlbwt, idx->sampled, idx->gbwt_fast_locate,
        src, TARGET_PATH_ID);

    print_result_header("path x -> path y", result);
    for (size_t i = 0; i < result.anchors.size(); ++i) {
        print_anchor(i, result.anchors[i]);
    }

    using S = panindexer::AnchorBuildResult::Status;
    ASSERT_EQ(result.status, S::Ok);
    ASSERT_EQ(result.target_path_length, 20u);  // y is 20 bases
    ASSERT_FALSE(result.anchors.empty());

    // Every emitted anchor's read range must be inside [0, 12), monotonically
    // increasing, and target base offsets must be inside the target.
    size_t last_read_end = 0;
    for (const auto& a : result.anchors) {
        EXPECT_GE(a.read_begin_offset, last_read_end);
        EXPECT_LE(a.read_end_offset, 12u);
        EXPECT_LE(a.path_offset_step_end, result.target_path_length);
        last_read_end = a.read_end_offset;
    }

    // The off-target source mapping (node 5 on x but not on y) must NOT appear
    // as a covered source mapping in any anchor.
    constexpr size_t NODE5_MAPPING_INDEX = 3;
    for (const auto& a : result.anchors) {
        EXPECT_FALSE(NODE5_MAPPING_INDEX >= a.source_mapping_begin &&
                     NODE5_MAPPING_INDEX <  a.source_mapping_end)
            << "node-5 mapping (index 3) should not be inside any anchor";
    }
}

// Empty source should produce EmptyAlignment status and no anchors.
TEST(SurjectAnchorBuilder, EmptySourceMappings) {
    auto idx = load_all();
    std::vector<panindexer::SourceMapping> empty;

    auto result = panindexer::build_surject_anchors_for_path(
        idx->gbz, idx->rlbwt, idx->sampled, idx->gbwt_fast_locate,
        empty, /*target_gbwt_path_id=*/1);

    print_result_header("empty source", result);
    EXPECT_EQ(result.status, panindexer::AnchorBuildResult::Status::EmptyAlignment);
    EXPECT_TRUE(result.anchors.empty());
}

// Out-of-range target path id should be reported cleanly.
TEST(SurjectAnchorBuilder, UnknownTargetPathId) {
    auto idx = load_all();
    auto src = make_source_mappings({{1, false, 4}});

    auto result = panindexer::build_surject_anchors_for_path(
        idx->gbz, idx->rlbwt, idx->sampled, idx->gbwt_fast_locate,
        src, /*target_gbwt_path_id=*/999);

    print_result_header("unknown target path id", result);
    EXPECT_EQ(result.status, panindexer::AnchorBuildResult::Status::UnknownPath);
}

// Source visits only nodes that the target does NOT visit → NoCommonNodes.
// Path y does not visit node 5 (length 1). Construct a source that only
// touches node 5.
TEST(SurjectAnchorBuilder, NoCommonNodesWithTarget) {
    auto idx = load_all();
    auto src = make_source_mappings({{5, false, 1}});

    auto result = panindexer::build_surject_anchors_for_path(
        idx->gbz, idx->rlbwt, idx->sampled, idx->gbwt_fast_locate,
        src, /*target_gbwt_path_id=*/1);

    print_result_header("source visits only node 5", result);
    EXPECT_EQ(result.status, panindexer::AnchorBuildResult::Status::NoCommonNodes);
    EXPECT_TRUE(result.anchors.empty());
}

// Source = path y itself, target = path y. Every source mapping should match
// the corresponding target step exactly; expect one anchor covering all 8
// mappings, or multiple if the algorithm picks a different visit at the
// repeat nodes — we just print and sanity-check.
TEST(SurjectAnchorBuilder, PathYOntoItself) {
    auto idx = load_all();

    // Path y = 1, 3, 4, 1, 2, 4, 6, 7 with lengths 4, 1, 3, 4, 1, 3, 1, 3.
    auto src = make_source_mappings({
        {1, false, 4},
        {3, false, 1},
        {4, false, 3},
        {1, false, 4},
        {2, false, 1},
        {4, false, 3},
        {6, false, 1},
        {7, false, 3},
    });

    auto result = panindexer::build_surject_anchors_for_path(
        idx->gbz, idx->rlbwt, idx->sampled, idx->gbwt_fast_locate,
        src, /*target_gbwt_path_id=*/1);

    print_result_header("path y -> path y (identity)", result);
    for (size_t i = 0; i < result.anchors.size(); ++i) {
        print_anchor(i, result.anchors[i]);
    }

    EXPECT_EQ(result.status, panindexer::AnchorBuildResult::Status::Ok);
    EXPECT_EQ(result.target_path_length, 20u);
    EXPECT_FALSE(result.anchors.empty());
}

// Reverse strand: source = the reverse-complement traversal of path x
// (nodes 7-, 5-, 4-, 2-, 1- in read order), target = path x. The builder must
// detect the reverse strand, normalize to the target's forward orientation,
// and produce forward-strand anchors on x that cover the SAME target span as
// surjecting forward x onto x.
TEST(SurjectAnchorBuilder, ReverseStrandPathXOntoX) {
    auto idx = load_all();

    // Forward x = 1+, 2+, 4+, 5+, 7+ (lengths 4,1,3,1,3 = 12 bp). A revcomp
    // read of x traverses those nodes in reverse order and flipped orientation.
    auto src_rev = make_source_mappings({
        {7, true, 3},
        {5, true, 1},
        {4, true, 3},
        {2, true, 1},
        {1, true, 4},
    });

    constexpr size_t TARGET_PATH_ID = 0;  // path x

    auto rev = panindexer::build_surject_anchors_for_path(
        idx->gbz, idx->rlbwt, idx->sampled, idx->gbwt_fast_locate,
        src_rev, TARGET_PATH_ID);

    print_result_header("revcomp(x) -> path x", rev);
    for (size_t i = 0; i < rev.anchors.size(); ++i) print_anchor(i, rev.anchors[i]);

    using S = panindexer::AnchorBuildResult::Status;
    ASSERT_EQ(rev.status, S::Ok);
    EXPECT_TRUE(rev.target_rev_strand);        // detected reverse strand
    EXPECT_EQ(rev.target_path_length, 12u);    // x is 12 bp
    ASSERT_FALSE(rev.anchors.empty());

    size_t rev_min_tb = SIZE_MAX, rev_max_tb = 0;
    for (const auto& a : rev.anchors) {
        EXPECT_LE(a.path_offset_step_begin, rev.target_path_length);
        EXPECT_LE(a.path_offset_step_end, rev.target_path_length);
        rev_min_tb = std::min(rev_min_tb, a.path_offset_step_begin);
        rev_max_tb = std::max(rev_max_tb, a.path_offset_step_end);
    }
    // x is colinear (no repeats): the revcomp read maps back onto all of x,
    // so the target span must start at base 0 (node 1).
    EXPECT_EQ(rev_min_tb, 0u);

    // Forward x -> x, for comparison: same target base span, forward strand.
    auto src_fwd = make_source_mappings({
        {1, false, 4}, {2, false, 1}, {4, false, 3}, {5, false, 1}, {7, false, 3},
    });
    auto fwd = panindexer::build_surject_anchors_for_path(
        idx->gbz, idx->rlbwt, idx->sampled, idx->gbwt_fast_locate,
        src_fwd, TARGET_PATH_ID);
    print_result_header("x -> path x (forward, for comparison)", fwd);
    for (size_t i = 0; i < fwd.anchors.size(); ++i) print_anchor(i, fwd.anchors[i]);

    ASSERT_EQ(fwd.status, S::Ok);
    EXPECT_FALSE(fwd.target_rev_strand);

    size_t fwd_min_tb = SIZE_MAX, fwd_max_tb = 0;
    for (const auto& a : fwd.anchors) {
        fwd_min_tb = std::min(fwd_min_tb, a.path_offset_step_begin);
        fwd_max_tb = std::max(fwd_max_tb, a.path_offset_step_end);
    }
    // The reverse and forward surjections of the same locus must cover the
    // same target base span.
    EXPECT_EQ(rev_min_tb, fwd_min_tb);
    EXPECT_EQ(rev_max_tb, fwd_max_tb);
}

}  // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
