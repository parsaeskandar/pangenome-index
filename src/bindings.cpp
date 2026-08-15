#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pangenome_server.hpp"

#include <iomanip>
#include <sstream>

namespace py = pybind11;

PYBIND11_MODULE(liftover_ext, m) {
    m.doc() = "Pangenome coordinate liftover extension";

    py::class_<TranslatedInterval>(m, "TranslatedInterval")
        .def(py::init<>())
        .def_readwrite("haplotype", &TranslatedInterval::haplotype)
        .def_readwrite("start",     &TranslatedInterval::start)
        .def_readwrite("end",       &TranslatedInterval::end)
        .def_readwrite("strand",    &TranslatedInterval::strand)
        .def("__repr__", [](const TranslatedInterval& ti) {
            return "<TranslatedInterval haplotype='" + ti.haplotype +
                   "' start=" + std::to_string(ti.start) +
                   " end=" + std::to_string(ti.end) +
                   " strand=" + std::string(1, ti.strand) + ">";
        });

    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const std::invalid_argument& e) {
            PyErr_SetString(PyExc_ValueError, e.what());
        }
    });

    py::class_<AnchorRecord>(m, "AnchorRecord")
        .def(py::init<>())
        .def_readwrite("source_mapping_begin",   &AnchorRecord::source_mapping_begin)
        .def_readwrite("source_mapping_end",     &AnchorRecord::source_mapping_end)
        .def_readwrite("read_begin_offset",      &AnchorRecord::read_begin_offset)
        .def_readwrite("read_end_offset",        &AnchorRecord::read_end_offset)
        .def_readwrite("path_offset_step_begin", &AnchorRecord::path_offset_step_begin)
        .def_readwrite("path_offset_step_end",   &AnchorRecord::path_offset_step_end)
        .def_readwrite("gbwt_edge_begin_node",   &AnchorRecord::gbwt_edge_begin_node)
        .def_readwrite("gbwt_edge_begin_offset", &AnchorRecord::gbwt_edge_begin_offset)
        .def_readwrite("gbwt_edge_end_node",     &AnchorRecord::gbwt_edge_end_node)
        .def_readwrite("gbwt_edge_end_offset",   &AnchorRecord::gbwt_edge_end_offset)
        .def("__repr__", [](const AnchorRecord& a) {
            return "<AnchorRecord read=[" +
                   std::to_string(a.read_begin_offset) + "," +
                   std::to_string(a.read_end_offset) + ") path=[" +
                   std::to_string(a.path_offset_step_begin) + "," +
                   std::to_string(a.path_offset_step_end) + "] src=[" +
                   std::to_string(a.source_mapping_begin) + "," +
                   std::to_string(a.source_mapping_end) + ")>";
        });

    py::class_<AnchorBuildPyResult>(m, "AnchorBuildResult")
        .def(py::init<>())
        .def_readwrite("status",             &AnchorBuildPyResult::status)
        .def_readwrite("anchors",            &AnchorBuildPyResult::anchors)
        .def_readwrite("target_path_length", &AnchorBuildPyResult::target_path_length)
        .def_readwrite("target_rev_strand",  &AnchorBuildPyResult::target_rev_strand)
        // find_sequences_for_tag LF-cost diagnostics (see AnchorBuildPyResult).
        .def_readwrite("find_seq_calls",     &AnchorBuildPyResult::find_seq_calls)
        .def_readwrite("find_seq_runs",      &AnchorBuildPyResult::find_seq_runs)
        .def_readwrite("find_seq_lf_steps",  &AnchorBuildPyResult::find_seq_lf_steps)
        .def_readwrite("find_seq_visits",    &AnchorBuildPyResult::find_seq_visits)
        .def_readwrite("last_run_nav_steps", &AnchorBuildPyResult::last_run_nav_steps)
        .def_readwrite("last_run_length",    &AnchorBuildPyResult::last_run_length)
        // target-path walk diagnostics
        .def_readwrite("walk_lf_steps",      &AnchorBuildPyResult::walk_lf_steps)
        .def_readwrite("walk_span",          &AnchorBuildPyResult::walk_span)
        .def_readwrite("first_anchor_base",  &AnchorBuildPyResult::first_anchor_base)
        .def_readwrite("last_anchor_base",   &AnchorBuildPyResult::last_anchor_base)
        // wall-clock attribution (ms)
        .def_readwrite("find_seq_ms",        &AnchorBuildPyResult::find_seq_ms)
        .def_readwrite("decompress_sa_ms",   &AnchorBuildPyResult::decompress_sa_ms)
        .def_readwrite("walk_ms",            &AnchorBuildPyResult::walk_ms)
        .def_readwrite("decompress_sa_calls",   &AnchorBuildPyResult::decompress_sa_calls)
        .def_readwrite("decompress_sa_entries", &AnchorBuildPyResult::decompress_sa_entries)
        .def_readwrite("n_target_subpaths",  &AnchorBuildPyResult::n_target_subpaths)
        .def_readwrite("n_source_mappings",  &AnchorBuildPyResult::n_source_mappings)
        .def("__repr__", [](const AnchorBuildPyResult& r) {
            return "<AnchorBuildResult status='" + r.status +
                   "' n_anchors=" + std::to_string(r.anchors.size()) +
                   " target_path_length=" + std::to_string(r.target_path_length) +
                   " find_seq_lf_steps=" + std::to_string(r.find_seq_lf_steps) + ">";
        });

    py::class_<HaplotypeCoverage>(m, "HaplotypeCoverage")
        .def(py::init<>())
        .def_readwrite("haplotype",  &HaplotypeCoverage::haplotype)
        .def_readwrite("covered_bp", &HaplotypeCoverage::covered_bp)
        .def_readwrite("coverage",   &HaplotypeCoverage::coverage)
        .def("__repr__", [](const HaplotypeCoverage& h) {
            std::ostringstream ss;
            ss << "<HaplotypeCoverage " << h.haplotype << " "
               << std::fixed << std::setprecision(1) << h.coverage
               << "% (" << h.covered_bp << " bp)>";
            return ss.str();
        });

    py::class_<Index>(m, "Index")
        .def(py::init<>())
        .def("load", &Index::load,
             py::arg("gbz_path"),
             py::arg("ri_path"),
             py::arg("tags_path"),
             py::arg("gbwt_ri_path"),
             py::arg("table1_path"),
             py::arg("table2_path"),
             "Load all index files into memory (call once at startup).")
        .def("translate", &Index::translate,
             py::arg("src_haplotype"),
             py::arg("start"),
             py::arg("end"),
             py::arg("tgt_haplotype"),
             "Translate coordinates from source to target haplotype.")
        .def("get_haplotype_names", &Index::get_haplotype_names,
             "Return list of valid haplotype names in the loaded index.")
        .def("haplotype_coverage", &Index::haplotype_coverage,
             py::arg("graph_alignment_gaf"),
             py::arg("min_coverage") = 0.0,
             py::arg("include_zero") = false,
             "Score every haplotype by the percentage of the alignment's aligned "
             "bases that lie on nodes it also visits. Returns HaplotypeCoverage "
             "records sorted by descending coverage. min_coverage=0 reports every "
             "haplotype sharing any node; include_zero also lists those sharing "
             "none, scored 0.")
        .def("translatable_haplotypes", &Index::translatable_haplotypes,
             py::arg("src_haplotype"),
             py::arg("start"),
             py::arg("end"),
             "For a source contig interval, return the sorted list of target "
             "haplotype names that have a homologous region overlapping it "
             "(a Table-2 overlap check; no coordinate trace).")
        .def("build_surject_anchors",
             [](const Index& self, const std::string& gaf_str,
                const std::string& target_haplotype) {
                 // Return the list of AnchorRecord directly so the smoke test
                 // can call len()/iterate on it. Status + path_length are
                 // available via build_surject_anchors_full() below.
                 AnchorBuildPyResult result =
                     self.build_surject_anchors(gaf_str, target_haplotype);
                 return result.anchors;
             },
             py::arg("graph_alignment_gaf"),
             py::arg("target_haplotype"),
             "Build surjection anchors from a graph GAF and a target haplotype name. "
             "Returns a list of AnchorRecord (empty if no anchors could be built).")
        .def("build_surject_anchors_full", &Index::build_surject_anchors,
             py::arg("graph_alignment_gaf"),
             py::arg("target_haplotype"),
             "Same as build_surject_anchors but returns the full "
             "AnchorBuildResult (status, anchors, target_path_length).");
}
