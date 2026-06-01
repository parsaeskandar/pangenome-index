#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pangenome_server.hpp"

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
        .def("__repr__", [](const AnchorBuildPyResult& r) {
            return "<AnchorBuildResult status='" + r.status +
                   "' n_anchors=" + std::to_string(r.anchors.size()) +
                   " target_path_length=" + std::to_string(r.target_path_length) + ">";
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
