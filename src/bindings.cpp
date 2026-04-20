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

    py::class_<GraphPathVisit>(m, "GraphPathVisit")
        .def(py::init<>())
        .def(py::init<int64_t, bool, int32_t>(),
             py::arg("node_id"),
             py::arg("is_reverse"),
             py::arg("from_length"))
        .def_readwrite("node_id", &GraphPathVisit::node_id)
        .def_readwrite("is_reverse", &GraphPathVisit::is_reverse)
        .def_readwrite("from_length", &GraphPathVisit::from_length);

    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const std::invalid_argument& e) {
            PyErr_SetString(PyExc_ValueError, e.what());
        }
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
        .def("accept_graph_path", &Index::accept_graph_path,
             py::arg("path"),
             "Pass graph path visits from a mapper (stub; reserved for haplotype resolution).");
}
