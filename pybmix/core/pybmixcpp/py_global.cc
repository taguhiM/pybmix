#include "py_global.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include <vector>
#include <pybind11/embed.h>
#include <memory>
#include <pybind11/pybind11.h>


namespace py = pybind11;
using namespace py::literals;

PYBIND11_MAKE_OPAQUE(std::vector<double>);

PYBIND11_EMBEDDED_MODULE(opaque_module, m)
{
py::bind_vector<std::vector<double>>(m, "VectorDouble");
}

py::scoped_interpreter guard{};

namespace py_global{
    py::module_ numpy = py::module_::import("numpy");
    py::module_ fun = py::module_::import("fun");
    py::module_ numpy_random = py::module_::import("numpy.random");
    py::object py_engine = numpy_random.attr("MT19937")();
    py::object py_gen = numpy_random.attr("Generator")(py_engine);
    py::object posterior_hypers_evaluator = fun.attr("compute_posterior_hypers");
    py::object like_lpdf_evaluator = fun.attr("like_lpdf");
    py::object marg_lpdf_evaluator = fun.attr("marg_lpdf");
    py::object initialize_state_evaluator = fun.attr("initialize_state");
    py::object draw_evaluator = fun.attr("draw");
    py::object update_summary_statistics_evaluator = fun.attr("update_summary_statistics");
    py::object clear_summary_statistics_evaluator = fun.attr("clear_summary_statistics");
};

py::module_ opaque_module = py::module_::import("opaque_module");