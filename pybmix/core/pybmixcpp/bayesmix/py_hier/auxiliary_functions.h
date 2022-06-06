#ifndef PYBMIX_AUXILIARY_FUNCTIONS_H
#define PYBMIX_AUXILIARY_FUNCTIONS_H
#include <random>
#include <list>
#include <Eigen/Dense>
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>
#include "pybind11/numpy.h"

#include <iostream>

namespace py = pybind11;
using namespace py::literals;

void synchronize_cpp_to_py_state(const std::mt19937 &cpp_gen,
                                 py::object &py_gen);

void synchronize_py_to_cpp_state(std::mt19937 &cpp_gen,
                                 const py::object &py_gen);

std::vector<double> list_to_vector(py::list &x);

std::list<Eigen::RowVectorXd> list_to_list_eigen(py::list &x);

#endif //PYBMIX_AUXILIARY_FUNCTIONS_H
