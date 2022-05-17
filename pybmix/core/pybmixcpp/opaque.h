#ifndef OPAQUE_
#define OPAQUE_

#include <pybind11/pybind11.h>
#include <vector>
#include <pybind11/embed.h>
#include <memory>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

PYBIND11_MAKE_OPAQUE(std::vector<double>);

#endif

