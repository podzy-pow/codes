#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>
#include "gf4.h"
#include "series.h"
#include "codes.h"
#include "search.h"

namespace py = pybind11;

using namespace cppcodes;

PYBIND11_MODULE(codeslib, m){
    m.doc() = "codeslib";

    py::class_<gf4>(m, "gf4")
        .def(py::init<>())
        .def(py::init<short>())
        .def(py::init<std::string>())
        .def("trace", &gf4::trace)
        .def(py::self + py::self)
        .def(py::self * py::self)
        .def("__repr__", &gf4::toString)
        .def("__eq__", [](const gf4& a, const gf4& b){ return a == b;})
        .def("__eq__", [](const gf4& a, const short v){ return a == v;})
        .def("__eq__", [](const gf4& a, const std::string& v){ return a == v;})
    ;

    py::class_<Series>(m, "Series")
        .def(py::init<>())
        .def(py::init<size_t, size_t>())
        .def(py::init<std::vector<gf4>, size_t>())
        .def(py::init<std::vector<gf4>>())
        .def(py::init<std::string, size_t>())
        .def(py::init<std::string>())
        .def("__getitem__", [](const Series& s, size_t i){ return s[i]; })
        .def("__setitem__", [](Series& s, size_t i, const gf4& v){ s[i] = v; })
        .def("strip", &Series::strip)
        .def("__repr__", &Series::toString)
        .def(py::self + py::self)
        .def(py::self * py::self)
        .def("strip", &Series::strip)
        .def("inverse", &Series::inverse)
        .def("conj", &Series::conj)
    ;

    py::class_<Code>(m, "Code")
        .def(py::init<size_t, size_t>())
        .def(py::init<size_t>())
        .def(py::init<std::vector<Series>, size_t, size_t>())
        .def(py::init<std::vector<Series>, size_t>())
        .def("__getitem__", [](const Code& c, size_t i){ return c[i]; })
        .def("__setitem__", [](Code& c, size_t i, const Series& v){ c[i] = v; })
        .def("add", &Code::add)
        .def("remove", &Code::remove)
        .def("validate", &Code::validate)
        .def("minDistance", &Code::minDistance)
        .def("isSelfOrthogonal", &Code::isSelfOrthogonal)
        .def_readonly("n", &Code::n)
        .def_readonly("k", &Code::k)
        .def_readwrite("generators", &Code::generators)
        .def("__repr__", &Code::toString)
    ;

    py::class_<SearchSelfOrthogonal>(m, "SearchSelfOrthogonal")
        .def(py::init<size_t, size_t>())
        .def_property_readonly("k", &SearchSelfOrthogonal::getK)
        .def_property_readonly("n", &SearchSelfOrthogonal::getN)
        .def_property_readonly("degree", &SearchSelfOrthogonal::getDegree)
        .def("find", &SearchSelfOrthogonal::find)
    ;

    py::class_<CodeGenerator>(m, "CodeGenerator")
        .def(py::init<size_t, size_t, size_t>())
        .def("next", &CodeGenerator::next)
    ;
}