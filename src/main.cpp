#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

#include "confpool.h"
#include "molproxy.h"

PYBIND11_MODULE(confpool, m) {
    py::class_<MolProxy>(m, "MolProxy")
        .def("__getitem__", &MolProxy::__getitem__)
        .def("__setitem__", &MolProxy::__setitem__)
        .def("__getattr__", &MolProxy::__getattr__)
        .def("__setattr__", &MolProxy::__setattr__)
        .def("l", &MolProxy::l)
        .def("v", &MolProxy::v)
        .def("z", &MolProxy::z)
        ;

    py::class_<Confpool>(m, "Confpool")
        .def(py::init<>())
        .def("__getitem__", &Confpool::__getitem__)
        .def("__setitem__", &Confpool::__setitem__)
        .def("__getattr__", &Confpool::__getattr__)
        .def("__setattr__", &Confpool::__setattr__)
        .def("__delitem__", &Confpool::__delitem__)
        .def("__len__", &Confpool::__len__)
        .def("include_from_file", &Confpool::include_from_file)
        .def("include_from_xyz", &Confpool::include_from_xyz)
        .def("include_subset", &Confpool::include_subset)
        .def("filter", &Confpool::filter)
        .def("count", &Confpool::count)
        .def("upper_cutoff", &Confpool::upper_cutoff)
        .def("lower_cutoff", &Confpool::lower_cutoff)
        .def("sort", &Confpool::sort)
        .def("rmsd_filter", &Confpool::rmsd_filter)
        .def("get_rmsd_matrix", &Confpool::get_rmsd_matrix)
        .def("save", &Confpool::save)
        .def("as_table", &Confpool::as_table)
        ;
}
