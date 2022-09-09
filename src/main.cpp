#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

#include "confpool.h"
#include "molproxy.h"

PYBIND11_MODULE(confpool, m) {
    py::class_<MolProxy>(m, "MolProxy")
        .def("__getitem__", &MolProxy::__getitem__)
        .def("__setitem__", &MolProxy::__setitem__)
        .def("descr", &MolProxy::descr)
        .def("xyz", &MolProxy::xyz)
        .def("l", &MolProxy::l)
        .def("v", &MolProxy::v)
        .def("z", &MolProxy::z)
        // .def("__setitem__", &Confpool::__setitem__)
        ;

    py::class_<Confpool>(m, "Confpool")
        .def(py::init<>())
        .def("__getitem__", &Confpool::__getitem__)
        .def("include_from_file", &Confpool::include_from_file)
        .def("key_from_description", &Confpool::key_from_description)
        .def("distance_to_key", &Confpool::distance_to_key)
        .def("vangle_to_key", &Confpool::vangle_to_key)
        .def("dihedral_to_key", &Confpool::dihedral_to_key)
        .def("filter", &Confpool::filter)
        .def("count", &Confpool::count)
        .def("upper_cutoff", &Confpool::upper_cutoff)
        .def("lower_cutoff", &Confpool::lower_cutoff)
        .def("sort", &Confpool::sort)
        .def("update_description", &Confpool::update_description)
        .def("rmsd_filter", &Confpool::rmsd_filter)
        // .def("get_structure", &Confpool::get_structure)
        .def("get_atom_symbols", &Confpool::get_atom_symbols)
        .def("size", &Confpool::size)
        .def("save", &Confpool::save)
        ;
}
