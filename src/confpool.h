#ifndef CONFPOOL_H
#define CONFPOOL_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <iostream>
#include <vector>
#include <unordered_map>
#include "utils.h"
// #include "molproxy.h"

namespace py = pybind11;

class MolProxy;

class Confpool {
    private:
        using CoordContainerType = Utils::Coordinates<Utils::BoostXyzContainer>;
        using SymVector = std::vector<std::string>;
        SymVector sym_;

        std::vector<CoordContainerType> coord_;
        std::vector<std::string> descr_;
        std::vector<MolProxy> proxies_;
        std::unordered_map<std::string, std::vector<double>> keys_;
        unsigned int natoms;
        void include(const std::string& filename);
        void resize();
        void remove_structure(const int& idx);
        void full_check() const;

    public:
        Confpool() : natoms(0) {}
        py::int_ size() { return py::cast(coord_.size()); }
        py::list get_atom_symbols() { return py::cast(sym_); }
        MolProxy __getitem__(const py::int_& idx);

        inline std::vector<double> key_access(const std::string& key) const { return keys_.at(key); }
        inline py::str descr_access(const int& i) const { return descr_[i]; }
        inline CoordContainerType& coord_access(const int& i) { return coord_[i]; }
        inline const int get_natoms() const noexcept { return natoms; }
        
        void delete_by_idx(const int& idx) { remove_structure(idx); resize(); } // create just one 'delete' which takes py::object
        void delete_by_proxy(const MolProxy& mol);

        void include_from_file(const py::str& py_filename);
        void filter(const py::function& py_parser);
        void upper_cutoff(const py::str& py_keyname, const py::float_& py_cutoff);
        void lower_cutoff(const py::str& py_keyname, const py::float_& py_cutoff);
        py::int_ count(const py::function& py_criterion);
        void update_description(const py::function& descr_f);

        void key_from_description(const py::str& py_keyname, const py::function& py_parser);
        void distance_to_key(const py::str& py_keyname, const py::int_& py_idxA, const py::int_& py_idxB);
        void vangle_to_key(const py::str& py_keyname, const py::int_& py_idxA, const py::int_& py_idxB, const py::int_& py_idxC);
        void dihedral_to_key(const py::str& py_keyname, const py::int_& py_idxA, const py::int_& py_idxB, const py::int_& py_idxC, const py::int_& py_idxD);

        void sort(const py::str& py_keyname);
        void rmsd_filter(const py::float_& py_rmsd_cutoff);
        void save(const py::str& py_filename);
};


#endif // CONFPOOL_H