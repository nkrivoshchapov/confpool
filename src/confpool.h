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
        std::unordered_map<std::string, std::vector<double>> keys_;
        std::vector<MolProxy> proxies_;
        unsigned int natoms;
        int include(const std::string& filename);
        void resize();
        void remove_structure(const int& idx);
        void remove_key(const std::string& keyname);
        void full_check() const;
        void modify_key(const std::string& key_name, const py::function& func);
        void modify_descr(const py::function& func);
        inline py::list key_to_list(const std::string& key) const;
        inline py::list get_atom_symbols() const noexcept { return py::cast(sym_); }

    public:
        Confpool() : natoms(0) {}
        py::object __getitem__(const py::object& key);
        void __setitem__(const py::object& key, const py::object& func);
        py::object __getattr__(const py::str& py_attr);
        void __setattr__(const py::str& py_attr, const py::object& func);
        void __delitem__(const py::object& py_key);

        py::int_ include_from_file(const py::str& py_filename);
        py::int_ filter(const py::function& py_parser);
        py::int_ upper_cutoff(const py::str& py_keyname, const py::float_& py_cutoff);
        py::int_ lower_cutoff(const py::str& py_keyname, const py::float_& py_cutoff);
        py::int_ rmsd_filter(const py::float_& py_rmsd_cutoff);

        py::int_ count(const py::function& py_criterion);

        void sort(const py::str& py_keyname, const py::kwargs& kwargs);
        void save(const py::str& py_filename) const;
        py::dict as_table() const;

        inline std::vector<double>& key_access(const std::string& key) { return keys_.at(key); }
        inline std::string& descr_access(const int& i) { return descr_[i]; }
        inline CoordContainerType& coord_access(const int& i) { return coord_[i]; }
        inline const int get_natoms() const noexcept { return natoms; }
        
        inline void prepare_key(const std::string& keyname) noexcept {
            if (keys_.find(keyname) == keys_.end())
                keys_[keyname] = std::vector<double>(coord_.size());
        }
};


#endif // CONFPOOL_H