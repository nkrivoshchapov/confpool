#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <iostream>
#include <vector>
#include <math.h>
#include <numeric>
#include <fstream>
#include <chrono>

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "fmt/core.h"
#include "fmt/args.h"

#include "confpool.h"
#include "molproxy.h"
#include "rmsd.h"
#include "utils.h"

namespace py = pybind11;


py::int_ Confpool::include_from_file(const py::str& py_filename) {
    const auto filename = py_filename.cast<std::string>();
    return include(filename);
}

int Confpool::include(const std::string& filename) {
    auto mylines = Utils::readlines(filename);
    // for (int i = 0; i < mylines.size(); ++i)
        // fmt::print("{} -- {}\n", i, mylines[i]);
    int cline = 0;
    int added_count = 0;
    while (cline < mylines.size()) {
        // std::cout << "Casting to int '" << mylines[cline] << "'" << "\n";
        boost::algorithm::trim(mylines[cline]);
        const unsigned int cur_natoms = boost::lexical_cast<int>(mylines[cline].c_str());
        if (natoms == 0) 
            natoms = cur_natoms;
        else if (natoms != cur_natoms) 
            throw std::runtime_error("Wrong numer of atoms");
        
        auto description = mylines[cline + 1];

        auto geom = CoordContainerType(natoms);
        SymVector atom_types;
        for (unsigned int i = 2; i < natoms + 2; ++i)
        {
            if (cline + i >= mylines.size())
                throw std::runtime_error(fmt::format("Unexpected number of atoms (expected {}, nlines={}). Check {}", natoms, mylines.size(), filename));
            
            std::vector<std::string> parts;
            boost::algorithm::trim(mylines[cline + i]);
            boost::split(parts, mylines[cline + i], boost::is_any_of(" "), boost::token_compress_on);
            if(parts.size() != 4)
                throw std::runtime_error("Unexpected number of parts in line. Check " + filename);
            
            atom_types.push_back(parts[0]);
            geom.set_atom(i - 2, {boost::lexical_cast<double>(parts[1].c_str()),
                                  boost::lexical_cast<double>(parts[2].c_str()),
                                  boost::lexical_cast<double>(parts[3].c_str())});
        }
        
        if (sym_.size() == 0)
            sym_ = atom_types;
        else if (sym_ != atom_types)
            throw std::runtime_error("Unexpected atom types. Check " + filename);
        
        coord_.push_back(geom);
        descr_.push_back(description);
        cline += 2 + natoms;
        added_count++;
    }
    resize();
    return added_count;
}

py::int_ Confpool::filter(const py::function& py_criterion) {
    full_check();

    unsigned int del_count = 0;
    for(int i = coord_.size() - 1; i >= 0; --i) {
        if (!py_criterion(proxies_[i]).cast<bool>()) {
            remove_structure(i);
            del_count += 1;
        }
    }
    resize();
    return del_count;
}

py::int_ Confpool::count(const py::function& py_criterion) {
    full_check();

    unsigned int res = 0;
    for(int i = coord_.size() - 1; i >= 0; --i)
        if (py_criterion(proxies_[i]).cast<bool>())
            res += 1;
    return res;
}

void Confpool::remove_key(const std::string& keyname) {
    keys_.erase(keyname);
}

void Confpool::remove_structure(const int& i) {
    auto c_it = coord_.begin();
    std::advance(c_it, i);
    coord_.erase(c_it);

    auto d_it = descr_.begin();
    std::advance(d_it, i);
    descr_.erase(d_it);

    for(auto& pair : keys_) { // pair = key, std::vector<double>
        auto it = pair.second.begin();
        std::advance(it, i);
        pair.second.erase(it);
    }
}

void Confpool::resize() {
    if (coord_.size() != descr_.size())
        throw std::runtime_error(fmt::format("Mismatch of container sizes (coord and descr): {} vs. {}", coord_.size(), descr_.size()));
    
    for(auto& pair : keys_) { // pair = key, std::vector<double>
        if (pair.second.size() > coord_.size())
            throw std::runtime_error(fmt::format("Mismatch of container sizes (key container > coord): {} < {}", pair.second.size(), coord_.size()));
        else if (pair.second.size() < coord_.size())
            pair.second.resize(coord_.size());
    }
    
    if (proxies_.size() != coord_.size()) {
        proxies_.resize(coord_.size());
        for (int i = 0; i < coord_.size(); ++i)
            proxies_[i] = MolProxy(this, i);
    }
}

py::int_ Confpool::upper_cutoff(const py::str& py_keyname, const py::float_& py_cutoff) {
    full_check();
    const auto cutoff = py_cutoff.cast<double>();
    if (cutoff <= 0.0)
        throw std::runtime_error(fmt::format("Cutoff value must be > 0. {} given.", cutoff));

    const auto& key_data = keys_[py_keyname.cast<std::string>()];

    auto minimal_value = key_data[0];
    for (const auto& current_val : key_data)
        if (current_val < minimal_value)
            minimal_value = current_val;
    
    unsigned int del_count = 0;
    for(int i = coord_.size() - 1; i >= 0; --i) {
        if (key_data[i] - minimal_value > cutoff) {
            remove_structure(i);
            del_count += 1;
        }
    }
    resize();
    return del_count;
}

py::int_ Confpool::lower_cutoff(const py::str& py_keyname, const py::float_& py_cutoff) {
    full_check();
    const auto cutoff = py_cutoff.cast<double>();
    if (cutoff <= 0.0)
        throw std::runtime_error(fmt::format("Cutoff value must be > 0. {} given.", cutoff));

    const auto& key_data = keys_[py_keyname.cast<std::string>()];

    auto maximal_value = key_data[0];
    for (const auto& current_val : key_data)
        if (maximal_value < current_val)
            maximal_value = current_val;
    
    unsigned int del_count = 0;
    for(int i = coord_.size() - 1; i >= 0; --i) {
        if (maximal_value - key_data[i] > cutoff) {
            remove_structure(i);
            del_count += 1;
        }
    }
    resize();
    return del_count;
}

void Confpool::save(const py::str& py_filename) const {
    full_check();
    const auto filename = py_filename.cast<std::string>();

    auto natoms_str = boost::lexical_cast<std::string>(natoms);
    std::vector<std::string> reslines;
    for(auto i = 0; i < coord_.size(); ++i) {
        reslines.push_back(natoms_str);
        reslines.push_back(descr_[i]);

        for(auto j = 0; j < natoms; ++j) {
            const auto* coords = coord_[i].get_atom_raw(j);
            reslines.push_back(fmt::format("{:>2}  {:12.8f}  {:12.8f}  {:12.8f}", sym_[j], coords[0], coords[1], coords[2]));
        }
    }

    auto joined = boost::algorithm::join(reslines, "\n");
    std::ofstream out(filename);
    out << joined << "\n\n\n";
    out.close();
}

py::dict Confpool::as_table() const {
    py::dict res;
    for(const auto& pair : keys_) { // pair = key, std::vector<double>
        res[py::cast(pair.first)] = py::list(py::cast(pair.second));
    }
    return res;
}

void Confpool::full_check() const {
    if (coord_.size() != descr_.size())
        throw std::runtime_error(fmt::format("Mismatch of container sizes (coord and descr): {} vs. {}", coord_.size(), descr_.size()));
    
    for(auto& pair : keys_) { // pair = key, std::vector<double>
        if (pair.second.size() != coord_.size())
            throw std::runtime_error(fmt::format("Mismatch of container sizes (key='{}' container != coord): {} != {}", pair.first, pair.second.size(), coord_.size()));
    }

    if (coord_.size() != proxies_.size())
        throw std::runtime_error(fmt::format("Mismatch of container sizes (coord and proxy container): {} vs. {}", coord_.size(), proxies_.size()));
    
    for (int i = 0; i < coord_.size(); ++i)
        if (proxies_[i].get_index() != i)
            throw std::runtime_error(fmt::format("MolProxy #{} has number {} (must be equal)", i, proxies_[i].get_index()));
}

void Confpool::sort(const py::str& py_keyname, const py::kwargs& kwargs) {
    full_check();
    
    bool ascend = true;
    if (kwargs.attr("__contains__")("ascending").cast<bool>())
        ascend = kwargs["ascending"].cast<bool>();

    const auto keyname = py_keyname.cast<std::string>();
    auto p = Utils::sort_permutation(keys_[keyname], ascend);

    Utils::apply_permutation_in_place(descr_, p);
    Utils::apply_permutation_in_place(coord_, p);
    for(auto& pair : keys_) // pair = key, std::vector<double>
        Utils::apply_permutation_in_place(pair.second, p);
}

// py::int_ Confpool::rmsd_filter(const py::float_& py_rmsd_cutoff) {
//     full_check();
//     using namespace std::chrono;
//     auto start_time = high_resolution_clock::now();

//     unsigned int del_count = 0;
//     const auto cutoff = py_rmsd_cutoff.cast<double>();
//     const auto atom_ints = Utils::generate_atom_ints(sym_);
//     RmsdCalculator rmsd(natoms, atom_ints);
//     int rmsd_calc_count = 0;
//     for (int i = coord_.size() - 1; i > 0; i--) {
//         auto& curgeom = coord_[i].to_boost_format();
//         for (int j = i - 1; j >= 0; j--) {
//             auto& testgeom = coord_[j].to_boost_format();
//             rmsd_calc_count++;
//             if (rmsd.calc(curgeom, testgeom) < cutoff) {
//                 remove_structure(i);
//                 del_count += 1;
//                 break;
//             }
//         }
//     }
//     auto stop_time = high_resolution_clock::now();
//     duration<double, std::milli> fp_ms = stop_time - start_time;
//     fmt::print("Time elapsed = {} ms\nPer RMSD calc = {} ms\n", fp_ms.count(), fp_ms.count() / rmsd_calc_count);

//     resize();
//     return del_count;
// }

py::dict Confpool::rmsd_filter(const py::float_& py_rmsd_cutoff) {
    full_check();
    using namespace std::chrono;
    auto start_time = high_resolution_clock::now();

    unsigned int del_count = 0;
    const auto cutoff = py_rmsd_cutoff.cast<double>();
    const auto atom_ints = Utils::generate_atom_ints(sym_);
    RmsdCalculator rmsd(natoms, atom_ints);
    double minimal_rmsd = 100.0; // FIX ME
    int min_pairA = -1, min_pairB = -1;
    int rmsd_calc_count = 0;
    for (int i = 0; i < coord_.size(); i++) {
        double cur_min_rmsd = minimal_rmsd;
        int cur_min_pair = -1;
        auto& curgeom = coord_[i].to_boost_format();
        bool removed = false;
        for (int j = i - 1; j >= 0; j--) {
            auto& testgeom = coord_[j].to_boost_format();
            auto rmsd_value = rmsd.calc(curgeom, testgeom);
            rmsd_calc_count++;
            if (rmsd_value < cutoff) {
                remove_structure(i);
                del_count += 1;
                removed = true;
                i--;
                break;
            } else if (rmsd_value < cur_min_rmsd) {
                cur_min_pair = j;
                cur_min_rmsd = rmsd_value;
            }
        }

        if (!removed && (cur_min_rmsd < minimal_rmsd)) {
            if (cur_min_pair == -1)
                throw std::runtime_error("A bug detected. 306");
            min_pairA = i;
            min_pairB = cur_min_pair;
            minimal_rmsd = cur_min_rmsd;
        }
    }
    auto stop_time = high_resolution_clock::now();
    duration<double, std::milli> fp_ms = stop_time - start_time;
    fmt::print("RMSD calc statistics:\nTime elapsed = {} ms\nPer RMSD calc = {} ms\n", fp_ms.count(), fp_ms.count() / rmsd_calc_count);

    if (min_pairA == -1)
        throw std::runtime_error("A bug detected. 314");
    if (min_pairB == -1)
        throw std::runtime_error("A bug detected. 316");
    py::dict res;
    res["DelCount"] = del_count;
    res["MinRMSD_pairA"] = min_pairA;
    res["MinRMSD_pairB"] = min_pairB;
    res["MinRMSD"] = minimal_rmsd;
    resize();
    return res;
}

py::object Confpool::__getitem__(const py::object& key) {
    if (py::isinstance<py::int_>(key))
        return py::cast(proxies_[key.cast<int>()]);
    else if (py::isinstance<py::str>(key))
        return key_to_list(key.cast<std::string>());
    else
        throw std::runtime_error(fmt::format("Expected either an integer (conformer index) or a string (a key). Got a {}", py::repr(key).cast<std::string>()));
}

void Confpool::__setitem__(const py::object& key, const py::object& func) {
    if (!py::isinstance<py::str>(key))
            throw std::runtime_error(fmt::format("Expected str as key. Got a {}", py::repr(key).cast<std::string>()));
    if (!py::isinstance<py::function>(func))
            throw std::runtime_error(fmt::format("Expected a function for keyvalue modification. Got a {}", py::repr(func).cast<std::string>()));
    modify_key(key.cast<std::string>(), func.cast<py::function>());
}

void Confpool::modify_key(const std::string& keyname, const py::function& func) {
    full_check();
    prepare_key(keyname);
    for(int i = 0; i < coord_.size(); ++i)
        keys_[keyname][i] = func(proxies_[i]).cast<double>();
}

void Confpool::modify_descr(const py::function& func) {
    full_check();
    for(int i = 0; i < coord_.size(); ++i)
        descr_[i] = func(proxies_[i]).cast<std::string>();
}

py::object Confpool::__getattr__(const py::str& py_attr) {
    const auto attr = py_attr.cast<std::string>();
    if (attr == "size")
        return py::cast(coord_.size());
    else if (attr == "atom_symbols")
        return get_atom_symbols();
    else
        throw std::runtime_error(fmt::format("Unknown attr {}", attr));
}

py::int_ Confpool::__len__() {
    return py::cast(coord_.size());
}

void Confpool::__setattr__(const py::str& py_attr, const py::object& value) {
    const auto attr = py_attr.cast<std::string>();
    if (attr == "descr") {
        if (!py::isinstance<py::function>(value))
            throw std::runtime_error(fmt::format("Expected a function for description modification. Got a {}", py::repr(value).cast<std::string>()));
        modify_descr(value.cast<py::function>());
    } else {
        throw std::runtime_error(fmt::format("Unknown attr {}", attr));
    }
}

void Confpool::__delitem__(const py::object& py_key) {
    if (py::isinstance<py::int_>(py_key)) {
        remove_structure(py_key.cast<int>());
        resize();
    } else if (py::isinstance<py::str>(py_key)) {
        remove_key(py_key.cast<std::string>());
    } else {
        throw std::runtime_error(fmt::format("Expected either keyname or structure index for deletion. Got {}", py::repr(py_key).cast<std::string>()));
    }
    full_check();
}

inline py::list Confpool::key_to_list(const std::string& key) const {
    py::list res;
    for (const auto& item : keys_.at(key))
        res.append(item);
    return res;
}

py::array_t<double> Confpool::get_rmsd_matrix() {
    using namespace std::chrono;
    auto start_time = high_resolution_clock::now();

    int nstructs = coord_.size();
    const size_t size = nstructs * nstructs;
    double *coord_array = new double[size];
    
    int rmsd_calc_count = 0;
    const auto atom_ints = Utils::generate_atom_ints(sym_);
    RmsdCalculator rmsd(natoms, atom_ints);
    for (size_t i = 0; i < nstructs; i++) {
        auto& curgeom = coord_[i].to_boost_format();
        for (size_t j = 0; j < i; j++) {
            auto& testgeom = coord_[j].to_boost_format();
            auto rmsd_value = rmsd.calc(curgeom, testgeom);
            coord_array[i * nstructs + j] = rmsd_value;
            coord_array[j * nstructs + i] = rmsd_value;
            rmsd_calc_count++;
        }
    }
    auto stop_time = high_resolution_clock::now();
    duration<double, std::milli> fp_ms = stop_time - start_time;
    fmt::print("RMSD matrix statistics:\nTime elapsed = {} ms\nPer RMSD calc = {} ms\n", fp_ms.count(), fp_ms.count() / rmsd_calc_count);
    
    py::capsule free_when_done(coord_array, [](void *f) {
        double *foo = reinterpret_cast<double *>(f);
        delete[] foo;
    });
    
    return py::array_t<double>(
        {nstructs, nstructs}, // shape
        {nstructs*8, 8}, // C-style contiguous strides for double
        coord_array, // the data pointer
        free_when_done);
}

void Confpool::include_from_xyz(const py::array_t<double>& xyz, const py::str& descr) {
    if (natoms == 0)
        throw std::runtime_error("The object is unintialized");
    
    py::buffer_info input_buf = xyz.request();
    if (input_buf.ndim != 2)
        throw std::runtime_error("numpy.ndarray dims must be 2!");

    if ((input_buf.shape[0] != natoms) || (input_buf.shape[1] != 3))
        throw std::runtime_error(fmt::format("Expected dimensions ({}, {}). Got ({}, {})", natoms, 3, input_buf.shape[0], input_buf.shape[1]));

    full_check();
    auto geom = CoordContainerType(natoms);
    double* ptr1 = (double*)input_buf.ptr;
    for (int i = 0; i < input_buf.shape[0]; i++)
    {
        geom.set_atom(i, {ptr1[i * input_buf.shape[1]],
                          ptr1[i * input_buf.shape[1] + 1],
                          ptr1[i * input_buf.shape[1] + 2]});
    }
    coord_.push_back(geom);
    descr_.push_back(descr.cast<std::string>());
    resize();
}

void Confpool::include_subset(Confpool& other, const py::list& py_idxs) {
    full_check();
    auto idxs = py_idxs.cast<std::vector<int>>();
    if (natoms == 0) {
        natoms = other.get_natoms();
        sym_ = other.symbols_access();
    }

    for(const auto& pair : other.key_full_access())
        prepare_key(pair.first);

    for (const auto& idx : idxs) {
        coord_.push_back(other.coord_access(idx));
        descr_.push_back(other.descr_access(idx));
    }
    
    for (const auto& idx : idxs)
        for(const auto& pair : other.key_full_access())
            keys_[pair.first].push_back(pair.second[idx]);
    resize();
}
