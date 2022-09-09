#include "molproxy.h"
#include "confpool.h"

#include "fmt/core.h"
#include "fmt/args.h"

namespace py = pybind11;

MolProxy::MolProxy(Confpool* base, const int& idx) :
    base_(base), idx_(idx)
{}

py::float_ MolProxy::__getitem__(const py::str& key)
{ return base_->key_access(key)[idx_]; }

void MolProxy::__setitem__(const py::str& key, const py::float_& value)
{ base_->key_access(key)[idx_] = value; }

py::str MolProxy::descr()
{ return base_->descr_access(idx_); }

py::float_ MolProxy::l(const py::int_& a, const py::int_& b)
{ return base_->coord_access(idx_).get_distance(a - 1, b - 1); }

py::float_ MolProxy::v(const py::int_& a, const py::int_& b, const py::int_& c)
{ return base_->coord_access(idx_).get_vangle(a - 1, b - 1, c - 1); }

py::float_ MolProxy::z(const py::int_& a, const py::int_& b, const py::int_& c, const py::int_& d)
{ return base_->coord_access(idx_).get_dihedral(a - 1, b - 1, c - 1, d - 1); }

py::array_t<double> MolProxy::xyz() {
    const size_t size = base_->get_natoms() * 3;
    double *coord_array = new double[size];
    const auto& coords = base_->coord_access(idx_).to_boost_format();
    for (size_t i = 0; i < base_->get_natoms(); i++)
        for (size_t j = 0; j < 3; j++)
            coord_array[i * 3 + j] = coords(i, j);

    py::capsule free_when_done(coord_array, [](void *f) {
        double *foo = reinterpret_cast<double *>(f);
        delete[] foo;
    });

    return py::array_t<double>(
        {static_cast<int>(base_->get_natoms()), 3}, // shape
        {3*8, 8}, // C-style contiguous strides for double
        coord_array, // the data pointer
        free_when_done);
}
