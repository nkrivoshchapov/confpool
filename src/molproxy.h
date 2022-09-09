#ifndef MOLPROXY_H
#define MOLPROXY_H

#include "confpool.h"

namespace py = pybind11;

class MolProxy {
    public:
        MolProxy() {}
        MolProxy(Confpool* base, const int& idx);
        py::float_ __getitem__(const py::str& key);
        void __setitem__(const py::str& key, const py::float_& value);
        py::str descr();
        py::array_t<double> xyz();

        py::float_ l(const py::int_& a, const py::int_& b);
        py::float_ v(const py::int_& a, const py::int_& b, const py::int_& c);
        py::float_ z(const py::int_& a, const py::int_& b, const py::int_& c, const py::int_& d);

        const int get_index() const { return idx_; }
    private:
        Confpool* base_;
        int idx_;
};


#endif // MOLPROXY_H
