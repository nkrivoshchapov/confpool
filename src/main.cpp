#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <vector>

#include <boost/type_index.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using boost::typeindex::type_id_with_cvr;
namespace py = pybind11;

void add(py::list x) {
    std::cout << "Argtype = " << type_id_with_cvr<decltype(x)>().pretty_name() << "\n";
    for(unsigned int i = 0; i < x.size(); ++i ){
        x[i] = 1;
        int xxx = x[i].cast<int>();
        std::cout << xxx << "\n";
    }
    x.append(10);
}

PYBIND11_MODULE(confpool, m) {
    m.def("add", &add);
    m.def("subtract", [](int i, int j) { return i - j; });
}
