#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j + 1;
}

namespace py = pybind11;

PYBIND11_MODULE(confpool, m) {
    m.def("add", &add);

    m.def("subtract", [](int i, int j) { return i - j; });
}
