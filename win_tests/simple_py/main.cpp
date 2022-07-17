#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <vector>
#include <math.h>
#include <numeric>
#include <fstream>

namespace py = pybind11;

py::object give_five() {
	return py::cast(5);
}

PYBIND11_MODULE(confpool, m) {
	m.def("give_five", &give_five);
}
