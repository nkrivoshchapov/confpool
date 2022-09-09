#include <iostream>
#include <fstream>
#include <math.h>
#include <numeric>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <Eigen/Dense>
#include <gsl/gsl_cblas.h>
#include "utils.h"

# define DOUBLE_THRESHOLD 0.9999999999
// #include <boost/type_index.hpp>
// using boost::typeindex::type_id_with_cvr;

namespace py = pybind11;


namespace Utils {
    std::vector<std::string> readlines(const std::string& filename) {
        std::ifstream infile(filename);
        std::string str;
        std::vector<std::string> reslines;

        if(!infile)
            throw std::runtime_error("File " + filename + " not found");

        while (std::getline(infile, str))
        {
            // if(str.size() > 0) {
			erase_remove_if(str, InvalidChar());
			reslines.push_back(str);
            // }
        }
        return reslines;
    }

    const std::vector<int> generate_atom_ints(std::vector<std::string> inp) {
        std::vector<int> res;
        res.reserve(inp.size());
        for(const auto& item : inp)
            res.push_back(NAMES_ELEMENT.at(item));
        return res;
    }

    std::string repr_matrix_buffer(const double* buf, const int& size1, const int& size2) {
        py::list res;
        for (unsigned int i = 0; i < size1; ++i) {
            py::list temp;
            for (unsigned int j = 0; j < size2; ++j)
                temp.append(buf[i * size2 + j]);
            res.append(temp);
        }
        return py::repr(res).cast<std::string>();
    }
}