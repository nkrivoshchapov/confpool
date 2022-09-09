#ifndef UTILS_H
#define UTILS_H

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

# define DOUBLE_THRESHOLD 0.9999999999
// #include <boost/type_index.hpp>
// using boost::typeindex::type_id_with_cvr;

namespace py = pybind11;

using BoostRow = boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double>>;

namespace Utils {
    const double H2KC = 627.509474063;
    const double DEG2RAD = 0.0174532925199432957692;
    const double RAD2DEG = 1 / DEG2RAD;

    template <typename C, typename P> 
    void erase_remove_if(C& c, P predicate) {
        c.erase(std::remove_if(c.begin(), c.end(), predicate), c.end());
    }

    std::vector<std::string> readlines(const std::string& filename);

    template<class T> double get_distance(const T&, const T&) noexcept;
    inline double get_distance_raw(const double* a_xyz, const double* b_xyz) noexcept {
        double res[3];
        cblas_dcopy(3, a_xyz, 1, &res[0], 1);
        cblas_daxpy(3, -1.0, b_xyz, 1, &res[0], 1);
        return cblas_dnrm2(3, &res[0], 1);
    }

    template<> 
    inline double get_distance<std::vector<double>>(const std::vector<double>& a_xyz, const std::vector<double>& b_xyz) noexcept {
        return get_distance_raw(&a_xyz[0], &b_xyz[0]);
    }
    
    template<> 
    inline double get_distance<BoostRow>(const BoostRow& a_xyz, const BoostRow& b_xyz) noexcept {
        return boost::numeric::ublas::norm_2(a_xyz - b_xyz);
    }

    template<class T> inline double get_vangle(const T&, const T&, const T&) noexcept;
    inline double get_vangle_raw(const double* a_xyz, const double* b_xyz, const double* c_xyz) noexcept {
        double dirA[3];
        double dirC[3];
        cblas_dcopy(3, a_xyz, 1, &dirA[0], 1);
        cblas_dcopy(3, c_xyz, 1, &dirC[0], 1);
        cblas_daxpy(3, -1.0, b_xyz, 1, &dirA[0], 1);
        cblas_daxpy(3, -1.0, b_xyz, 1, &dirC[0], 1);
        
        double normA = cblas_dnrm2(3, &dirA[0], 1);
        double normC = cblas_dnrm2(3, &dirC[0], 1);
        cblas_dscal(3, 1 / normA, &dirA[0], 1);
        cblas_dscal(3, 1 / normC, &dirC[0], 1);
        return acos(cblas_ddot(3, &dirA[0], 1, &dirC[0], 1)) * RAD2DEG;
    }

    template<>
    inline double get_vangle<std::vector<double>>(const std::vector<double>& a_xyz, const std::vector<double>& b_xyz, const std::vector<double>& c_xyz) noexcept {
        return get_vangle_raw(&a_xyz[0], &b_xyz[0], &c_xyz[0]);
    }

    // template<>
    // inline double get_vangle<BoostRow>(const BoostRow& a_xyz, const BoostRow& b_xyz, const BoostRow& c_xyz) noexcept {
    //     return 0.0;
    // }

    template<class T> inline double get_dihedral(const T&, const T&, const T&, const T&) noexcept;
    inline double get_dihedral_raw(const double* a_xyz, const double* b_xyz, const double* c_xyz, const double* d_xyz) noexcept {
        using Eigen::Vector3d;
        using Eigen::Matrix3d;
        Eigen::Map<const Eigen::Vector3d> r_a(a_xyz);
        Eigen::Map<const Eigen::Vector3d> r_b(b_xyz);
        Eigen::Map<const Eigen::Vector3d> r_c(c_xyz);
        Eigen::Map<const Eigen::Vector3d> r_d(d_xyz);

        Vector3d fr1_side = r_a - r_b;
        Vector3d fr1_mid = r_c - r_b;
        Vector3d fr2_mid = -fr1_mid;
        Vector3d fr2_side = r_d - r_c;
        fr1_side -= (fr1_side.dot(fr1_mid) / fr1_mid.squaredNorm()) * fr1_mid;
        fr2_side -= (fr2_side.dot(fr2_mid) / fr2_mid.squaredNorm()) * fr2_mid;
        fr1_side.normalize();
        fr2_side.normalize();

        auto dotprod = fr1_side.dot(fr2_side);
        if (dotprod >= 1.0)
            dotprod = DOUBLE_THRESHOLD;
        else if (dotprod <= -1.0)
            dotprod = -DOUBLE_THRESHOLD;
        auto ang = acos(dotprod);

        Matrix3d mymatr;
        mymatr << fr1_side, fr1_mid, fr2_side;
        if (mymatr.determinant() < 0)
            ang = -ang;
        return ang * RAD2DEG;
    }
    
    template<>
    inline double get_dihedral<std::vector<double>>(const std::vector<double>& a_xyz, const std::vector<double>& b_xyz, const std::vector<double>& c_xyz, const std::vector<double>& d_xyz) noexcept {
        return get_dihedral_raw(&a_xyz[0], &b_xyz[0], &c_xyz[0], &d_xyz[0]);
    }
    
    // template<>
    // inline double get_dihedral<BoostRow>(const BoostRow& a_xyz, const BoostRow& b_xyz, const BoostRow& c_xyz, const BoostRow& d_xyz) noexcept {
    //     return 0.0;
    // }

    template<class MType>
    void make_centered(MType& xyzs) {
        using boost::numeric::ublas::matrix_column;
        using boost::numeric::ublas::sum;
        auto size = xyzs.size1();
        for(int i = 0; i < 3; ++i) {
            matrix_column<MType> c(xyzs, i);
            double mean = sum(c)/size;
            for(auto& coord : c)
                coord -= mean;
        }
    }

    template<class MType>
    void fill_xyz_matrix(MType& xyz_matr, const py::list& py_xyzs) {
        for(int i = 0; i < xyz_matr.size1(); ++i) {
            auto xyz = py_xyzs[i].cast< std::array<double, 3> >();
            xyz_matr(i, 0) = xyz[0];
            xyz_matr(i, 1) = xyz[1];
            xyz_matr(i, 2) = xyz[2];
        }
    }

    template<class MType>
    inline double det3x3(const MType& m) {
        return m(0, 0)*m(1, 1)*m(2, 2)+
               m(1, 0)*m(2, 1)*m(0, 2)+
               m(2, 0)*m(1, 2)*m(0, 1)-
               m(2, 0)*m(1, 1)*m(0, 2)-
               m(1, 0)*m(0, 1)*m(2, 2)-
               m(2, 1)*m(1, 2)*m(0, 0);
    }

    template <typename T>
    std::vector<std::size_t> sort_permutation(
        const std::vector<T>& vec)
    {
        std::vector<std::size_t> p(vec.size());
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return vec[i] < vec[j]; });
        return p;
    }

    template <typename T>
    void apply_permutation_in_place(
        std::vector<T>& vec,
        const std::vector<std::size_t>& p)
    {
        std::vector<bool> done(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i)
        {
            if (done[i])
            {
                continue;
            }
            done[i] = true;
            std::size_t prev_j = i;
            std::size_t j = p[i];
            while (i != j)
            {
                std::swap(vec[prev_j], vec[j]);
                done[j] = true;
                prev_j = j;
                j = p[j];
            }
        }
    }

    class VectorXyzContainer {
        private:
            using CoordMatrix = std::vector<std::vector<double>>;
            unsigned int nrows;
            CoordMatrix xyz_;
        
        public:
            static constexpr int DIM = 3;

            VectorXyzContainer() {}
            VectorXyzContainer(decltype(nrows) natoms) : nrows(natoms), 
                                                         xyz_(CoordMatrix(natoms, std::vector<double>(DIM))) {}

            inline const std::vector<double>& get_atom(const int& index) const noexcept
            { return xyz_[index]; }

            inline const double* get_atom_raw(const int& index) const noexcept
            { return &(xyz_[index][0]); }

            inline void set_atom(const int& index, const std::vector<double>& new_coords) noexcept
            { xyz_[index] = new_coords; }

            inline double get_distance(const int& a_idx, const int& b_idx) const noexcept
            { return Utils::get_distance(get_atom(a_idx), get_atom(b_idx)); }
            
            inline double get_vangle(const int& a_idx, const int& b_idx, const int& c_idx) const noexcept
            { return Utils::get_vangle(get_atom(a_idx), get_atom(b_idx), get_atom(c_idx)); }
                        
            inline double get_dihedral(const int& a_idx, const int& b_idx, const int& c_idx, const int& d_idx) const noexcept
            { return Utils::get_dihedral(get_atom(a_idx), get_atom(b_idx), get_atom(c_idx), get_atom(d_idx)); }

            inline std::vector<std::vector<double>>::iterator begin() noexcept { return xyz_.begin(); }
            inline std::vector<std::vector<double>>::iterator end() noexcept { return xyz_.end(); }
            inline std::vector<std::vector<double>>::const_iterator cbegin() const noexcept { return xyz_.cbegin(); }
            inline std::vector<std::vector<double>>::const_iterator cend() const noexcept { return xyz_.cend(); }

            std::string repr() const noexcept {
                py::list res;
                for (auto it = xyz_.cbegin(); it < xyz_.cend(); it++) {
                    py::list temp;
                    for (unsigned int j = 0; j < 3; ++j)
                        temp.append(py::cast((*it)[j]));
                    res.append(temp);
                }
                return py::repr(res).cast<std::string>();
            }
    };
    
    class BoostXyzContainer {
        private:
            using CoordMatrix = boost::numeric::ublas::matrix<double>;
            using MatrixRow = boost::numeric::ublas::matrix_row<CoordMatrix>;
            CoordMatrix xyz_;
        
        public:
            static constexpr int DIM = 3;

            BoostXyzContainer() {}
            BoostXyzContainer(unsigned int natoms) : xyz_(CoordMatrix(natoms, DIM)) { }

            inline const MatrixRow get_atom(const int& index) noexcept
            { return MatrixRow(xyz_, index); }

            inline const double* get_atom_raw(const int& index) const noexcept
            { return &xyz_(index, 0); }

            void set_atom(const int& index, std::initializer_list<double> new_coords)
            { std::copy(new_coords.begin(), new_coords.end(), &xyz_(index, 0)); }

            inline double get_distance(const int& a_idx, const int& b_idx) noexcept
            { return Utils::get_distance(get_atom(a_idx), get_atom(b_idx)); }
            
            inline double get_vangle(const int& a_idx, const int& b_idx, const int& c_idx) const noexcept
            { return Utils::get_vangle_raw(get_atom_raw(a_idx), get_atom_raw(b_idx), get_atom_raw(c_idx)); }
                        
            inline double get_dihedral(const int& a_idx, const int& b_idx, const int& c_idx, const int& d_idx) const noexcept
            { return Utils::get_dihedral_raw(get_atom_raw(a_idx), get_atom_raw(b_idx), get_atom_raw(c_idx), get_atom_raw(d_idx)); }
            
            inline CoordMatrix& to_boost_format() noexcept 
            { return xyz_; }

            // inline std::vector<std::vector<double>>::iterator begin() noexcept { return xyz_.begin(); }
            // inline std::vector<std::vector<double>>::iterator end() noexcept { return xyz_.end(); }
            // inline std::vector<std::vector<double>>::const_iterator cbegin() const noexcept { return xyz_.cbegin(); }
            // inline std::vector<std::vector<double>>::const_iterator cend() const noexcept { return xyz_.cend(); }
    };

    
    template<class C>
    class Coordinates {
        using ContainerType = C;
        private:
            ContainerType xyz_;

        public:
            Coordinates() {}
            Coordinates(const int &natoms) : xyz_(ContainerType(natoms)) {}

            inline decltype(auto) get_atom(const int & index) noexcept
            { return xyz_.get_atom(index); }
            
            inline const double* get_atom_raw(const int & index) const noexcept
            { return xyz_.get_atom_raw(index); }
            
            void set_atom(const int & index, std::initializer_list<double> new_coords)
            { xyz_.set_atom(index, new_coords); }

            inline double get_distance(const int& a_idx, const int& b_idx) noexcept
            { return xyz_.get_distance(a_idx, b_idx); }
            
            inline double get_vangle(const int& a_idx, const int& b_idx, const int& c_idx) noexcept
            { return xyz_.get_vangle(a_idx, b_idx, c_idx); }
            
            inline double get_dihedral(const int& a_idx, const int& b_idx, const int& c_idx, const int& d_idx) noexcept
            { return xyz_.get_dihedral(a_idx, b_idx, c_idx, d_idx); }
            
            inline boost::numeric::ublas::matrix<double>& to_boost_format() noexcept 
            { return xyz_.to_boost_format(); }

            // inline decltype(auto) begin() noexcept { return xyz_.begin(); }
            // inline decltype(auto) end() noexcept { return xyz_.end(); }
            // inline decltype(auto) cbegin() const noexcept { return xyz_.cbegin(); }
            // inline decltype(auto) cend() const noexcept { return xyz_.cend(); }

            inline std::string repr() const noexcept { return xyz_.repr(); }
    };

    struct InvalidChar {
        bool operator()(char c) const {
            return !isprint(static_cast<unsigned char>(c));
        }
    };

    // template <typename MType> std::string repr_matrix(const MType& m);
    // template <typename MType> std::string repr_eigen_matrix(const MType& m);
    // template <typename MType> std::string repr_gsl_matrix(const MType* m);
    // template <typename MType> std::string repr_gsl_vector(const MType* m);
    // template <typename VType> std::string repr_vector(const VType& v);
    // template <typename Type> std::string repr_sorted(const Type& v);
    template <typename MType>
    std::string repr_matrix(const MType& m) {
        py::list res;
        for (unsigned int i = 0; i < m.size1 (); ++i) {
            py::list temp;
            for (unsigned int j = 0; j < m.size2 (); ++j)
                temp.append(py::cast(m(i, j)));
            res.append(temp);
        }
        return py::repr(res).cast<std::string>();
    }
    
    template <typename MType>
    std::string repr_eigen_matrix(const MType& m) {
        py::list res;
        for (unsigned int i = 0; i < m.rows (); ++i) {
            py::list temp;
            for (unsigned int j = 0; j < m.cols (); ++j)
                temp.append(py::cast(m(i, j)));
            res.append(temp);
        }
        return py::repr(res).cast<std::string>();
    }
    
    template <typename MType>
    std::string repr_gsl_matrix(const MType* m) {
        py::list res;
        for (unsigned int i = 0; i < m->size1; ++i) {
            py::list temp;
            for (unsigned int j = 0; j < m->size2; ++j)
                temp.append(py::cast(gsl_matrix_get(m, i, j)));
            res.append(temp);
        }
        return py::repr(res).cast<std::string>();
    }
    
    template <typename MType>
    std::string repr_gsl_vector(const MType* m) {
        py::list res;
        for (unsigned int i = 0; i < m->size; ++i)
            res.append(py::cast(gsl_vector_get(m, i)));
        return py::repr(res).cast<std::string>();
    }
    
    template <typename VType>
    std::string repr_vector(const VType& v) {
        py::list res;
        for (unsigned int i = 0; i < v.size(); ++ i)
            res.append(v[i]);
        return py::repr(res).cast<std::string>();
    }
    
    template <typename Type>
    std::string repr_sorted(const Type& v) {
        py::list res;
        for (const auto& item : v)
            res.append(item);
        res.attr("sort")();
        return py::repr(res).cast<std::string>();
    }
    std::string repr_matrix_buffer(const double* buf, const int& size1, const int& size2);


    const std::unordered_map<std::string, int> NAMES_ELEMENT = {
        {"H", 1},
        {"He", 2},
        {"Li", 3},
        {"Be", 4},
        {"B", 5},
        {"C", 6},
        {"N", 7},
        {"O", 8},
        {"F", 9},
        {"Ne", 10},
        {"Na", 11},
        {"Mg", 12},
        {"Al", 13},
        {"Si", 14},
        {"P", 15},
        {"S", 16},
        {"Cl", 17},
        {"Ar", 18},
        {"K", 19},
        {"Ca", 20},
        {"Sc", 21},
        {"Ti", 22},
        {"V", 23},
        {"Cr", 24},
        {"Mn", 25},
        {"Fe", 26},
        {"Co", 27},
        {"Ni", 28},
        {"Cu", 29},
        {"Zn", 30},
        {"Ga", 31},
        {"Ge", 32},
        {"As", 33},
        {"Se", 34},
        {"Br", 35},
        {"Kr", 36},
        {"Rb", 37},
        {"Sr", 38},
        {"Y", 39},
        {"Zr", 40},
        {"Nb", 41},
        {"Mo", 42},
        {"Tc", 43},
        {"Ru", 44},
        {"Rh", 45},
        {"Pd", 46},
        {"Ag", 47},
        {"Cd", 48},
        {"In", 49},
        {"Sn", 50},
        {"Sb", 51},
        {"Te", 52},
        {"I", 53},
        {"Xe", 54},
        {"Cs", 55},
        {"Ba", 56},
        {"La", 57},
        {"Ce", 58},
        {"Pr", 59},
        {"Nd", 60},
        {"Pm", 61},
        {"Sm", 62},
        {"Eu", 63},
        {"Gd", 64},
        {"Tb", 65},
        {"Dy", 66},
        {"Ho", 67},
        {"Er", 68},
        {"Tm", 69},
        {"Yb", 70},
        {"Lu", 71},
        {"Hf", 72},
        {"Ta", 73},
        {"W", 74},
        {"Re", 75},
        {"Os", 76},
        {"Ir", 77},
        {"Pt", 78},
        {"Au", 79},
        {"Hg", 80},
        {"Tl", 81},
        {"Pb", 82},
        {"Bi", 83},
        {"Po", 84},
        {"At", 85},
        {"Rn", 86},
        {"Fr", 87},
        {"Ra", 88},
        {"Ac", 89},
        {"Th", 90},
        {"Pa", 91},
        {"U", 92},
        {"Np", 93},
        {"Pu", 94},
        {"Am", 95},
        {"Cm", 96},
        {"Bk", 97},
        {"Cf", 98},
        {"Es", 99},
        {"Fm", 100},
        {"Md", 101},
        {"No", 102},
        {"Lr", 103},
        {"Rf", 104},
        {"Db", 105},
        {"Sg", 106},
        {"Bh", 107},
        {"Hs", 108},
        {"Mt", 109},
        {"Ds", 110},
        {"Rg", 111},
        {"Cn", 112},
        {"Uuq", 114},
        {"Uuh", 116}
    };

    const std::vector<int> generate_atom_ints(std::vector<std::string> inp);
}

#endif // UTILS_H
