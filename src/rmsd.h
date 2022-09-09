#ifndef RMSD_H
#define RMSD_H

#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "fmt/core.h"
#include "fmt/args.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "hungarian.h"
#include "utils.h"

namespace py = pybind11;

using Matrix = boost::numeric::ublas::matrix<double>;
template <std::size_t N, std::size_t M> using FMatrix = boost::numeric::ublas::fixed_matrix<double, N, M>;
template <std::size_t N> using FVector = boost::numeric::ublas::fixed_vector<double, N>;
using ElemVector = boost::numeric::ublas::vector<int>;
using XYZVector = boost::numeric::ublas::fixed_vector<int, 3>;
using MatrixSlice = boost::numeric::ublas::matrix_slice<Matrix>;

using IndArray = boost::numeric::ublas::indirect_array<boost::numeric::ublas::vector<int>>;
using Permutation = boost::numeric::ublas::permutation_matrix<int>;
using MatrixView = boost::numeric::ublas::matrix_indirect<Matrix, IndArray>;
template <class T> using VectorView = boost::numeric::ublas::vector_indirect<boost::numeric::ublas::vector<T>, IndArray>;

class RmsdCalculator {
    private:
        Matrix p_coord, q_coord;
        ElemVector atom_ints;
        std::unordered_set<int> unique_atoms;
        IndArray coord_dir;
        std::vector<IndArray> element_views;

    public:
        RmsdCalculator(const int& natoms, const std::vector<int>& src_atom_ints) {
            p_coord = Matrix(natoms, 3);
            q_coord = Matrix(natoms, 3);

            atom_ints = ElemVector(src_atom_ints.size());
            std::copy(src_atom_ints.cbegin(), src_atom_ints.cend(), &atom_ints[0]);
            
            unique_atoms = std::unordered_set<int>(atom_ints.begin(), atom_ints.end());

            using namespace boost::numeric::ublas;
            vector<int> xyz_view(3);
            xyz_view(0) = 0;
            xyz_view(1) = 1;
            xyz_view(2) = 2;
            coord_dir = IndArray(&(*xyz_view.begin()), &(*xyz_view.begin()) + xyz_view.size());
            element_views.reserve(unique_atoms.size());
            for(const auto& element : unique_atoms) {
                auto elem_count = std::count(atom_ints.begin(), atom_ints.end(), element);
                vector<int> atom_view(elem_count);
                auto idx = 0;
                for (int i = 0; i < atom_ints.size(); ++i)
                    if (atom_ints(i) == element)
                        atom_view(idx++) = i;
                IndArray atom_dir(&(*atom_view.begin()), &(*atom_view.begin()) + atom_view.size());
                element_views.push_back(atom_dir);
            }
        }

        const double calc(Matrix& p_all, Matrix& q_all) {
            using namespace boost::numeric::ublas;
            const auto n_atoms = p_all.size1();

            // std::cout << "P = " << Utils::repr_matrix(p_coord) << std::endl;
            // std::cout << "P = " << Utils::repr_matrix(q_coord) << std::endl;
            for(int i = 0; i < n_atoms; ++i) {
                for(int j = 0; j < 3; ++j) {
                    p_coord(i, j) = p_all(i, j);
                    q_coord(i, j) = q_all(i, j);
                }
            }

            Utils::make_centered(p_coord);
            Utils::make_centered(q_coord);

            Permutation view_reorder(n_atoms);
            for(auto& item : view_reorder)
                item = -1;

            for(const auto& element_view : element_views) {
                MatrixView A_coord(p_coord, element_view, coord_dir);
                MatrixView B_coord(q_coord, element_view, coord_dir);

                const int n_elem = element_view.size();
                Matrix distances (n_elem, n_elem);
                for(int i = 0; i < n_elem; ++i) {
                    matrix_row<MatrixView> a_xyz(A_coord, i);
                    for(int j = 0; j < n_elem; ++j) {
                        matrix_row<MatrixView> b_xyz(B_coord, j);
                        distances(i, j) = norm_2(b_xyz - a_xyz);
                    }
                }

                std::vector<std::vector<double>> costMatrix;
                for(int i = 0; i < n_elem; ++i) {
                    std::vector<double> temp;
                    matrix_row<MatrixView> a_xyz(A_coord, i);
                    for(int j = 0; j < n_elem; ++j)
                        temp.push_back(distances(i, j));
                    costMatrix.push_back(temp);
                }

                HungarianAlgorithm HungAlgo;
                std::vector<int> assignment_temp;
                HungAlgo.Solve(costMatrix, assignment_temp);
                VectorView<int> lhs_view(view_reorder, element_view);
                ElemVector rhs_vec(n_elem);
                for(int i = 0; i < n_elem; i++)
                    rhs_vec(i) = element_view(i);
                IndArray rhs_map(&(*assignment_temp.begin()), &(*assignment_temp.begin()) + assignment_temp.size());
                VectorView<int> rhs_view(rhs_vec, rhs_map);
                for(int i = 0; i < n_elem; i++)
                    lhs_view(i) = rhs_view(i);
            }

            Matrix q_temp(q_coord);
            IndArray q_map(&(*view_reorder.begin()), &(*view_reorder.begin()) + view_reorder.size());
            MatrixView q_view(q_temp, q_map, coord_dir);
            for(int i = 0; i < n_atoms; i++)
                for(int j = 0; j < 3; j++)
                    q_coord(i, j) = q_view(i, j);

            FMatrix<3, 3> C, V, W;
            FVector<3> S;
            noalias(C) = prod(trans(p_coord), q_coord);
            
            Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> C_ei(&C(0, 0));
            Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(C_ei);
            std::copy(svd.singularValues().cbegin(), svd.singularValues().cend(), S.begin());
            
            auto U_ei = svd.matrixU();
            using boost::typeindex::type_id_with_cvr;
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = 0; j < 3; ++j)
                    V(i, j) = U_ei(i, j);
            auto V_ei = svd.matrixV();
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = 0; j < 3; ++j)
                    W(i, j) = V_ei(j, i);

            if ((Utils::det3x3(V) * Utils::det3x3(W)) < 0.0) {
                for(int i = 0; i < 3; i++)
                    V(i, 2) = -V(i, 2);
            }
            FMatrix<3, 3> U;
            noalias(U) = prod(V, W);
            p_coord = prod(p_coord, U);

            Matrix diff(p_coord);
            noalias(diff) -= q_coord;
            for (auto i = &diff(0, 0); i < &diff(n_atoms - 1, 2) + 1; ++i)
                *i *= *i; // Nice
            
            double res = sqrt(std::accumulate(&diff(0, 0), &diff(n_atoms - 1, 2) + 1, 0.0) / n_atoms);
            // std::cout << "RMSD = " << res << std::endl;
            return res;
        }
};

#endif // RMSD_H
