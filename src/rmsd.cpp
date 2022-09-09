#include <math.h>
#include <numeric>

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/type_index.hpp>

// #include <gsl/gsl_blas.h>
// #include <gsl/gsl_math.h>
// #include <gsl/gsl_eigen.h>
// #include <gsl/gsl_linalg.h>

#include "utils.h"
#include "rmsdpair.h"
#include "Hungarian.h"
#include "svd3.h"
#include <Eigen/Dense>

using Utils::repr_matrix;
using Utils::repr_vector;
using Utils::repr_sorted;
using Utils::repr_eigen_matrix;
using Permutation = boost::numeric::ublas::permutation_matrix<int>;


RmsdPair::RmsdPair(const py::list& py_xyzA, const py::list& py_xyzB, const py::list& py_atomsints, const py::bool_& py_hungarian_only) {
    p_all = Matrix(py_xyzA.size(), 3);
    q_all = Matrix(py_xyzB.size(), 3);
    Utils::fill_xyz_matrix(p_all, py_xyzA);
    Utils::fill_xyz_matrix(q_all, py_xyzB);

    hungarian_only = py_hungarian_only.cast<bool>();

    atom_ints = ElemVector(py_atomsints.size());
    for(int i = 0; i < py_atomsints.size(); ++i)
        atom_ints(i) = py_atomsints[i].cast<int>();
    
    unique_atoms = std::unordered_set<int>(atom_ints.begin(), atom_ints.end());

    using namespace boost::numeric::ublas;
    vector<int> xyz_view(3);
    xyz_view(0) = 0;
    xyz_view(1) = 1;
    xyz_view(2) = 2;
    coord_dir = IndArray(&(*xyz_view.begin()), &(*xyz_view.begin()) + xyz_view.size());
    for(const auto& element : unique_atoms) {
        auto elem_count = std::count(atom_ints.begin(), atom_ints.end(), element);
        vector<int> atom_view(elem_count);
        auto idx = 0;
        for (int i = 0; i < atom_ints.size(); ++i)
            if (atom_ints(i) == element)
                atom_view(idx++) = i;
        IndArray atom_dir(&(*atom_view.begin()), &(*atom_view.begin()) + atom_view.size());
        element_views.push_back(atom_dir);

        #ifdef RMSD_LOG
            atom_types.push_back(element);
        #endif
    }

    #ifdef RMSD_LOG
        log("[CTOR] NEW CALL");
        log(fmt::format("INIT:: xyzA = {}", repr_matrix(p_all)));
        log(fmt::format("INIT:: xyzB = {}", repr_matrix(q_all)));
        log(fmt::format("INIT:: sym = {}", repr_vector(atom_ints)));
    #endif
}

py::float_ RmsdPair::calc_rmsd() const {
    using namespace boost::numeric::ublas;
    const auto n_atoms = p_all.size1();
    if (q_all.size1() != n_atoms)
        throw std::runtime_error(fmt::format("Mismatch in number of atoms: {} vs. {}", p_all.size1(), q_all.size1()));

    Matrix p_coord(n_atoms, 3), q_coord(n_atoms, 3);
    p_coord = p_all;
    q_coord = q_all;
    
    if(!hungarian_only) {
        Utils::make_centered(p_coord);
        Utils::make_centered(q_coord);
        
        #ifdef RMSD_LOG
            log_rstart("center");
            log(fmt::format(" IN:: xyzA = {}", repr_matrix(p_all)));
            log(fmt::format(" IN:: xyzB = {}", repr_matrix(q_all)));
            log(fmt::format(" IN:: sym = {}", repr_vector(atom_ints)));
            log(fmt::format(" CHECK:: xyzA_coord = {}", repr_matrix(p_coord)));
            log(fmt::format(" CHECK:: xyzB_coord = {}", repr_matrix(q_coord)));
            log_rend("center");
        #endif
    }

    Permutation view_reorder(n_atoms);
    for(auto& item : view_reorder)
        item = -1;
    
    #ifdef RMSD_LOG
        log_rstart("hungarian_prep");
        log(fmt::format(" IN:: xyzA = {}", repr_matrix(p_all)));
        log(fmt::format(" IN:: xyzB = {}", repr_matrix(q_all)));
        log(fmt::format(" IN:: sym = {}", repr_vector(atom_ints)));
        log(fmt::format(" CHECK:: unique_atoms = {}", repr_sorted(unique_atoms)));
        log(fmt::format(" CHECK:: view_reorder = {}", repr_vector(view_reorder)));
        log_rend("hungarian_prep");
    #endif
    
    #ifdef RMSD_LOG
        int j = 0;
    #endif

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
        // VectorView<int> view_reorder_local(view_reorder, element_view);
        // for(int i = 0; i < n_elem; i++)
        //     view_reorder_local(i) = element_view(i);
        // Utils::reorder_non_destructive(assignment_temp.begin(), assignment_temp.end(), view_reorder_local.begin());

        #ifdef RMSD_LOG
            log_rstart("hungarian_step");
            log(fmt::format(" IN:: xyzA = {}", repr_matrix(p_all)));
            log(fmt::format(" IN:: xyzB = {}", repr_matrix(q_all)));
            log(fmt::format(" IN:: sym = {}", repr_vector(atom_ints)));
            log(fmt::format(" IN:: atom = {}", atom_types[j]));
            log(fmt::format(" CHECK:: atom_idx = {}", repr_vector(element_view)));
            log(fmt::format(" CHECK:: A_coord = {}", repr_matrix(A_coord)));
            log(fmt::format(" CHECK:: B_coord = {}", repr_matrix(B_coord)));
            log(fmt::format(" CHECK:: distances = {}", repr_matrix(distances)));
            log(fmt::format(" CHECK:: indices_b = {}", repr_vector(assignment_temp)));
            log(fmt::format(" CHECK:: atom_idx[indices_b] = {}", repr_vector(rhs_view)));
            // log(fmt::format(" CHECK:: view_reorder = {}", repr_vector(view_reorder)));
            log_rend("hungarian_step");
            j++;
        #endif
    }

    Matrix q_temp(q_coord);
    IndArray q_map(&(*view_reorder.begin()), &(*view_reorder.begin()) + view_reorder.size());
    MatrixView q_view(q_temp, q_map, coord_dir);
    for(int i = 0; i < n_atoms; i++)
        for(int j = 0; j < 3; j++)
            q_coord(i, j) = q_view(i, j);

    #ifdef RMSD_LOG
        log_rstart("hungarian_done");
        log(fmt::format(" IN:: xyzA = {}", repr_matrix(p_all)));
        log(fmt::format(" IN:: xyzB = {}", repr_matrix(q_all)));
        log(fmt::format(" IN:: sym = {}", repr_vector(atom_ints)));
        log(fmt::format(" CHECK:: view_reorder = {}", repr_vector(view_reorder)));
        log(fmt::format(" CHECK:: q_coord = {}", repr_matrix(q_coord)));
        log_rend("hungarian_done");
    #endif

    if(!hungarian_only) {
        FMatrix<3, 3> C, V, W;
        FVector<3> S;
        noalias(C) = prod(trans(p_coord), q_coord);
        // double   a11 = C(0, 0), a12 = C(0, 1), a13 = C(0, 2), 
        //          a21 = C(1, 0), a22 = C(1, 1), a23 = C(1, 2), 
        //          a31 = C(2, 0), a32 = C(2, 1), a33 = C(2, 2);
        // double   u11, u12, u13, 
        //          u21, u22, u23, 
        //          u31, u32, u33;
        // double   s11, s12, s13, 
        //          s21, s22, s23, 
        //          s31, s32, s33;
        // double   v11, v12, v13, 
        //          v21, v22, v23, 
        //          v31, v32, v33;
        // svd(a11, a12, a13, a21, a22, a23, a31, a32, a33,
        //     u11, u12, u13, u21, u22, u23, u31, u32, u33,
        //     s11, s12, s13, s21, s22, s23, s31, s32, s33,
        //     v11, v12, v13, v21, v22, v23, v31, v32, v33);
        
        // W(0, 0) = -v11; W(0, 1) = -v21; W(0, 2) = -v31;
        // W(1, 0) = v12; W(1, 1) = v22; W(1, 2) =   v32;
        // W(2, 0) =  v12; W(2, 1) = v23; W(2, 2) = -v33;
        // // fmt::print(" TRASH:: W = {}", repr_matrix(W));

        // V(0, 0) = -u11; V(0, 1) = u12; V(0, 2) = -u13;
        // V(1, 0) = -u21; V(1, 1) = u22; V(1, 2) = -u23;
        // V(2, 0) = -u31; V(2, 1) = u32; V(2, 2) = -u33;
        
        // S(0) = s11; S(1) = s22; S(2) = s33;
        
        Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> C_ei(&C(0, 0));
        fmt::print("\n\n EIGEN:: C_ei = {}\n", repr_eigen_matrix(C_ei));
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
            
        // fmt::print(" EIGEN:: W = {}\n", repr_matrix(W));
        // fmt::print(" EIGEN:: V = {}\n", repr_matrix(V));
        // fmt::print(" EIGEN:: S = {}\n", repr_vector(S));

        // gsl_matrix* c_m = gsl_matrix_alloc(3, 3);
        // gsl_matrix* v_m = gsl_matrix_alloc(3, 3);
        // gsl_vector* s_v = gsl_vector_alloc(3);
        // gsl_vector* work_v = gsl_vector_alloc(3);

        // for (unsigned int i = 0; i < 3; ++i)
        //     for (unsigned int j = 0; j < 3; ++j)
        //         gsl_matrix_set(c_m, i, j, C(i, j));

        // fmt::print(" GSL:: C = {}\n", Utils::repr_gsl_matrix(c_m));
        // gsl_linalg_SV_decomp(c_m, v_m, s_v, work_v);
        // fmt::print(" GSL:: V = {}\n", Utils::repr_gsl_matrix(v_m));
        // fmt::print(" GSL:: W = {}\n", Utils::repr_gsl_matrix(c_m));
        // fmt::print(" GSL:: S = {}\n", Utils::repr_gsl_vector(s_v));

        // for (unsigned int i = 0; i < 3; ++i)
        //     for (unsigned int j = 0; j < 3; ++j)
        //         W(i, j) = gsl_matrix_get(c_m, i, j);
        // W(1, 2) = -W(1, 2);
        // W(2, 0) = -W(2, 0);
        // W(2, 2) = -W(2, 2);

        // for (unsigned int i = 0; i < 3; ++i)
        //     for (unsigned int j = 0; j < 3; ++j)
        //         V(i, j) = gsl_matrix_get(v_m, i, j);
        // for (unsigned int i = 0; i < 3; ++i)
        //     V(i, 2) = -V(i, 2);
        
        // for (unsigned int i = 0; i < 3; ++i)
        //     S(i) = gsl_vector_get(s_v, i);

        // gsl_matrix_free(c_m);
        // gsl_matrix_free(v_m);
        // gsl_vector_free(s_v);
        // gsl_vector_free(work_v);
        // throw std::runtime_error(fmt::format("BAKA"));

        #ifdef RMSD_LOG
            log_rstart("opt_A");
            log(fmt::format(" IN:: xyzA = {}", repr_matrix(p_all)));
            log(fmt::format(" IN:: xyzB = {}", repr_matrix(q_all)));
            log(fmt::format(" IN:: sym = {}", repr_vector(atom_ints)));
            log(fmt::format(" CHECK:: C = {}", repr_matrix(C)));
            // log(fmt::format(" CHECK:: V = {}", repr_matrix(V)));
            log(fmt::format(" CHECK:: S = {}", repr_vector(S)));
            // log(fmt::format(" CHECK:: W = {}", repr_matrix(W)));
            // log(fmt::format(" CHECK:: det(V) = {}", Utils::det3x3(V)));
            // log(fmt::format(" CHECK:: det(W) = {}", Utils::det3x3(W)));
            log_rend("opt_A");
        #endif

        if ((Utils::det3x3(V) * Utils::det3x3(W)) < 0.0) {
            for(int i = 0; i < 3; i++)
                V(i, 2) = -V(i, 2);
        }
        FMatrix<3, 3> U;
        noalias(U) = prod(V, W);
        p_coord = prod(p_coord, U);

        #ifdef RMSD_LOG
            log_rstart("opt_done");
            log(fmt::format(" IN:: xyzA = {}", repr_matrix(p_all)));
            log(fmt::format(" IN:: xyzB = {}", repr_matrix(q_all)));
            log(fmt::format(" IN:: sym = {}", repr_vector(atom_ints)));
            // log(fmt::format(" CHECK:: S = {}", repr_vector(S)));
            // log(fmt::format(" CHECK:: V = {}", repr_matrix(V)));
            log(fmt::format(" CHECK:: U = {}", repr_matrix(U)));
            log(fmt::format(" CHECK:: p_coord = {}", repr_matrix(p_coord)));
            log_rend("opt_done");
        #endif
    }

    Matrix diff(p_coord);
    noalias(diff) -= q_coord;
    for (auto i = &diff(0, 0); i < &diff(n_atoms - 1, 2) + 1; ++i)
        *i *= *i;
    
    double res = sqrt(std::accumulate(&diff(0, 0), &diff(n_atoms - 1, 2) + 1, 0.0) / n_atoms);
    #ifdef RMSD_LOG
        log_rstart("rmsd_done");
        log(fmt::format(" IN:: xyzA = {}", repr_matrix(p_all)));
        log(fmt::format(" IN:: xyzB = {}", repr_matrix(q_all)));
        log(fmt::format(" IN:: sym = {}", repr_vector(atom_ints)));
        // log(fmt::format(" CHECK:: diff = {}", repr_matrix(diff)));
        log(fmt::format(" CHECK:: RMSD = {}", res));
        log_rend("rmsd_done");
    #endif
    return res;
}
