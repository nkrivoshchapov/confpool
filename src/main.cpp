#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <iostream>
#include <vector>
#include <math.h>
#include <numeric>

#include <boost/type_index.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "fmt/core.h"
#include "fmt/args.h"

#include <gsl/gsl_cblas.h>
#include <Eigen/Dense>

using boost::typeindex::type_id_with_cvr;
namespace py = pybind11;

# define DOUBLE_THRESHOLD 0.9999999999

namespace Utils {
    const double H2KC = 627.509474063;
    const double DEG2RAD = 0.0174532925199432957692;
    const double RAD2DEG = 1 / DEG2RAD;

    class VectorXyzContainer
    {
        private:
            typedef std::vector<std::vector<double>> CoordMatrix;
            unsigned int nrows;
            CoordMatrix xyz_;
        
        public:
            static constexpr int DIM = 3;

            VectorXyzContainer() {}
            VectorXyzContainer(decltype(nrows) natoms) : nrows(natoms), 
                                                         xyz_(CoordMatrix(natoms, std::vector<double>(DIM))) {}

            const std::vector<double> & get_atom(const int & index) const
            { return xyz_[index]; }

            void set_atom(const int & index, const std::vector<double> & new_coords)
            { xyz_[index] = new_coords; }

            inline std::vector<std::vector<double>>::iterator begin() noexcept { return xyz_.begin(); }
            inline std::vector<std::vector<double>>::iterator end() noexcept { return xyz_.end(); }
    };

    
    class Coordinates
    {
        private:
            VectorXyzContainer xyz_;

        public:
            Coordinates() {}
            Coordinates(const int &natoms) : xyz_(VectorXyzContainer(natoms)) {}

            decltype(auto) get_atom(const int & index) const
            { return xyz_.get_atom(index); }

            void set_atom(const int & index, std::initializer_list<double> new_coords)
            { xyz_.set_atom(index, new_coords); }

            inline decltype(auto) begin() noexcept { return xyz_.begin(); }
            inline decltype(auto) end() noexcept { return xyz_.end(); } // Use cbegin/cend ?

            void printout() // How to make it const???
            {
                for (auto xyz : xyz_)
                {
                    std::cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << "\n";
                }
            }
    };


    struct InvalidChar
    {
        bool operator()(char c) const {
            return !isprint(static_cast<unsigned char>(c));
        }
    };


    template <typename C, typename P> 
    void erase_remove_if(C& c, P predicate) {
        c.erase(std::remove_if(c.begin(), c.end(), predicate), c.end());
    }


    std::vector<std::string> readlines(const std::string& filename)
    {
        std::ifstream infile(filename);
        std::string str;
        std::vector<std::string> reslines;

        if(!infile)
            throw std::runtime_error("File " + filename + " not found");

        while (std::getline(infile, str))
        {
            if(str.size() > 0) {
                erase_remove_if(str, InvalidChar());
                reslines.push_back(str);
            }
        }
        return reslines;
    }

    double get_distance(const std::vector<double>& a_xyz, const std::vector<double>& b_xyz) {
        double res[3];
        cblas_dcopy(3, &a_xyz[0], 1, &res[0], 1);
        cblas_daxpy(3, -1.0, &b_xyz[0], 1, &res[0], 1);
        return cblas_dnrm2(3, &res[0], 1);
    }

    double get_vangle(const std::vector<double>& a_xyz, const std::vector<double>& b_xyz, const std::vector<double>& c_xyz) {
        double dirA[3];
        double dirC[3];
        cblas_dcopy(3, &a_xyz[0], 1, &dirA[0], 1);
        cblas_dcopy(3, &c_xyz[0], 1, &dirC[0], 1);
        cblas_daxpy(3, -1.0, &b_xyz[0], 1, &dirA[0], 1);
        cblas_daxpy(3, -1.0, &b_xyz[0], 1, &dirC[0], 1);
        
        double normA = cblas_dnrm2(3, &dirA[0], 1);
        double normC = cblas_dnrm2(3, &dirC[0], 1);
        cblas_dscal(3, 1 / normA, &dirA[0], 1);
        cblas_dscal(3, 1 / normC, &dirC[0], 1);
        return acos(cblas_ddot(3, &dirA[0], 1, &dirC[0], 1)) * RAD2DEG;
    }

    double get_dihedral(const std::vector<double>& a_xyz, const std::vector<double>& b_xyz, const std::vector<double>& c_xyz, const std::vector<double>& d_xyz) {
        using Eigen::Vector3d;
        using Eigen::Matrix3d;
        Eigen::Map<const Eigen::Vector3d> r_a(a_xyz.data());
        Eigen::Map<const Eigen::Vector3d> r_b(b_xyz.data());
        Eigen::Map<const Eigen::Vector3d> r_c(c_xyz.data());
        Eigen::Map<const Eigen::Vector3d> r_d(d_xyz.data());

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

    template< typename order_iterator, typename value_iterator >
    void reorder_non_destructive( order_iterator order_begin, order_iterator order_end, value_iterator v )  {   
        typedef typename std::iterator_traits< value_iterator >::value_type value_t;
        typedef typename std::iterator_traits< order_iterator >::value_type index_t;
        typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;
        
        diff_t remaining = order_end - 1 - order_begin;
        for ( index_t s = index_t(), d; remaining > 0; ++ s ) {
            for ( d = order_begin[s]; d > s; d = order_begin[d] ) ;
            if ( d == s ) {
                -- remaining;
                value_t temp = v[s];
                while ( d = order_begin[d], d != s ) {
                    std::swap( temp, v[d] );
                    -- remaining;
                }
                v[s] = temp;
            }
        }
    }

    template< typename order_iterator, typename value_iterator >
    void reorder_destructive( order_iterator order_begin, order_iterator order_end, value_iterator v )  {
        typedef typename std::iterator_traits< value_iterator >::value_type value_t;
        typedef typename std::iterator_traits< order_iterator >::value_type index_t;
        typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;
        
        diff_t remaining = order_end - 1 - order_begin;
        for ( index_t s = index_t(); remaining > 0; ++ s ) {
            index_t d = order_begin[s];
            if ( d == (diff_t) -1 ) continue;
            -- remaining;
            value_t temp = v[s];
            for ( index_t d2; d != s; d = d2 ) {
                std::swap( temp, v[d] );
                std::swap( order_begin[d], d2 = (diff_t) -1 );
                -- remaining;
            }
            v[s] = temp;
        }
    }
}


class Confpool {
    private:
        typedef std::vector<Utils::Coordinates> CoordVector;
        typedef std::vector<std::string> SymVector;
        typedef std::vector<double> EnerVector;
        typedef std::vector<std::string> DescrVector;
        CoordVector coord_;
        SymVector sym_;
        EnerVector ener_;
        DescrVector descr_;
        unsigned int natoms;

    public:
        Confpool() : natoms(0) {}

        void include_from_file(py::str& py_filename, const py::kwargs& kwargs) {
            const auto filename = py_filename.cast<std::string>();
            // std::cout << "The type of filename = "<< type_id_with_cvr<decltype(filename)>().pretty_name() << "\n";

            py::object energy_f = py::none();
            if (kwargs.attr("__contains__")("energy").cast<bool>())
                energy_f = kwargs["energy"].cast<py::object>();
            
            if (energy_f.is(py::none())) {
                if (ener_.size() > 0)
                    throw std::runtime_error("Cannot read without energies when structures with energies have already been included");
                _include(filename, false, energy_f);
            } else {
                if ((ener_.size() == 0) && (coord_.size() > 0))
                    throw std::runtime_error("Cannot read energies when structures without energies have already been included");
                _include(filename, true, energy_f);
            }
        }

        void _include(const std::string& filename, bool get_energy, const py::object& energy_obj) {
            auto ener_function = py::function();
            if (get_energy)
                ener_function = energy_obj.cast<py::function>();
            
            auto mylines = Utils::readlines(filename);
            int cline = 0;
            double energy = 0.0;
            while (cline < mylines.size()) {
                // std::cout << "Casting to int '" << mylines[cline] << "'" << "\n";
                boost::algorithm::trim(mylines[cline]);
                const unsigned int cur_natoms = boost::lexical_cast<int>(mylines[cline].c_str());
                if (natoms == 0) 
                    natoms = cur_natoms;
                else if (natoms != cur_natoms) 
                    throw std::runtime_error("Wrong numer of atoms");
                
                auto description = mylines[cline + 1];
                if (get_energy) {
                    energy = ener_function(py::cast(description)).cast<double>();
                }

                auto geom = Utils::Coordinates(natoms);
                SymVector atom_types;
                for (unsigned int i = 2; i < natoms + 2; ++i)
                {
                    if (cline + i >= mylines.size())
                        throw std::runtime_error("Unexpected number of atoms. Check " + filename);
                    
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
                if (get_energy)
                    ener_.push_back(energy);
                cline += 2 + natoms;
            }
        }

        void update_description(const py::function& descr_f) {
            if ((ener_.size() != 0) && (ener_.size() != coord_.size()))
                throw std::runtime_error(fmt::format("Energy list size is {} but the structures list size is different ({})", ener_.size(), coord_.size()));
            
            for(auto i = 0; i < coord_.size(); ++i) {
                py::object ener = py::none();
                if (ener_.size() == coord_.size())
                    ener = py::cast(ener_[i]);
                descr_[i] = descr_f(ener, py::cast(descr_[i])).cast<std::string>();
            }
        }

        void energy_filter(const py::float_& py_maxener, const py::kwargs& kwargs) {
            if (ener_.size() != coord_.size())
                throw std::runtime_error(fmt::format("Energy list size is {} but the structures list size is different ({})", ener_.size(), coord_.size()));

            const auto maxener = py_maxener.cast<double>();

            std::string etype = "kcal/mol";
            py::object energy_f = py::none();
            if (kwargs.attr("__contains__")("etype").cast<bool>())
                etype = kwargs["etype"].cast<std::string>();
            
            double mult;
            if (etype == "kcal/mol")
                mult = Utils::H2KC;
            else if ((etype == "a.u.") || (etype == "hartree"))
                mult = 1.0;
            else
                throw std::runtime_error("Unknown energy type - " + etype);
            
            double minener = ener_[0];
            for (const auto& myener : ener_)
                if (myener < minener)
                    minener = minener;

            unsigned int del_count = 0;
            for(int i = coord_.size() - 1; i >= 0; --i) {
                if ((ener_[i] - minener) * mult > maxener) {
                    remove_structure(i);
                    del_count += 1;
                }
            }
            fmt::print("Deleted {} structures\n", del_count);
        }

        void distance_filter(const py::int_& py_a_idx, const py::int_& py_b_idx, const py::function& dist_condition) {
            const int a_idx = py_a_idx.cast<int>() - 1;
            const int b_idx = py_b_idx.cast<int>() - 1;

            unsigned int del_count = 0;
            for(int i = coord_.size() - 1; i >= 0; --i) {
                const auto& a_xyz = coord_[i].get_atom(a_idx);
                const auto& b_xyz = coord_[i].get_atom(b_idx);
                auto py_dist = py::cast(Utils::get_distance(a_xyz, b_xyz));
                // descr_[i] = fmt::format("Dist = {}", Utils::get_distance(a_xyz, b_xyz));
                if (!(dist_condition(py_dist).cast<bool>())) {
                    remove_structure(i);
                    del_count += 1;
                }
            }
            fmt::print("Deleted {} structures\n", del_count);
        }

        void valence_filter(const py::int_& py_a_idx, const py::int_& py_b_idx, const py::int_& py_c_idx, const py::function& vangle_condition) {
            const int a_idx = py_a_idx.cast<int>() - 1;
            const int b_idx = py_b_idx.cast<int>() - 1;
            const int c_idx = py_c_idx.cast<int>() - 1;

            unsigned int del_count = 0;
            for(int i = coord_.size() - 1; i >= 0; --i) {
                const auto& a_xyz = coord_[i].get_atom(a_idx);
                const auto& b_xyz = coord_[i].get_atom(b_idx);
                const auto& c_xyz = coord_[i].get_atom(c_idx);
                auto py_dist = py::cast(Utils::get_vangle(a_xyz, b_xyz, c_xyz));
                // descr_[i] = fmt::format("Vangle = {}", Utils::get_vangle(a_xyz, b_xyz, c_xyz));
                if (!(vangle_condition(py_dist).cast<bool>())) {
                    remove_structure(i);
                    del_count += 1;
                }
            }
            fmt::print("Deleted {} structures\n", del_count);
        }

        void dihedral_filter(const py::int_& py_a_idx, const py::int_& py_b_idx, const py::int_& py_c_idx, const py::int_& py_d_idx, const py::function& dihedral_condition) {
            const int a_idx = py_a_idx.cast<int>() - 1;
            const int b_idx = py_b_idx.cast<int>() - 1;
            const int c_idx = py_c_idx.cast<int>() - 1;
            const int d_idx = py_d_idx.cast<int>() - 1;

            unsigned int del_count = 0;
            for(int i = coord_.size() - 1; i >= 0; --i) {
                const auto& a_xyz = coord_[i].get_atom(a_idx);
                const auto& b_xyz = coord_[i].get_atom(b_idx);
                const auto& c_xyz = coord_[i].get_atom(c_idx);
                const auto& d_xyz = coord_[i].get_atom(d_idx);
                auto py_dist = py::cast(Utils::get_dihedral(a_xyz, b_xyz, c_xyz, d_xyz));
                // descr_[i] = fmt::format("Dihedral = {}", Utils::get_dihedral(a_xyz, b_xyz, c_xyz, d_xyz));
                if (!(dihedral_condition(py_dist).cast<bool>())) {
                    remove_structure(i);
                    del_count += 1;
                }
            }
            fmt::print("Deleted {} structures\n", del_count);
        }

        py::int_ energy_count(const py::float_& py_maxener, const py::kwargs& kwargs) {
            if (ener_.size() != coord_.size())
                throw std::runtime_error(fmt::format("Energy list size is {} but the structures list size is different ({})", ener_.size(), coord_.size()));

            const auto maxener = py_maxener.cast<double>();

            std::string etype = "kcal/mol";
            py::object energy_f = py::none();
            if (kwargs.attr("__contains__")("etype").cast<bool>())
                etype = kwargs["etype"].cast<std::string>();
            
            double mult;
            if (etype == "kcal/mol")
                mult = Utils::H2KC;
            else if ((etype == "a.u.") || (etype == "hartree"))
                mult = 1.0;
            else
                throw std::runtime_error("Unknown energy type - " + etype);
            
            double minener = ener_[0];
            for (const auto& myener : ener_)
                if (myener < minener)
                    minener = minener;

            unsigned int count = 0;
            for(int i = coord_.size() - 1; i >= 0; --i) {
                if ((ener_[i] - minener) * mult > maxener)
                    count += 1;
            }
            return py::cast(count);
        }


        py::int_ distance_count(const py::int_& py_a_idx, const py::int_& py_b_idx, const py::function& dist_condition) {
            const int a_idx = py_a_idx.cast<int>() - 1;
            const int b_idx = py_b_idx.cast<int>() - 1;

            unsigned int count = 0;
            for(int i = coord_.size() - 1; i >= 0; --i) {
                const auto& a_xyz = coord_[i].get_atom(a_idx);
                const auto& b_xyz = coord_[i].get_atom(b_idx);
                auto py_dist = py::cast(Utils::get_distance(a_xyz, b_xyz));
                if (!(dist_condition(py_dist).cast<bool>()))
                    count += 1;
            }
            return py::cast(count);
        }

        py::int_ valence_count(const py::int_& py_a_idx, const py::int_& py_b_idx, const py::int_& py_c_idx, const py::function& vangle_condition) {
            const int a_idx = py_a_idx.cast<int>() - 1;
            const int b_idx = py_b_idx.cast<int>() - 1;
            const int c_idx = py_c_idx.cast<int>() - 1;

            unsigned int count = 0;
            for(int i = coord_.size() - 1; i >= 0; --i) {
                const auto& a_xyz = coord_[i].get_atom(a_idx);
                const auto& b_xyz = coord_[i].get_atom(b_idx);
                const auto& c_xyz = coord_[i].get_atom(c_idx);
                auto py_dist = py::cast(Utils::get_vangle(a_xyz, b_xyz, c_xyz));
                if (!(vangle_condition(py_dist).cast<bool>()))
                    count += 1;
            }
            return py::cast(count);
        }

        py::int_ dihedral_count(const py::int_& py_a_idx, const py::int_& py_b_idx, const py::int_& py_c_idx, const py::int_& py_d_idx, const py::function& dihedral_condition) {
            const int a_idx = py_a_idx.cast<int>() - 1;
            const int b_idx = py_b_idx.cast<int>() - 1;
            const int c_idx = py_c_idx.cast<int>() - 1;
            const int d_idx = py_d_idx.cast<int>() - 1;

            unsigned int count = 0;
            for(int i = coord_.size() - 1; i >= 0; --i) {
                const auto& a_xyz = coord_[i].get_atom(a_idx);
                const auto& b_xyz = coord_[i].get_atom(b_idx);
                const auto& c_xyz = coord_[i].get_atom(c_idx);
                const auto& d_xyz = coord_[i].get_atom(d_idx);
                auto py_dist = py::cast(Utils::get_dihedral(a_xyz, b_xyz, c_xyz, d_xyz));
                if (!(dihedral_condition(py_dist).cast<bool>()))
                    count += 1;
            }
            return py::cast(count);
        }

        void sort() {
            if (ener_.size() != coord_.size())
                throw std::runtime_error(fmt::format("Energy list size is {} but the structures list size is different ({})", ener_.size(), coord_.size()));

            std::vector<int> indices(coord_.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(),
                        [&](int A, int B) -> bool {
                                return ener_[A] < ener_[B];
                        });
            Utils::reorder_non_destructive(indices.cbegin(), indices.cend(), ener_.begin());
            Utils::reorder_non_destructive(indices.cbegin(), indices.cend(), descr_.begin());
            Utils::reorder_destructive(indices.begin(), indices.end(), coord_.begin()); // I am speeeed____
        }

        void remove_structure(const int& i) {
            if (ener_.size() > 0) {
                auto e_it = ener_.begin();
                std::advance(e_it, i);
                ener_.erase(e_it);
            }

            auto c_it = coord_.begin();
            std::advance(c_it, i);
            coord_.erase(c_it);

            auto d_it = descr_.begin();
            std::advance(d_it, i);
            descr_.erase(d_it);
        }

        void save(const py::str& py_filename) {
            const auto filename = py_filename.cast<std::string>();

            auto natoms_str = boost::lexical_cast<std::string>(natoms);
            std::vector<std::string> reslines;
            for(auto i = 0; i < coord_.size(); ++i) {
                reslines.push_back(natoms_str);
                reslines.push_back(descr_[i]);

                for(auto j = 0; j < natoms; ++j) {
                    const auto& coords = coord_[i].get_atom(j);
                    // std::cout << "The type of coords = "<< type_id_with_cvr<decltype(coords)>().pretty_name() << "\n";
                    reslines.push_back(fmt::format("{:>2}  {:12.8f}  {:12.8f}  {:12.8f}", sym_[j], coords[0], coords[1], coords[2]));
                }
            }

            auto joined = boost::algorithm::join(reslines, "\n");
            std::ofstream out(filename);
            out << joined << "\n\n\n";
            out.close();
        }

        py::int_ size() {
            return py::cast(coord_.size());
        }

        py::dict get_structure(const py::int_& py_index) {
            const auto index = py_index.cast<int>();
            if (index < 0) {
                throw std::runtime_error(fmt::format("Given index ({}) is less than zero", index));
            }
            else if (index >= coord_.size()) {
                throw std::runtime_error(fmt::format("Given index ({}) is bigger than the maximum list element ({})", index, coord_.size()));
            }

            py::dict res;

            if (index < ener_.size())
                res["energy"] = py::cast(ener_[index]);
            res["descr"] = py::cast(descr_[index]);

            const size_t size = natoms * 3;
            double *coord_array = new double[size];
            for (size_t i = 0; i < natoms; i++) {
                const auto& coords = coord_[index].get_atom(i);
                coord_array[i * 3 + 0] = coords[0];
                coord_array[i * 3 + 1] = coords[1];
                coord_array[i * 3 + 2] = coords[2];
            }

            py::capsule free_when_done(coord_array, [](void *f) {
                double *foo = reinterpret_cast<double *>(f);
                delete[] foo;
            });

            res["xyz"] = py::array_t<double>(
                {static_cast<int>(natoms), 3}, // shape
                {3*8, 8}, // C-style contiguous strides for double
                coord_array, // the data pointer
                free_when_done);
            return res;
        }

        py::list get_atom_symbols() {
            return py::cast(sym_);
        }
};


PYBIND11_MODULE(confpool, m) {
    py::class_<Confpool>(m, "Confpool")
        .def(py::init<>())
        .def("include_from_file", &Confpool::include_from_file)
        .def("update_description", &Confpool::update_description)
        .def("energy_filter", &Confpool::energy_filter)
        .def("distance_filter", &Confpool::distance_filter)
        .def("valence_filter", &Confpool::valence_filter)
        .def("dihedral_filter", &Confpool::dihedral_filter)
        .def("energy_count", &Confpool::energy_count)
        .def("distance_count", &Confpool::distance_count)
        .def("valence_count", &Confpool::valence_count)
        .def("dihedral_count", &Confpool::dihedral_count)
        .def("get_structure", &Confpool::get_structure)
        .def("get_atom_symbols", &Confpool::get_atom_symbols)
        .def("size", &Confpool::size)
        .def("sort", &Confpool::sort)
        .def("save", &Confpool::save);
}
