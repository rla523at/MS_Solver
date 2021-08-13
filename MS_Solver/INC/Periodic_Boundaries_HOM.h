#pragma once
#include "Grid_Builder.h"
#include "Reconstruction_Method.h"
#include "Numerical_Flux_Function.h"


//HOM이면 공통으로 사용하는 variable & method
template<typename Reconstruction_Method, typename Numerical_Flux_Function>
class Periodic_Boundaries_HOM
{
private:
    static constexpr ushort space_dimension_    = Reconstruction_Method::space_dimension();
    static constexpr ushort num_basis_          = Reconstruction_Method::num_basis();
    static constexpr ushort num_equation_       = Numerical_Flux_Function::num_equation();

    using This_                 = Periodic_Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>;
    using Space_Vector_         = Euclidean_Vector<space_dimension_>;
    using Solution_Coefficient_ = Matrix<num_equation_, num_basis_>;
    using Residual_             = Matrix<num_equation_, num_basis_>;

protected:
    inline static std::vector<std::pair<uint, uint>> oc_nc_index_pairs_;
    inline static std::vector<std::pair<Dynamic_Matrix, Dynamic_Matrix>> oc_nc_side_basis_qnode_pairs_;
    inline static std::vector<std::vector<Space_Vector_>> set_of_normals_;
    inline static std::vector<std::pair<Dynamic_Matrix, Dynamic_Matrix>> oc_nc_side_basis_weight_pairs_;

private:
    Periodic_Boundaries_HOM(void) = delete;

public:
    static void initialize(Grid<space_dimension_>&& grid);
    static void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients);
};


// template definition part
template<typename Reconstruction_Method, typename Numerical_Flux_Function>
void Periodic_Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>::initialize(Grid<space_dimension_>&& grid) {
    SET_TIME_POINT;

    This_::oc_nc_index_pairs_ = std::move(grid.connectivity.periodic_boundary_oc_nc_index_pairs);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& periodic_boundary_element_pairs = grid.elements.periodic_boundary_element_pairs;

    const auto num_periodic_pair = periodic_boundary_element_pairs.size();
    This_::oc_nc_side_basis_qnode_pairs_.reserve(num_periodic_pair);
    This_::set_of_normals_.reserve(num_periodic_pair);
    This_::oc_nc_side_basis_weight_pairs_.reserve(num_periodic_pair);

    constexpr auto integrand_order = 2 * Reconstruction_Method::solution_order() + 1;

    for (uint i = 0; i < num_periodic_pair; ++i) {
        const auto& [oc_side_element, nc_side_element] = periodic_boundary_element_pairs[i];
        const auto& oc_side_geometry = oc_side_element.geometry_;
        const auto& nc_side_geometry = nc_side_element.geometry_;


        const auto [oc_index, nc_index] = This_::oc_nc_index_pairs_[i];
        const auto& oc_element = cell_elements[oc_index];

        const auto& oc_side_quadrature_rule = oc_side_geometry.get_quadrature_rule(integrand_order);
        const auto& oc_side_qnodes = oc_side_quadrature_rule.points;
        const auto& oc_side_qweights = oc_side_quadrature_rule.weights;

        const auto& nc_side_quadrature_rule = nc_side_geometry.get_quadrature_rule(integrand_order);
        const auto& nc_side_qnodes = nc_side_quadrature_rule.points;
        const auto& nc_side_qweights = nc_side_quadrature_rule.weights;


        const auto num_qnode = oc_side_qnodes.size();

        auto oc_side_basis_qnode = Reconstruction_Method::calculate_basis_nodes(oc_index, oc_side_qnodes);
        auto nc_side_basis_qnode = Reconstruction_Method::calculate_basis_nodes(nc_index, nc_side_qnodes);
        This_::oc_nc_side_basis_qnode_pairs_.push_back({ std::move(oc_side_basis_qnode), std::move(nc_side_basis_qnode) });

        std::vector<Space_Vector_> normals(num_qnode);
        Dynamic_Matrix oc_side_basis_weight(num_qnode, This_::num_basis_);
        Dynamic_Matrix nc_side_basis_weight(num_qnode, This_::num_basis_);

        for (ushort q = 0; q < num_qnode; ++q) {
            normals[q] = oc_side_element.normalized_normal_vector(oc_element, oc_side_qnodes[q]);
            //oc_side_basis_weight.change_row(q, Reconstruction_Method::calculate_basis_node(oc_index, oc_side_qnodes[q]) * oc_side_qweights[q] * -1);    //performance test해서 더 빠른거로 결정하기
            oc_side_basis_weight.change_row(q, Reconstruction_Method::calculate_basis_node(oc_index, oc_side_qnodes[q]) * oc_side_qweights[q]);    
            nc_side_basis_weight.change_row(q, Reconstruction_Method::calculate_basis_node(nc_index, nc_side_qnodes[q]) * nc_side_qweights[q]);         //수정 예정
        }

        This_::set_of_normals_.push_back(std::move(normals));
        This_::oc_nc_side_basis_weight_pairs_.push_back({ std::move(oc_side_basis_weight), std::move(nc_side_basis_weight) });
    }

    Log::content_ << std::left << std::setw(50) << "@ Periodic Boundaries HOM precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template<typename Reconstruction_Method, typename Numerical_Flux_Function>
void Periodic_Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) {
    const auto num_periodic_pairs = This_::oc_nc_index_pairs_.size();

    for (uint i = 0; i < num_periodic_pairs; ++i) {
        const auto [oc_index, nc_index] = This_::oc_nc_index_pairs_[i];

        const auto& oc_solution_coefficient = solution_coefficients[oc_index];
        const auto& nc_solution_coefficient = solution_coefficients[nc_index];

        const auto& [oc_side_basis_qnode, nc_side_basis_qnode] = This_::oc_nc_side_basis_qnode_pairs_[i];
        const auto oc_side_cvariables = oc_solution_coefficient * oc_side_basis_qnode;
        const auto nc_side_cvariables = nc_solution_coefficient * nc_side_basis_qnode;

        const auto [num_equation, num_qnode] = oc_side_cvariables.size();
        const auto& normals = This_::set_of_normals_[i];

        Dynamic_Matrix numerical_flux_quadrature(This_::num_equation_, num_qnode);
        for (ushort q = 0; q < num_qnode; ++q) {
            const auto oc_side_cvariable = oc_side_cvariables.column<This_::num_equation_>(q);
            const auto nc_side_cvariable = nc_side_cvariables.column<This_::num_equation_>(q);

            numerical_flux_quadrature.change_column(q, NUMERICAL_FLUX_FUNCTION::calculate(oc_side_cvariable, nc_side_cvariable, normals[q]));
        }


        This_::Residual_ owner_side_delta_rhs, neighbor_side_delta_rhs;
        const auto& [oc_side_basis_weight, nc_side_basis_weight] = This_::oc_nc_side_basis_weight_pairs_[i];
        ms::gemm(numerical_flux_quadrature, oc_side_basis_weight, owner_side_delta_rhs);
        ms::gemm(numerical_flux_quadrature, nc_side_basis_weight, neighbor_side_delta_rhs);

        RHS[oc_index] -= owner_side_delta_rhs;
        RHS[nc_index] += neighbor_side_delta_rhs;
    }
}

//
////HOM이면 공통으로 사용하는 variable
//template <ushort space_dimension>
//class Periodic_Boundaries_FVM_Base
//{
//private:
//    using This_ = Periodic_Boundaries_FVM_Base<space_dimension>;
//    using Space_Vector_ = Euclidean_Vector<space_dimension>;
//
//protected:
//    inline static std::vector<Space_Vector_> normals_;
//    inline static std::vector<std::pair<uint, uint>> oc_nc_index_pairs_;
//    inline static std::vector<double> areas_;
//
//private:
//    Periodic_Boundaries_FVM_Base(void) = delete;
//
//public:
//    static void initialize(Grid<space_dimension>&& grid);
//};
//
//
////FVM이고 Constant Reconstruction이면 사용하는 variable & method
//template<typename Numerical_Flux_Function>
//class Periodic_Boundaries_FVM_Constant : public Periodic_Boundaries_FVM_Base<Numerical_Flux_Function::space_dimension()>
//{
//private:
//    static constexpr ushort space_dimension_ = Numerical_Flux_Function::space_dimension();
//    static constexpr ushort num_equation_ = Numerical_Flux_Function::num_equation();
//
//    using This_ = Periodic_Boundaries_FVM_Constant<Numerical_Flux_Function>;
//    using Space_Vector_ = Euclidean_Vector<This_::space_dimension_>;
//    using Solution_ = Euclidean_Vector<This_::num_equation_>;
//    using Residual_ = Euclidean_Vector< This_::num_equation_>;
//
//private:
//    Periodic_Boundaries_FVM_Constant(void) = delete;
//
//public:
//    static void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions);
//};
//
//
////FVM이고 Linear Reconstruction이면 공통으로 사용하는 variable & method
//template<typename Reconstruction_Method, typename Numerical_Flux_Function>
//class Periodic_Boundaries_FVM_Linear : public Periodic_Boundaries_FVM_Base<Numerical_Flux_Function::space_dimension()>
//{
//private:
//    static constexpr ushort space_dimension_ = Numerical_Flux_Function::space_dimension();
//    static constexpr ushort num_equation_ = Numerical_Flux_Function::num_equation();
//
//    using This_ = Periodic_Boundaries_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>;
//    using Base_ = Periodic_Boundaries_FVM_Base<This_::space_dimension_>;
//    using Space_Vector_ = Euclidean_Vector<This_::space_dimension_>;
//    using Solution_ = Euclidean_Vector<This_::num_equation_>;
//    using Residual_ = Euclidean_Vector<This_::num_equation_>;
//
//protected:
//    inline static std::vector<std::pair<Space_Vector_, Space_Vector_>> oc_nc_to_oc_nc_side_face_vector_pairs_;
//
//private:
//    Periodic_Boundaries_FVM_Linear(void) = delete;
//
//public:
//    static void initialize(Grid<This_::space_dimension_>&& grid);
//    static void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions);
//};
//
//
//
////template definition part
//template <ushort space_dimension>
//void Periodic_Boundaries_FVM_Base<space_dimension>::initialize(Grid<space_dimension>&& grid) {
//    SET_TIME_POINT;
//
//    This_::oc_nc_index_pairs_ = std::move(grid.connectivity.periodic_boundary_oc_nc_index_pairs);
//
//    const auto& cell_elements = grid.elements.cell_elements;
//    const auto& pbdry_element_pairs = grid.elements.periodic_boundary_element_pairs;
//
//    const auto num_pbdry_pair = pbdry_element_pairs.size();
//    This_::areas_.reserve(num_pbdry_pair);
//    This_::normals_.reserve(num_pbdry_pair);
//
//    for (uint i = 0; i < num_pbdry_pair; ++i) {
//        const auto& [oc_side_element, nc_side_element] = pbdry_element_pairs[i];
//
//        This_::areas_.push_back(oc_side_element.geometry_.volume());
//
//        const auto [oc_index, nc_index] = This_::oc_nc_index_pairs_[i];
//        const auto& oc_element = cell_elements[oc_index];
//
//        const auto oc_side_center = oc_side_element.geometry_.center_node();
//
//        This_::normals_.push_back(oc_side_element.normalized_normal_vector(oc_element, oc_side_center));
//    }
//
//    Log::content_ << std::left << std::setw(50) << "@ Periodic boundaries FVM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
//    Log::print();
//}
//
//template<typename Numerical_Flux_Function>
//void Periodic_Boundaries_FVM_Constant<Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) {
//    const auto numerical_fluxes = Numerical_Flux_Function::calculate(solutions, This_::normals_, This_::oc_nc_index_pairs_);
//
//    const auto num_pbdry_pair = This_::normals_.size();
//    for (size_t i = 0; i < num_pbdry_pair; ++i) {
//        const auto [oc_index, nc_index] = This_::oc_nc_index_pairs_[i];
//        const auto delta_RHS = This_::areas_[i] * numerical_fluxes[i];
//        RHS[oc_index] -= delta_RHS;
//        RHS[nc_index] += delta_RHS;
//    }
//}
//
//template<typename Reconstruction_Method, typename Numerical_Flux_Function>
//void Periodic_Boundaries_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>::initialize(Grid<This_::space_dimension_>&& grid) {
//    SET_TIME_POINT;
//
//    Base_::initialize(std::move(grid));
//
//    const auto& cell_elements = grid.elements.cell_elements;
//    const auto& pbdry_element_pairs = grid.elements.periodic_boundary_element_pairs;
//
//    const auto num_pbdry_pair = pbdry_element_pairs.size();
//    This_::oc_nc_to_oc_nc_side_face_vector_pairs_.reserve(num_pbdry_pair);
//
//    for (size_t i = 0; i < num_pbdry_pair; ++i) {
//        const auto& [oc_index, nc_index] = This_::oc_nc_index_pairs_[i];
//        const auto& oc_geometry = cell_elements[oc_index].geometry_;
//        const auto& nc_geometry = cell_elements[nc_index].geometry_;
//
//        const auto& [oc_side_element, nc_side_element] = pbdry_element_pairs[i];
//        const auto& oc_side_geometry = oc_side_element.geometry_;
//        const auto& nc_side_geometry = nc_side_element.geometry_;
//
//        const auto oc_center = oc_geometry.center_node();
//        const auto nc_center = nc_geometry.center_node();
//        const auto oc_side_center = oc_side_geometry.center_node();
//        const auto nc_side_center = nc_side_geometry.center_node();
//
//        const auto oc_to_oc_side_vector = oc_side_center - oc_center;
//        const auto nc_to_nc_side_vector = nc_side_center - nc_center;
//
//        This_::oc_nc_to_oc_nc_side_face_vector_pairs_.push_back(std::make_pair(oc_to_oc_side_vector, nc_to_nc_side_vector));
//    }
//
//    Log::content_ << std::left << std::setw(50) << "@ Periodic boundaries FVM linear precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
//    Log::print();
//}
//
//template<typename Reconstruction_Method, typename Numerical_Flux_Function>
//void Periodic_Boundaries_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) {
//    const auto& solution_gradients = Reconstruction_Method::get_solution_gradients();
//
//    const auto num_pbdry_pair = Base_::normals_.size();
//
//    for (size_t i = 0; i < num_pbdry_pair; ++i) {
//        const auto [oc_index, nc_index] = Base_::oc_nc_index_pairs_[i];
//        const auto& oc_solution = solutions[oc_index];
//        const auto& nc_solution = solutions[nc_index];
//
//        const auto& oc_solution_gradient = solution_gradients[oc_index];
//        const auto& nc_solution_gradient = solution_gradients[nc_index];
//
//        const auto& [oc_to_oc_side_face_vector, nc_to_nc_side_face_vector] = This_::oc_nc_to_oc_nc_side_face_vector_pairs_[i];
//
//        const auto oc_side_solution = oc_solution + oc_solution_gradient * oc_to_oc_side_face_vector;
//        const auto nc_side_solution = nc_solution + nc_solution_gradient * nc_to_nc_side_face_vector;
//        const auto& pbdry_normal = Base_::normals_[i];
//
//        const auto numerical_flux = Numerical_Flux_Function::calculate(oc_side_solution, nc_side_solution, pbdry_normal);
//        const auto delta_RHS = Base_::areas_[i] * numerical_flux;
//        RHS[oc_index] -= delta_RHS;
//        RHS[nc_index] += delta_RHS;
//    }
//}