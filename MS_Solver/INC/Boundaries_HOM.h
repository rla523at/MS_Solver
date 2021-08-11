#pragma once
#include "Boundary_Flux_Function.h"
#include "Grid_Builder.h"
#include "Reconstruction_Method.h"


//HOM이면 공통으로 사용하는 variable
template <typename Governing_Equation, typename Reconstruction_Method>
class Boundaries_HOM_Base
{
private:
    static constexpr ushort space_dimension_    = Governing_Equation::space_dimension();
    static constexpr ushort solution_order_     = Reconstruction_Method::solution_order();

    using This_         = Boundaries_HOM_Base<Governing_Equation, Reconstruction_Method>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;

private:
    Boundaries_HOM_Base(void) = delete;

protected:
    inline static std::vector<std::unique_ptr<Boundary_Flux_Function<Governing_Equation>>> boundary_flux_functions_;

    inline static std::vector<uint> oc_indexes_;
    inline static std::vector<Dynamic_Matrix> set_of_oc_side_basis_qnodes_;
    inline static std::vector<std::vector<Space_Vector_>> set_of_normals_;

    inline static std::vector<Dynamic_Matrix> basis_weights_;

public:
    static void initialize(Grid<space_dimension_>&& grid);
};


////FVM이고 Constant Reconstruction이면 공통으로 사용하는 variable & Method
//template <typename Governing_Equation>
//class Boundaries_FVM_Constant : public Boundaries_HOM_Base<Governing_Equation>
//{
//private:
//    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
//    static constexpr size_t num_equation_ = Governing_Equation::num_equation();
//
//    using Base_ = Boundaries_HOM_Base<Governing_Equation>;
//    using Solution_ = typename Governing_Equation::Solution_;
//    using Boundary_Flux_ = Euclidean_Vector<num_equation_>;
//
//private:
//    Boundaries_FVM_Constant(void) = delete;
//
//public:
//    static void calculate_RHS(std::vector<Boundary_Flux_>& RHS, const std::vector<Solution_>& solutions);
//};
//
//
////FVM이고 Constant Reconstruction이면 공통으로 사용하는 variable & Method
//template <typename Governing_Equation, typename Reconstruction_Method>
//class Boundaries_FVM_Linear : public Boundaries_HOM_Base<Governing_Equation>
//{
//private:
//    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
//    static constexpr size_t num_equation_ = Governing_Equation::num_equation();
//
//    using Base_ = Boundaries_HOM_Base<Governing_Equation>;
//    using This_ = Boundaries_FVM_Linear<Governing_Equation, Reconstruction_Method>;
//    using Space_Vector_ = Euclidean_Vector <space_dimension_>;
//    using Solution_ = Euclidean_Vector<num_equation_>;
//    using Boundary_Flux_ = Euclidean_Vector<num_equation_>;
//
//private:
//    inline static std::vector<Space_Vector_> oc_to_boundary_vectors_;
//
//private:
//    Boundaries_FVM_Linear(void) = delete;
//
//public:
//    static void initialize(Grid<space_dimension_>&& grid);
//    static void calculate_RHS(std::vector<Boundary_Flux_>& RHS, const std::vector<Solution_>& solutions);
//};
//


//template definition
template <typename Governing_Equation, typename Reconstruction_Method>
void Boundaries_HOM_Base<Governing_Equation, Reconstruction_Method>::initialize(Grid<space_dimension_>&& grid) {
    SET_TIME_POINT;

    This_::oc_indexes_ = std::move(grid.connectivity.boundary_oc_indexes);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& boundary_elements = grid.elements.boundary_elements;

    const auto num_boundaries = boundary_elements.size();
    This_::boundary_flux_functions_.reserve(num_boundaries);
    This_::set_of_oc_side_basis_qnodes_.reserve(num_boundaries);
    This_::set_of_normals_.reserve(num_boundaries);
    This_::basis_weights_.reserve(num_boundaries);

    
    constexpr auto integrand_order = 2 * This_::solution_order_ + 1;
    for (uint i = 0; i < num_boundaries; ++i) {
        const auto& boundary_element = boundary_elements[i];
        const auto& boundary_geometry = boundary_element.geometry_;

        This_::boundary_flux_functions_.push_back(Boundary_Flux_Function_Factory<Governing_Equation>::make(boundary_element.type()));

        const auto oc_index = This_::oc_indexes_[i];
        const auto& oc_element = cell_elements[oc_index];

        const auto& quadrature_rule = boundary_geometry.get_quadrature_rule(integrand_order);
        const auto& qnodes = quadrature_rule.points;
        const auto num_qnode = qnodes.size();

        This_::set_of_oc_side_basis_qnodes_.push_back(Reconstruction_Method::calculate_basis_nodes(oc_index, qnodes));

        std::vector<Space_Vector_> normals(num_qnode);
        
        for (ushort q = 0; q < num_qnode; ++q) 
            normals[q] = boundary_element.normalized_normal_vector(oc_element, qnodes[q]);
        
        This_::set_of_normals_.push_back(std::move(normals));


    }

    Log::content_ << std::left << std::setw(50) << "@ Boundaries FVM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}
//
//template <typename Governing_Equation>
//void Boundaries_FVM_Constant<Governing_Equation>::calculate_RHS(std::vector<Boundary_Flux_>& RHS, const std::vector<Solution_>& solutions) {
//    const auto num_boundary = Base_::normals_.size();
//
//    for (size_t i = 0; i < num_boundary; ++i) {
//        const auto& boundary_flux_function = Base_::boundary_flux_functions_.at(i);
//
//        const auto oc_index = Base_::oc_indexes_[i];
//        const auto& normal = Base_::normals_[i];
//
//        const auto boundary_flux = boundary_flux_function->calculate(solutions[oc_index], normal);
//        const auto delta_RHS = Base_::areas_[i] * boundary_flux;
//
//        RHS[oc_index] -= delta_RHS;
//    }
//}
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//void Boundaries_FVM_Linear<Governing_Equation, Reconstruction_Method>::initialize(Grid<space_dimension_>&& grid) {
//    SET_TIME_POINT;
//
//    Base_::initialize(std::move(grid));
//
//    const auto num_boundary = Base_::normals_.size();
//    This_::oc_to_boundary_vectors_.reserve(num_boundary);
//
//    const auto& cell_elements = grid.elements.cell_elements;
//    const auto& boundary_elements = grid.elements.boundary_elements;
//
//    for (size_t i = 0; i < num_boundary; ++i) {
//        const auto oc_index = Base_::oc_indexes_[i];
//        const auto& oc_geometry = cell_elements[oc_index].geometry_;
//        const auto& boundary_geometry = boundary_elements[i].geometry_;
//
//        const auto oc_center = oc_geometry.center_node();
//        const auto boundary_center = boundary_geometry.center_node();
//
//        const auto oc_to_face_vector = boundary_center - oc_center;
//
//        This_::oc_to_boundary_vectors_.push_back(oc_to_face_vector);
//    }
//
//    Log::content_ << std::left << std::setw(50) << "@ Boundaries FVM linear precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
//    Log::print();
//};
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//void Boundaries_FVM_Linear<Governing_Equation, Reconstruction_Method>::calculate_RHS(std::vector<Boundary_Flux_>& RHS, const std::vector<Solution_>& solutions) {
//    const auto& solution_gradients = Reconstruction_Method::get_solution_gradients();
//
//    const auto num_boundary = Base_::normals_.size();
//
//    for (size_t i = 0; i < num_boundary; ++i) {
//        const auto& boundary_flux_function = Base_::boundary_flux_functions_.at(i);
//        const auto& normal = Base_::normals_.at(i);
//
//        const auto oc_index = Base_::oc_indexes_.at(i);
//
//        const auto& oc_solution = solutions[oc_index];
//        const auto& oc_solution_gradient = solution_gradients[oc_index];
//        const auto& oc_to_face_vector = This_::oc_to_boundary_vectors_[i];
//
//        const auto oc_side_solution = oc_solution + oc_solution_gradient * oc_to_face_vector;
//
//        const auto boundary_flux = boundary_flux_function->calculate(oc_side_solution, normal);
//
//        const auto delta_RHS = Base_::areas_[i] * boundary_flux;
//
//        RHS[oc_index] -= delta_RHS;
//    }
//}