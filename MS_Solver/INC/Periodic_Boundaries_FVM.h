#pragma once
#include "Grid_Builder.h"
#include "Reconstruction_Method.h"


//FVM이면 공통으로 사용하는 variable
template <ushort space_dimension>
class Periodic_Boundaries_FVM_Base
{
    using Space_Vector_ = Euclidean_Vector<space_dimension>;

protected:
    size_t num_pbdry_pair_ = 0;
    std::vector<Space_Vector_> normals_;
    std::vector<std::pair<uint, uint>> oc_nc_index_pairs_;
    std::vector<double> areas_;

public:
    Periodic_Boundaries_FVM_Base(Grid<space_dimension>&& grid);
};


//FVM이고 Constant Reconstruction이면 사용하는 variable & method
template <ushort space_dimension>
class Periodic_Boundaries_FVM_Constant : public Periodic_Boundaries_FVM_Base<space_dimension>
{
private:
    using Space_Vector_ = Euclidean_Vector<space_dimension>;

protected:
    std::vector<std::pair<Space_Vector_, Space_Vector_>> oc_nc_to_oc_nc_side_face_vector_pairs_;

public:
    Periodic_Boundaries_FVM_Constant(Grid<space_dimension>&& grid) : Periodic_Boundaries_FVM_Base<space_dimension>(std::move(grid)) {};

    template<typename Numerical_Flux_Function, typename Residual, typename Solution>
    void calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const;
};


//FVM이고 Linear Reconstruction이면 공통으로 사용하는 variable & method
template <ushort space_dimension>
class Periodic_Boundaries_FVM_Linear : public Periodic_Boundaries_FVM_Base<space_dimension>
{
private:
    using Space_Vector_ = Euclidean_Vector<space_dimension>;

protected:
    std::vector<std::pair<Space_Vector_, Space_Vector_>> oc_nc_to_oc_nc_side_face_vector_pairs_;

public:
    Periodic_Boundaries_FVM_Linear(Grid<space_dimension>&& grid);

    template<typename Numerical_Flux_Function, ushort num_equation>
    void calculate_RHS(std::vector<Euclidean_Vector<num_equation>>& RHS, const Linear_Reconstructed_Solution<num_equation, space_dimension>& linear_reconstructed_solution) const;
};



//template definition part
template <ushort space_dimension>
Periodic_Boundaries_FVM_Base<space_dimension>::Periodic_Boundaries_FVM_Base(Grid<space_dimension>&& grid) {
    SET_TIME_POINT;

    this->oc_nc_index_pairs_ = std::move(grid.connectivity.periodic_boundary_oc_nc_index_pairs);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& periodic_boundary_element_pairs = grid.elements.periodic_boundary_element_pairs;

    this->num_pbdry_pair_ = periodic_boundary_element_pairs.size();
    this->areas_.reserve(this->num_pbdry_pair_);
    this->normals_.reserve(this->num_pbdry_pair_);

    for (uint i = 0; i < this->num_pbdry_pair_; ++i) {
        const auto& [oc_side_element, nc_side_element] = periodic_boundary_element_pairs[i];

        this->areas_.push_back(oc_side_element.geometry_.volume());

        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& oc_element = cell_elements[oc_index];

        const auto oc_side_center = oc_side_element.geometry_.center_node();

        this->normals_.push_back(oc_side_element.normalized_normal_vector(oc_element, oc_side_center));
    }

    Log::content_ << std::left << std::setw(50) << "@ Periodic boundaries FVM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <ushort space_dimension>
template<typename Numerical_Flux_Function, typename Residual, typename Solution>
void Periodic_Boundaries_FVM_Constant<space_dimension>::calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const {
    const auto numerical_fluxes = Numerical_Flux_Function::calculate(solutions, this->normals_, this->oc_nc_index_pairs_);
    for (size_t i = 0; i < this->num_pbdry_pair_; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto delta_RHS = this->areas_[i] * numerical_fluxes[i];
        RHS[oc_index] -= delta_RHS;
        RHS[nc_index] += delta_RHS;
    }
}

template <ushort space_dimension>
Periodic_Boundaries_FVM_Linear<space_dimension>::Periodic_Boundaries_FVM_Linear(Grid<space_dimension>&& grid) :Periodic_Boundaries_FVM_Base<space_dimension>(std::move(grid)) {
    SET_TIME_POINT;

    this->oc_nc_to_oc_nc_side_face_vector_pairs_.reserve(this->num_pbdry_pair_);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& pbdry_element_pairs = grid.elements.periodic_boundary_element_pairs;
    for (size_t i = 0; i < this->num_pbdry_pair_; ++i) {
        const auto& [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& oc_geometry = cell_elements[oc_index].geometry_;
        const auto& nc_geometry = cell_elements[nc_index].geometry_;

        const auto& [oc_side_element, nc_side_element] = pbdry_element_pairs[i];
        const auto& oc_side_geometry = oc_side_element.geometry_;
        const auto& nc_side_geometry = nc_side_element.geometry_;

        const auto oc_center = oc_geometry.center_node();
        const auto nc_center = nc_geometry.center_node();
        const auto oc_side_center = oc_side_geometry.center_node();
        const auto nc_side_center = nc_side_geometry.center_node();

        const auto oc_to_oc_side_face_vector = oc_side_center - oc_center;
        const auto nc_to_nc_side_face_vector = nc_side_center - nc_center;

        this->oc_nc_to_oc_nc_side_face_vector_pairs_.push_back(std::make_pair(oc_to_oc_side_face_vector, nc_to_nc_side_face_vector));
    }

    Log::content_ << std::left << std::setw(50) << "@ Periodic boundaries FVM linear precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <ushort space_dimension>
template <typename Numerical_Flux_Function, ushort num_equation>
void Periodic_Boundaries_FVM_Linear<space_dimension>::calculate_RHS(std::vector<Euclidean_Vector<num_equation>>& RHS, const Linear_Reconstructed_Solution<num_equation, space_dimension>& linear_reconstructed_solution) const {
    const auto& solutions = linear_reconstructed_solution.solutions;
    const auto& solution_gradients = linear_reconstructed_solution.solution_gradients;

    for (size_t i = 0; i < this->num_pbdry_pair_; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& oc_solution = solutions[oc_index];
        const auto& nc_solution = solutions[nc_index];

        const auto& oc_solution_gradient = solution_gradients[oc_index];
        const auto& nc_solution_gradient = solution_gradients[nc_index];

        const auto& [oc_to_oc_side_face_vector, nc_to_nc_side_face_vector] = this->oc_nc_to_oc_nc_side_face_vector_pairs_[i];

        const auto oc_side_solution = oc_solution + oc_solution_gradient * oc_to_oc_side_face_vector;
        const auto nc_side_solution = nc_solution + nc_solution_gradient * nc_to_nc_side_face_vector;
        const auto& pbdry_normal = this->normals_[i];

        const auto numerical_flux = Numerical_Flux_Function::calculate(oc_side_solution, nc_side_solution, pbdry_normal);
        const auto delta_RHS = this->areas_[i] * numerical_flux;
        RHS[oc_index] -= delta_RHS;
        RHS[nc_index] += delta_RHS;
    }
}