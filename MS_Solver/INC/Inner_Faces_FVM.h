#pragma once
#include "Grid.h"
#include "Reconstruction_Method_FVM.h"


//FVM이면 공통으로 사용하는 variable
template <ushort space_dimension>
class Inner_Faces_FVM_Base
{
private:
    using Space_Vector_ = Euclidean_Vector<space_dimension>;

protected:
    std::vector<Space_Vector_> normals_;
    std::vector<std::pair<uint, uint>> oc_nc_index_pairs_;
    std::vector<double> areas_;

public:
    Inner_Faces_FVM_Base(const Grid<space_dimension>& grid);
};


//FVM이고 Constant Reconstruction이면 공통으로 사용하는 variable & Method
template<typename Numerical_Flux_Function>
class Inner_Faces_FVM_Constant : public Inner_Faces_FVM_Base<Numerical_Flux_Function::space_dimension()>
{
private:
    static constexpr ushort space_dimension_ = Numerical_Flux_Function::space_dimension();
    static constexpr ushort num_equation_ = Numerical_Flux_Function::num_equation();

    using Space_Vector_ = Euclidean_Vector<space_dimension_>;
    using Solution_ = Euclidean_Vector<num_equation_>;
    using Residual_ = Euclidean_Vector< num_equation_>;

public:
    Inner_Faces_FVM_Constant(const Grid<space_dimension_>& grid) : Inner_Faces_FVM_Base<space_dimension_>(grid) {};

    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const;
};


//FVM이고 Linear Reconstruction이면 공통으로 사용하는 variable & Method
template <typename Reconstruction_Method, typename Numerical_Flux_Function>
class Inner_Faces_FVM_Linear : public Inner_Faces_FVM_Base<Numerical_Flux_Function::space_dimension()>
{
private:
    static constexpr ushort space_dimension_    = Numerical_Flux_Function::space_dimension();
    static constexpr ushort num_equation_       = Numerical_Flux_Function::num_equation();

    using This_ = Inner_Faces_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>;
    using Space_Vector_ = Euclidean_Vector<This_::space_dimension_>;
    using Solution_ = Euclidean_Vector<This_::num_equation_>;
    using Residual_ = Euclidean_Vector<This_::num_equation_>;

private:
    std::vector<std::pair<Space_Vector_, Space_Vector_>> oc_nc_to_face_vector_pairs_;
    const Reconstruction_Method& reconstruction_method_;

public:
    Inner_Faces_FVM_Linear(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method);

public:
    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const;
};


// template definition part
template <ushort space_dimension>
Inner_Faces_FVM_Base<space_dimension>::Inner_Faces_FVM_Base(const Grid<space_dimension>& grid) {
    SET_TIME_POINT;

    this->oc_nc_index_pairs_ = grid.inner_face_oc_nc_index_pairs();

    const auto& grid_elements = grid.get_grid_elements();
    const auto& cell_elements = grid_elements.cell_elements;
    const auto& inner_face_elements = grid_elements.inner_face_elements;

    const auto num_inner_face = inner_face_elements.size();
    this->areas_.reserve(num_inner_face);
    this->normals_.reserve(num_inner_face);

    for (uint i = 0; i < num_inner_face; ++i) {
        const auto& inner_face_element = inner_face_elements[i];

        this->areas_.push_back(inner_face_element.geometry_.volume());

        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& oc_element = cell_elements[oc_index];

        const auto inner_face_center = inner_face_element.geometry_.center_node();

        this->normals_.push_back(inner_face_element.normalized_normal_vector(oc_element, inner_face_center));
    }

    Log::content_ << std::left << std::setw(50) << "@ Inner faces FVM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template<typename Numerical_Flux_Function>
void Inner_Faces_FVM_Constant<Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const {
    const auto num_inner_face = this->normals_.size();

    const auto numerical_fluxes = Numerical_Flux_Function::calculate(solutions, this->normals_, this->oc_nc_index_pairs_);

    for (size_t i = 0; i < num_inner_face; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto delta_RHS = this->areas_[i] * numerical_fluxes[i];
        RHS[oc_index] -= delta_RHS;
        RHS[nc_index] += delta_RHS;
    }
}

template<typename Reconstruction_Method, typename Numerical_Flux_Function>
Inner_Faces_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>::Inner_Faces_FVM_Linear(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
    : Inner_Faces_FVM_Base<space_dimension_>(grid), reconstruction_method_(reconstruction_method) {
    SET_TIME_POINT;

    const auto num_inner_face = this->normals_.size();
    this->oc_nc_to_face_vector_pairs_.reserve(num_inner_face);

    const auto& grid_elements = grid.get_grid_elements();
    const auto& cell_elements = grid_elements.cell_elements;
    const auto& inner_face_elements = grid_elements.inner_face_elements;

    for (size_t i = 0; i < num_inner_face; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];

        const auto& oc_geometry = cell_elements[oc_index].geometry_;
        const auto& nc_geometry = cell_elements[nc_index].geometry_;
        const auto& inner_face_geometry = inner_face_elements[i].geometry_;

        const auto oc_center = oc_geometry.center_node();
        const auto nc_center = nc_geometry.center_node();
        const auto inner_face_center = inner_face_geometry.center_node();

        const auto oc_to_face_vector = inner_face_center - oc_center;
        const auto nc_to_face_vector = inner_face_center - nc_center;

        this->oc_nc_to_face_vector_pairs_.push_back(std::make_pair(oc_to_face_vector, nc_to_face_vector));
    }

    Log::content_ << std::left << std::setw(50) << "@ Inner faces FVM linear precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
};


template<typename Reconstruction_Method, typename Numerical_Flux_Function>
void Inner_Faces_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const {
    const auto& solution_gradients = this->reconstruction_method_.get_solution_gradients();

    const auto num_inner_face = this->normals_.size();

    for (size_t i = 0; i < num_inner_face; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& oc_solution = solutions[oc_index];
        const auto& nc_solution = solutions[nc_index];

        const auto& oc_solution_gradient = solution_gradients[oc_index];
        const auto& nc_solution_gradient = solution_gradients[nc_index];

        const auto& [oc_to_face_vector, nc_to_face_vector] = this->oc_nc_to_face_vector_pairs_[i];

        const auto oc_side_solution = oc_solution + oc_solution_gradient * oc_to_face_vector;
        const auto nc_side_solution = nc_solution + nc_solution_gradient * nc_to_face_vector;
        const auto& inner_face_normal = this->normals_[i];

        const auto numerical_flux = Numerical_Flux_Function::calculate(oc_side_solution, nc_side_solution, inner_face_normal);
        const auto delta_RHS = this->areas_[i] * numerical_flux;
        RHS[oc_index] -= delta_RHS;
        RHS[nc_index] += delta_RHS;
    }
}