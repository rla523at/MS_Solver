#pragma once
#include "Grid_Builder.h"


template <ushort num_equation_, ushort space_dimension_>
class Least_Square_Base
{
private:    
    using Solution_             = Euclidean_Vector<num_equation_>;
    using Solution_gradinet_    = Matrix<num_equation_, space_dimension_>;

protected:
    std::vector<std::vector<size_t>> near_cell_indexes_set_;
    std::vector<Dynamic_Matrix> least_square_matrixes_;

public:
    std::vector<Solution_gradinet_> calculate_solution_gradients(const std::vector<Solution_>& solutions) const;

private:
    std::vector<Dynamic_Matrix> calculate_solution_delta_matrixes(const std::vector<Solution_>& solutions) const;

public:
    static constexpr ushort space_dimension(void) { return space_dimension_; };
    static constexpr ushort num_equation(void) { return num_equation_; };
};


template <ushort num_equation_, ushort space_dimension_>
class Vertex_Least_Square : public Least_Square_Base<num_equation_, space_dimension_>
{
public:
    Vertex_Least_Square(const Grid<space_dimension_>& grid);

    static std::string name(void) { return "Vertex_Least_Square"; };
};


template <ushort num_equation_, ushort space_dimension_>
class Face_Least_Square : public Least_Square_Base<num_equation_, space_dimension_>
{
public:
    Face_Least_Square(const Grid<space_dimension_>& grid);

    static std::string name(void) { return "Face_Least_Square"; };
};



//template definition part
template <ushort num_equation_, ushort space_dimension_>
std::vector<Matrix<num_equation_, space_dimension_>> Least_Square_Base<num_equation_, space_dimension_>::calculate_solution_gradients(const std::vector<Solution_>& solutions)  const {
    const auto num_solution = solutions.size();

    std::vector<Matrix<num_equation_, space_dimension_>> solution_gradients(num_solution);

    const auto solution_delta_matrixes = this->calculate_solution_delta_matrixes(solutions);
    for (size_t i = 0; i < num_solution; ++i)
        ms::gemm(solution_delta_matrixes[i], this->least_square_matrixes_[i], solution_gradients[i]);

    return solution_gradients;
}

template <ushort num_equation_, ushort space_dimension_>
std::vector<Dynamic_Matrix> Least_Square_Base<num_equation_, space_dimension_>::calculate_solution_delta_matrixes(const std::vector<Solution_>& solutions) const {
    const auto num_solution = solutions.size();

    std::vector<Dynamic_Matrix> solution_delta_matrixes;
    solution_delta_matrixes.reserve(num_solution);

    for (size_t i = 0; i < num_solution; ++i) {
        const auto& near_cell_indexes = this->near_cell_indexes_set_.at(i);
        const auto num_near_cell = near_cell_indexes.size();

        Dynamic_Matrix solution_delta_matrix(num_equation_, num_near_cell);
        for (size_t j = 0; j < num_near_cell; ++j) {
            const auto solution_delta = solutions[near_cell_indexes[j]] - solutions[i];
            solution_delta_matrix.change_column(j, solution_delta);
        }

        solution_delta_matrixes.push_back(std::move(solution_delta_matrix));
    }

    return solution_delta_matrixes;
}

template <ushort num_equation_, ushort space_dimension_>
Vertex_Least_Square<num_equation_, space_dimension_>::Vertex_Least_Square(const Grid<space_dimension_>& grid) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    this->near_cell_indexes_set_.reserve(num_cell);
    this->least_square_matrixes_.reserve(num_cell);

    auto set_of_vertex_share_cell_indexes = grid.calculate_set_of_vertex_share_cell_indexes();

    for (size_t i = 0; i < num_cell; ++i) {
        auto& vertex_share_cell_indexes = set_of_vertex_share_cell_indexes[i];
        const auto num_vertex_share_cell = vertex_share_cell_indexes.size();

        const auto& element = cell_elements[i];
        const auto& geometry = element.geometry_;

        Dynamic_Matrix center_to_center_matrix(space_dimension_, num_vertex_share_cell);

        const auto this_center = geometry.center_node();

        for (size_t j = 0; j < num_vertex_share_cell; ++j) {
            const auto& neighbor_geometry = cell_elements[vertex_share_cell_indexes[j]].geometry_;
            const auto neighbor_center = neighbor_geometry.center_node();
            const auto center_to_center = neighbor_center - this_center;

            center_to_center_matrix.change_column(j, center_to_center);
        }

        const auto& Rc = center_to_center_matrix;
        auto RcT = Rc.transpose();
        auto least_square_matrix = RcT * (Rc * RcT).be_inverse();

        this->near_cell_indexes_set_.push_back(std::move(vertex_share_cell_indexes));
        this->least_square_matrixes_.push_back(std::move(least_square_matrix));
    }

    Log::content_ << std::left << std::setw(50) << "@ Vertex Least Sqaure precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}



template <ushort num_equation_, ushort space_dimension_>
Face_Least_Square<num_equation_, space_dimension_>::Face_Least_Square(const Grid<space_dimension_>& grid) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    this->near_cell_indexes_set_.reserve(num_cell);
    this->least_square_matrixes_.reserve(num_cell);

    auto set_of_face_share_cell_indexes = grid.calculate_set_of_face_share_cell_indexes();

    for (size_t i = 0; i < num_cell; ++i) {        
        auto& face_share_cell_indexes = set_of_face_share_cell_indexes[i];         
        const auto num_face_share_cell = face_share_cell_indexes.size();

        const auto target_cell_element = cell_elements[i];
        const auto target_cell_center = target_cell_element.geometry_.center_node();

        Dynamic_Matrix center_to_center_matrix(space_dimension_, num_face_share_cell);
        for (size_t j = 0; j < num_face_share_cell; ++j) {
            const auto& neighbor_geometry = cell_elements[face_share_cell_indexes[j]].geometry_;
            const auto neighbor_center = neighbor_geometry.center_node();
            const auto center_to_center = neighbor_center - target_cell_center;

            center_to_center_matrix.change_column(j, center_to_center);
        }

        const auto& Rc = center_to_center_matrix;
        auto RcT = Rc.transpose();
        auto least_square_matrix = RcT * (Rc * RcT).be_inverse();

        this->near_cell_indexes_set_.push_back(std::move(face_share_cell_indexes));
        this->least_square_matrixes_.push_back(std::move(least_square_matrix));
    }

    Log::content_ << std::left << std::setw(50) << "@ Face Least Sqaure precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}