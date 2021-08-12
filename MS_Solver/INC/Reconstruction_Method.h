#pragma once
#include <type_traits>
#include <string>

#include "Gradient_Method.h"
#include "PostAI.h"

class RM {};	// Reconstruction Method


class Constant_Reconstruction : public RM 
{
private:
    Constant_Reconstruction(void) = delete;

public:
    template <size_t space_dimension>
    static void initialize(const Grid<space_dimension>& grid) {}; //because of same form in main
    static std::string name(void) { return "Constant_Reconstruction"; };
};


template <typename Gradient_Method>
class Linear_Reconstruction : public RM
{
    //Gradient_Method<num_equation, space_dimension>
private:
    static constexpr size_t num_equation_       = Gradient_Method::num_equation();
    static constexpr size_t space_dimension_    = Gradient_Method::space_dimension();

    using This_                 = Linear_Reconstruction<Gradient_Method>;
    using Solution_             = Euclidean_Vector<num_equation_>;
    using Solution_Gradient_    = Matrix<num_equation_, space_dimension_>;

private:
    inline static std::vector<Solution_Gradient_> solution_gradients_;

private:
    Linear_Reconstruction(void) = delete;

public:
    static void initialize(const Grid<space_dimension_>& grid) { Gradient_Method::initialize(grid); };
    static void reconstruct(const std::vector<Solution_>& solutions);
    static std::string name(void) { return "Linear_Reconstruction_" + Gradient_Method::name(); };
    static const std::vector<Solution_Gradient_>& get_solution_gradients(void) { return solution_gradients_; };
};


template <typename Gradient_Method>
class MLP_u1 : public RM
{
private:
    static constexpr size_t num_equation_       = Gradient_Method::num_equation();
    static constexpr size_t space_dimension_    = Gradient_Method::space_dimension();

    using This_                 = MLP_u1<Gradient_Method>;
    using Solution_             = Euclidean_Vector<num_equation_>;
    using Solution_Gradient_    = Matrix<num_equation_, space_dimension_>;

private:
    MLP_u1(void) = delete;

protected:
    inline static std::vector<Solution_Gradient_> solution_gradients_;
	inline static std::vector<std::vector<uint>> vnode_indexes_set_;
    inline static std::vector<Dynamic_Matrix> center_to_vertex_matrixes_;
    inline static std::unordered_map<uint, std::set<uint>> vnode_index_to_share_cell_indexes_;

public:
    static void initialize(Grid<space_dimension_>&& grid);
    static void reconstruct(const std::vector<Solution_>& solutions);
    static std::string name(void) { return "MLP_u1_" + Gradient_Method::name(); };
    static const std::vector<Solution_Gradient_>& get_solution_gradients(void) { return solution_gradients_; };

protected:
    static auto calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution_>& solutions);
    static double limit(const double vertex_solution_delta, const double center_solution, const double min_solution, const double max_solution);
};


struct ANN_Model
{
    std::vector<Dynamic_Matrix> weights;
    std::vector<Dynamic_Euclidean_Vector> biases;
};


class ReLU {
public:
    double operator()(const double value) {
        return (value < 0) ? 0 : value;
    };
};

class HardTanh {
public:
    double operator()(const double value) {
        return (value > 1) ? 1 : (value < 0) ? 0 : value;
    };
};


template <typename Gradient_Method>
class ANN_limiter : public RM
{
private:
    static constexpr ushort num_equation_ = Gradient_Method::num_equation();
    static constexpr ushort space_dimension_ = Gradient_Method::space_dimension();

    using This_ = ANN_limiter<Gradient_Method>;
    using Solution_ = Euclidean_Vector<num_equation_>;
    using Solution_Gradient_ = Matrix<num_equation_, space_dimension_>;

private:
    inline static std::vector<Solution_Gradient_> solution_gradients_;
    inline static std::vector<std::vector<size_t>> set_of_vertex_share_cell_indexes;
    inline static std::vector<std::vector<size_t>> set_of_face_share_cell_indexes;

private:
    ANN_limiter(void) = delete;

public:
    static void initialize(const Grid<space_dimension_>& grid);
    static void reconstruct(const std::vector<Solution_>& solutions);
    static std::string name(void) { return "ANN_Reconstruction_" + Gradient_Method::name(); };
    static const std::vector<Solution_Gradient_>& get_solution_gradients(void) { return solution_gradients_; };

private:
    static bool is_constant_region(const std::vector<Solution_>& solutions, const size_t target_cell_index);
    static std::vector<size_t> ordering_function(const std::vector<Solution_>& solutions, const size_t target_cell_index);

    static void limit(Dynamic_Euclidean_Vector& feature);
    static ANN_Model read_model(void);
};


class HOM_Reconstruction : public RM {};

template <ushort space_dimension_, ushort solution_order_>
class Polynomial_Reconstruction : public HOM_Reconstruction
{
private:
    using This_         = Polynomial_Reconstruction<space_dimension_, solution_order_>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;
        
    static constexpr ushort num_basis_ = ms::combination_with_repetition(1 + space_dimension_, solution_order_);

private:
    inline static std::vector<Vector_Function<Polynomial<space_dimension_>, num_basis_>> set_of_basis_functions_;

private:
    Polynomial_Reconstruction(void) = delete;

public:
    static void initialize(const Grid<space_dimension_>& grid);
    static auto calculate_set_of_transposed_gradient_basis(void);
    static auto calculate_basis_node(const uint cell_index, const Space_Vector_& node);
    static Dynamic_Matrix calculate_basis_nodes(const uint cell_index, const std::vector<Space_Vector_>& nodes);
    static double calculate_P0_basis_value(const uint cell_index, const Space_Vector_& center_node);
    static std::vector<Dynamic_Matrix> calculate_set_of_basis_nodes(const std::vector<std::vector<Space_Vector_>>& set_of_nodes);
    static std::vector<double> calculate_P0_basis_values(const std::vector<Space_Vector_>& center_nodes);
    static constexpr ushort num_basis(void) { return This_::num_basis_; };
    static constexpr ushort solution_order(void) { return solution_order_; };
    static constexpr ushort space_dimension(void) { return space_dimension_; };
};


namespace ms {
	template <typename T>
	inline constexpr bool is_reconsturction_method = std::is_base_of_v<RM, T>;

    template <typename T>
    inline constexpr bool is_default_reconstruction = std::is_same_v<Constant_Reconstruction, T>; /*|| std::is_same_v<Polynomial_Reconstruction, T>;*/
}


//template definition part
template <typename Gradient_Method>
void Linear_Reconstruction<Gradient_Method>::reconstruct(const std::vector<Solution_>& solutions) {
    const auto num_cell = solutions.size();
    const auto solution_gradients_temp = Gradient_Method::calculate_solution_gradients(solutions);

    //dynamic matrix to matrix
    This_::solution_gradients_.reserve(num_cell);

    for (const auto& solution_gradient : solution_gradients_temp)
        solution_gradients_.push_back(solution_gradient);
}

template <typename Gradient_Method>
void MLP_u1<Gradient_Method>::initialize(Grid<space_dimension_>&& grid) {
    SET_TIME_POINT;
    Gradient_Method::initialize(grid);

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    This_::vnode_indexes_set_.reserve(num_cell);
    This_::center_to_vertex_matrixes_.reserve(num_cell);
    This_::solution_gradients_.resize(num_cell);

    //vnode index to share cell indexes
    This_::vnode_index_to_share_cell_indexes_ = std::move(grid.connectivity.vnode_index_to_share_cell_indexes);

    for (size_t i = 0; i < num_cell; ++i) {
        const auto& element = cell_elements[i];
        const auto& geometry = cell_elements[i].geometry_;

        // vnode indexes set
        This_::vnode_indexes_set_.push_back(element.vertex_node_indexes());

        //center to vertex matrix
        const auto center_node = geometry.center_node();
        const auto vertex_nodes = geometry.vertex_nodes();
        const auto num_vertex = vertex_nodes.size();

        Dynamic_Matrix center_to_vertex_matrix(space_dimension_, num_vertex);
        for (size_t i = 0; i < num_vertex; ++i) {
            const auto center_to_vertex = vertex_nodes[i] - center_node;
            center_to_vertex_matrix.change_column(i, center_to_vertex);
        }

        This_::center_to_vertex_matrixes_.push_back(std::move(center_to_vertex_matrix));
    }

    Log::content_ << std::left << std::setw(50) << "@ MLP u1 precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}


template <typename Gradient_Method>
void MLP_u1<Gradient_Method>::reconstruct(const std::vector<Solution_>& solutions) {
    const auto solution_gradients = Gradient_Method::calculate_solution_gradients(solutions);
    const auto vnode_index_to_min_max_solution = This_::calculate_vertex_node_index_to_min_max_solution(solutions);

    Post_AI_Data::record_solution_datas(solutions, solution_gradients);

    const auto num_cell = solutions.size();

    for (uint i = 0; i < num_cell; ++i) {
        const auto& gradient = solution_gradients[i];
        const auto vertex_solution_delta_matrix = gradient * This_::center_to_vertex_matrixes_[i];

        std::array<double, num_equation_> limiting_values;
        limiting_values.fill(1);

        const auto& vnode_indexes = This_::vnode_indexes_set_[i];
        const auto num_vertex = vnode_indexes.size();

        for (ushort j = 0; j < num_vertex; ++j) {
            const auto vnode_index = vnode_indexes[j];
            const auto& [min_solution, max_solution] = vnode_index_to_min_max_solution.at(vnode_index);

            for (ushort e = 0; e < num_equation_; ++e) {
                const auto limiting_value = This_::limit(vertex_solution_delta_matrix.at(e, j), solutions[i].at(e), min_solution.at(e), max_solution.at(e));
                limiting_values[e] = min(limiting_values[e], limiting_value);
            }
        }

        Post_AI_Data::record_limiting_value(i, limiting_values);

        //Matrix<num_equation_, num_equation_> limiting_value_matrix = limiting_values;
        const auto limiting_value_matrix = Matrix<num_equation_, num_equation_>::diagonal_matrix(limiting_values);

        This_::solution_gradients_[i] = limiting_value_matrix * gradient;
    }

    Post_AI_Data::post();
}



template <typename Gradient_Method>
auto MLP_u1<Gradient_Method>::calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution_>& solutions) {
    const auto num_vnode = This_::vnode_index_to_share_cell_indexes_.size();

    std::unordered_map<uint, std::pair<Solution_, Solution_>> vnode_index_to_min_max_solution;
    vnode_index_to_min_max_solution.reserve(num_vnode);

    for (const auto& [vnode_index, share_cell_indexes] : This_::vnode_index_to_share_cell_indexes_) {
        const auto num_share_cell = share_cell_indexes.size();
        std::array<std::vector<double>, num_equation_> equation_wise_solutions;

        for (ushort i = 0; i < num_equation_; ++i)
            equation_wise_solutions[i].reserve(num_share_cell);

        for (const auto cell_index : share_cell_indexes) {
            for (ushort i = 0; i < num_equation_; ++i)
                equation_wise_solutions[i].push_back(solutions[cell_index].at(i));
        }

        std::array<double, num_equation_> min_solution;
        std::array<double, num_equation_> max_solution;
        for (ushort i = 0; i < num_equation_; ++i) {
            min_solution[i] = *std::min_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
            max_solution[i] = *std::max_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
        }

        Solution_ min_sol = min_solution;
        Solution_ max_sol = max_solution;

        vnode_index_to_min_max_solution.emplace(vnode_index, std::make_pair(min_sol, max_sol));
    }

    return vnode_index_to_min_max_solution;
}


template <typename Gradient_Method>
double MLP_u1<Gradient_Method>::limit(const double vertex_solution_delta, const double center_solution, const double min_solution, const double max_solution) {
    if (vertex_solution_delta < 0)
        return min((min_solution - center_solution) / vertex_solution_delta, 1);
    else
        return min((max_solution - center_solution) / vertex_solution_delta, 1);
}


template <typename Gradient_Method>
void ANN_limiter<Gradient_Method>::initialize(const Grid<space_dimension_>& grid) {
    Gradient_Method::initialize(grid);

    This_::set_of_face_share_cell_indexes = grid.calculate_set_of_face_share_cell_indexes();
    This_::set_of_vertex_share_cell_indexes = grid.calculate_set_of_vertex_share_cell_indexes();
}

template <typename Gradient_Method>
void ANN_limiter<Gradient_Method>::reconstruct(const std::vector<Solution_>& solutions) {
    static const auto num_solution = solutions.size();

    This_::solution_gradients_ = Gradient_Method::calculate_solution_gradients(solutions);

    for (size_t i = 0; i < num_solution; ++i) {
        if (This_::is_constant_region(solutions, i))
            continue;

        const auto ordered_indexes = This_::ordering_function(solutions, i);

        const auto num_vertex_share_cell = This_::set_of_vertex_share_cell_indexes.at(i).size();
        std::vector<double> input_values(num_vertex_share_cell * 3);

        for (size_t j = 0; j < num_equation_; ++j) {
            for (size_t k = 0; k < num_vertex_share_cell; ++k) {
                const auto index = ordered_indexes[k];

                const auto& solution = solutions[index];
                const auto& solution_gradient = solution_gradients_[index];

                const auto solution_start_index = 0;
                const auto solution_gradient_x_start_index = num_vertex_share_cell;
                const auto solution_gradient_y_start_index = 2 * num_vertex_share_cell;

                input_values[solution_start_index + k] = solution.at(j);
                input_values[solution_gradient_x_start_index + k] = solution_gradient.at(j, 0);
                input_values[solution_gradient_y_start_index + k] = solution_gradient.at(j, 1);
            }
        }

        Dynamic_Euclidean_Vector input = std::move(input_values);
        This_::limit(input);

        This_::solution_gradients_[i] *= input[0]; //temporal code
    }
}


template <typename Gradient_Method>
bool ANN_limiter<Gradient_Method>::is_constant_region(const std::vector<Solution_>& solutions, const size_t target_cell_index) {
    const auto& vertex_share_cell_index_set = This_::set_of_vertex_share_cell_indexes.at(target_cell_index);
    const auto num_vertex_share_cell = vertex_share_cell_index_set.size();

    std::vector<double> vertex_share_cell_solutions;
    vertex_share_cell_solutions.reserve(num_vertex_share_cell);

    for (const auto vertex_share_cell_index : vertex_share_cell_index_set)
        vertex_share_cell_solutions.push_back(solutions[vertex_share_cell_index][0]);

    const auto min_solution = *std::min_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
    const auto max_solution = *std::max_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
    const auto solution_diff = max_solution - min_solution;

    return solution_diff < 0.01;
}


template <typename Gradient_Method>
std::vector<size_t> ANN_limiter<Gradient_Method>::ordering_function(const std::vector<Solution_>& solutions, const size_t target_cell_index) {
    const auto& vnode_share_cell_indexes = This_::set_of_vertex_share_cell_indexes.at(target_cell_index);

    const auto num_vnode_share_cell = vnode_share_cell_indexes.size();
    std::vector<size_t> ordered_indexes;
    ordered_indexes.reserve(num_vnode_share_cell);

    ordered_indexes.push_back(target_cell_index);

    while (ordered_indexes.size() != num_vnode_share_cell) {
        const auto& face_share_cell_indexes = This_::set_of_face_share_cell_indexes.at(ordered_indexes.back());

        std::vector<size_t> face_share_cell_indexes_in_chunk;
        std::set_intersection(vnode_share_cell_indexes.begin(), vnode_share_cell_indexes.end(), face_share_cell_indexes.begin(), face_share_cell_indexes.end(), std::back_inserter(face_share_cell_indexes_in_chunk));

        std::vector<size_t> candidate_cell_indexes;
        std::set<size_t> ordered_index_set(ordered_indexes.begin(), ordered_indexes.end());
        std::set_difference(face_share_cell_indexes_in_chunk.begin(), face_share_cell_indexes_in_chunk.end(), ordered_index_set.begin(), ordered_index_set.end(), std::back_inserter(candidate_cell_indexes));

        if (1 < candidate_cell_indexes.size()) {
            std::vector<double> temporary_solutions;
            temporary_solutions.reserve(candidate_cell_indexes.size());

            for (const auto candidate_cell_index : candidate_cell_indexes)
                temporary_solutions.push_back(solutions[candidate_cell_index].at(0));

            const auto max_solution_iter = std::max_element(temporary_solutions.begin(), temporary_solutions.end());
            const auto pos = max_solution_iter - temporary_solutions.begin();

            const auto index = *std::next(candidate_cell_indexes.begin(), pos);
            ordered_indexes.push_back(index);
        }
        else
            ordered_indexes.push_back(candidate_cell_indexes.front());
    }

    return ordered_indexes;
}

template <typename Gradient_Method>
void ANN_limiter<Gradient_Method>::limit(Dynamic_Euclidean_Vector& feature) {
    static const auto model = This_::read_model();
    static const auto num_layer = model.weights.size();
    static const ReLU activation_function;
    static const HardTanh output_function;

    for (ushort i = 0; i < num_layer; ++i) {
        Dynamic_Euclidean_Vector new_feature = model.biases[i];

        ms::gemvpv(model.weights[i], feature, new_feature);
        new_feature.apply(activation_function);
        feature = std::move(new_feature);
    }
    feature.apply(output_function);
}

template <typename Gradient_Method>
ANN_Model ANN_limiter<Gradient_Method>::read_model() {
    std::ifstream file("RSC/model.bin", std::ios::binary);
    dynamic_require(file.is_open(), "Model file should be open");

    int num_layer;
    file.read((char*)&num_layer, sizeof(int));

    ANN_Model model;
    model.weights.reserve(num_layer);
    model.biases.reserve(num_layer);

    for (int i = 0; i < num_layer; i++) {
        int num_row;
        int num_column;

        file.read((char*)&num_row, sizeof(int));
        file.read((char*)&num_column, sizeof(int));

        std::vector<double> weight(num_row * num_column);
        std::vector<double> bias(num_row);

        file.read((char*)&weight[0], sizeof(double) * num_row * num_column);
        file.read((char*)&bias[0], sizeof(double) * num_row);

        model.weights.push_back({ static_cast<size_t>(num_row), static_cast<size_t>(num_column), std::move(weight) });
        model.biases.push_back({ std::move(bias) });
    }

    return model;
}

template <ushort space_dimension_, ushort solution_order_>
void Polynomial_Reconstruction<space_dimension_, solution_order_>::initialize(const Grid<space_dimension_>& grid) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    This_::set_of_basis_functions_.reserve(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        const auto& cell_element = cell_elements[i];
        const auto& cell_geometry = cell_element.geometry_;

        This_::set_of_basis_functions_.push_back(cell_geometry.orthonormal_basis_function<solution_order_>());
    }    

    Log::content_ << std::left << std::setw(50) << "@ Polynomial reconstruction precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <ushort space_dimension_, ushort solution_order_>
auto Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_set_of_transposed_gradient_basis(void) {
    const auto num_cell = This_::set_of_basis_functions_.size();
    
    std::vector<Matrix_Function<Polynomial<space_dimension_>, space_dimension_, This_::num_basis_>> transposed_basis_Jacobians(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        auto& transposed_jacobian_basis = transposed_basis_Jacobians[i];
        const auto& basis_function = This_::set_of_basis_functions_[i];

        for (ushort j = 0; j < This_::num_basis_; ++j)
            transposed_jacobian_basis.change_column(j, basis_function[j].gradient());
    }

    return transposed_basis_Jacobians;
}

template <ushort space_dimension_, ushort solution_order_>
auto Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_basis_node(const uint cell_index, const Space_Vector_& node) {
    const auto& basis_functions = This_::set_of_basis_functions_[cell_index];

    std::array<double, This_::num_basis_> values;
    for (size_t i = 0; i < This_::num_basis_; ++i)
        values[i] = basis_functions[i](node);

    Euclidean_Vector<This_::num_basis_> basis_value = values;
    return basis_value;
}


template <ushort space_dimension_, ushort solution_order_>
Dynamic_Matrix Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_basis_nodes(const uint cell_index, const std::vector<Space_Vector_>& nodes) {
    const auto& basis_function = This_::set_of_basis_functions_[cell_index];

    const auto num_node = nodes.size();

    Dynamic_Matrix basis_nodes(This_::num_basis_, num_node);
    for (size_t j = 0; j < num_node; ++j)
        basis_nodes.change_column(j, basis_function(nodes[j]));

    return basis_nodes;
}

template <ushort space_dimension_, ushort solution_order_>
double Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_P0_basis_value(const uint cell_index, const Space_Vector_& center_node) {
    const auto& basis_function = This_::set_of_basis_functions_[cell_index];
    const auto& P0_basis_function = basis_function[0];

    return P0_basis_function(center_node);
}

template <ushort space_dimension_, ushort solution_order_>
std::vector<Dynamic_Matrix> Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_set_of_basis_nodes(const std::vector<std::vector<Space_Vector_>>& set_of_nodes) {
    const auto num_cell = This_::set_of_basis_functions_.size();

    dynamic_require(set_of_nodes.size() == num_cell, "set size should be same with num cell");

    std::vector<Dynamic_Matrix> set_of_basis_nodes;
    set_of_basis_nodes.reserve(num_cell);

    for (uint i = 0; i < num_cell; ++i) 
        set_of_basis_nodes.push_back(This_::calculate_basis_nodes(i, set_of_nodes[i]));
    

    return set_of_basis_nodes;
}

template <ushort space_dimension_, ushort solution_order_>
std::vector<double> Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_P0_basis_values(const std::vector<Space_Vector_>& center_nodes) {
    const auto num_cell = center_nodes.size();
    std::vector<double> P0_basis_values(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        const auto& basis_function = This_::set_of_basis_functions_[i];
        const auto& P0_basis_function = basis_function[0];
        const auto& cell_center = center_nodes[i];

        P0_basis_values[i] = P0_basis_function(cell_center);
    }

    return P0_basis_values;
}





//template <typename Gradient_Method>
//AI_limiter<Gradient_Method>::AI_limiter(const Grid<space_dimension_>& grid) : gradient_method_(grid) {
//    const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;
//    const auto& cell_elements = grid.elements.cell_elements;
//    this->num_data_ = cell_elements.size();
//
//    this->vertex_share_cell_indexes_set_.reserve(this->num_data_);
//    this->ai_limiter_text_set_.resize(this->num_data_);
//    this->target_cell_indexes_.reserve(this->num_data_);
//
//    const auto face_share_cell_indexes_set = this->calculate_face_share_cell_indexes_set(grid);
//    //const auto vnodes_coordinate_string_set = this->calculate_vertex_nodes_coordinate_string_set(grid);
//
//    for (size_t i = 0; i < this->num_data_; ++i) {
//        const auto& cell_element = cell_elements[i];
//        const auto& cell_geometry = cell_element.geometry_;
//
//        //vertex share cell indexes temp
//        std::set<size_t> vertex_share_cell_indexes_temp;
//
//        const auto vnode_indexes = cell_element.vertex_node_indexes();
//        for (const auto& vnode_index : vnode_indexes) {
//            const auto& vnode_share_cell_indexes = vnode_index_to_share_cell_indexes.at(vnode_index);
//            vertex_share_cell_indexes_temp.insert(vnode_share_cell_indexes.begin(), vnode_share_cell_indexes.end());
//        }
//
//        //chunk edge connectivities //quad3에서는 제대로 작동하지 않는 algorithm
//        std::set<std::set<size_t>> face_share_cell_index_pairs;
//
//        for (const auto chunk_cell_index : vertex_share_cell_indexes_temp) {
//            const auto& face_share_cell_indexes = face_share_cell_indexes_set.at(chunk_cell_index);
//
//            std::vector<size_t> face_share_cell_indexes_in_chunk;
//            std::set_intersection(vertex_share_cell_indexes_temp.begin(), vertex_share_cell_indexes_temp.end(), face_share_cell_indexes.begin(), face_share_cell_indexes.end(), std::back_inserter(face_share_cell_indexes_in_chunk));
//
//            for (const auto face_share_cell_index_in_chunk : face_share_cell_indexes_in_chunk)
//                face_share_cell_index_pairs.insert({ chunk_cell_index, face_share_cell_index_in_chunk });
//        }
//
//        //vertex_share_cell_indexes
//        vertex_share_cell_indexes_temp.erase(i);
//
//        std::vector<size_t> vertex_share_cell_indexes;
//        vertex_share_cell_indexes.push_back(i);
//        vertex_share_cell_indexes.insert(vertex_share_cell_indexes.end(), vertex_share_cell_indexes_temp.begin(), vertex_share_cell_indexes_temp.end());
//
//        // header string
//        this->ai_limiter_text_set_[i] << "#########################";
//
//        // node number string
//        const auto num_node = vertex_share_cell_indexes.size();
//        this->ai_limiter_text_set_[i] << "@nodeNumber\n" + std::to_string(num_node);
//
//        // node index order string
//        std::string node_index_order_string = "@nodeIndexOrder\n";
//        for (const auto& node_index : vertex_share_cell_indexes)
//            node_index_order_string += std::to_string(node_index) + "\t";
//        this->ai_limiter_text_set_[i] << std::move(node_index_order_string);
//
//        //edge number string
//        const auto num_edge = face_share_cell_index_pairs.size();
//        this->ai_limiter_text_set_[i] << "@edgeNumber\n" + std::to_string(num_edge);
//
//        //connectivity string
//        std::string node_connectivity_string = "@connectivity\n";
//        for (const auto& chunk_edge_connectivity : face_share_cell_index_pairs) {
//            for (const auto& node_index : chunk_edge_connectivity)
//                node_connectivity_string += std::to_string(node_index) + "\t";
//            node_connectivity_string += "\n";
//        }
//        node_connectivity_string.pop_back();
//        this->ai_limiter_text_set_[i] << std::move(node_connectivity_string);
//
//        ////cell coords string
//        //std::string cell_coords_string = "@cellCoords\n";
//        //for (const auto vertex_share_cell_index : vertex_share_cell_indexes)
//        //    cell_coords_string += vnodes_coordinate_string_set[vertex_share_cell_index];
//        //cell_coords_string.pop_back();
//
//        //this->ai_limiter_text_set_[i] << std::move(cell_coords_string);
//
//        //vertex_share_cell_indexes
//        this->vertex_share_cell_indexes_set_.push_back(std::move(vertex_share_cell_indexes));
//    }
//}
//
//template <typename Gradient_Method>
//auto AI_limiter<Gradient_Method>::reconstruct_solutions(const std::vector<Solution_>& solutions) {
//    auto solution_gradients = this->gradient_method_.calculate_solution_gradients(solutions);
//
//    this->record_solution_datas(solutions, solution_gradients);
//    auto ai_limiter_str = this->make_ai_limiter_str();
//
//
//    //dynamic matrix to matrix
//    std::vector<Matrix<num_equation_, space_dimension_>> limited_solution_gradient;
//    //limited_solution_gradient.reserve(num_cell);
//
//    for (const auto& solution_gradient : solution_gradients)
//        limited_solution_gradient.push_back(solution_gradient);
//
//    return Linear_Reconstructed_Solution<num_equation_, space_dimension_>{ solutions, limited_solution_gradient };
//}
//
//template <typename Gradient_Method>
//auto AI_limiter<Gradient_Method>::calculate_face_share_cell_indexes_set(const Grid<space_dimension_>& grid) const {
//    const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;
//    const auto& cell_elements = grid.elements.cell_elements;
//    const auto num_cell = cell_elements.size();
//
//    //face share cell indexes set
//    std::vector<std::set<size_t>> face_share_cell_indexes_set;
//    face_share_cell_indexes_set.reserve(num_cell);
//
//    for (size_t i = 0; i < num_cell; ++i) {
//        const auto& element = cell_elements[i];
//        const auto& geometry = cell_elements[i].geometry_;
//
//        const auto face_vnode_indexes_set = element.face_vertex_node_indexes_set();
//        const auto num_face = face_vnode_indexes_set.size();
//
//        std::set<size_t> face_share_cell_indexes;
//
//        for (const auto& face_vnode_indexes : face_vnode_indexes_set) {
//            std::vector<size_t> this_face_share_cell_indexes;
//
//            const auto num_face_vnode = face_vnode_indexes.size();
//
//            const auto& set_0 = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[0]);
//            const auto& set_1 = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[1]);
//            std::set_intersection(set_0.begin(), set_0.end(), set_1.begin(), set_1.end(), std::back_inserter(this_face_share_cell_indexes));
//
//            if (2 < num_face_vnode) {
//                std::vector<size_t> buffer;
//                for (size_t i = 2; i < num_face_vnode; ++i) {
//                    const auto& set_i = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[i]);
//
//                    buffer.clear();
//                    std::set_intersection(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), set_i.begin(), set_i.end(), std::back_inserter(buffer));
//                    std::swap(this_face_share_cell_indexes, buffer);
//                }
//            }
//
//            const auto my_index_pos_iter = std::find(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), i);
//            dynamic_require(my_index_pos_iter != this_face_share_cell_indexes.end(), "my index should be included in this face share cell indexes");
//
//            this_face_share_cell_indexes.erase(my_index_pos_iter);
//            dynamic_require(this_face_share_cell_indexes.size() == 1, "face share cell should be unique");
//
//            face_share_cell_indexes.insert(this_face_share_cell_indexes.front());
//        }
//
//        face_share_cell_indexes_set.push_back(std::move(face_share_cell_indexes));
//    }
//
//    return face_share_cell_indexes_set;
//}
//
//
//template <typename Gradient_Method>
//auto AI_limiter<Gradient_Method>::calculate_vertex_nodes_coordinate_string_set(const Grid<space_dimension_>& grid) const {
//    const auto& cell_elements = grid.elements.cell_elements;
//    const auto num_cell = cell_elements.size();
//
//    std::vector<std::string> vnodes_coordinate_string_set;
//    vnodes_coordinate_string_set.reserve(num_cell);
//
//    for (size_t i = 0; i < num_cell; ++i) {
//        const auto& element = cell_elements[i];
//        const auto& geometry = element.geometry_;
//
//        const auto vnodes = geometry.vertex_nodes();
//        const auto num_vnode = vnodes.size();
//        std::string vnodes_coordinate_string = std::to_string(num_vnode) + "\n";
//
//        for (const auto& vnode : vnodes) {
//            for (size_t i = 0; i < space_dimension_; ++i)
//                vnodes_coordinate_string += ms::double_to_str_sp(vnode[i]) + "\t";
//
//            vnodes_coordinate_string += "\n";
//        }
//
//        vnodes_coordinate_string_set.push_back(std::move(vnodes_coordinate_string));
//    }
//
//    return vnodes_coordinate_string_set;
//}
//
//
//template <typename Gradient_Method>
//void AI_limiter<Gradient_Method>::record_solution_datas(const std::vector<Solution_>& solutions, const std::vector<Dynamic_Matrix_>& solution_gradients) {
//    dynamic_require(num_data_ == solutions.size(), "number of solution should be same with number of data");
//    dynamic_require(num_data_ == solution_gradients.size(), "number of solution gradient should be same with number of data");
//
//    const auto solution_strings = this->convert_to_solution_strings(solutions);
//    const auto solution_gradient_strings = this->convert_to_solution_gradient_strings(solution_gradients);
//
//    std::string cell_average_string;
//    std::string cell_gradient_string;
//    for (size_t i = 0; i < num_data_; ++i) {
//        const auto& vertex_share_cell_indexes = vertex_share_cell_indexes_set_[i];
//        const auto num_vertex_share_cell = vertex_share_cell_indexes.size();
//
//
//        std::vector<double> vertex_share_cell_solutions(num_vertex_share_cell);
//        for (size_t i = 0; i < num_vertex_share_cell; ++i)
//            vertex_share_cell_solutions[i] = solutions[vertex_share_cell_indexes[i]][0];
//
//        const auto min_solution = *std::min_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
//        const auto max_solution = *std::max_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
//        const auto solution_diff = max_solution - min_solution;
//
//        if (solution_diff < 0.01)
//            continue;
//
//        this->target_cell_indexes_.push_back(i);
//
//        cell_average_string = "@cellAverage\n";
//        cell_gradient_string = "@cellGradient\n";
//        for (const auto vertex_share_cell_index : vertex_share_cell_indexes) {
//            cell_average_string += solution_strings[vertex_share_cell_index] + "\n";
//            cell_gradient_string += solution_gradient_strings[vertex_share_cell_index] + "\n";
//        }
//        cell_average_string.pop_back();
//        cell_gradient_string.pop_back();
//
//        this->ai_limiter_text_set_[i] << std::move(cell_average_string) << std::move(cell_gradient_string);
//    }
//}
//
//template <typename Gradient_Method>
//auto AI_limiter<Gradient_Method>::convert_to_solution_strings(const std::vector<Solution_>& solutions) const {
//    const auto num_solution = solutions.size();
//
//    std::vector<std::string> solution_strings;
//    solution_strings.reserve(num_solution);
//
//    std::string solution_string;
//    for (size_t i = 0; i < num_solution; ++i) {
//
//        const auto& solution = solutions[i];
//        for (size_t j = 0; j < num_equation_; ++j)
//            solution_string += ms::double_to_str_sp(solution[j]) + "\t";
//
//        solution_strings.push_back(std::move(solution_string));
//    }
//
//    return solution_strings;
//}
//
//template <typename Gradient_Method>
//auto AI_limiter<Gradient_Method>::convert_to_solution_gradient_strings(const std::vector<Dynamic_Matrix_>& solution_gradients) const {
//    const auto num_solution = solution_gradients.size();
//    const auto [num_equation, space_dimension] = solution_gradients.front().size();
//
//    std::vector<std::string> solution_gradient_strings;
//    solution_gradient_strings.reserve(num_solution);
//
//    std::string solution_gradient_string;
//    for (size_t i = 0; i < num_solution; ++i) {
//
//        const auto& solution_gradient = solution_gradients[i];
//        for (size_t j = 0; j < num_equation; ++j)
//            for (size_t k = 0; k < space_dimension; ++k)
//                solution_gradient_string += ms::double_to_str_sp(solution_gradient.at(j, k)) + "\t";
//
//        solution_gradient_strings.push_back(std::move(solution_gradient_string));
//    }
//
//    return solution_gradient_strings;
//}
//
//template <typename Gradient_Method>
//auto AI_limiter<Gradient_Method>::make_ai_limiter_str(void) {
//    std::string ai_limiter_str;
//
//    constexpr size_t num_solution_str = 2;
//    for (const auto target_cell_index : this->target_cell_indexes_) {
//        auto& ai_limiter_text = this->ai_limiter_text_set_.at(target_cell_index);
//
//        for (const auto& str : ai_limiter_text)
//            ai_limiter_str += str + "\n";
//
//        ai_limiter_text.erase(ai_limiter_text.end() - num_solution_str, ai_limiter_text.end());
//    }
//
//    this->target_cell_indexes_.clear();
//
//    return ai_limiter_str;
//}
