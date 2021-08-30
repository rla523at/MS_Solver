#pragma once
#include "Spatial_Discrete_Method.h"
#include "Gradient_Method.h"
#include "PostAI.h"
#include "Tecplot.h"


class RM {};


class Constant_Reconstruction : public RM 
{
public:
    template <size_t space_dimension>
    Constant_Reconstruction(const Grid<space_dimension>& grid) {};
    static std::string name(void) { return "Constant_Reconstruction"; };
};


template <typename Gradient_Method>
class Linear_Reconstruction : public RM
{
private:
    static constexpr size_t num_equation_       = Gradient_Method::num_equation();
    static constexpr size_t space_dimension_    = Gradient_Method::space_dimension();

    using Solution_             = Euclidean_Vector<num_equation_>;
    using Solution_Gradient_    = Matrix<num_equation_, space_dimension_>;

private:
    Gradient_Method gradient_method_;
    std::vector<Solution_Gradient_> solution_gradients_;

public:
    Linear_Reconstruction(const Grid<space_dimension_>& grid) : gradient_method_(grid) {};
    
    void reconstruct(const std::vector<Solution_>& solutions);
        
    const std::vector<Solution_Gradient_>& get_solution_gradients(void) const { return solution_gradients_; };

public:
    static std::string name(void) { return "Linear_Reconstruction_" + Gradient_Method::name(); };
};


class MLP_u1_Limiting_Strategy
{
private:
    MLP_u1_Limiting_Strategy(void) = delete;

public:
    static double calculate_limiting_value(const double P1_mode_solution, const double P0_solution, const double allowable_min, const double allowable_max) {
        if (P1_mode_solution == 0)
            return 1.0;

        if (P1_mode_solution < 0)
            return (std::min)((allowable_min - P0_solution) / P1_mode_solution, 1.0);
        else
            return (std::min)((allowable_max - P0_solution) / P1_mode_solution, 1.0);
    };
};


template <typename Gradient_Method>
class MLP_u1 : public RM
{
private:
    static constexpr size_t num_equation_ = Gradient_Method::num_equation();
    static constexpr size_t space_dimension_ = Gradient_Method::space_dimension();

    using Solution_ = Euclidean_Vector<num_equation_>;
    using Solution_Gradient_ = Matrix<num_equation_, space_dimension_>;

private:
    Gradient_Method gradient_method_;
    std::unordered_map<uint, std::set<uint>> vnode_index_to_share_cell_indexes_;
    std::vector<std::vector<uint>> vnode_indexes_set_;
    std::vector<Dynamic_Matrix> center_to_vertex_matrixes_;
    std::vector<Solution_Gradient_> solution_gradients_;

public:
    MLP_u1(Grid<space_dimension_>&& grid);

    void reconstruct(const std::vector<Solution_>& solutions);
        
public:
    const std::vector<Solution_Gradient_>& get_solution_gradients(void) const { return this->solution_gradients_; };

protected:
    auto calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution_>& solutions) const;

public:
    static std::string name(void) { return "MLP_u1_" + Gradient_Method::name(); };
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

    using Solution_             = Euclidean_Vector<num_equation_>;
    using Solution_Gradient_    = Matrix<num_equation_, space_dimension_>;

private:
    Gradient_Method gradient_method_;
    std::vector<std::vector<size_t>> set_of_vertex_share_cell_indexes_;
    std::vector<std::vector<size_t>> set_of_face_share_cell_indexes_;
    std::vector<Solution_Gradient_> solution_gradients_;

    std::vector<size_t> ordered_indexes_;   //�߰�

public:
    ANN_limiter(const Grid<space_dimension_>& grid);

public:
    void reconstruct(const std::vector<Solution_>& solutions);    
    const std::vector<Solution_Gradient_>& get_solution_gradients(void) const { return solution_gradients_; };

private:
    bool is_constant_region(const std::vector<Solution_>& solutions, const size_t target_cell_index) const;
<<<<<<< HEAD
    std::vector<size_t> ordering_function_using_solutions(const std::vector<Solution_>& solutions, const size_t target_cell_index) const;
    std::vector<size_t> ordering_function_using_cell_indexes(const size_t target_cell_index) const;

=======
    std::vector<size_t> ordering_function(const std::vector<Solution_>& solutions, const size_t target_cell_index) const;
>>>>>>> ms/dev/DFM
    void limit(Dynamic_Euclidean_Vector& feature) const;
    ANN_Model read_model(void) const;

public:
    static std::string name(void) { return "ANN_Reconstruction_" + Gradient_Method::name(); };
};





namespace ms {
	template <typename T>
	inline constexpr bool is_reconsturction_method = std::is_base_of_v<RM, T>;

    template <typename T>
    inline constexpr bool is_constant_reconustruction = std::is_same_v<Constant_Reconstruction, T>;


    template <typename Spatial_Discrete_Method, typename Reconstruction_Method, typename = void>
    inline constexpr bool is_default_reconstruction;
        
    template <typename Spatial_Discrete_Method, typename Reconstruction_Method>
    inline constexpr bool is_default_reconstruction<typename Spatial_Discrete_Method, typename Reconstruction_Method, std::enable_if_t<std::is_same_v<FVM, Spatial_Discrete_Method>>>
        = ms::is_constant_reconustruction<Reconstruction_Method>;
}


//template definition part
template <typename Gradient_Method>
void Linear_Reconstruction<Gradient_Method>::reconstruct(const std::vector<Solution_>& solutions) {
    this->solution_gradients_ = this->gradient_method_.calculate_solution_gradients(solutions);
}

template <typename Gradient_Method>
MLP_u1<Gradient_Method>::MLP_u1(Grid<space_dimension_>&& grid) : gradient_method_(grid) {
    SET_TIME_POINT;

    //vnode index to share cell indexes
    this->vnode_index_to_share_cell_indexes_ = std::move(grid.connectivity.vnode_index_to_share_cell_index_set);

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    this->vnode_indexes_set_.reserve(num_cell);
    this->center_to_vertex_matrixes_.reserve(num_cell);
    this->solution_gradients_.resize(num_cell);

    for (size_t i = 0; i < num_cell; ++i) {
        const auto& element = cell_elements[i];
        const auto& geometry = element.geometry_;

        // vnode indexes set
        this->vnode_indexes_set_.push_back(element.vertex_node_indexes());

        //center to vertex matrix
        const auto center_node = geometry.center_node();
        const auto vertex_nodes = geometry.vertex_nodes();
        const auto num_vertex = vertex_nodes.size();

        Dynamic_Matrix center_to_vertex_matrix(space_dimension_, num_vertex);
        for (size_t i = 0; i < num_vertex; ++i) {
            const auto center_to_vertex = vertex_nodes[i] - center_node;
            center_to_vertex_matrix.change_column(i, center_to_vertex);
        }

        this->center_to_vertex_matrixes_.push_back(std::move(center_to_vertex_matrix));
    }

    Log::content_ << std::left << std::setw(50) << "@ MLP u1 precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}


template <typename Gradient_Method>
void MLP_u1<Gradient_Method>::reconstruct(const std::vector<Solution_>& solutions) {
    const auto solution_gradients = this->gradient_method_.calculate_solution_gradients(solutions);
    const auto vnode_index_to_min_max_solution = this->calculate_vertex_node_index_to_min_max_solution(solutions);
    const auto num_cell = solutions.size();

    for (uint i = 0; i < num_cell; ++i) {
        const auto& gradient = solution_gradients[i];
        const auto P1_mode_solution_vnodes = gradient * this->center_to_vertex_matrixes_[i];

        std::array<double, num_equation_> limiting_values;
        limiting_values.fill(1);

        const auto& vnode_indexes = this->vnode_indexes_set_[i];
        const auto num_vertex = vnode_indexes.size();

        for (ushort j = 0; j < num_vertex; ++j) {
            const auto vnode_index = vnode_indexes[j];
            const auto& [min_solution, max_solution] = vnode_index_to_min_max_solution.at(vnode_index);

            for (ushort e = 0; e < num_equation_; ++e) {
                const auto limiting_value = MLP_u1_Limiting_Strategy::calculate_limiting_value(P1_mode_solution_vnodes.at(e, j), solutions[i].at(e), min_solution.at(e), max_solution.at(e));
                limiting_values[e] = min(limiting_values[e], limiting_value);
            }
        }
                
        const Matrix limiting_value_matrix = limiting_values;
        this->solution_gradients_[i] = limiting_value_matrix * gradient;
    }
};


template <typename Gradient_Method>
auto MLP_u1<Gradient_Method>::calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution_>& solutions) const {
    const auto num_vnode = this->vnode_index_to_share_cell_indexes_.size();

    std::unordered_map<uint, std::pair<Solution_, Solution_>> vnode_index_to_min_max_solution;
    vnode_index_to_min_max_solution.reserve(num_vnode);

    for (const auto& [vnode_index, share_cell_indexes] : this->vnode_index_to_share_cell_indexes_) {
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
ANN_limiter<Gradient_Method>::ANN_limiter(const Grid<space_dimension_>& grid) :gradient_method_(grid) {
    this->set_of_face_share_cell_indexes_    = grid.calculate_set_of_face_share_cell_indexes();
    this->set_of_vertex_share_cell_indexes_  = grid.calculate_set_of_vertex_share_cell_indexes();

    //sorting
    for (auto& face_share_cell_indexes : this->set_of_face_share_cell_indexes_)
        std::sort(face_share_cell_indexes.begin(), face_share_cell_indexes.end());
    for (auto& vertex_share_cell_indexes : this->set_of_vertex_share_cell_indexes_)
        std::sort(vertex_share_cell_indexes.begin(), vertex_share_cell_indexes.end());

}

template <typename Gradient_Method>
void ANN_limiter<Gradient_Method>::reconstruct(const std::vector<Solution_>& solutions) {
    this->solution_gradients_ = this->gradient_method_.calculate_solution_gradients(solutions);

    const auto num_solution = solutions.size();    

    for (size_t i = 0; i < num_solution; ++i) {

        if (this->is_constant_region(solutions, i)) {
            //this->solution_gradients_[i] *= 0.0; //temporal code
            continue;
        }
            
        //const auto ordered_indexes = this->ordering_function_using_solutions(solutions, i);
        const auto ordered_indexes = this->ordering_function_using_cell_indexes(i);
        const auto num_ordered_index = ordered_indexes.size();

        const auto num_input_values = 3 * num_ordered_index; //temporal code
        std::vector<double> input_values(num_input_values);

        for (size_t j = 0; j < num_equation_; ++j) {
            for (size_t k = 0; k < num_ordered_index; ++k) {
                const auto index = ordered_indexes[k];

                const auto& solution = solutions[index];
                const auto& solution_gradient = solution_gradients_[index];

                const auto solution_start_index = 0;
                const auto solution_gradient_x_start_index = num_ordered_index;
                const auto solution_gradient_y_start_index = 2 * num_ordered_index;

                input_values[solution_start_index + k] = solution.at(j);
                input_values[solution_gradient_x_start_index + k] = solution_gradient.at(j, 0);
                input_values[solution_gradient_y_start_index + k] = solution_gradient.at(j, 1);
            }
        }

        Dynamic_Euclidean_Vector input = std::move(input_values);
        this->limit(input);
        this->solution_gradients_[i] *= input[0]; //temporal code                
    }
}


template <typename Gradient_Method>
bool ANN_limiter<Gradient_Method>::is_constant_region(const std::vector<Solution_>& solutions, const size_t target_cell_index) const {
    const auto& vertex_share_cell_index_set = this->set_of_vertex_share_cell_indexes_.at(target_cell_index);
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
std::vector<size_t> ANN_limiter<Gradient_Method>::ordering_function_using_solutions(const std::vector<Solution_>& solutions, const size_t target_cell_index) const {
    const auto& vnode_share_cell_indexes = this->set_of_vertex_share_cell_indexes_.at(target_cell_index);

    const auto num_vnode_share_cell = vnode_share_cell_indexes.size();
    const auto num_ordered_indexes = num_vnode_share_cell + 1; // include target cell

    std::vector<size_t> ordered_indexes;
    ordered_indexes.reserve(num_ordered_indexes);

    ordered_indexes.push_back(target_cell_index);

    while (ordered_indexes.size() != num_ordered_indexes) {
        const auto& face_share_cell_indexes = this->set_of_face_share_cell_indexes_.at(ordered_indexes.back());

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

            std::reverse(temporary_solutions.begin(), temporary_solutions.end());
            std::reverse(candidate_cell_indexes.begin(), candidate_cell_indexes.end());

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
std::vector<size_t> ANN_limiter<Gradient_Method>::ordering_function_using_cell_indexes(const size_t target_cell_index) const {
    const auto& vnode_share_cell_indexes = this->set_of_vertex_share_cell_indexes_.at(target_cell_index);

    const auto num_vnode_share_cell = vnode_share_cell_indexes.size();
    const auto num_ordered_indexes = num_vnode_share_cell + 1; // include target cell

    std::vector<size_t> ordered_indexes;
    ordered_indexes.reserve(num_ordered_indexes);

    ordered_indexes.push_back(target_cell_index);

    while (ordered_indexes.size() != num_ordered_indexes) {
        const auto& face_share_cell_indexes = this->set_of_face_share_cell_indexes_.at(ordered_indexes.back());

        std::vector<size_t> face_share_cell_indexes_in_chunk;
        std::set_intersection(vnode_share_cell_indexes.begin(), vnode_share_cell_indexes.end(), face_share_cell_indexes.begin(), face_share_cell_indexes.end(), std::back_inserter(face_share_cell_indexes_in_chunk));

        std::vector<size_t> candidate_cell_indexes;
        std::set<size_t> ordered_index_set(ordered_indexes.begin(), ordered_indexes.end());
        std::set_difference(face_share_cell_indexes_in_chunk.begin(), face_share_cell_indexes_in_chunk.end(), ordered_index_set.begin(), ordered_index_set.end(), std::back_inserter(candidate_cell_indexes));

        if (1 < candidate_cell_indexes.size()) {

            std::reverse(candidate_cell_indexes.begin(), candidate_cell_indexes.end());

            const auto max_index_iter = std::max_element(candidate_cell_indexes.begin(), candidate_cell_indexes.end());
            const auto pos = max_index_iter - candidate_cell_indexes.begin();

            const auto index = *std::next(candidate_cell_indexes.begin(), pos);
            ordered_indexes.push_back(index);
        }
        else
            ordered_indexes.push_back(candidate_cell_indexes.front());
    }

    return ordered_indexes;
}

template <typename Gradient_Method>
void ANN_limiter<Gradient_Method>::limit(Dynamic_Euclidean_Vector& feature) const {
    static const auto model = this->read_model();
    static const auto num_layer = model.weights.size();
    static const ReLU activation_function;
    static const HardTanh output_function;

    //hidden layer
    for (ushort i = 0; i < num_layer - 1; ++i) {
        auto new_feature = model.biases[i];

        ms::gemvpv(model.weights[i], feature, new_feature);
        new_feature.apply(activation_function);
        feature = std::move(new_feature);
    }

    //output layer
    auto new_feature = model.biases[num_layer - 1];

    ms::gemvpv(model.weights[num_layer - 1], feature, new_feature);
    new_feature.apply(output_function);
    feature = std::move(new_feature);
}

template <typename Gradient_Method>
ANN_Model ANN_limiter<Gradient_Method>::read_model(void) const {
    //std::ifstream file("RSC/model.bin", std::ios::binary);
    std::ifstream file("RSC/case17.bin", std::ios::binary);

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