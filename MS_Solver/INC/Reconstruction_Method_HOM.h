#pragma once
#include "Debugger.h" //debug
#include "Reconstruction_Method_FVM.h"
#include "Polynomial.h"
#include "Setting.h"

class HOM_Reconstruction : public RM {};

template <ushort space_dimension_, ushort solution_order_>
class Polynomial_Reconstruction : public HOM_Reconstruction
{
protected:
    using This_         = Polynomial_Reconstruction<space_dimension_, solution_order_>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;

    static constexpr ushort num_basis_ = ms::combination_with_repetition(1 + space_dimension_, solution_order_);

public:
    static constexpr ushort num_basis(void) { return num_basis_; };
    static constexpr ushort solution_order(void) { return solution_order_; };
    static constexpr ushort space_dimension(void) { return space_dimension_; };
    static std::string name(void) { return "Polynomial_Reconstruction_P" + std::to_string(solution_order_); };

protected:
    std::vector<Vector_Function<Polynomial<space_dimension_>, num_basis_>> basis_vector_functions_;

public:
    Polynomial_Reconstruction(const Grid<space_dimension_>& grid);

public:
    auto calculate_set_of_transposed_gradient_basis(void) const;
    auto calculate_basis_node(const uint cell_index, const Space_Vector_& node) const;
    Dynamic_Matrix basis_nodes(const uint cell_index, const std::vector<Space_Vector_>& nodes) const;    
    double calculate_P0_basis_value(const uint cell_index) const;
    std::vector<double> calculate_P0_basis_values(void) const;

public:
    template<ushort num_equation>
    auto solution_function(const uint cell_index, const Matrix<num_equation, num_basis_>& coefficient) const;
};


class P1_Projected_MLP_Condition
{
private:
    P1_Projected_MLP_Condition(void) = delete;

public:
    static bool is_satisfy(const double P1_projected_value, const double allowable_min, const double allowable_max) {
        return allowable_min <= P1_projected_value && P1_projected_value <= allowable_max;
    }
};


class MLP_Smooth_Extrema_Detector
{
private:
    MLP_Smooth_Extrema_Detector(void) = delete;

public:
    static bool is_smooth_extrema(const double solution, const double higher_mode_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) {
        if (P1_mode_solution > 0 && higher_mode_solution < 0 && solution > allowable_min)
            return true;
        else if (P1_mode_solution < 0 && higher_mode_solution > 0 && solution < allowable_max)
            return true;
        else
            return false;
    }
};


class Constant_Region_Detector
{
private:
    Constant_Region_Detector(void) = delete;

public:
    static bool is_constant(const double solution, const double P0_solution, const double volume) {
        const auto constant_criterion = (std::max)(1.0E-3 * std::abs(P0_solution), volume);
        return std::abs(solution - P0_solution) <= constant_criterion;
    }
};


template <ushort space_dimension_, ushort solution_order_>
class hMLP_Base : public Polynomial_Reconstruction<space_dimension_, solution_order_>
{
private:
    static_require(1 <= solution_order_, "solution order should be greater than 1 for hMLP reconstruction");

    using This_ = hMLP_Base<space_dimension_, solution_order_>;

protected:
    static constexpr ushort criterion_variable_index_ = 0;    

protected:
    const std::unordered_map<uint, std::set<uint>>& vnode_index_to_share_cell_index_set_;
    std::array<ushort, solution_order_ + 1> num_Pn_projection_basis_;
    std::vector<std::vector<uint>> set_of_vnode_indexes_;
    std::vector<double> volumes_;
    std::vector<double> P0_basis_values_;

protected:
    hMLP_Base(const Grid<space_dimension_>& grid);

protected:
    auto Pn_projection_matrix(const ushort Pn) const;
    auto Pn_mode_matrix(const ushort Pn) const;
    auto limiting_matrix(const double limiting_value) const;
};

template <ushort space_dimension_, ushort solution_order_>
class hMLP_Reconstruction : public hMLP_Base<space_dimension_, solution_order_>
{
private:
    using This_ = hMLP_Reconstruction<space_dimension_, solution_order_>;

public:
    static std::string name(void) { return "hMLP_Reconstruction_P" + std::to_string(solution_order_); };

protected:
    std::vector<Dynamic_Matrix> set_of_basis_vnodes_;
    std::vector<Dynamic_Matrix> set_of_P1_projected_basis_vnodes_;     
    
public:
    hMLP_Reconstruction(const Grid<space_dimension_>& grid);

public:
    template <ushort num_equation>
    void reconstruct(std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const;

protected:
    template <ushort num_equation>
    auto calculate_P0_criterion_values(const std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const;        
    auto calculate_vertex_node_index_to_allowable_min_max_criterion_value(const std::vector<double>& P0_criterion_values) const;
};


enum class BD_Type {
    standard, noBD, improved,
    no_typeI, typeI_1, typeI_2, typeI_3, typeI_4,
    no_typeII, typeI_1_without_typeII, typeI_4_without_typeII, typeI_1_cc, typeI_1_dc,
    noBD_dc
};


template <ushort space_dimension_, ushort solution_order_>
class hMLP_BD_Reconstruction : public hMLP_Base<space_dimension_, solution_order_>
{
private:
    static constexpr ushort num_P1_projected_basis_ = ms::combination_with_repetition(1 + space_dimension_, 1);
    
    using This_ = hMLP_BD_Reconstruction<space_dimension_, solution_order_>;

private:
    std::unordered_map<uint, std::set<uint>> pbdry_vnode_index_to_matched_vnode_index_set_;
    std::vector<std::vector<uint>> set_of_face_share_cell_indexes_;

    std::vector<Dynamic_Matrix> set_of_basis_vnodes_;
    std::vector<Dynamic_Matrix> set_of_simplex_P1_projected_basis_vnodes_;
    std::vector<Dynamic_Matrix> set_of_simplex_P0_projected_basis_vnodes_;

    std::vector<double> face_volumes_;
    std::vector<double> face_characteristic_lengths_;
    std::vector<std::pair<uint, uint>> face_oc_nc_index_pairs_;
    std::vector<std::pair<Dynamic_Matrix, Dynamic_Matrix>> face_oc_nc_side_basis_jump_qnodes_pairs_;    
    std::vector<Dynamic_Euclidean_Vector> set_of_jump_qweights_;

public:
    hMLP_BD_Reconstruction(const Grid<space_dimension_>& grid);

public:
    template <ushort num_equation>
    void reconstruct(std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const;

private:
    template <ushort num_equation>
    auto calculate_set_of_vertex_node_index_to_simplex_P0_criterion_value(const std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const;
    auto calculate_vertex_node_index_to_allowable_min_max_criterion_value(const std::vector<std::map<uint, double>>& cell_index_to_vnode_index_to_simplex_P0_criterion_values) const;
    template <ushort num_equation>
    std::vector<ushort> calculate_set_of_num_troubled_boundary(const std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const;
    template <ushort projection_order>
    auto calculate_simplex_Pn_projection_basis_vector_function(const uint cell_index, const Geometry<space_dimension_>& sub_simplex_geometry) const;
    bool is_typeI_subcell_oscillation(const ushort num_trouble_boundaries) const;
    bool is_typeII_subcell_oscillation(const ushort num_trouble_boundary) const;   
    
    template <ushort num_equation>
    std::vector<Euclidean_Vector<num_equation>> calculate_P0_solutions(const std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const;

    template <ushort num_equation>
    std::vector<bool> is_shock(const std::vector<Euclidean_Vector<num_equation>>& P0_solutions) const;//temporary code

    template <ushort num_equation>
    std::vector<bool> is_discontinuity(const std::vector<Euclidean_Vector<num_equation>>& P0_solutions) const;//temporary code

public:
    static std::string name(void) { return "hMLP_BD_Reconstruction_P" + std::to_string(solution_order_); };
};


namespace ms {
    template <typename T>
    inline constexpr bool is_HOM_reconsturction_method = std::is_base_of_v<HOM_Reconstruction, T>;

    template <typename T>
    inline constexpr bool is_polynomial_reconustruction = std::is_same_v<Polynomial_Reconstruction<T::space_dimension(), T::solution_order()>, T>;

    template <typename Spatial_Discrete_Method, typename Reconstruction_Method>
    inline constexpr bool is_default_reconstruction<typename Spatial_Discrete_Method, typename Reconstruction_Method, std::enable_if_t<std::is_same_v<HOM, Spatial_Discrete_Method>>>
        = ms::is_polynomial_reconustruction<Reconstruction_Method>;

    template <typename T, typename... Args>
    std::vector<T>& merge(std::vector<T>& vec1, std::vector<T>&& vec2, Args&&... args) {
        static_require((... && std::is_same_v<Args, std::vector<T>>), "every arguments should be vector of same type");

        vec1.reserve(vec1.size() + vec2.size());
        vec1.insert(vec1.end(), std::make_move_iterator(vec2.begin()), std::make_move_iterator(vec2.end()));

        if constexpr (sizeof...(Args) == 0) 
            return vec1;
        else 
            return ms::merge(vec1, std::move(args)...);
    }

    inline std::string to_string(const BD_Type bd_type) {
        switch (bd_type)
        {
        case BD_Type::noBD: return "noBD";
        case BD_Type::no_typeI: return "no_typeI";
        case BD_Type::no_typeII: return "no_typeII";
        case BD_Type::standard: return "standard";
        case BD_Type::improved: return "improved";
        case BD_Type::typeI_1: return "typeI_1";
        case BD_Type::typeI_2: return "typeI_2";
        case BD_Type::typeI_3: return "typeI_3";
        case BD_Type::typeI_4: return "typeI_4";
        case BD_Type::typeI_1_without_typeII: return "typeI_1_without_typeII";
        case BD_Type::typeI_4_without_typeII: return "typeI_4_without_typeII";
        case BD_Type::typeI_1_cc: return "typeI_1_cc";
        case BD_Type::typeI_1_dc: return "typeI_1_dc";
        case BD_Type::noBD_dc: return "noBD_dc";
        default:
            throw std::runtime_error("wrong bd_type");
            return "";
        }
    }
}


//template definition part
template <ushort space_dimension_, ushort solution_order_>
Polynomial_Reconstruction<space_dimension_, solution_order_>::Polynomial_Reconstruction(const Grid<space_dimension_>& grid) {
    SET_TIME_POINT;

    this->basis_vector_functions_ = grid.cell_basis_vector_functions<solution_order_>();

    Log::content_ << std::left << std::setw(50) << "@ Polynomial reconstruction precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <ushort space_dimension_, ushort solution_order_>
auto Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_set_of_transposed_gradient_basis(void) const {
    const auto num_cell = this->basis_vector_functions_.size();
    std::vector<Matrix_Function<Polynomial<space_dimension_>, space_dimension_, this->num_basis_>> set_of_transposed_gradient_basis(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        auto& transposed_gradient_basis = set_of_transposed_gradient_basis[i];
        const auto& basis_function = this->basis_vector_functions_[i];

        for (ushort j = 0; j < this->num_basis_; ++j)
            transposed_gradient_basis.change_column(j, basis_function[j].gradient());
    }

    return set_of_transposed_gradient_basis;
}

template <ushort space_dimension_, ushort solution_order_>
auto Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_basis_node(const uint cell_index, const Space_Vector_& node) const {
    const auto& basis_functions = this->basis_vector_functions_[cell_index];
    return basis_functions(node);
}

template <ushort space_dimension_, ushort solution_order_>
Dynamic_Matrix Polynomial_Reconstruction<space_dimension_, solution_order_>::basis_nodes(const uint cell_index, const std::vector<Space_Vector_>& nodes) const {
    const auto& basis_functions = this->basis_vector_functions_[cell_index];

    const auto num_node = nodes.size();
    Dynamic_Matrix basis_nodes(this->num_basis_, num_node);

    for (ushort j = 0; j < num_node; ++j)
        basis_nodes.change_column(j, basis_functions(nodes[j]));

    return basis_nodes;
}

template <ushort space_dimension_, ushort solution_order_>
double Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_P0_basis_value(const uint cell_index) const {
    const auto& basis_function = this->basis_vector_functions_[cell_index];
    const auto& P0_basis_function = basis_function[0];

    return P0_basis_function.to_constant();
}

template <ushort space_dimension_, ushort solution_order_>
std::vector<double> Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_P0_basis_values(void) const {
    const auto num_cell = this->basis_vector_functions_.size();

    std::vector<double> P0_basis_values(num_cell);
    for (uint i = 0; i < num_cell; ++i)
        P0_basis_values[i] = this->calculate_P0_basis_value(i);

    return P0_basis_values;
} 

template <ushort space_dimension_, ushort solution_order_>
template<ushort num_equation>
auto Polynomial_Reconstruction<space_dimension_, solution_order_>::solution_function(const uint cell_index, const Matrix<num_equation, num_basis_>& coefficient) const {
    dynamic_require(cell_index < basis_vector_functions_.size(), "cell index should be less than number of cell");
    return coefficient * this->basis_vector_functions_[cell_index];
}

template <ushort space_dimension_, ushort solution_order_>
hMLP_Base<space_dimension_, solution_order_>::hMLP_Base(const Grid<space_dimension_>& grid)
    : Polynomial_Reconstruction<space_dimension_, solution_order_>(grid), vnode_index_to_share_cell_index_set_(grid.get_vnode_index_to_share_cell_index_set_consider_pbdry()) {
    SET_TIME_POINT;

    for (ushort i = 0; i <= solution_order_; ++i)
        this->num_Pn_projection_basis_[i] = ms::combination_with_repetition(1 + space_dimension_, i);

    this->set_of_vnode_indexes_ = grid.cell_set_of_vnode_indexes();
    this->volumes_ = grid.cell_volumes();
    this->P0_basis_values_ = this->calculate_P0_basis_values();

    Log::content_ << std::left << std::setw(50) << "@ hMLP Base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <ushort space_dimension_, ushort solution_order_>
auto hMLP_Base<space_dimension_, solution_order_>::Pn_projection_matrix(const ushort Pn) const {
    dynamic_require(Pn <= solution_order_, "Projection order should be less then solution order");

    const auto num_Pn_projection_basis = this->num_Pn_projection_basis_[Pn];

    std::array<double, This_::num_basis_> projection_value = { 0 };
    for (ushort j = 0; j < num_Pn_projection_basis; ++j)
        projection_value[j] = 1.0; //preserve Pn_projection_basis

    return Matrix(projection_value);
}

template <ushort space_dimension_, ushort solution_order_>
auto hMLP_Base<space_dimension_, solution_order_>::Pn_mode_matrix(const ushort Pn) const {
    dynamic_require(Pn <= solution_order_, "mode should be less then solution order");

    const auto num_Pnm1_projection_basis = this->num_Pn_projection_basis_[Pn - 1];
    const auto num_Pn_projection_basis = this->num_Pn_projection_basis_[Pn];

    std::array<double, This_::num_basis_> projection_value = { 0 };
    for (ushort j = num_Pnm1_projection_basis; j < num_Pn_projection_basis; ++j)
        projection_value[j] = 1.0; //preserve Pn_mode_basis

    return Matrix(projection_value);
}

template <ushort space_dimension_, ushort solution_order_>
auto hMLP_Base<space_dimension_, solution_order_>::limiting_matrix(const double limiting_value) const {
    std::array<double, This_::num_basis_> limiting_values;
    limiting_values.fill(limiting_value);
    limiting_values[0] = 1.0; // keep P0 coefficient
        
    return Matrix(limiting_values);
}

template <ushort space_dimension_, ushort solution_order_>
hMLP_Reconstruction<space_dimension_, solution_order_>::hMLP_Reconstruction(const Grid<space_dimension_>& grid) 
    : hMLP_Base<space_dimension_, solution_order_>(grid){
    SET_TIME_POINT;

    constexpr ushort P1 = 1;
    const auto set_of_vnodes = grid.cell_set_of_vnodes();

    const auto num_cell = this->volumes_.size();
    this->set_of_basis_vnodes_.reserve(num_cell);
    this->set_of_P1_projected_basis_vnodes_.reserve(num_cell);
    this->P0_basis_values_.reserve(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        auto basis_vnodes = this->basis_nodes(i, set_of_vnodes[i]);
        auto P1_projected_basis_vnodes = this->Pn_projection_matrix(P1) * basis_vnodes;
        const auto P0_basis_value = P1_projected_basis_vnodes.at(0, 0);

        this->set_of_basis_vnodes_.push_back(std::move(basis_vnodes));
        this->set_of_P1_projected_basis_vnodes_.push_back(std::move(P1_projected_basis_vnodes));
        this->P0_basis_values_.push_back(P0_basis_value);
    }

    Log::content_ << std::left << std::setw(50) << "@ hMLP precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <ushort space_dimension_, ushort solution_order_>
template <ushort num_equation>
void hMLP_Reconstruction<space_dimension_, solution_order_>::reconstruct(std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const {
    const auto P0_criterion_values = this->calculate_P0_criterion_values(solution_coefficients);
    const auto vnode_index_to_allowable_min_max_criterion_value = this->calculate_vertex_node_index_to_allowable_min_max_criterion_value(P0_criterion_values);

    const auto num_cell = solution_coefficients.size();

    for (uint i = 0; i < num_cell; ++i) {
        auto temporal_solution_order = solution_order_;              

        const auto& vnode_indexes = this->set_of_vnode_indexes_[i];
        const auto num_vnode = vnode_indexes.size();

        const auto P0_criterion_value = P0_criterion_values[i];
        const auto volume = this->volumes_[i];

        auto& solution_coefficient = solution_coefficients[i];

        while (true) {
            bool is_normal = true;

            auto solution_vnodes = solution_coefficient * this->set_of_basis_vnodes_[i];
            auto P1_projected_solution_vnodes = solution_coefficient * this->set_of_P1_projected_basis_vnodes_[i];
            
            for (ushort j = 0; j < num_vnode; ++j) {
                const auto vnode_index = vnode_indexes[j];
                const auto [allowable_min, allowable_max] = vnode_index_to_allowable_min_max_criterion_value.at(vnode_index);

                const auto criterion_value = solution_vnodes.at(This_::criterion_variable_index_, j);
                const auto P1_projected_criterion_value = P1_projected_solution_vnodes.at(This_::criterion_variable_index_, j);

                const auto higher_mode_criterion_value = criterion_value - P1_projected_criterion_value;
                const auto P1_mode_criterion_value = P1_projected_criterion_value - P0_criterion_value;

                if (!Constant_Region_Detector::is_constant(criterion_value, P0_criterion_value, volume) &&
                    !P1_Projected_MLP_Condition::is_satisfy(P1_projected_criterion_value, allowable_min, allowable_max) &&
                    !MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, higher_mode_criterion_value, P1_mode_criterion_value, allowable_min, allowable_max)) {
                    is_normal = false;
                    break;
                }
            }

            if (is_normal)
                break;

            if (temporal_solution_order <= 2) {
                temporal_solution_order = 1;
                solution_coefficient *= this->Pn_projection_matrix(temporal_solution_order);
                P1_projected_solution_vnodes = solution_coefficient * this->set_of_P1_projected_basis_vnodes_[i];

                double limiting_value = 1.0;
                for (ushort j = 0; j < num_vnode; ++j) {
                    const auto vnode_index = vnode_indexes[j];
                    const auto [allowable_min, allowable_max] = vnode_index_to_allowable_min_max_criterion_value.at(vnode_index);
    
                    const auto P1_projected_criterion_value = P1_projected_solution_vnodes.at(This_::criterion_variable_index_, j);
                    const auto P1_mode_criterion_value = P1_projected_criterion_value - P0_criterion_value;
    
                    limiting_value = (std::min)(limiting_value, MLP_u1_Limiting_Strategy::limiter_function(P1_mode_criterion_value, P0_criterion_value, allowable_min, allowable_max));
                }
    
                solution_coefficient *= this->limiting_matrix(limiting_value);
                break;
            }
            else
            solution_coefficient *= this->Pn_projection_matrix(--temporal_solution_order); //limiting highest mode       

        }
    }
}

template <ushort space_dimension_, ushort solution_order_>
template <ushort num_equation>
auto hMLP_Reconstruction<space_dimension_, solution_order_>::calculate_P0_criterion_values(const std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const {
    const auto num_cell = solution_coefficients.size();
    std::vector<double> P0_criterion_values(num_cell);

    for (uint i = 0; i < num_cell; ++i) 
        P0_criterion_values[i] = solution_coefficients[i].at(This_::criterion_variable_index_, 0) * this->P0_basis_values_[i];          

    return P0_criterion_values;
}

template <ushort space_dimension_, ushort solution_order_>
auto hMLP_Reconstruction<space_dimension_, solution_order_>::calculate_vertex_node_index_to_allowable_min_max_criterion_value(const std::vector<double>& P0_criterion_values) const {
    const auto num_vnode = this->vnode_index_to_share_cell_index_set_.size();

    std::unordered_map<uint, std::pair<double, double>> vnode_index_to_allowable_min_max_criterion_value;
    vnode_index_to_allowable_min_max_criterion_value.reserve(num_vnode);

    for (const auto& [vnode_index, share_cell_indexes] : this->vnode_index_to_share_cell_index_set_) {
        const auto num_share_cell = share_cell_indexes.size();
        std::vector<double> criterion_variables;
        criterion_variables.reserve(num_share_cell);

        for (const auto cell_index : share_cell_indexes)
            criterion_variables.push_back(P0_criterion_values[cell_index]);

        const auto min_criterion_value = *std::min_element(criterion_variables.begin(), criterion_variables.end());
        const auto max_criterion_value = *std::max_element(criterion_variables.begin(), criterion_variables.end());

        vnode_index_to_allowable_min_max_criterion_value.emplace(vnode_index, std::make_pair(min_criterion_value, max_criterion_value));
    }

    return vnode_index_to_allowable_min_max_criterion_value;
}


template <ushort space_dimension_, ushort solution_order_>
hMLP_BD_Reconstruction<space_dimension_, solution_order_>::hMLP_BD_Reconstruction(const Grid<space_dimension_>& grid) 
    : hMLP_Base<space_dimension_, solution_order_>(grid) {
    SET_TIME_POINT;

    this->pbdry_vnode_index_to_matched_vnode_index_set_ = grid.pbdry_vnode_index_to_matched_pbdry_vnode_index_set();
    this->set_of_face_share_cell_indexes_ = grid.set_of_face_share_cell_indexes_consider_pbdry();

    const auto set_of_vnodes = grid.cell_set_of_vnodes();
    const auto simplex_flags = grid.cell_simplex_flags();
    const auto set_of_sub_simplex_geometries = grid.cell_set_of_sub_simplex_geometries(simplex_flags);

    const auto num_cell = set_of_vnodes.size();
    this->set_of_basis_vnodes_.reserve(num_cell);
    this->set_of_simplex_P1_projected_basis_vnodes_.reserve(num_cell);
    this->set_of_simplex_P0_projected_basis_vnodes_.reserve(num_cell);

    constexpr ushort P0 = 0;
    constexpr ushort P1 = 1;
    for (uint i = 0; i < num_cell; ++i) {
        const auto& vnodes = set_of_vnodes[i];
        auto basis_vnodes = this->basis_nodes(i, vnodes);
                
        if (simplex_flags[i]) {
            auto P0_projected_basis_vnodes = this->Pn_projection_matrix(P0) * basis_vnodes;
            auto P1_projected_basis_vnodes = this->Pn_projection_matrix(P1) * basis_vnodes;

            this->set_of_simplex_P0_projected_basis_vnodes_.push_back(std::move(P0_projected_basis_vnodes));
            this->set_of_simplex_P1_projected_basis_vnodes_.push_back(std::move(P1_projected_basis_vnodes));            
        }
        else {
            const auto& sub_simplex_geometries = set_of_sub_simplex_geometries[i];

            const auto num_vnode = vnodes.size();
            Dynamic_Matrix simplex_P0_projected_basis_vnodes(This_::num_basis_, num_vnode);
            Dynamic_Matrix simplex_P1_projected_basis_vnodes(This_::num_basis_, num_vnode);

            for (ushort j = 0; j < num_vnode; ++j) {
                const auto simplex_P0_projection_basis_vector_function = this->calculate_simplex_Pn_projection_basis_vector_function<P0>(i, sub_simplex_geometries[j]);
                const auto simplex_P1_projection_basis_vector_function = this->calculate_simplex_Pn_projection_basis_vector_function<P1>(i, sub_simplex_geometries[j]);

                const auto& vnode = vnodes[j];

                simplex_P0_projected_basis_vnodes.change_column(j, simplex_P0_projection_basis_vector_function(vnode));
                simplex_P1_projected_basis_vnodes.change_column(j, simplex_P1_projection_basis_vector_function(vnode));                    
            }

            this->set_of_simplex_P0_projected_basis_vnodes_.push_back(std::move(simplex_P0_projected_basis_vnodes));
            this->set_of_simplex_P1_projected_basis_vnodes_.push_back(std::move(simplex_P1_projected_basis_vnodes));
        }

        this->set_of_basis_vnodes_.push_back(std::move(basis_vnodes));
    }

    ms::merge(this->face_volumes_, grid.inner_face_volumes(), grid.periodic_boundary_volumes());
    ms::merge(this->face_oc_nc_index_pairs_, grid.inner_face_oc_nc_index_pairs(), grid.periodic_boundary_oc_nc_index_pairs());

    const auto num_face_without_bdry = face_volumes_.size();
    this->face_characteristic_lengths_.reserve(num_face_without_bdry);
    this->face_oc_nc_side_basis_jump_qnodes_pairs_.reserve(num_face_without_bdry);
    this->set_of_jump_qweights_.reserve(num_face_without_bdry);

    const auto inner_face_quadrature_rules = grid.inner_face_quadrature_rules(solution_order_);
    const auto pbdry_quadrature_rule_pairs = grid.periodic_boundary_quadrature_rule_pairs(solution_order_);

    const auto num_inner_face = inner_face_quadrature_rules.size();
    const auto num_pbdry_pair = pbdry_quadrature_rule_pairs.size();

    for (uint i = 0; i < num_inner_face; ++i) {
        const auto charactersitic_length = std::pow(this->face_volumes_[i], 1.0 / (space_dimension_ - 1));
        this->face_characteristic_lengths_.push_back(charactersitic_length);

        const auto [oc_index, nc_index] = this->face_oc_nc_index_pairs_[i];
        const auto& quadrature_rule = inner_face_quadrature_rules[i];

        this->face_oc_nc_side_basis_jump_qnodes_pairs_.push_back({ this->basis_nodes(oc_index, quadrature_rule.points), this->basis_nodes(nc_index, quadrature_rule.points) });
        this->set_of_jump_qweights_.push_back(quadrature_rule.weights);
    }

    for (uint i = 0; i < num_pbdry_pair; ++i) {
        const auto characteristic_length = std::pow(this->face_volumes_[num_inner_face + i], 1.0 / (space_dimension_ - 1));
        this->face_characteristic_lengths_.push_back(characteristic_length);

        const auto [oc_index, nc_index] = this->face_oc_nc_index_pairs_[num_inner_face + i];
        const auto& [oc_side_quadrature_rule, nc_side_quadrature_rule] = pbdry_quadrature_rule_pairs[i];

        this->face_oc_nc_side_basis_jump_qnodes_pairs_.push_back({ this->basis_nodes(oc_index, oc_side_quadrature_rule.points), this->basis_nodes(nc_index, nc_side_quadrature_rule.points) });
        this->set_of_jump_qweights_.push_back(oc_side_quadrature_rule.weights); // oc side qweights == nc side qweights
    }

    Log::content_ << std::left << std::setw(50) << "@ hMLP_BD precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <ushort space_dimension_, ushort solution_order_>
template <ushort num_equation>
void hMLP_BD_Reconstruction<space_dimension_, solution_order_>::reconstruct(std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const {
    const auto set_of_vnode_index_to_simplex_P0_criterion_value = this->calculate_set_of_vertex_node_index_to_simplex_P0_criterion_value(solution_coefficients);
    const auto vnode_index_to_allowable_min_max_criterion_value = this->calculate_vertex_node_index_to_allowable_min_max_criterion_value(set_of_vnode_index_to_simplex_P0_criterion_value);
    const auto set_of_num_troubled_boundary = this->calculate_set_of_num_troubled_boundary(solution_coefficients);     

    //temporary
    const auto P0_solutions = this->calculate_P0_solutions(solution_coefficients);
    //const auto shock_flags = this->is_shock(P0_solutions);
    //const auto discontinuity_flags = this->is_discontinuity(P0_solutions);
    //temporary

    const auto num_cell = solution_coefficients.size();

    ////post
    //const auto before_limiting = solution_coefficients;

    //const auto P1_projection_matrix = this->Pn_projection_matrix(1);

    //auto P1_projection_coefficients = solution_coefficients;
    //for (auto& P1_projection_coefficient : P1_projection_coefficients)
    //    P1_projection_coefficient *= P1_projection_matrix;

    //std::vector<ushort> typeI_flags(num_cell);
    //std::vector<ushort> type_II_flag(num_cell);    
    ////post

    for (uint i = 0; i < num_cell; ++i) {
        auto temporal_solution_order = solution_order_;

        const auto& vnode_indexes = this->set_of_vnode_indexes_[i];
        const auto num_vnode = vnode_indexes.size();

        const auto& vnode_index_to_simplex_P0_criterion_value = set_of_vnode_index_to_simplex_P0_criterion_value[i];

        const auto volume = this->volumes_[i];
        const auto num_troubled_boundary = set_of_num_troubled_boundary[i];
        auto& solution_coefficient = solution_coefficients[i];

        while (true) {
            bool is_normal_cell = true;

            auto solution_vnodes = solution_coefficient * this->set_of_basis_vnodes_[i];
            auto simplex_P1_projected_solution_vnodes = solution_coefficient * this->set_of_simplex_P1_projected_basis_vnodes_[i];

            for (ushort j = 0; j < num_vnode; ++j) {
                const auto vnode_index = vnode_indexes[j];
                const auto [allowable_min, allowable_max] = vnode_index_to_allowable_min_max_criterion_value.at(vnode_index);

                const auto simplex_P0_criterion_value = vnode_index_to_simplex_P0_criterion_value.at(vnode_index);

                const auto criterion_value = solution_vnodes.at(This_::criterion_variable_index_, j);
                const auto simplex_P1_projected_criterion_value = simplex_P1_projected_solution_vnodes.at(This_::criterion_variable_index_, j);

                const auto simplex_higher_mode_criterion_value = criterion_value - simplex_P1_projected_criterion_value;
                const auto simplex_P1_mode_criterion_value = simplex_P1_projected_criterion_value - simplex_P0_criterion_value;
 
                //if constexpr (__hMLP_BD_TYPE__ == BD_Type::noBD) {
                //    if (!Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume) &&
                //        !P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max) &&
                //        !MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                if constexpr (__hMLP_BD_TYPE__ == BD_Type::standard) {
                    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                        continue;

                    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {

                        //if (this->is_typeI_subcell_oscillation(num_troubled_boundary) && shock_flags[i]) {
                        if (this->is_typeI_subcell_oscillation(num_troubled_boundary)) {
                            temporal_solution_order = 1;
                            is_normal_cell = false;
                            break;
                        }
                        continue;
                    }

                    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max) ||
                        this->is_typeII_subcell_oscillation(num_troubled_boundary)) {

                        is_normal_cell = false;
                        break;
                    }
                }
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::improved) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        //if (this->is_typeI_subcell_oscillation(num_troubled_boundary) && shock_flags[i]) {
                //        if (this->is_typeI_subcell_oscillation(num_troubled_boundary)) {
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max)) {
                //        if (this->is_typeII_subcell_oscillation(num_troubled_boundary) && discontinuity_flags[i]) {
                //            is_normal_cell = false;
                //            break;
                //        }
                //    }
                //    else {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::no_typeI) {
                //    ////post
                //    //if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max) && shock_flags[i])
                //    //    type_I_flag[i] = 1;
                //    ////post

                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max))
                //        continue;
                //    
                //    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max) ||
                //        this->is_typeII_subcell_oscillation(num_troubled_boundary)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::no_typeII) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        if (this->is_typeI_subcell_oscillation(num_troubled_boundary) && shock_flags[i]) { //temporary
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::typeI_1) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        if (num_troubled_boundary >= 1 && shock_flags[i]) { //temporary
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max) ||
                //        this->is_typeII_subcell_oscillation(num_troubled_boundary)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::typeI_2) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        if (num_troubled_boundary >= 2 && shock_flags[i]) { //temporary
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max) ||
                //        this->is_typeII_subcell_oscillation(num_troubled_boundary)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::typeI_3) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        if (num_troubled_boundary >= 3 && shock_flags[i]) { //temporary
                //            //typeI_flags[i] = 1; //post
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max) ||
                //        this->is_typeII_subcell_oscillation(num_troubled_boundary)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::typeI_4) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        if (num_troubled_boundary >= 4 && shock_flags[i]) { //temporary
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max) ||
                //        this->is_typeII_subcell_oscillation(num_troubled_boundary)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::typeI_1_without_typeII) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        //if (num_troubled_boundary >= 1) {
                //        if (num_troubled_boundary >= 1 && shock_flags[i]) {//temporary
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::typeI_4_without_typeII) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        //if (num_troubled_boundary >= 4 ) {
                //        if (num_troubled_boundary >= 4 && shock_flags[i]) {//temporary
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::typeI_1_cc) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        if (num_troubled_boundary >= 1 && shock_flags[i]) {//temporary
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max)) {
                //        if (discontinuity_flags[i]) {
                //            is_normal_cell = false;
                //            break;
                //        }
                //    }
                //    else {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::typeI_1_dc) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        if (num_troubled_boundary >= 1 && discontinuity_flags[i]) {//temporary
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (!MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max)) {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                //else if constexpr (__hMLP_BD_TYPE__ == BD_Type::noBD_dc) {
                //    if (Constant_Region_Detector::is_constant(criterion_value, simplex_P0_criterion_value, volume))
                //        continue;

                //    if (P1_Projected_MLP_Condition::is_satisfy(simplex_P1_projected_criterion_value, allowable_min, allowable_max)) {
                //        if (shock_flags[i]) {
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //        continue;
                //    }

                //    if (MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, simplex_higher_mode_criterion_value, simplex_P1_mode_criterion_value, allowable_min, allowable_max)) {
                //        if (discontinuity_flags[i]) {
                //            temporal_solution_order = 1;
                //            is_normal_cell = false;
                //            break;
                //        }
                //    }
                //    else {
                //        is_normal_cell = false;
                //        break;
                //    }
                //}
                else
                    throw std::runtime_error("wrong bd type");
            }

            if (is_normal_cell)
                break;

            //test
            std::cout << "activate limiter on " << i << "cell\n";
            std::exit(523);
            //test

            if (temporal_solution_order <= 2) {
                temporal_solution_order = 1;
                solution_coefficient *= this->Pn_projection_matrix(temporal_solution_order);
                simplex_P1_projected_solution_vnodes = solution_coefficient * this->set_of_simplex_P1_projected_basis_vnodes_[i];

                double limiting_value = 1.0;

                for (ushort j = 0; j < num_vnode; ++j) {
                    const auto vnode_index = vnode_indexes[j];
                    const auto [allowable_min, allowable_max] = vnode_index_to_allowable_min_max_criterion_value.at(vnode_index);

                    const auto simplex_P0_criterion_value = vnode_index_to_simplex_P0_criterion_value.at(vnode_index);

                    const auto simplex_P1_projected_criterion_value = simplex_P1_projected_solution_vnodes.at(This_::criterion_variable_index_, j);
                    const auto simplex_P1_mode_criterion_value = simplex_P1_projected_criterion_value - simplex_P0_criterion_value;

                    limiting_value = (std::min)(limiting_value, MLP_u1_Limiting_Strategy::limiter_function(simplex_P1_mode_criterion_value, simplex_P0_criterion_value, allowable_min, allowable_max));
               }

                solution_coefficient *= this->limiting_matrix(limiting_value);
                break;
            }
            else
                solution_coefficient *= this->Pn_projection_matrix(--temporal_solution_order); //limiting highest mode
        }
    }

    //debug
    //if (Debugger::conditions_[0])
    //    std::exit(523);
    //debug

    //post
    
    //auto limited_P1_projection_coefficients = solution_coefficients;
    //for (auto& P1_projection_coefficient : limited_P1_projection_coefficients)
    //    P1_projection_coefficient *= P1_projection_matrix;

    if (Tecplot::post_condition_ == true) {
        Tecplot::record_cell_indexes();
        //Tecplot::record_cell_variables("num_troubled_boundary", set_of_num_troubled_boundary);
        //Tecplot::record_cell_variables("TypeI_flag", typeI_flags);
        //Tecplot::record_cell_variables("TypeII_flag", type_II_flag);
        //Tecplot::post_solution(before_limiting, "before");
        //Tecplot::post_solution(P1_projection_coefficients, "P1_projection");
        Tecplot::post_solution(solution_coefficients);
        //Tecplot::post_solution(limited_P1_projection_coefficients, "Limited_P1_projection");
    }

    //post
}

template <ushort space_dimension_, ushort solution_order_>
template <ushort num_equation>
auto hMLP_BD_Reconstruction<space_dimension_, solution_order_>::calculate_set_of_vertex_node_index_to_simplex_P0_criterion_value(const std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const {
    const auto num_cell = solution_coefficients.size();
    const auto num_vnode = this->vnode_index_to_share_cell_index_set_.size();

    std::vector<std::map<uint, double>> vnode_index_to_simplex_P0_criterion_values;
    vnode_index_to_simplex_P0_criterion_values.reserve(num_cell);

    for (uint i = 0; i < num_cell; ++i) {   
        std::map<uint, double> vnode_index_to_simplex_P0_criterion_value;

        const auto simplex_P0_solution_vnodes = solution_coefficients[i] * this->set_of_simplex_P0_projected_basis_vnodes_[i];
        const auto simplex_P0_criterion_value_vnodes = simplex_P0_solution_vnodes.row(This_::criterion_variable_index_);

        const auto& vnode_indexes = this->set_of_vnode_indexes_[i];
        const auto num_vnode = vnode_indexes.size();

        for (ushort j = 0; j < num_vnode; ++j)
            vnode_index_to_simplex_P0_criterion_value.emplace(vnode_indexes[j], simplex_P0_criterion_value_vnodes[j]);

        vnode_index_to_simplex_P0_criterion_values.push_back(std::move(vnode_index_to_simplex_P0_criterion_value));
    }

    return vnode_index_to_simplex_P0_criterion_values;
}

template <ushort space_dimension_, ushort solution_order_>
auto hMLP_BD_Reconstruction<space_dimension_, solution_order_>::calculate_vertex_node_index_to_allowable_min_max_criterion_value(const std::vector<std::map<uint, double>>& set_of_vnode_index_to_simplex_P0_criterion_value) const {
    const auto num_vnode = this->vnode_index_to_share_cell_index_set_.size();

    //gathering P0 criterion values at vertex
    std::unordered_map<uint, std::vector<double>> vnode_index_to_simplex_P0_criterion_values;
    vnode_index_to_simplex_P0_criterion_values.reserve(num_vnode);

    for (const auto& vnode_index_to_simplex_P0_criterion_value : set_of_vnode_index_to_simplex_P0_criterion_value) {
        for (const auto& [vnode_index, simplex_P0_criterion_value] : vnode_index_to_simplex_P0_criterion_value) {
            if (!vnode_index_to_simplex_P0_criterion_values.contains(vnode_index))
                vnode_index_to_simplex_P0_criterion_values.emplace(vnode_index, std::vector<double>());

            vnode_index_to_simplex_P0_criterion_values.at(vnode_index).push_back(simplex_P0_criterion_value);
        }
    }

    //consider pbdry
    std::unordered_set<uint> considered_pbdry_node_index_set;

    for (const auto& [vnode_index, matched_vnode_index_set] : this->pbdry_vnode_index_to_matched_vnode_index_set_) {
        if (considered_pbdry_node_index_set.contains(vnode_index))
            continue;
        
        //gathering
        std::vector<double> gathered_values = vnode_index_to_simplex_P0_criterion_values.at(vnode_index);
        for (const auto matched_vnode_index : matched_vnode_index_set) {
            const auto& values = vnode_index_to_simplex_P0_criterion_values.at(matched_vnode_index);
            gathered_values.insert(gathered_values.end(), values.begin(), values.end());             
        }

        //copy
        for (const auto matched_vnode_index : matched_vnode_index_set) 
            vnode_index_to_simplex_P0_criterion_values.at(matched_vnode_index) = gathered_values;

        vnode_index_to_simplex_P0_criterion_values.at(vnode_index) = std::move(gathered_values);

        //prevent repetition
        considered_pbdry_node_index_set.insert(matched_vnode_index_set.begin(), matched_vnode_index_set.end());
        considered_pbdry_node_index_set.insert(vnode_index);
    }

    std::unordered_map<uint, std::pair<double, double>> vnode_index_to_allowable_min_max_criterion_value;
    vnode_index_to_allowable_min_max_criterion_value.reserve(num_vnode);

    for (const auto& [vnode_index, simplex_P0_criterion_values] : vnode_index_to_simplex_P0_criterion_values) {
        const auto min_criterion_value = *std::min_element(simplex_P0_criterion_values.begin(), simplex_P0_criterion_values.end());
        const auto max_criterion_value = *std::max_element(simplex_P0_criterion_values.begin(), simplex_P0_criterion_values.end());

        vnode_index_to_allowable_min_max_criterion_value.emplace(vnode_index, std::make_pair(min_criterion_value, max_criterion_value));
    }

    return vnode_index_to_allowable_min_max_criterion_value;
}

template <ushort space_dimension_, ushort solution_order_>
template <ushort num_equation>
std::vector<ushort> hMLP_BD_Reconstruction<space_dimension_, solution_order_>::calculate_set_of_num_troubled_boundary(const std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const {
    const auto num_cell = solution_coefficients.size();

    std::vector<ushort> set_of_num_trouble_boundary(num_cell);

    const auto num_face = this->face_characteristic_lengths_.size();
    for (uint i = 0; i < num_face; ++i) {
        //calculate jump criterion
        const auto characteristic_length = this->face_characteristic_lengths_[i];
        const auto threshold_value = std::pow(characteristic_length, 0.5 * (solution_order_ + 1));

        //calculate jump
        const auto [oc_index, nc_index] = this->face_oc_nc_index_pairs_[i];
        const auto& [oc_side_basis_jump_qnodes, nc_side_basis_jump_qnodes] = this->face_oc_nc_side_basis_jump_qnodes_pairs_[i];

        const auto oc_side_solution_jump_qnodes = solution_coefficients[oc_index] * oc_side_basis_jump_qnodes;
        const auto nc_side_solution_jump_qnodes = solution_coefficients[nc_index] * nc_side_basis_jump_qnodes;

        const auto oc_side_criterion_solution_jump_qnodes = oc_side_solution_jump_qnodes.row(This_::criterion_variable_index_);
        const auto nc_side_criterion_solution_jump_qnodes = nc_side_solution_jump_qnodes.row(This_::criterion_variable_index_);
        auto criterion_solution_diff_jump_qnodes = oc_side_criterion_solution_jump_qnodes - nc_side_criterion_solution_jump_qnodes;

        criterion_solution_diff_jump_qnodes.be_absolute();

        const auto jump = criterion_solution_diff_jump_qnodes.inner_product(this->set_of_jump_qweights_[i]);
        const auto smooth_boundary_indicator = jump / this->face_volumes_[i];

        //compare and reflect
        if (smooth_boundary_indicator > threshold_value) {
            set_of_num_trouble_boundary[oc_index]++;
            set_of_num_trouble_boundary[nc_index]++;
        }
    }

    return set_of_num_trouble_boundary;
}

template <ushort space_dimension_, ushort solution_order_>
template <ushort projection_order>
auto hMLP_BD_Reconstruction<space_dimension_, solution_order_>::calculate_simplex_Pn_projection_basis_vector_function(const uint cell_index, const Geometry<space_dimension_>& sub_simplex_geometry) const {        
    constexpr auto Pk = projection_order;
    const auto& basis_vector_function = this->basis_vector_functions_[cell_index];
    const auto Pk_simplex_basis_vector_function = sub_simplex_geometry.orthonormal_basis_vector_function<Pk>();
        
    std::array<Polynomial<space_dimension_>, This_::num_basis_> simplex_Pn_projection_basis_functions = { 0 };

    for (ushort i = 0; i < This_::num_basis_; ++i) 
        for (ushort j = 0; j < this->num_Pn_projection_basis_[Pk]; ++j) 
            simplex_Pn_projection_basis_functions[i] += ms::inner_product(basis_vector_function[i], Pk_simplex_basis_vector_function[j], sub_simplex_geometry) * Pk_simplex_basis_vector_function[j];

    Vector_Function<Polynomial<space_dimension_>, This_::num_basis_> simplex_Pn_projection_basis_vector_function = simplex_Pn_projection_basis_functions;
    return simplex_Pn_projection_basis_vector_function;
}


template <ushort space_dimension_, ushort solution_order_>
bool hMLP_BD_Reconstruction<space_dimension_, solution_order_>::is_typeI_subcell_oscillation(const ushort num_troubled_boundary) const {
    constexpr ushort typeI_threshold_value = 2;
    return typeI_threshold_value <= num_troubled_boundary;
}

template <ushort space_dimension_, ushort solution_order_>
bool hMLP_BD_Reconstruction<space_dimension_, solution_order_>::is_typeII_subcell_oscillation(const ushort num_troubled_boundary) const {
    constexpr ushort typeI_threshold_value = 1;
    return typeI_threshold_value <= num_troubled_boundary;
}

template <ushort space_dimension_, ushort solution_order_>
template <ushort num_equation>
std::vector<Euclidean_Vector<num_equation>> hMLP_BD_Reconstruction<space_dimension_, solution_order_>::calculate_P0_solutions(const std::vector<Matrix<num_equation, This_::num_basis_>>& solution_coefficients) const {
    const auto num_cell = solution_coefficients.size();    
    std::vector<Euclidean_Vector<num_equation>> P0_solutions(num_cell);

    for (uint i = 0; i < num_cell; ++i)
        P0_solutions[i] = solution_coefficients[i].column(0) * this->P0_basis_values_[i];

    return P0_solutions;
}


template <ushort space_dimension_, ushort solution_order_>
template <ushort num_equation>
std::vector<bool> hMLP_BD_Reconstruction<space_dimension_, solution_order_>::is_shock(const std::vector<Euclidean_Vector<num_equation>>& P0_solutions) const {    
    const auto pressures = Euler<space_dimension_>::pressures(P0_solutions); //temporary code

    const auto num_cell = P0_solutions.size();
    std::vector<bool> is_shock(num_cell, false);

    for (uint i = 0; i < num_cell; ++i) {
        const auto& face_share_cell_indexes = this->set_of_face_share_cell_indexes_[i];

        const auto my_pressure = pressures[i];
        for (const auto face_share_cell_index : face_share_cell_indexes) {
            const auto other_pressure = pressures[face_share_cell_index];

            if (std::abs(my_pressure - other_pressure) >= 0.1 * (std::min)(my_pressure,other_pressure)) {
                is_shock[i] = true;
                break;
            }
        }
    }

    return is_shock;
}


template <ushort space_dimension_, ushort solution_order_>
template <ushort num_equation>
std::vector<bool> hMLP_BD_Reconstruction<space_dimension_, solution_order_>::is_discontinuity(const std::vector<Euclidean_Vector<num_equation>>& P0_solutions) const {

    const auto num_cell = P0_solutions.size();
    std::vector<bool> discontinuity_flags(num_cell, false);

    for (uint i = 0; i < num_cell; ++i) {
        const auto& face_share_cell_indexes = this->set_of_face_share_cell_indexes_[i];

        const auto my_density = P0_solutions[i][0];
        for (const auto face_share_cell_index : face_share_cell_indexes) {
            const auto other_density = P0_solutions[face_share_cell_index][0];

            if (std::abs(my_density - other_density) >= 0.1 * (std::min)(my_density, other_density)) {
                discontinuity_flags[i] = true;
                break;
            }
        }
    }

    return discontinuity_flags;
}