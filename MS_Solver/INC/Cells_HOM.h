#pragma once
#include "Governing_Equation.h"
#include "Grid_Builder.h"
#include "Reconstruction_Method_HOM.h"
#include "Solution_Scaling_Method.h"

//HOM이면 공통으로 사용하는 variable & method
template <typename Governing_Equation, typename Reconstruction_Method>
class Cells_HOM
{
private:
    static constexpr ushort space_dimension_    = Governing_Equation::space_dimension();
    static constexpr ushort num_equation_       = Governing_Equation::num_equation();
    static constexpr ushort num_basis_          = Reconstruction_Method::num_basis();

    using This_                 = Cells_HOM<Governing_Equation, Reconstruction_Method>;
    using Space_Vector_         = Euclidean_Vector<space_dimension_>;
    using Solution_Coefficient_ = Matrix<num_equation_, num_basis_>;
    using Residual_             = Matrix<num_equation_, num_basis_>;

public:
    using Discretized_Solution_ = Matrix<num_equation_, num_basis_>;

protected:    
    const Reconstruction_Method& reconstruction_method_;

    std::vector<double> volumes_;
    std::vector<std::array<double, space_dimension_>> projected_volumes_;
    std::vector<const Quadrature_Rule<space_dimension_>*> quadrature_rule_ptrs_;
    std::vector<Dynamic_Matrix> set_of_basis_qnodes_;
    std::vector<Dynamic_Matrix> gradient_basis_weights_;
    std::vector<double> P0_basis_values_;

public:
    Cells_HOM(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method);

public:
    double calculate_time_step(const std::vector<Solution_Coefficient_>& solution_coefficients, const double cfl) const;
    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) const;

    template <typename Initial_Condition>
    auto calculate_initial_solutions(void) const;

    template <typename Initial_Condition>
    void estimate_error(const std::vector<Solution_Coefficient_>& solution_coefficients, const double time) const;

public:
    void initialize_scaling_method(void) const;
};


//template definition
template <typename Governing_Equation, typename Reconstruction_Method>
Cells_HOM<Governing_Equation, Reconstruction_Method>::Cells_HOM(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
    : reconstruction_method_(reconstruction_method) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    this->volumes_.reserve(num_cell);
    this->projected_volumes_.reserve(num_cell);
    this->quadrature_rule_ptrs_.reserve(num_cell);
    this->set_of_basis_qnodes_.reserve(num_cell); // 2 * solution_order
    this->gradient_basis_weights_.reserve(num_cell);
    this->P0_basis_values_.reserve(num_cell);

    const auto set_of_transposed_gradient_basis = this->reconstruction_method_.calculate_set_of_transposed_gradient_basis();
    constexpr auto solution_order = Reconstruction_Method::solution_order();
    constexpr auto integrand_order = 2 * solution_order;

    for (uint i = 0; i < num_cell; ++i) {
        const auto& geometry = cell_elements[i].geometry_;

        this->volumes_.push_back(geometry.volume());
        this->projected_volumes_.push_back(geometry.projected_volume());
        this->P0_basis_values_.push_back(this->reconstruction_method_.calculate_P0_basis_value(i, geometry.center_node()));

        const auto& quadrature_rule = geometry.get_quadrature_rule(integrand_order);
        this->quadrature_rule_ptrs_.push_back(&quadrature_rule);
        this->set_of_basis_qnodes_.push_back(this->reconstruction_method_.calculate_basis_nodes(i, quadrature_rule.points));
        
        const auto& transposed_gradient_basis = set_of_transposed_gradient_basis[i];        
        const auto num_quadrature_point = quadrature_rule.points.size();
        
        Dynamic_Matrix gradient_basis_weight(num_quadrature_point * this->space_dimension_, this->num_basis_);

        for (ushort q = 0; q < num_quadrature_point; ++q) {
            const auto& quadrature_point = quadrature_rule.points[q];
            const auto quadrature_weight = quadrature_rule.weights[q];
            
            const auto part_of_gradient_basis_weight = transposed_gradient_basis(quadrature_point) * quadrature_weight;

            gradient_basis_weight.change_rows(q * this->space_dimension_, part_of_gradient_basis_weight);
        }

        this->gradient_basis_weights_.push_back(std::move(gradient_basis_weight));
    }

    Log::content_ << std::left << std::setw(50) << "@ Cells HOM precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
};


template <typename Governing_Equation, typename Reconstruction_Method>
double Cells_HOM<Governing_Equation, Reconstruction_Method>::calculate_time_step(const std::vector<Solution_Coefficient_>& solution_coefficients, const double cfl) const {
    const auto num_cell = solution_coefficients.size();

    std::vector<Euclidean_Vector<This_::num_equation_>> P0_solutions(num_cell);
    for (size_t i = 0; i < num_cell; ++i) {
        const auto P0_coefficient = solution_coefficients[i].column(0);
        P0_solutions[i] = P0_coefficient * this->P0_basis_values_[i];
    }

    const auto coordinate_projected_maximum_lambdas = Governing_Equation::calculate_coordinate_projected_maximum_lambdas(P0_solutions);
    
    std::vector<double> local_time_step(num_cell);
    for (size_t i = 0; i < num_cell; ++i) {
        const auto [y_projected_volume, x_projected_volume] = this->projected_volumes_[i];
        const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda] = coordinate_projected_maximum_lambdas[i];

        const auto x_radii = y_projected_volume * x_projeced_maximum_lambda;
        const auto y_radii = x_projected_volume * y_projeced_maximum_lambda;

        dynamic_require(std::isfinite(x_radii) && std::isfinite(y_radii), "time step should be finite number");
        local_time_step[i] = cfl * this->volumes_[i] / (x_radii + y_radii);        
    }

    static const auto solution_order = Reconstruction_Method::solution_order();
    static const auto c = static_cast<double>(1.0 / (2.0 * solution_order + 1.0));
    return *std::min_element(local_time_step.begin(), local_time_step.end()) * c;
}



 
template <typename Governing_Equation, typename Reconstruction_Method>
void Cells_HOM<Governing_Equation, Reconstruction_Method>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) const {
    const auto num_solution = solution_coefficients.size();
        
    for (uint i = 0; i < num_solution; ++i) {
        const auto& solution_coefficient = solution_coefficients[i];
        const auto& basis_quadrature_nodes = This_::set_of_basis_qnodes_[i];
        const auto solution_quadrature_nodes = solution_coefficient * basis_quadrature_nodes;

        const auto [num_eq, num_quadrature_node] = solution_quadrature_nodes.size();
        Dynamic_Matrix flux_quadrature_points(num_eq, This_::space_dimension_ * num_quadrature_node);

        for (size_t j = 0; j < num_quadrature_node; ++j) {
            const auto physical_flux = Governing_Equation::physical_flux(solution_quadrature_nodes.column<This_::num_equation_>(j));
            flux_quadrature_points.change_columns(j * This_::space_dimension_, physical_flux);
        }

        Residual_ delta_rhs;
        ms::gemm(flux_quadrature_points, This_::gradient_basis_weights_[i], delta_rhs);

        RHS[i] += delta_rhs;

        //if ((Debugger::conditions_[0] || Debugger::conditions_[1]) && (i == 1200 || i == 1202)) {
        //    std::cout << "cell_index " << i << "\n";
        //    std::cout << "Cell_RHS\n" << delta_rhs << "\n\n"; //debug
        //}
    }

    //if (Debugger::conditions_[0])//debug
    //    Debugger::conditions_[0] = false;

    //if (Debugger::conditions_[1])//debug
    //    exit(523);
}


template <typename Governing_Equation, typename Reconstruction_Method>
template <typename Initial_Condition>
auto Cells_HOM<Governing_Equation, Reconstruction_Method>::calculate_initial_solutions(void) const {
    const auto num_cell = this->quadrature_rule_ptrs_.size();
    std::vector<Matrix<This_::num_equation_, This_::num_basis_>> initial_solution_coefficients(num_cell);
    
    for (uint i = 0; i < num_cell; ++i) {
        const auto& qnodes = this->quadrature_rule_ptrs_[i]->points;
        const auto& qweights = this->quadrature_rule_ptrs_[i]->weights;

        const auto num_qnode = qnodes.size();

        Dynamic_Matrix initial_solution_qnodes(This_::num_equation_, num_qnode);
        Dynamic_Matrix basis_weight(num_qnode, This_::num_basis_);

        for (ushort q = 0; q < num_qnode; ++q) {
            initial_solution_qnodes.change_column(q, Initial_Condition::calculate_solution(qnodes[q]));
            basis_weight.change_row(q, this->reconstruction_method_.calculate_basis_node(i, qnodes[q]) * qweights[q]);
        }

        ms::gemm(initial_solution_qnodes, basis_weight, initial_solution_coefficients[i]);
    }

    return initial_solution_coefficients;
}

template <typename Governing_Equation, typename Reconstruction_Method>
template <typename Initial_Condition>
void Cells_HOM<Governing_Equation, Reconstruction_Method>::estimate_error(const std::vector<Solution_Coefficient_>& solution_coefficients, const double time) const {
    Log::content_ << "================================================================================\n";
    Log::content_ << "\t\t\t\t Error Anlysis\n";
    Log::content_ << "================================================================================\n";
    
    if constexpr (ms::is_Linear_Advection_2D<Governing_Equation>) {
        ////Version 1
        //const auto num_cell = solution_coefficients.size();

        //double global_L1_error = 0.0;
        //double global_L2_error = 0.0;
        //double global_Linf_error = 0.0;

        //for (size_t i = 0; i < num_cell; ++i) {
        //    const auto& qnodes = this->quadrature_rule_ptrs_[i]->points;
        //    const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Governing_Equation>(qnodes, time);
        //    const auto computed_solutions = solution_coefficients[i] * this->set_of_basis_qnodes_[i];

        //    const auto& qweights = this->quadrature_rule_ptrs_[i]->weights;

        //    const auto num_qnode = qnodes.size();

        //    double local_error = 0.0;
        //    double volume = 0.0;
        //    for (size_t q = 0; q < num_qnode; ++q) {
        //        local_error += (exact_solutions[q] - computed_solutions.column<This_::num_equation_>(q)).L1_norm() * qweights[q];
        //        volume += qweights[q];
        //    }
        //    local_error = local_error / volume;

        //    global_L1_error += local_error;
        //    global_L2_error += local_error * local_error;
        //    global_Linf_error = max(global_Linf_error, local_error);

        //}

        //global_L1_error = global_L1_error / num_cell;
        //global_L2_error = global_L2_error / num_cell;

        //global_L2_error = std::sqrt(global_L2_error);

        //Log::content_ << "L1 error \t\tL2 error \t\tLinf error \n";
        //Log::content_ << ms::double_to_string(global_L1_error) << "\t" << ms::double_to_string(global_L2_error) << "\t" << ms::double_to_string(global_Linf_error) << "\n\n";


        //Version 2
        const auto num_cell = solution_coefficients.size();

        double global_L1_error = 0.0;
        double global_L2_error = 0.0;
        double global_Linf_error = 0.0;

        for (size_t i = 0; i < num_cell; ++i) {
            const auto& qnodes = this->quadrature_rule_ptrs_[i]->points;
            const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Governing_Equation>(qnodes, time);
            const auto computed_solutions = solution_coefficients[i] * this->set_of_basis_qnodes_[i];

            const auto& qweights = this->quadrature_rule_ptrs_[i]->weights;

            const auto num_qnode = qnodes.size();

            double local_L1_error = 0.0;
            double local_L2_error = 0.0;
            double local_Linf_error = 0.0;

            double volume = 0.0;
            for (size_t q = 0; q < num_qnode; ++q) {
                const auto local_diff = (exact_solutions[q] - computed_solutions.column<This_::num_equation_>(q)).L1_norm();
                local_L1_error += local_diff * qweights[q];
                local_L2_error += local_diff * local_diff * qweights[q];
                local_Linf_error = max(local_Linf_error, local_diff);

                volume += qweights[q];
            }

            local_L1_error = local_L1_error / volume;
            local_L2_error = std::sqrt(local_L2_error / volume);


            global_L1_error += local_L1_error;
            global_L2_error += local_L2_error;
            global_Linf_error += local_Linf_error;
        }

        global_L1_error = global_L1_error / num_cell;
        global_L2_error = global_L2_error / num_cell;
        global_Linf_error = global_Linf_error / num_cell;

        Log::content_ << "L1 error \t\tL2 error \t\tLinf error \n";
        Log::content_ << ms::double_to_string(global_L1_error) << "\t" << ms::double_to_string(global_L2_error) << "\t" << ms::double_to_string(global_Linf_error);
    }
    else
        Log::content_ << Governing_Equation::name() << " does not provide error analysis result.\n\n";

    Log::print();
}

template <typename Governing_Equation, typename Reconstruction_Method>
void Cells_HOM<Governing_Equation, Reconstruction_Method>::initialize_scaling_method(void) const {
    Solution_Scaler::record_cell_basis_qnodes(this->set_of_basis_qnodes_);
}