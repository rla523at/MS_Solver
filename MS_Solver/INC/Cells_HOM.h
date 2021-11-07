#pragma once
#include "Discretized_Solution.h"

class Cells
{
public:
	virtual void calculate_RHS(double* RHS, const Discretized_Solution& discretized_solution) abstract;
};


class Cells_HOM : public Cells
{
public:
	void calculate_RHS(double* RHS, const Discretized_Solution& discretized_solution) override
	{
        for (uint i = 0; i < this->num_cells_; ++i) {
            const auto solution_at_quadrature_nodes = discretized_solution.calculate_solution_at_quadrature_nodes(i);
            
            Matrix flux_quadrature_points(num_eq, This_::space_dimension_ * num_quadrature_node);

            for (size_t j = 0; j < num_quadrature_node; ++j) {
                const auto physical_flux = this->governing_equation_->calculate_physical_flux(solution_at_quadrature_nodes[j]);
                flux_quadrature_points.change_columns(j * This_::space_dimension_, physical_flux);
            }

            ms::gemm(flux_quadrature_points, This_::qweights_gradient_basis_[i], delta_rhs);

            RHS[i] += delta_rhs;
        }
    }

private:
    size_t num_cells_;
    std::unique_ptr<Governing_Equation> governing_equation_;
};




//#include "Governing_Equation.h"
//#include "Grid.h"
//#include "Reconstruction_Method_HOM.h"
//#include "Solution_Scaler.h"
//
////HOM이면 공통으로 사용하는 variable & method
//template <typename Governing_Equation, typename Reconstruction_Method>
//class Cells_HOM
//{
//private:
//    static constexpr ushort space_dimension_    = Governing_Equation::space_dimension();
//    static constexpr ushort num_equation_       = Governing_Equation::num_equation();
//    static constexpr ushort num_basis_          = Reconstruction_Method::num_basis();
//
//    using This_                 = Cells_HOM<Governing_Equation, Reconstruction_Method>;
//    using Space_Vector_         = Euclidean_Vector<space_dimension_>;
//    using Solution_Coefficient_ = Static_Matrix<num_equation_, num_basis_>;
//    using Residual_             = Static_Matrix<num_equation_, num_basis_>;
//
//public:
//    using Discretized_Solution_ = Static_Matrix<num_equation_, num_basis_>;
//
//protected:    
//    const Reconstruction_Method& reconstruction_method_;
//
//    std::vector<double> volumes_;
//    std::vector<std::array<double, space_dimension_>> projected_volumes_;
//    std::vector<Quadrature_Rule<space_dimension_>> quadrature_rules_;
//    std::vector<Matrix> set_of_basis_qnodes_;
//    std::vector<Matrix> qweights_gradient_basis_;
//    std::vector<double> P0_basis_values_;
//
//public:
//    Cells_HOM(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method);
//
//public:
//    double calculate_time_step(const std::vector<Solution_Coefficient_>& solution_coefficients, const double cfl) const;
//    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) const;
//    void estimate_error(const std::vector<Solution_Coefficient_>& solution_coefficients, const double time) const;
//
//    template <typename Initial_Condition>
//    auto calculate_initial_solutions(void) const;
//
//public:
//    void initialize_scaling_method(void) const;
//};
//
//
////template definition
//template <typename Governing_Equation, typename Reconstruction_Method>
//Cells_HOM<Governing_Equation, Reconstruction_Method>::Cells_HOM(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
//    : reconstruction_method_(reconstruction_method) {
//    SET_TIME_POINT;
//
//    const auto set_of_transposed_gradient_basis = this->reconstruction_method_.calculate_set_of_transposed_gradient_basis();
//    constexpr auto solution_degree = Reconstruction_Method::solution_order();
//    constexpr auto integrand_degree = 2 * solution_degree;
//
//    this->volumes_ = grid.cell_volumes();
//    this->projected_volumes_ = grid.cell_projected_volumes();
//    this->quadrature_rules_ = grid.cell_quadrature_rules(integrand_degree);
//
//    const auto num_cell = this->volumes_.size();
//    this->P0_basis_values_.reserve(num_cell);
//    this->set_of_basis_qnodes_.reserve(num_cell);
//    this->qweights_gradient_basis_.reserve(num_cell);
//
//    for (uint i = 0; i < num_cell; ++i) {
//        const auto& transposed_gradient_basis = set_of_transposed_gradient_basis[i];
//        const auto& qnodes = this->quadrature_rules_[i].points;
//        const auto& qweights = this->quadrature_rules_[i].weights;
//        const auto num_qnode = qnodes.size();
//
//        Matrix qweight_gradient_basis(num_qnode * this->space_dimension_, this->num_basis_);
//
//        for (ushort q = 0; q < num_qnode; ++q) {
//            const auto part_of_gradient_basis_weight = transposed_gradient_basis(qnodes[q]) * qweights[q];
//            qweight_gradient_basis.change_rows(q * this->space_dimension_, part_of_gradient_basis_weight);
//        }
//
//        this->P0_basis_values_.push_back(this->reconstruction_method_.calculate_P0_basis_value(i));
//        this->set_of_basis_qnodes_.push_back(this->reconstruction_method_.basis_nodes(i, qnodes));
//        this->qweights_gradient_basis_.push_back(std::move(qweight_gradient_basis));
//    }
//
//    Log::content_ << std::left << std::setw(50) << "@ Cells HOM precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
//    Log::print();
//};
//
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//double Cells_HOM<Governing_Equation, Reconstruction_Method>::calculate_time_step(const std::vector<Solution_Coefficient_>& solution_coefficients, const double cfl) const {
//    const auto num_cell = solution_coefficients.size();
//
//    std::vector<Euclidean_Vector<This_::num_equation_>> P0_solutions(num_cell);
//    for (size_t i = 0; i < num_cell; ++i) {
//        const auto P0_coefficient = solution_coefficients[i].column(0);
//        P0_solutions[i] = P0_coefficient * this->P0_basis_values_[i];
//    }
//
//    const auto coordinate_projected_maximum_lambdas = Governing_Equation::calculate_coordinate_projected_maximum_lambdas(P0_solutions);
//    
//    std::vector<double> local_time_step(num_cell);
//    for (size_t i = 0; i < num_cell; ++i) {
//        if constexpr (This_::space_dimension_ == 2) {
//            const auto [y_projected_volume, x_projected_volume] = this->projected_volumes_[i];
//            const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda] = coordinate_projected_maximum_lambdas[i];
//
//            const auto x_radii = y_projected_volume * x_projeced_maximum_lambda;
//            const auto y_radii = x_projected_volume * y_projeced_maximum_lambda;
//
//            local_time_step[i] = cfl * this->volumes_[i] / (x_radii + y_radii);
//        }
//        else if (This_::space_dimension_ == 3) {
//            const auto [yz_projected_volume, xz_projected_volume, xy_projected_volume] = this->projected_volumes_[i];
//            const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda, z_projeced_maximum_lambda] = coordinate_projected_maximum_lambdas[i];
//
//            const auto x_radii = yz_projected_volume * x_projeced_maximum_lambda;
//            const auto y_radii = xz_projected_volume * y_projeced_maximum_lambda;
//            const auto z_radii = xy_projected_volume * z_projeced_maximum_lambda;
//
//            local_time_step[i] = cfl * this->volumes_[i] / (x_radii + y_radii + z_radii);
//        }
//        else
//            throw std::runtime_error("not supproted order");
//    }
//
//    constexpr auto solution_order = Reconstruction_Method::solution_order();
//    constexpr auto c = static_cast<double>(1.0 / (2.0 * solution_order + 1.0));
//    return *std::min_element(local_time_step.begin(), local_time_step.end()) * c;
//}
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//void Cells_HOM<Governing_Equation, Reconstruction_Method>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) const {
//    const auto num_solution = solution_coefficients.size();
//        
//    for (uint i = 0; i < num_solution; ++i) {
//        const auto& solution_coefficient = solution_coefficients[i];
//        const auto& basis_qnodes = This_::set_of_basis_qnodes_[i];
//        const auto solution_qnodes = solution_coefficient * basis_qnodes;
//
//        const auto [num_eq, num_quadrature_node] = solution_qnodes.size();
//        Matrix flux_quadrature_points(num_eq, This_::space_dimension_ * num_quadrature_node);
//
//        for (size_t j = 0; j < num_quadrature_node; ++j) {
//            const auto physical_flux = Governing_Equation::physical_flux(solution_qnodes.column<This_::num_equation_>(j));
//            flux_quadrature_points.change_columns(j * This_::space_dimension_, physical_flux);
//        }
//
//        Residual_ delta_rhs;
//        ms::gemm(flux_quadrature_points, This_::qweights_gradient_basis_[i], delta_rhs);
//
//        RHS[i] += delta_rhs;
//    }
//}
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//template <typename Initial_Condition>
//auto Cells_HOM<Governing_Equation, Reconstruction_Method>::calculate_initial_solutions(void) const {
//    const auto num_cell = this->quadrature_rules_.size();
//    std::vector<Static_Matrix<This_::num_equation_, This_::num_basis_>> initial_solution_coefficients(num_cell);
//    
//    for (uint i = 0; i < num_cell; ++i) {
//        const auto& qnodes = this->quadrature_rules_[i].points;
//        const auto& qweights = this->quadrature_rules_[i].weights;
//
//        const auto num_qnode = qnodes.size();
//
//        Matrix initial_solution_qnodes(This_::num_equation_, num_qnode);
//        Matrix basis_weight(num_qnode, This_::num_basis_);
//
//        for (ushort q = 0; q < num_qnode; ++q) {
//            initial_solution_qnodes.change_column(q, Initial_Condition::calculate_solution(qnodes[q]));
//            basis_weight.change_row(q, this->reconstruction_method_.calculate_basis_node(i, qnodes[q]) * qweights[q]);
//        }
//
//        ms::gemm(initial_solution_qnodes, basis_weight, initial_solution_coefficients[i]);
//    }
//
//    return initial_solution_coefficients;
//}
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//void Cells_HOM<Governing_Equation, Reconstruction_Method>::estimate_error(const std::vector<Solution_Coefficient_>& solution_coefficients, const double time) const {
//    Log::content_ << "================================================================================\n";
//    Log::content_ << "\t\t\t\t Error Anlysis\n";
//    Log::content_ << "================================================================================\n";
//
//    const auto num_cell = solution_coefficients.size();
//    
//    double arithmetic_mean_L1_error = 0.0;
//    double arithmetic_mean_L2_error = 0.0;
//    double arithmetic_mean_Linf_error = 0.0;
//
//    for (size_t i = 0; i < num_cell; ++i) {
//        const auto& qnodes = this->quadrature_rules_[i].points;
//        const auto& qweights = this->quadrature_rules_[i].weights;
//        const auto num_qnode = qnodes.size();
//
//        const auto exact_solutions = Sine_Wave<space_dimension_>::calculate_exact_solutions(qnodes, time);
//        const auto computed_solutions = solution_coefficients[i] * this->set_of_basis_qnodes_[i];
//
//        double local_L1_error = 0.0;
//        double local_L2_error = 0.0;
//        double local_Linf_error = 0.0;
//
//        double volume = 0.0;
//        for (size_t q = 0; q < num_qnode; ++q) {
//            const auto local_diff = (exact_solutions[q] - computed_solutions.column<This_::num_equation_>(q)).L1_norm();
//            local_L1_error += local_diff * qweights[q];
//            local_L2_error += local_diff * local_diff * qweights[q];
//            local_Linf_error = (std::max)(local_Linf_error, local_diff);
//
//            volume += qweights[q];
//        }
//
//        local_L1_error = local_L1_error / volume;
//        local_L2_error = std::sqrt(local_L2_error / volume);
//
//        arithmetic_mean_L1_error += local_L1_error;
//        arithmetic_mean_L2_error += local_L2_error;
//        arithmetic_mean_Linf_error += local_Linf_error;
//    }
//
//    arithmetic_mean_L1_error = arithmetic_mean_L1_error / num_cell;
//    arithmetic_mean_L2_error = arithmetic_mean_L2_error / num_cell;
//    arithmetic_mean_Linf_error = arithmetic_mean_Linf_error / num_cell;
//
//    Log::content_ << std::left << std::setprecision(16);
//    Log::content_ << std::setw(25) << "L1 error" << std::setw(25) << "L2 error" << "Linf error \n";
//    Log::content_ << std::setw(25) << arithmetic_mean_L1_error << std::setw(25) << arithmetic_mean_L2_error << arithmetic_mean_Linf_error << "\n";
//
//    std::string error_str = ms::double_to_string(arithmetic_mean_L1_error) + " " + ms::double_to_string(arithmetic_mean_L2_error) + " " + ms::double_to_string(arithmetic_mean_Linf_error) + "\n";
//    Log::write_error_text(error_str);
//}
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//void Cells_HOM<Governing_Equation, Reconstruction_Method>::initialize_scaling_method(void) const {
//    Solution_Scaler<This_::space_dimension_>::record_cell_basis_qnodes(this->set_of_basis_qnodes_);
//}