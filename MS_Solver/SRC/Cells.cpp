#include "../INC/Cells.h"

Cells::Cells(const std::shared_ptr<Governing_Equation>& governing_equation, std::unique_ptr<Time_Step_Calculator>&& time_step_calculator, const Grid& grid)
    :governing_equation_(governing_equation),
     time_step_calculator_(std::move(time_step_calculator))
{
    this->num_cells_ = grid.num_cells();
}

Cells_DG::Cells_DG(const std::shared_ptr<Governing_Equation>& governing_equation, std::unique_ptr<Time_Step_Calculator>&& time_step_calculator, const Grid& grid, Discrete_Solution_DG& discrete_solution)
    : Cells(governing_equation, std::move(time_step_calculator), grid)
{
    Profiler::set_time_point();

    const auto space_dimension = this->governing_equation_->space_dimension();
    this->set_of_QWs_gradient_basis_m_.resize(this->num_cells_);    

    //for construct optimization
    const auto num_equations = this->governing_equation_->num_equations();
    const auto num_solutions = discrete_solution.num_solutions();

    this->set_of_num_QPs_.resize(this->num_cells_);
    this->set_of_flux_QPs_m_.resize(this->num_cells_);
    //this->solution_v_at_QPs_.resize(this->max_num_QPs, Euclidean_Vector(num_solutions));
    this->solution_v_at_QPs_.fill(Euclidean_Vector(num_solutions));
    this->physical_flux_ = Matrix(num_equations, space_dimension);
    //

    //for precalculation
    std::vector<Quadrature_Rule> quadrature_rules(this->num_cells_);
    //

    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto solution_degree = discrete_solution.solution_degree(cell_index);
        const auto integrand_degree = solution_degree * 2;
        const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

        const auto& QPs = quadrature_rule.points;
        const auto& QWs = quadrature_rule.weights;
        const auto num_QPs = QPs.size();

        const auto num_basis = discrete_solution.num_basis(cell_index);
        Matrix QWs_gradient_basis_m(num_QPs * space_dimension, num_basis);

        const auto transposed_gradient_basis = discrete_solution.calculate_tranposed_gradient_basis(cell_index);
        for (int q = 0; q < num_QPs; ++q)
        {
            const auto start_row_index = q * space_dimension;
            const auto part_of_QWs_gradient_basis_m = transposed_gradient_basis(QPs[q]) * QWs[q];
            QWs_gradient_basis_m.change_rows(start_row_index, part_of_QWs_gradient_basis_m);
        }

        this->set_of_QWs_gradient_basis_m_[cell_index] = std::move(QWs_gradient_basis_m);

        //for construct optimization
        this->set_of_num_QPs_[cell_index] = static_cast<ushort>(num_QPs);
        this->set_of_flux_QPs_m_[cell_index] = Matrix(num_equations, space_dimension * num_QPs);
        //

        //for precalculation
        quadrature_rules[cell_index] = quadrature_rule;
        //
    }

    //for precalculation
    discrete_solution.precalculate_cell_P0_basis_values();
    discrete_solution.precalcualte_cell_QPs_basis_values(quadrature_rules);
    //

    LOG << std::left << std::setw(50) << "@ Cells DG precalculation" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << LOG.print_;
}

double Cells_DG::calculate_time_step(const Discrete_Solution_DG& discrete_solution) const
{
    const auto P0_solutions = discrete_solution.calculate_P0_solutions();
    const auto time_step = this->time_step_calculator_->calculate(P0_solutions, *this->governing_equation_);

    const auto maximum_solution_degree = static_cast<double>(discrete_solution.maximum_solution_degree());
    const auto relaxation_factor = 1.0 / (2 * maximum_solution_degree + 1);

    return time_step * relaxation_factor;
}

void Cells_DG::calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_solution) const
{
    const auto num_equations = this->governing_equation_->num_equations();
    const auto space_dimension = this->governing_equation_->space_dimension();

    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {  
        discrete_solution.calculate_solution_at_cell_QPs(this->solution_v_at_QPs_.data(), cell_index);

        const auto num_QPs = this->set_of_num_QPs_[cell_index];
        auto& flux_QPs_m = this->set_of_flux_QPs_m_[cell_index];

        for (int j = 0; j < num_QPs; ++j)
        {
            this->governing_equation_->calculate_physical_flux(this->physical_flux_, this->solution_v_at_QPs_[j]);

            const auto start_column_index = j * space_dimension;
            flux_QPs_m.change_columns(start_column_index, this->physical_flux_);
        }

        const auto num_values = discrete_solution.num_values(cell_index);
        this->reset_residual_values(num_values);

        ms::gemm(flux_QPs_m, this->set_of_QWs_gradient_basis_m_[cell_index], this->residual_values_.data());
        residual.update_rhs(cell_index, this->residual_values_.data());
    }
}

void Cells_DG::reset_residual_values(const uint num_values) const
{
    std::fill(this->residual_values_.begin(), this->residual_values_.begin() + num_values, 0.0);
}