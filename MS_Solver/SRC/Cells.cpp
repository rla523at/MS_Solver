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
    const auto num_equations = this->governing_equation_->num_equations();
    const auto space_dimension = this->governing_equation_->space_dimension();

    this->set_of_QWs_gradient_basis_m_.resize(this->num_cells_);

    std::vector<Quadrature_Rule> quadrature_rules(this->num_cells_);

    for (int cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto solution_degree = discrete_solution.solution_degree(cell_index);
        const auto integrand_degree = solution_degree * 2;
        auto quadrature_rule = grid.cell_quadrature_rule(cell_index, solution_degree);

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
        quadrature_rules[cell_index] = std::move(quadrature_rule);
    }

    discrete_solution.precalculate_cell_P0_basis_values();
    discrete_solution.precalcualte_cell_QPs_basis_values(quadrature_rules);
}

double Cells_DG::calculate_time_step(const Discrete_Solution_DG& discrete_solution) const
{
    const auto P0_solutions = discrete_solution.calculate_P0_solutions();
    return this->time_step_calculator_->calculate(P0_solutions, *this->governing_equation_);
}

void Cells_DG::calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_solution) const
{
    const auto num_equations = this->governing_equation_->num_equations();
    const auto space_dimension = this->governing_equation_->space_dimension();

    //this routine can be changed using sparse matrix
    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto solution_at_QPs = discrete_solution.calculate_solution_at_cell_QPs(cell_index);
        const auto num_QPs = solution_at_QPs.size();

        Matrix flux_QPs_m(num_equations, space_dimension * num_QPs);

        for (int j = 0; j < num_QPs; ++j)
        {
            const auto start_column_index = j * space_dimension;
            const auto physical_flux = this->governing_equation_->calculate_physical_flux(solution_at_QPs[j]);
            flux_QPs_m.change_columns(start_column_index, physical_flux);
        }

        const auto delta_rhs = flux_QPs_m * this->set_of_QWs_gradient_basis_m_[cell_index];
        residual.update_rhs(cell_index, delta_rhs);
    }
}































//Cells::Cells(const Configuration& configuration, const Grid& grid)
//{    
//    this->governing_equation_ = Governing_Equation_Factory::make(configuration);
//    this->time_step_calculator_ = Time_Step_Calculator_Factory::make(configuration, grid);
//    this->num_cells_ = grid.num_cells();
//}
//
//Cells_DG::Cells_DG(const Configuration& configuration, const Grid& grid)
//    :   Cells(configuration, grid),
//        discrete_solution_(configuration, grid, *this->governing_equation_)
//{
//    const auto num_equations = this->governing_equation_->num_equations();
//    const auto space_dimension = this->governing_equation_->space_dimension();
//
//    this->P0_basis_values_.resize(this->num_cells_);
//    this->set_of_basis_QPs_m_.resize(this->num_cells_);
//    this->set_of_QWs_gradient_basis_m_.resize(this->num_cells_);
//    for (int cell_index = 0; cell_index < this->num_cells_; ++cell_index)
//    {        
//        const auto solution_degree = this->discrete_solution_.solution_degree(cell_index);
//        const auto integrand_degree = solution_degree * 2;
//        const auto quadrature_rule = grid.cell_quadrature_rule(cell_index, solution_degree);
//
//        const auto& QPs = quadrature_rule.points;
//        const auto& QWs = quadrature_rule.weights;
//        const auto num_QPs = QPs.size();
//
//        const auto num_basis = this->discrete_solution_.num_basis(cell_index);
//        Matrix QWs_gradient_basis_m(num_QPs * space_dimension, num_basis);
//
//        const auto transposed_gradient_basis = this->discrete_solution_.calculate_tranposed_gradient_basis(cell_index);
//        for (int q = 0; q < num_QPs; ++q)
//        {
//            const auto part_of_QWs_gradient_basis_m = transposed_gradient_basis(QPs[q]) * QWs[q];
//            QWs_gradient_basis_m.change_rows(q * space_dimension, part_of_QWs_gradient_basis_m);
//        }
//
//        this->P0_basis_values_[cell_index] = this->discrete_solution_.calculate_P0_basis_value(cell_index);
//        this->set_of_basis_QPs_m_[cell_index] = this->discrete_solution_.calculate_basis_points_m(cell_index, QPs);
//        this->set_of_QWs_gradient_basis_m_[cell_index] = std::move(QWs_gradient_basis_m);        
//    }
//}
//
//
//void Cells_DG::update_solution(Euclidean_Vector&& updated_solution_v)
//{
//    this->discrete_solution_.update_solution(std::move(updated_solution_v));
//}
//
//double Cells_DG::calculate_time_step(void) const
//{
//    const auto P0_solutions = this->discrete_solution_.calculate_P0_solutions(this->P0_basis_values_);
//    return this->time_step_calculator_->calculate(P0_solutions, *this->governing_equation_);
//}
//
//void Cells_DG::calculate_RHS(double* rhs) const
//{
//    const auto num_equations = this->governing_equation_->num_equations();
//    const auto space_dimension = this->governing_equation_->space_dimension();
//
//    //this routine can be changed using sparse matrix
//    for (uint i = 0; i < this->num_cells_; ++i)     
//    {
//        const auto solution_at_QPs = this->discrete_solution_.calculate_solution_at_points(i, this->set_of_basis_QPs_m_[i]);
//        const auto num_QPs = solution_at_QPs.size();
//
//        Matrix flux_QPs_m(num_equations, space_dimension * num_QPs);
//
//        for (int j = 0; j < num_QPs; ++j)
//        {
//            const auto physical_flux = this->governing_equation_->calculate_physical_flux(solution_at_QPs[j]);
//            flux_QPs_m.change_columns(j * space_dimension, physical_flux);
//        }
//
//        const auto delta_rhs = flux_QPs_m * this->set_of_QWs_gradient_basis_m_[i];
//        this->update_rhs(i, rhs, delta_rhs);
//    }
//}
//
//const Euclidean_Vector& Cells_DG::get_solution_vector(void) const
//{
//    return this->discrete_solution_.get_solution_vector();
//}
//
//void Cells_DG::update_rhs(const uint cell_index, double* rhs_ptr, const Matrix& delta_rhs) const
//{
//    const auto n = static_cast<MKL_INT>(delta_rhs.num_values());
//    const auto a = 1.0;
//    const auto incx = 1;
//    const auto incy = 1;
//
//    auto icell_rhs_ptr = rhs_ptr + this->discrete_solution_.coefficient_start_index(cell_index);
//    cblas_daxpy(n, a, delta_rhs.data(), incx, icell_rhs_ptr, incy);
//}
//
////std::vector<ushort> Cells_DG::calculate_integrand_degrees(const std::vector<ushort>& solution_degrees) const
////{
////    std::vector<ushort> integrand_degrees(this->num_cells_);
////
////    for (int i = 0; i < this->num_cells_; ++i)
////        integrand_degrees[i] = 2 * solution_degrees[i];
////
////    return integrand_degrees;
////}
//
//size_t Cells_DG::num_solution_values(void) const
//{
//    return this->discrete_solution_.num_values();
//}
//
