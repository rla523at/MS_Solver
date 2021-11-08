#include "../INC/Cells.h"

Cells::Cells(const Configuration& configuration, const Grid& grid)
{    
    this->governing_equation_ = Governing_Equation_Factory::make(configuration);
    this->time_step_calculator_ = Time_Step_Calculator_Factory::make(configuration, grid);
    this->num_cells_ = grid.num_cells();
}

Cells_HOM::Cells_HOM(const Configuration& configuration, const Grid& grid)
    :   Cells(configuration, grid),
        discrete_solution_(configuration, *this->governing_equation_, grid)    
{
    this->num_equations_ = this->governing_equation_->num_equations();
    this->space_dimension_ = this->governing_equation_->space_dimension();

    const auto& solution_degrees = this->discrete_solution_.get_solution_degrees();

    const auto quadrature_rules = grid.cell_quadrature_rules(this->discrete_solution_.degree() * 2);

}


void Cells_HOM::update_solution(Euclidean_Vector&& updated_solution_v)
{
    this->discrete_solution_.update_solution(std::move(updated_solution_v));
}

double Cells_HOM::calculate_time_step(void) const
{
    const auto P0_solutions = this->discrete_solution_.calculate_P0_solutions(this->P0_basis_values_);
    return this->time_step_calculator_->calculate(P0_solutions, *this->governing_equation_);
}

void Cells_HOM::calculate_RHS(double* rhs) const
{
    //this routine can be changed using sparse matrix
    for (uint i = 0; i < this->num_cells_; ++i)     
    {
        const auto solution_at_QPs = this->discrete_solution_.calculate_solution_at_points(i, this->set_of_basis_QPs_m_[i]);
        const auto num_QPs = solution_at_QPs.size();

        Matrix flux_QPs_m(this->num_equations_, this->space_dimension_ * num_QPs);

        for (int j = 0; j < num_QPs; ++j)
        {
            const auto physical_flux = this->governing_equation_->calculate_physical_flux(solution_at_QPs[j]);
            flux_QPs_m.change_columns(j * this->space_dimension_, physical_flux);
        }

        const auto delta_rhs = flux_QPs_m * this->set_of_QWs_gradient_basis_m_[i];
        this->update_rhs(i, rhs, delta_rhs);
    }
}

const Euclidean_Vector& Cells_HOM::get_solution_vector(void) const
{
    return this->discrete_solution_.get_solution_vector();
}

void Cells_HOM::update_rhs(const uint cell_index, double* rhs_ptr, const Matrix& delta_rhs) const
{
    const auto n = static_cast<MKL_INT>(delta_rhs.num_values());
    const auto a = 1.0;
    const auto incx = 1;
    const auto incy = 1;

    auto icell_rhs_ptr = rhs_ptr + this->discrete_solution_.coefficient_start_index(cell_index);
    cblas_daxpy(n, a, delta_rhs.data(), incx, icell_rhs_ptr, incy);
}

size_t Cells_HOM::num_solution_values(void) const
{
    return this->discrete_solution_.num_values();
}