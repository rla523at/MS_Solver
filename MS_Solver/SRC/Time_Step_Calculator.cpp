#include "../INC/Time_Step_Calculator.h"

CFL::CFL(const double cfl, const Grid& grid)
    :cfl_(cfl)
{
    this->cell_volumes_ = grid.cell_volumes();
    this->cell_projected_volumes_ = grid.cell_projected_volumes();
}

double CFL::calculate(const Discrete_Solution & discrete_solution, const Governing_Equation & governing_equation) const
{
    const auto space_dimension = governing_equation.space_dimension();

    const auto P0_solutions = discrete_solution.calculate_P0_solutions();
    const auto coordinate_projected_maximum_lambdas = governing_equation.calculate_coordinate_projected_maximum_lambda(P0_solutions);

    const auto num_cell = P0_solutions.size();
    std::vector<double> local_time_step(num_cell);

    for (size_t i = 0; i < num_cell; ++i)
    {
        double radii = 0;
        for (ushort j = 0; j < space_dimension; ++j)
        {
            radii += this->cell_projected_volumes_[i][j] * coordinate_projected_maximum_lambdas[i][j];
        }

        local_time_step[i] = this->cfl_ * this->cell_volumes_[i] / radii;
    }

    return *std::min_element(local_time_step.begin(), local_time_step.end());
}

double Constant_Dt::calculate(const Discrete_Solution& discrete_solution, const Governing_Equation& governing_equation) const
{
    return this->constant_dt_;
}