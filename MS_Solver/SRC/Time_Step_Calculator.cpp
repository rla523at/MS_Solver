#include "../INC/Time_Step_Calculator.h"

CFL::CFL(const double cfl, const Grid& grid)
    :cfl_(cfl)
{
    this->cell_index_to_volume_reciprocal_table_ = grid.cell_index_to_volume_table();
    this->cell_index_to_projected_volumes_table_ = grid.cell_index_to_projected_volumes_table();
}

double CFL::calculate(const std::vector<Euclidean_Vector>& P0_solutions, const Governing_Equation& governing_equation) const
{
    const auto space_dimension = governing_equation.space_dimension();

    const auto cell_index_to_coordinate_projected_maximum_lambdas_table = governing_equation.calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(P0_solutions);

    const auto num_cell = P0_solutions.size();
    std::vector<double> allowable_local_time_steps(num_cell);

    for (size_t cell_index = 0; cell_index < num_cell; ++cell_index)
    {
        const auto& projected_volumes = this->cell_index_to_projected_volumes_table_[cell_index];
        const auto& coordinate_projected_maximum_lambdas = cell_index_to_coordinate_projected_maximum_lambdas_table[cell_index];

        double radii = 0;
        for (ushort j = 0; j < space_dimension; ++j)
        {
            radii += projected_volumes[j] * coordinate_projected_maximum_lambdas[j];
        }

        allowable_local_time_steps[cell_index] = this->cfl_ * this->cell_index_to_volume_reciprocal_table_[cell_index] / radii;
    }

    return *std::min_element(allowable_local_time_steps.begin(), allowable_local_time_steps.end());
}

double Constant_Dt::calculate(const std::vector<Euclidean_Vector>& P0_solutions, const Governing_Equation& governing_equation) const
{
    return this->constant_dt_;
}

std::unique_ptr<Time_Step_Calculator> Time_Step_Calculator_Factory::make_unique(const Configuration& configuration, const Grid& grid)
{
    const auto time_step_method = configuration.get_time_step_method();
    if (ms::contains_icase(time_step_method, "CFL"))
    {
        const auto cfl_number = configuration.CFL_number();
        return std::make_unique<CFL>(cfl_number, grid);
    }
    else if (ms::contains_icase(time_step_method, "Constant", "dt"))
    {
        const auto constant_dt = configuration.constant_dt();
        return std::make_unique<Constant_Dt>(constant_dt);
    }
    else
    {
        EXCEPTION("time step calculator type in configuration file is not supported");
        return nullptr;
    }
}