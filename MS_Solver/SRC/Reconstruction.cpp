#include "../INC/Reconstruction.h"

hMLP_Reconstruction_DG::hMLP_Reconstruction_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution)
	: stability_criterion_(grid, discrete_solution)
	, indicator_(grid, discrete_solution) 
{
	this->num_cells_ = static_cast<uint>(grid.num_cells());
	this->set_of_P1_projected_criterion_values_at_verticies_.resize(this->num_cells_);
};


void hMLP_Reconstruction_DG::reconstruct(Discrete_Solution_DG& discrete_solution)
{
	this->precalculate(discrete_solution);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		auto solution_degree = discrete_solution.solution_degree(cell_index);

		this->indicator_.set_precalculated_result(&this->set_of_P1_projected_criterion_values_at_verticies_[cell_index]);
		this->limiter_.set_precalculated_result(&this->set_of_P1_projected_criterion_values_at_verticies_[cell_index]);

		while (true)
		{
			const auto cell_type = this->indicator_.indicate_opt(discrete_solution, cell_index, this->stability_criterion_);

			if (cell_type == cell_type::trouble)
			{
				if (solution_degree <= 2)
				{
					const auto limiting_value = this->limiter_.limiter_function_opt(discrete_solution,cell_index, this->stability_criterion_);

					solution_degree = 1;
					discrete_solution.project_to_Pn_space(cell_index, solution_degree);
					discrete_solution.limit_slope(cell_index, limiting_value);
					break;
				}
				else
				{
					discrete_solution.project_to_Pn_space(cell_index, --solution_degree);
				}
			}
			else
			{
				break;
			}
		}
	}
};

void hMLP_Reconstruction_DG::precalculate(const Discrete_Solution_DG& discrete_solution)
{
	const auto criterion_equation_index = this->stability_criterion_.get_criterion_equation_index();

	for (uint i = 0; i < this->num_cells_; ++i)
	{
		this->set_of_P1_projected_criterion_values_at_verticies_[i] = discrete_solution.calculate_P1_projected_nth_solution_at_vertices(i, criterion_equation_index);
	}

	this->stability_criterion_.caclulate(discrete_solution);
}

std::unique_ptr<Reconstruction_DG> Reconstruction_DG_Factory::make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution)
{
	const auto& reconstruction_scheme = configuration.get_reconstruction_scheme();

	if (ms::contains_icase(reconstruction_scheme, "no"))
	{
		return std::make_unique<No_Reconstruction_DG>();
	}
	else if (ms::contains_icase(reconstruction_scheme, "hMLP"))
	{
		return std::make_unique<hMLP_Reconstruction_DG>(grid, discrete_solution);
	}
	else
	{
		EXCEPTION("reconstruction shceme in configuration file is not supported");
		return nullptr;
	}
}