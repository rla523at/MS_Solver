#include "../INC/Reconstruction.h"

hMLP_Reconstruction_DG::hMLP_Reconstruction_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution)
	: stability_criterion_(grid, discrete_solution)
	, MLP_indicator_(grid, discrete_solution) 
	, discontinuity_indicator_(grid, discrete_solution, 0) // test
{
	this->num_cells_ = grid.num_cells();

	const auto set_of_num_vertices = grid.cell_set_of_num_vertices();
	
	this->set_of_P1_projected_criterion_values_at_verticies_.resize(this->num_cells_);
	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		this->set_of_P1_projected_criterion_values_at_verticies_[cell_index].resize(set_of_num_vertices[cell_index]);
	}
};

//hMLP_Reconstruction_DG::hMLP_Reconstruction_DG(MLP_Criterion&& stability_criterion, MLP_Indicator&& MLP_indicator, MLP_u1_Limiter&& MLP_u1_limiter, const Grid& grid)
//	: stability_criterion_(std::move(stability_criterion))
//	, MLP_indicator_(std::move(MLP_indicator))
//	, MLP_u1_limiter_(std::move(MLP_u1_limiter))
//{
//	this->num_cells_ = grid.num_cells();
//
//	const auto set_of_num_vertices = grid.cell_set_of_num_vertices();
//
//	this->set_of_P1_projected_criterion_values_at_verticies_.resize(this->num_cells_);
//	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
//	{
//		this->set_of_P1_projected_criterion_values_at_verticies_[cell_index].resize(set_of_num_vertices[cell_index]);
//	}
//}


void hMLP_Reconstruction_DG::reconstruct(Discrete_Solution_DG& discrete_solution)
{
	//test
	this->discontinuity_indicator_.precalculate(discrete_solution);
	const auto& discontinuity_factor = this->discontinuity_indicator_.get_discontinuity_factor();
	Post_Processor::record_solution();
	Post_Processor::record_variables("discontinuity_factor", discontinuity_factor);
	Post_Processor::post_solution();
	//test

	this->precalculate(discrete_solution);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		auto solution_degree = discrete_solution.solution_degree(cell_index);

		this->MLP_indicator_.set_precalculated_result(&this->set_of_P1_projected_criterion_values_at_verticies_[cell_index]);
		this->MLP_u1_limiter_.set_precalculated_result(&this->set_of_P1_projected_criterion_values_at_verticies_[cell_index]);

		while (true)
		{
			const auto cell_type = this->MLP_indicator_.indicate_opt(discrete_solution, cell_index, this->stability_criterion_);				

			if (cell_type == cell_type::trouble)
			{
				if (solution_degree <= 2)
				{
					solution_degree = 1;
					discrete_solution.project_to_Pn_space(cell_index, solution_degree);

					const auto limiting_value = this->MLP_u1_limiter_.limiter_function_opt(discrete_solution,cell_index, this->stability_criterion_);
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

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		discrete_solution.calculate_P1_projected_nth_solution_at_vertices(this->set_of_P1_projected_criterion_values_at_verticies_[cell_index].data(), cell_index, criterion_equation_index);
	}

	this->stability_criterion_.precaclulate(discrete_solution);
}

hMLP_BD_Reconstruction_DG::hMLP_BD_Reconstruction_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::string& governing_equation_name)
	: stability_criterion_(grid, discrete_solution)
	, MLP_indicator_(grid,discrete_solution)
	, subcell_oscillation_indicator(grid,discrete_solution)
	, shock_indicator_(grid,discrete_solution,governing_equation_name)
	, discontinuity_indicator_(grid,discrete_solution, 0)
{
	this->num_cells_ = grid.num_cells();

	const auto set_of_num_vertices = grid.cell_set_of_num_vertices();

	this->set_of_P1_projected_criterion_values_at_verticies_.resize(this->num_cells_);
	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		this->set_of_P1_projected_criterion_values_at_verticies_[cell_index].resize(set_of_num_vertices[cell_index]);
	}
}

//
//hMLP_BD_Reconstruction_DG::hMLP_BD_Reconstruction_DG(Simplex_Decomposed_MLP_Criterion&& stability_criterion, MLP_Indicator&& MLP_indicator, MLP_u1_Limiter&& MLP_u1_limiter,
//	Subcell_Oscillation_Indicator&& boundary_indicator, Shock_Indicator&& shock_indicator, const Grid& grid)
//	: stability_criterion_(std::move(stability_criterion))
//	, MLP_indicator_(std::move(MLP_indicator))
//	, MLP_u1_limiter_(std::move(MLP_u1_limiter))
//	, subcell_oscillation_indicator(std::move(boundary_indicator))
//	, shock_indicator_(std::move(shock_indicator))
//{
//	this->num_cells_ = grid.num_cells();
//
//	const auto set_of_num_vertices = grid.cell_set_of_num_vertices();
//
//	this->set_of_P1_projected_criterion_values_at_verticies_.resize(this->num_cells_);
//	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
//	{
//		this->set_of_P1_projected_criterion_values_at_verticies_[cell_index].resize(set_of_num_vertices[cell_index]);
//	}
//}

void hMLP_BD_Reconstruction_DG::reconstruct(Discrete_Solution_DG& discrete_solution)
{
	//test
	this->discontinuity_indicator_.precalculate(discrete_solution);
	const auto& discontinuity_factor = this->discontinuity_indicator_.get_discontinuity_factor();
	Post_Processor::record_solution();
	Post_Processor::record_variables("discontinuity_factor", discontinuity_factor);
	Post_Processor::post_solution();
	//test

	this->precalculate(discrete_solution);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		this->is_end_limting_ = false;

		auto solution_degree = discrete_solution.solution_degree(cell_index);

		this->MLP_indicator_.set_precalculated_result(&this->set_of_P1_projected_criterion_values_at_verticies_[cell_index]);
		this->MLP_u1_limiter_.set_precalculated_result(&this->set_of_P1_projected_criterion_values_at_verticies_[cell_index]);

		while (true)
		{
			const auto cell_type = this->MLP_indicator_.indicate_opt(discrete_solution, cell_index, this->stability_criterion_);
			
			switch (cell_type)
			{
			case cell_type::normal:
			{
				if (this->subcell_oscillation_indicator.is_typeI_cell(cell_index) && this->shock_indicator_.is_shock(cell_index))
				{
					solution_degree = 1;
					this->limit_solution(discrete_solution, cell_index, solution_degree);
				}
				else
				{
					this->is_end_limting_ = true;
				}

				break;
			}
			case cell_type::smooth_extrema:
			{
				if (this->subcell_oscillation_indicator.is_typeII_cell(cell_index))
				{
					this->limit_solution(discrete_solution, cell_index, solution_degree);
				}
				else
				{
					this->is_end_limting_ = true;
				}

				break;
			}
			case cell_type::trouble:
			{
				this->limit_solution(discrete_solution, cell_index, solution_degree);
				break;
			}
			default:
				EXCEPTION("wrong cell type");
				break;
			}

			if (this->is_end_limting_)
			{
				break;
			}
		}
	}
};

void hMLP_BD_Reconstruction_DG::precalculate(const Discrete_Solution_DG& discrete_solution)
{
	const auto criterion_equation_index = this->stability_criterion_.get_criterion_equation_index();

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		discrete_solution.calculate_simplex_P1_projected_nth_solution_at_vertices(this->set_of_P1_projected_criterion_values_at_verticies_[cell_index].data(), cell_index, criterion_equation_index);
	}

	this->stability_criterion_.precaclulate(discrete_solution);
	this->subcell_oscillation_indicator.precalculate(discrete_solution);
	this->shock_indicator_.precalculate(discrete_solution);
}

void hMLP_BD_Reconstruction_DG::limit_solution(Discrete_Solution_DG& discrete_solution, const uint cell_index, ushort& solution_degree) const
{
	if (solution_degree <= 2)
	{
		solution_degree = 1;
		discrete_solution.project_to_Pn_space(cell_index, solution_degree);

		const auto limiting_value = this->MLP_u1_limiter_.limiter_function_opt(discrete_solution, cell_index, this->stability_criterion_);
		discrete_solution.limit_slope(cell_index, limiting_value);
		this->is_end_limting_ = true;
	}
	else
	{
		discrete_solution.project_to_Pn_space(cell_index, --solution_degree);
	}
}

void Test_Reconstuction_DG::reconstruct(Discrete_Solution_DG& discrete_solution)
{
	this->discontinuity_indicator_.precalculate(discrete_solution);
	const auto& discontinuity_factor = this->discontinuity_indicator_.get_discontinuity_factor();
	Post_Processor::record_solution();
	Post_Processor::record_variables("discontinuity_factor", discontinuity_factor);
	Post_Processor::post_solution();
};

std::unique_ptr<Reconstruction_DG> Reconstruction_DG_Factory::make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution)
{
	const auto& reconstruction_scheme = configuration.get_reconstruction_scheme();

	if (ms::compare_icase(reconstruction_scheme, "no"))
	{
		return std::make_unique<No_Reconstruction_DG>();
	}
	else if (ms::compare_icase(reconstruction_scheme, "hMLP"))
	{
		//MLP_Criterion MLP_criterion(grid, discrete_solution);
		//MLP_Indicator MLP_indicator(grid, discrete_solution);
		//MLP_u1_Limiter MLP_u1_limiter;

		//return std::make_unique<hMLP_Reconstruction_DG>(std::move(MLP_criterion), std::move(MLP_indicator), std::move(MLP_u1_limiter), grid);
		return std::make_unique<hMLP_Reconstruction_DG>(grid, discrete_solution);
	}
	else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD"))
	{
		const auto& governing_equation_name = configuration.get_governing_equation();

		//Simplex_Decomposed_MLP_Criterion simplex_decomposed_MLP_criterion(grid, discrete_solution);
		//MLP_Indicator MLP_indicator(grid, discrete_solution);
		//MLP_u1_Limiter MLP_u1_limiter;
		//Subcell_Oscillation_Indicator subcell_oscillation_indicator(grid, discrete_solution);
		//Shock_Indicator shock_indicator(grid, discrete_solution, governing_equation_name);

		//return std::make_unique<hMLP_BD_Reconstruction_DG>(std::move(simplex_decomposed_MLP_criterion), std::move(MLP_indicator),
		//	std::move(MLP_u1_limiter), std::move(subcell_oscillation_indicator), std::move(shock_indicator), grid);

		return std::make_unique<hMLP_BD_Reconstruction_DG>(grid, discrete_solution, governing_equation_name);
	}
	else if (ms::compare_icase(reconstruction_scheme, "Test"))
	{
		Discontinuity_Indicator discontinuity_indicator(grid, discrete_solution, 0);
		return std::make_unique<Test_Reconstuction_DG>(std::move(discontinuity_indicator));
	}
	else
	{
		EXCEPTION("reconstruction shceme in configuration file is not supported");
		return nullptr;
	}
}