#include "../INC/Reconstruction_Impl.h"

std::unique_ptr<Discontinuity_Indicating_Function> Discontinuity_Indicator_Factory::make_unique(const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
{
	if (ms::compare_icase(type_name, "Always_Flase"))
	{
		return std::make_unique<Always_False_Discontinuity_Indicator>();
	}
	else if (ms::compare_icase(type_name, "Always_True"))
	{
		return std::make_unique<Always_True_Discontinuity_Indicator>();
	}
	else if (ms::compare_icase(type_name, "Scaled_Difference"))
	{
		const auto num_cells = grid.num_cells();
		constexpr auto allowable_diff = 0.1; //scaled Difference가 10%보다 크면 discontinuity라 보겠다.

		std::vector<double> cell_index_to_threshold_value_table(num_cells, allowable_diff);
		auto measuring_function = std::make_unique<Scaled_Difference>(grid, discrete_solution, criterion_solution_index);

		return std::make_unique<Greater_is_Discontinuity>(cell_index_to_threshold_value_table, std::move(measuring_function));
	}
	else if (ms::compare_icase(type_name, "Extrapolation_Difference"))
	{
		auto cell_index_to_characteristic_length_table = grid.cell_index_to_characteristic_length_table();//h보다 크면 discontinuity로 보겠다.
		auto measuring_function = std::make_unique<Extrapolation_Difference>(grid, discrete_solution, criterion_solution_index); 

		return std::make_unique<Greater_is_Discontinuity>(cell_index_to_characteristic_length_table, std::move(measuring_function));
	}
	else if (ms::compare_icase(type_name, "Difference_Of_Extrapolation_Difference"))
	{
		auto cell_index_to_characteristic_length_table = grid.cell_index_to_characteristic_length_table();//h보다 크면 discontinuity로 보겠다.
		auto measuring_function = std::make_unique<Difference_Of_Extrapolatiion_Difference>(grid, discrete_solution, criterion_solution_index);

		return std::make_unique<Greater_is_Discontinuity>(cell_index_to_characteristic_length_table, std::move(measuring_function));
	}
	else if (ms::compare_icase(type_name, "Max_Divergence_of_Velocity"))
	{
		const auto num_cells = grid.num_cells();
		constexpr auto zero = 0; // 0 보다 작으면 discontinuity로 보겠다.

		std::vector<double> cell_index_to_threshold_value_table(num_cells, zero);
		auto measuring_function = std::make_unique<Max_Divergence_of_Velocity>(grid, discrete_solution, criterion_solution_index);

		return std::make_unique<Less_is_Discontinuity>(cell_index_to_threshold_value_table, std::move(measuring_function));
	}
	else
	{
		EXCEPTION(type_name + " is not supported discontinuity indicator type");
		return nullptr;
	}
}

std::unique_ptr<Discontinuity_Indicating_Function> Discontinuity_Indicator_Factory::make_always_false(void)
{
	return std::make_unique<Always_False_Discontinuity_Indicator>();
}

std::unique_ptr<Discontinuity_Indicating_Function> Discontinuity_Indicator_Factory::make_always_true(void)
{
	return std::make_unique<Always_True_Discontinuity_Indicator>();
}

std::unique_ptr<Discontinuity_Indicating_Function> Shock_Indicator_Factory::make_unique(const std::string& governing_equation_name, const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
{
	if (ms::compare_icase(governing_equation_name, "Linear_Advection"))
	{
		return Discontinuity_Indicator_Factory::make_always_false();
	}
	else if (ms::compare_icase(governing_equation_name, "Burgers"))
	{
		constexpr auto criterion_solution_index = 0;
		return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, criterion_solution_index);
	}
	else if (ms::compare_icase(governing_equation_name, "Euler"))
	{
		const auto space_dimension = grid.space_dimension();

		if (space_dimension == 2)
		{
			return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, Euler_2D::pressure_index());
		}
		else if (space_dimension == 3)
		{
			return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, Euler_3D::pressure_index());
		}
		else
		{
			EXCEPTION("not supported space_dimension");
			return nullptr;
		}
	}
	else
	{
		EXCEPTION("not supported governing equation");
		return nullptr;
	}
}

std::unique_ptr<Discontinuity_Indicating_Function> Contact_Indicator_Factory::make_unique(const std::string& type_name, const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
{
	constexpr auto solution_index = 0;
	return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, solution_index);
}

std::unique_ptr<Indicator> Indicator_Factory::make_hMLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
{
	return std::make_unique<hMLP_Indicator>(grid, discrete_solution, criterion_equation_index);
}

std::unique_ptr<Indicator> Indicator_Factory::make_hMLP_BD_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
{
	auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Scaled_Difference", grid, discrete_solution);
	auto contact_indicator = Discontinuity_Indicator_Factory::make_always_true();
	return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
}

std::unique_ptr<Indicator> Indicator_Factory::make_hMLP_BD_Off_TypeII_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
{
	auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Scaled_Difference", grid, discrete_solution);
	auto contact_indicator = Discontinuity_Indicator_Factory::make_always_false();
	return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
}

std::unique_ptr<Indicator> Indicator_Factory::make_Improved_hMLP_BD1_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
{
	auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Scaled_Difference", grid, discrete_solution);
	auto contact_indicator = Contact_Indicator_Factory::make_unique("Scaled_Difference", governing_equation_name, grid, discrete_solution);
	return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
}

std::unique_ptr<Indicator> Indicator_Factory::make_Improved_hMLP_BD2_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
{
	auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Scaled_Difference", grid, discrete_solution);
	auto contact_indicator = Contact_Indicator_Factory::make_unique("Extrapolation_Difference", governing_equation_name, grid, discrete_solution);
	return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
}

std::unique_ptr<Indicator> Indicator_Factory::make_Improved_hMLP_BD3_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
{
	auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Scaled_Difference", grid, discrete_solution);
	auto contact_indicator = Contact_Indicator_Factory::make_unique("Difference_of_Extrapolation_Difference", governing_equation_name, grid, discrete_solution);
	return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
}

std::unique_ptr<Indicator> Indicator_Factory::make_Improved_hMLP_BD4_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
{
	auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Scaled_Difference", grid, discrete_solution);
	auto contact_indicator = Contact_Indicator_Factory::make_unique("Max_Divergence_of_Velocity", governing_equation_name, grid, discrete_solution);
	return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
}

std::unique_ptr<Reconstruction_DG> Reconstruction_DG_Factory::make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution)
{
	const auto& reconstruction_scheme = configuration.get_reconstruction_scheme();

	if (ms::compare_icase(reconstruction_scheme, "no"))
	{
		return std::make_unique<No_Reconstruction_DG>();
	}
	else if (ms::compare_icase(reconstruction_scheme, "hMLP"))
	{
		const auto criterion_equation_index = 0;
		auto stability_criterion = std::make_unique<MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
		auto indicator = Indicator_Factory::make_hMLP_Indicator(grid, discrete_solution, criterion_equation_index);
		auto limiter = std::make_unique <hMLP_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
	}
	else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD"))
	{
		const auto criterion_equation_index = 0;
		const auto& governing_equation_name = configuration.get_governing_equation();

		auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
		auto indicator = Indicator_Factory::make_hMLP_BD_Indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);
		auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
	}
	else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_Off_TypeII"))
	{
		const auto criterion_equation_index = 0;
		const auto& governing_equation_name = configuration.get_governing_equation();

		auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
		auto indicator = Indicator_Factory::make_hMLP_BD_Off_TypeII_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);
		auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
	}
	else if (ms::compare_icase(reconstruction_scheme, "Improved_hMLP_BD1"))
	{
		const auto criterion_equation_index = 0;
		const auto& governing_equation_name = configuration.get_governing_equation();

		auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
		auto indicator = Indicator_Factory::make_Improved_hMLP_BD1_Indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);
		auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
	}
	else if (ms::compare_icase(reconstruction_scheme, "Improved_hMLP_BD2"))
	{
		const auto criterion_equation_index = 0;
		const auto& governing_equation_name = configuration.get_governing_equation();

		auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
		auto indicator = Indicator_Factory::make_Improved_hMLP_BD2_Indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);
		auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
	}
	else if (ms::compare_icase(reconstruction_scheme, "Improved_hMLP_BD3"))
	{
		const auto criterion_equation_index = 0;
		const auto& governing_equation_name = configuration.get_governing_equation();

		auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
		auto indicator = Indicator_Factory::make_Improved_hMLP_BD3_Indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);
		auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
	}
	else if (ms::compare_icase(reconstruction_scheme, "Improved_hMLP_BD4"))
	{
		const auto criterion_equation_index = 0;
		const auto& governing_equation_name = configuration.get_governing_equation();

		auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
		auto indicator = Indicator_Factory::make_Improved_hMLP_BD4_Indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);
		auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
	}
	else
	{
		EXCEPTION("reconstruction shceme in configuration file is not supported");
		return nullptr;
	}
}