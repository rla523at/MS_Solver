#include "../INC/Reconstruction_Impl.h"

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
		auto indicator = std::make_unique<hMLP_Indicator>(grid, discrete_solution, criterion_equation_index);
		auto limiter = std::make_unique <hMLP_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
	}
	else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD"))
	{
		const auto& governing_equation_name = configuration.get_governing_equation();
		const auto criterion_equation_index = 0;
		auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
		auto indicator = std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, governing_equation_name);
		auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));

		//Simplex_Decomposed_MLP_Criterion simplex_decomposed_MLP_criterion(grid, discrete_solution);
		//MLP_Indicator MLP_indicator(grid, discrete_solution);
		//MLP_u1_Limiter MLP_u1_limiter;
		//Subcell_Oscillation_Indicator subcell_oscillation_indicator(grid, discrete_solution);
		//Shock_Indicator shock_indicator(grid, discrete_solution, governing_equation_name);

		//return std::make_unique<hMLP_BD_Reconstruction_DG>(std::move(simplex_decomposed_MLP_criterion), std::move(MLP_indicator),
		//	std::move(MLP_u1_limiter), std::move(subcell_oscillation_indicator), std::move(shock_indicator), grid);

		//return std::make_unique<hMLP_BD_Reconstruction_DG>(grid, discrete_solution, governing_equation_name);
	}
	//else if (ms::compare_icase(reconstruction_scheme, "Test"))
	//{
	//	//Discontinuity_Indicator discontinuity_indicator(grid, discrete_solution, 0);
	//	//return std::make_unique<Test_Reconstuction_DG>(std::move(discontinuity_indicator));
	//}
	else
	{
		EXCEPTION("reconstruction shceme in configuration file is not supported");
		return nullptr;
	}
}