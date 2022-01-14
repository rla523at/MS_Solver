#pragma once
#include "Cell_Indicator_Factory.h"
#include "Limiter_Impl.h"
#include "Reconstruction.h"
#include "Stability_Criterion_Impl.h"

class Reconstruction_DG_Factory//static class
{
public://Query
	static std::unique_ptr<Reconstruction_DG> make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Governing_Equation>& governing_equation)
	{
		const auto& reconstruction_scheme = configuration.get_reconstruction_scheme();

		if (ms::compare_icase(reconstruction_scheme, "Not_Use"))
		{
			return std::make_unique<No_Reconstruction_DG>();
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP"))
		{
			const auto criterion_equation_index = 0;
			auto stability_criterion = std::make_unique<MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_Indicator(grid, discrete_solution, criterion_equation_index);
			auto limiter = std::make_unique <hMLP_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD"))
		{
			const auto criterion_equation_index = 0;

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_11_1_F"))
		{
			const auto criterion_equation_index = 0;
			const auto& governing_equation_name = configuration.get_governing_equation();

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_11_1_F_indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_11_1_1"))
		{
			const auto criterion_equation_index = 0;
			const auto& governing_equation_name = configuration.get_governing_equation();

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_11_1_1_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_11_1_21"))
		{
			const auto criterion_equation_index = 0;
			const auto& governing_equation_name = configuration.get_governing_equation();

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_11_1_21_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_11_1_22"))
		{
			const auto criterion_equation_index = 0;
			const auto& governing_equation_name = configuration.get_governing_equation();

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_11_1_22_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_11_1_3"))
		{
			const auto criterion_equation_index = 0;
			const auto& governing_equation_name = configuration.get_governing_equation();

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_11_1_3_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_11_1_41"))
		{
			const auto criterion_equation_index = 0;
			const auto& governing_equation_name = configuration.get_governing_equation();

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_11_1_41_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_11_1_42"))
		{
			const auto criterion_equation_index = 0;
			const auto& governing_equation_name = configuration.get_governing_equation();

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_11_1_42_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_12_1_T"))
		{
			const auto criterion_equation_index = 0;
			const auto& governing_equation_name = configuration.get_governing_equation();

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_12_1_T_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_12_1_42"))
		{
		const auto criterion_equation_index = 0;
		const auto& governing_equation_name = configuration.get_governing_equation();

		auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
		auto indicator = Cell_Indicator_Factory::make_hMLP_BD_12_1_42_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
		auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

		return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else if (ms::compare_icase(reconstruction_scheme, "hMLP_BD_II_22_1_T"))
		{
			const auto criterion_equation_index = 0;
			const auto& governing_equation_name = configuration.get_governing_equation();

			auto stability_criterion = std::make_unique<Simplex_Decomposed_MLP_Criterion>(grid, discrete_solution, criterion_equation_index);
			auto indicator = Cell_Indicator_Factory::make_hMLP_BD_22_1_T_Indicator(grid, discrete_solution, criterion_equation_index, governing_equation);
			auto limiter = std::make_unique <hMLP_BD_Limiter>(grid);

			return std::make_unique<Hierarchical_Limiting_DG>(std::move(stability_criterion), std::move(indicator), std::move(limiter));
		}
		else
		{
			EXCEPTION("reconstruction shceme in configuration file is not supported");
			return nullptr;
		}
	}

private:
    Reconstruction_DG_Factory(void) = delete;
};