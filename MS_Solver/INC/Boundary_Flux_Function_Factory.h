#pragma once
#include "Boundary_Flux_Function_Impl.h"
#include "Grid.h"

class Boundary_Flux_Function_Factory//static class
{
public:
	static std::vector<std::unique_ptr<Boundary_Flux_Function>> make_bdry_flux_functions(const Grid& grid, const std::shared_ptr<Governing_Equation>& governing_equation, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function)
	{
		const auto num_bdrys = grid.num_boundaries();

		std::vector<std::unique_ptr<Boundary_Flux_Function>> bdry_flux_functions(num_bdrys);

		for (uint i = 0; i < num_bdrys; ++i)
		{
			const auto bdry_type = grid.boundary_type(i);

			switch (bdry_type)
			{
			case ElementType::initial_constant_BC:			
				bdry_flux_functions[i] = make_initial_constant_BC(numerical_flux_function);
				break;
			case ElementType::slip_wall:
				bdry_flux_functions[i] = make_slip_wall_BC(governing_equation, numerical_flux_function);
				break;
			default:
				EXCEPTION("not supported boundary type");
				bdry_flux_functions[i] = nullptr;
			}
		}
		
		return bdry_flux_functions;
	}

private:
	static std::unique_ptr<Boundary_Flux_Function> make_initial_constant_BC(const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function)
	{
		return std::make_unique<Initial_Constant_BC>(numerical_flux_function);
	};

	static std::unique_ptr<Boundary_Flux_Function> make_slip_wall_BC(const std::shared_ptr<Governing_Equation>& governing_equation, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function)
	{
		return std::make_unique<Slip_Wall_BC>(governing_equation, numerical_flux_function);
	};

private:
	Boundary_Flux_Function_Factory(void) = delete;
};