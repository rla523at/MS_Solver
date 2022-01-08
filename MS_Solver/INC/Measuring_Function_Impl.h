#pragma once
#include "Measuring_Function.h"

class Default_Average_Solution_Jump_Measurer : public Average_Solution_Jump_Measurer
{
public:
	Default_Average_Solution_Jump_Measurer(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
		:Average_Solution_Jump_Measurer(grid, discrete_solution, criterion_solution_index) {};

private:
	double calculate_scail_factor(const Discrete_Solution_DG& discrete_solution, const uint inner_face_index) const override { return 1.0; };
};

class Scaled_Average_Solution_Jump_Measurer : public Average_Solution_Jump_Measurer
{
public:
	Scaled_Average_Solution_Jump_Measurer(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
		:Average_Solution_Jump_Measurer(grid, discrete_solution, criterion_solution_index) 
	{
		//precalculation
		discrete_solution.precalculate_cell_P0_basis_values();
	};

private:
	double calculate_scail_factor(const Discrete_Solution_DG& discrete_solution, const uint inner_face_index) const override
	{  
		const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[inner_face_index];

		const auto oc_avg_sol = discrete_solution.calculate_P0_nth_solution(oc_index, this->criterion_solution_index_);
		const auto nc_avg_sol = discrete_solution.calculate_P0_nth_solution(nc_index, this->criterion_solution_index_);
		
		return 0.5 * (oc_avg_sol + nc_avg_sol);
	};
};
