#pragma once
#include "Measuring_Function.h"

// INT_{w_j} (q - q_j) / ||w_j|| ==> average solution jump
class Face_Jump_Measurer_Type1 : public Face_Jump_Measurer
{
public:
	Face_Jump_Measurer_Type1(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
		:Face_Jump_Measurer(grid, discrete_solution, criterion_solution_index) {};

private:
	double calculate_scail_factor(const Discrete_Solution_DG& discrete_solution, const uint inner_face_index) const override { return 1.0; };
};

// INT_{w_j} (q - q_j) / ( 0.5 * (q + q_j) * ||w_j|| ) ==> scaled average solution jump
class Face_Jump_Measurer_Type2 : public Face_Jump_Measurer
{
public:
	Face_Jump_Measurer_Type2(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
		:Face_Jump_Measurer(grid, discrete_solution, criterion_solution_index) 
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
