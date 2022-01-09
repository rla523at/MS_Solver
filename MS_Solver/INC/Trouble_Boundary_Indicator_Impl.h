#pragma once
#include "Indicator.h"
#include "Measuring_Function.h"

class Troubeld_Boudary_Indicator_with_Jump_Measurer_Base : public Trouble_Boundary_Indicator
{
public:
	Troubeld_Boudary_Indicator_with_Jump_Measurer_Base(const Grid& grid, std::unique_ptr<Face_Jump_Measurer>&& face_jump_measurer)
		: num_infcs_(grid.num_inner_faces())
		, infc_index_to_characteristic_length_table_(grid.inner_face_index_to_characteristic_length_table())
		, infc_index_to_oc_nc_index_pair_table_(grid.inner_face_index_to_oc_nc_index_pair_table())
		, face_jump_measurer_(std::move(face_jump_measurer))
	{
		this->cell_index_to_num_troubled_boundaries_table_.resize(grid.num_cells());
	}

public:
	void check(const Discrete_Solution_DG& discrete_solution)
	{
		std::fill(this->cell_index_to_num_troubled_boundaries_table_.begin(), this->cell_index_to_num_troubled_boundaries_table_.end(), 0);

		const auto infc_index_to_sol_jump_table = this->face_jump_measurer_->measure_inner_face_index_to_solution_jump_table(discrete_solution);

		for (uint infc_index = 0; infc_index < this->num_infcs_; ++infc_index)
		{
			const auto threshold_value = this->calculate_threshold_value(discrete_solution, infc_index);
			const auto sol_jump = infc_index_to_sol_jump_table[infc_index];
			
			if (threshold_value < sol_jump)
			{
				const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];

				this->cell_index_to_num_troubled_boundaries_table_[oc_index]++;
				this->cell_index_to_num_troubled_boundaries_table_[nc_index]++;
			}
		}
	}

private:
	virtual double calculate_threshold_value(const Discrete_Solution_DG& discrete_solution, const uint inner_face_index) const abstract;

protected:
	std::vector<double> infc_index_to_characteristic_length_table_;
	std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;

private:
	uint num_infcs_ = 0;
	std::unique_ptr<Face_Jump_Measurer> face_jump_measurer_;
};

// if h ^ (k+1/ 2) <= face jump ==> troubled boundary
// h := face characteristic legnth
class Trouble_Boundary_Indicator_Type1 : public Troubeld_Boudary_Indicator_with_Jump_Measurer_Base
{
public:
	Trouble_Boundary_Indicator_Type1(const Grid& grid, std::unique_ptr<Face_Jump_Measurer>&& face_jump_measurer)
		: Troubeld_Boudary_Indicator_with_Jump_Measurer_Base(grid, std::move(face_jump_measurer)) {};

private:
	double calculate_threshold_value(const Discrete_Solution_DG& discrete_solution, const uint inner_face_index) const override
	{
		const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[inner_face_index];

		const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
		const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
		const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

		const auto characteristic_length = this->infc_index_to_characteristic_length_table_[inner_face_index];
		return std::pow(characteristic_length, 0.5 * (max_solution_degree + 1));
	}
};

// if h <= face jump ==> troubled boundary
// h := face characteristic legnth
class Trouble_Boundary_Indicator_Type2 : public Troubeld_Boudary_Indicator_with_Jump_Measurer_Base
{
public:
	Trouble_Boundary_Indicator_Type2(const Grid& grid, std::unique_ptr<Face_Jump_Measurer>&& face_jump_measurer)
		: Troubeld_Boudary_Indicator_with_Jump_Measurer_Base(grid, std::move(face_jump_measurer)) {};

private:
	double calculate_threshold_value(const Discrete_Solution_DG& discrete_solution, const uint inner_face_index) const override
	{
		return this->infc_index_to_characteristic_length_table_[inner_face_index];
	}
};