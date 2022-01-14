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
	void check(const Discrete_Solution_DG& discrete_solution);

private:
	virtual double calculate_threshold_value(const Discrete_Solution_DG& discrete_solution, const uint inner_face_index) const abstract;

protected:
	std::vector<double> infc_index_to_characteristic_length_table_;
	std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;

private:
	uint num_infcs_ = 0;
	std::unique_ptr<Face_Jump_Measurer> face_jump_measurer_;
};

// if h ^ (k + 1 / 2) <= face jump ==> troubled boundary
// h := face characteristic legnth
class Trouble_Boundary_Indicator_Type1 : public Troubeld_Boudary_Indicator_with_Jump_Measurer_Base
{
public:
	Trouble_Boundary_Indicator_Type1(const Grid& grid, std::unique_ptr<Face_Jump_Measurer>&& face_jump_measurer)
		: Troubeld_Boudary_Indicator_with_Jump_Measurer_Base(grid, std::move(face_jump_measurer)) {};

private:
	double calculate_threshold_value(const Discrete_Solution_DG& discrete_solution, const uint inner_face_index) const override;
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