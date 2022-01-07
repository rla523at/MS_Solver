#include "../INC/Indicator.h"

Subcell_Oscillation_Indicator::Subcell_Oscillation_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
	:criterion_equation_index_(criterion_equation_index)
{
	const auto space_dimension = grid.space_dimension();

	this->set_of_num_troubled_boundaries_.resize(grid.num_cells());

	const auto num_inner_faces = grid.num_inner_faces();

	this->inner_face_volumes_.resize(num_inner_faces);
	this->inner_face_oc_nc_index_pairs_.resize(num_inner_faces);
	this->inner_face_characteristic_lengths_.resize(num_inner_faces);
	this->set_of_jump_QWs_v_.resize(num_inner_faces);

	std::vector<uint> oc_indexes(num_inner_faces);
	std::vector<uint> nc_indexes(num_inner_faces);
	std::vector<std::vector<Euclidean_Vector>> set_of_ocs_QPs(num_inner_faces);
	std::vector<std::vector<Euclidean_Vector>> set_of_ncs_QPs(num_inner_faces);

	for (uint i = 0; i < num_inner_faces; ++i)
	{
		this->inner_face_volumes_[i] = grid.inner_face_volume(i);
		this->inner_face_characteristic_lengths_[i] = std::pow(this->inner_face_volumes_[i], 1.0 / (space_dimension - 1));
		this->inner_face_oc_nc_index_pairs_[i] = grid.inner_face_oc_nc_index_pair(i);

		const auto [oc_index, nc_index] = this->inner_face_oc_nc_index_pairs_[i];
		const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
		const auto nc_solution_degree = discrete_solution.solution_degree(oc_index);
		const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

		const auto& [ocs_quadrature_rule, ncs_quadrature_rule] = grid.inner_face_quadrature_rule(i, max_solution_degree);
		this->set_of_jump_QWs_v_[i] = ocs_quadrature_rule.weights;

		//for precalculation
		oc_indexes[i] = oc_index;
		nc_indexes[i] = nc_index;
		set_of_ocs_QPs[i] = ocs_quadrature_rule.points;
		set_of_ncs_QPs[i] = ncs_quadrature_rule.points;
	}

	//precalculation
	discrete_solution.precalculate_infs_ocs_jump_QPs_basis_values(oc_indexes, set_of_ocs_QPs);
	discrete_solution.precalculate_infs_ncs_jump_QPs_basis_values(nc_indexes, set_of_ncs_QPs);
}

void Subcell_Oscillation_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution)
{
	std::fill(this->set_of_num_troubled_boundaries_.begin(), this->set_of_num_troubled_boundaries_.end(), 0);

	const auto num_infcs = this->inner_face_characteristic_lengths_.size();

	for (uint infc_index = 0; infc_index < num_infcs; ++infc_index)
	{
		const auto [oc_index, nc_index] = this->inner_face_oc_nc_index_pairs_[infc_index];

		//calculate smooth_boundary_indicator
		discrete_solution.calculate_nth_solution_at_infc_ocs_jump_QPs(this->value_at_ocs_jump_QPs_.data(), infc_index, oc_index, this->criterion_equation_index_);
		discrete_solution.calculate_nth_solution_at_infc_ncs_jump_QPs(this->value_at_ncs_jump_QPs_.data(), infc_index, nc_index, this->criterion_equation_index_);

		const auto num_QPs = static_cast<int>(this->set_of_jump_QWs_v_[infc_index].size());
		ms::BLAS::x_minus_y(num_QPs, this->value_at_ocs_jump_QPs_.data(), this->value_at_ncs_jump_QPs_.data(), this->value_diff_at_jump_QPs_.data());
		ms::BLAS::abs_x(num_QPs, this->value_diff_at_jump_QPs_.data());

		const auto jump = ms::BLAS::x_dot_y(num_QPs, this->value_diff_at_jump_QPs_.data(), this->set_of_jump_QWs_v_[infc_index].data());
		const auto smooth_boundary_indicator = jump / this->inner_face_volumes_[infc_index];

		//calculate threshold
		const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
		const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
		const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

		const auto characteristic_length = this->inner_face_characteristic_lengths_[infc_index];
		const auto threshold_value = std::pow(characteristic_length, 0.5 * (max_solution_degree + 1));
		//const auto threshold_value = characteristic_length;

		if (threshold_value < smooth_boundary_indicator)
		{
			this->set_of_num_troubled_boundaries_[oc_index]++;
			this->set_of_num_troubled_boundaries_[nc_index]++;
		}
	}
}

bool Subcell_Oscillation_Indicator::is_typeI_cell(const uint cell_index) const
{
	constexpr ushort typeI_threshold_number = 2;
	return typeI_threshold_number <= this->set_of_num_troubled_boundaries_[cell_index];
}

bool Subcell_Oscillation_Indicator::is_typeII_cell(const uint cell_index) const
{
	constexpr ushort typeII_threshold_number = 1;
	return typeII_threshold_number <= this->set_of_num_troubled_boundaries_[cell_index];
}