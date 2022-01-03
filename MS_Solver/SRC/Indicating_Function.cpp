#include "../INC/Indicating_Function.h"

MLP_Indicator::MLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_index)
	:criterion_equation_index_(criterion_index)
{
	this->num_cells_ = grid.num_cells();
	this->cell_index_to_volume_table_ = grid.cell_index_to_volume_table();
	this->cell_index_to_num_vertices_table_ = grid.cell_index_to_num_vertices_table();
	this->cell_index_to_P1_projected_value_at_vertices_table_.resize(this->num_cells_);

	const auto cell_index_to_num_vertices_table = grid.cell_index_to_num_vertices_table();
	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{		
		this->cell_index_to_P1_projected_value_at_vertices_table_[cell_index].resize(cell_index_to_num_vertices_table[cell_index]);
	}

	//precalculate
	const auto set_of_verticies = grid.cell_set_of_verticies();
	discrete_solution.precalculate_cell_vertices_basis_values(set_of_verticies);
}

void MLP_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution)
{
	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		discrete_solution.calculate_P1_projected_nth_solution_at_vertices(this->cell_index_to_P1_projected_value_at_vertices_table_[cell_index].data(), cell_index, this->criterion_equation_index_);
	}
}

//void MLP_Indicator::set_precalculated_result(const std::vector<double>* set_of_P1_projected_value_at_vertices_ptr)
//{
//	this->set_of_P1_projected_value_at_vertices_ptr_ = set_of_P1_projected_value_at_vertices_ptr;
//}
//
//cell_type MLP_Indicator::indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion& criterion) const
//{
//	discrete_solution.calculate_P1_projected_nth_solution_at_vertices(this->P1_projected_value_at_vertices_.data(), cell_index, this->criterion_equation_index_);
//	discrete_solution.calculate_nth_solution_at_vertices(this->value_at_vertices_.data(), cell_index, this->criterion_equation_index_);
//
//	return this->check_cell_type(cell_index, this->P1_projected_value_at_vertices_.data(), this->value_at_vertices_.data(), criterion);
//}
//
//cell_type MLP_Indicator::indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const Simplex_Decomposed_MLP_Criterion& criterion) const
//{
//	const auto criterion_equation_index = criterion.get_criterion_equation_index();
//
//	discrete_solution.calculate_simplex_P1_projected_nth_solution_at_vertices(this->P1_projected_value_at_vertices_.data(), cell_index, criterion_equation_index);
//	discrete_solution.calculate_nth_solution_at_vertices(this->value_at_vertices_.data(), cell_index, criterion_equation_index);
//
//	return this->check_cell_type(cell_index, this->P1_projected_value_at_vertices_.data(), this->value_at_vertices_.data(), criterion);
//}

cell_type MLP_Indicator::indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& criterion) const
{
	const auto criterion_equation_index = criterion.get_criterion_equation_index();

	discrete_solution.calculate_nth_solution_at_vertices(this->value_at_vertices_.data(), cell_index, criterion_equation_index);

	return this->check_cell_type(cell_index, this->value_at_vertices_.data(), criterion);
}

cell_type MLP_Indicator::check_cell_type(const uint cell_index, const double* value_at_vertices, const MLP_Criterion_Base& criterion) const
{
	const auto& P1_projected_value_at_vertices = this->cell_index_to_P1_projected_value_at_vertices_table_[cell_index]; 
	const auto& P0_value_at_vertices = criterion.get_P0_value_at_vertices(cell_index);
	const auto& allowable_min_maxs = criterion.get_allowable_min_max_value_at_vertices(cell_index);

	bool has_smooth_extrma = false;

	for (ushort j = 0; j < this->cell_index_to_num_vertices_table_[cell_index]; ++j)
	{
		const auto P0_value = P0_value_at_vertices[j];
		const auto [allowable_min, allowable_max] = allowable_min_maxs[j];

		const auto value = value_at_vertices[j];
		const auto P1_projected_value = P1_projected_value_at_vertices[j];

		const auto higher_mode_value = value - P1_projected_value;
		const auto P1_mode_value = P1_projected_value - P0_value;

		if (!this->is_constant(value, P0_value, cell_index) && !this->is_satisfy_MLP_condition(P1_projected_value, allowable_min, allowable_max))
		{
			if (this->is_smooth_extrema(value, higher_mode_value, P1_mode_value, allowable_min, allowable_max))
			{
				has_smooth_extrma = true;
			}
			else
			{
				return cell_type::trouble;
			}
		}
	}

	if (has_smooth_extrma)
	{
		return cell_type::smooth_extrema;
	}
	else
	{
		return cell_type::normal;
	}
}

bool MLP_Indicator::is_constant(const double value, const double P0_value, const uint cell_index) const
{
	const auto constant_criterion = (std::max)(1.0E-3 * std::abs(P0_value), this->cell_index_to_volume_table_[cell_index]);
	return std::abs(value - P0_value) <= constant_criterion;
}

bool MLP_Indicator::is_satisfy_MLP_condition(const double P1_projected_value, const double allowable_min, const double allowable_max) const
{
	return allowable_min <= P1_projected_value && P1_projected_value <= allowable_max;
}

bool MLP_Indicator::is_smooth_extrema(const double value, const double higher_mode_value, const double P1_mode_value, const double allowable_min, const double allowable_max) const
{
	if (P1_mode_value > 0 && higher_mode_value < 0 && value > allowable_min)
	{
		return true;
	}
	else if (P1_mode_value < 0 && higher_mode_value > 0 && value < allowable_max)
	{
		return true;
	}
	else
	{
		return false;
	}
}

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

		if (smooth_boundary_indicator > threshold_value)
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

Shock_Indicator::Shock_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::string& governing_equation_name)
{
	this->set_criterion_solution_index(grid.space_dimension(), governing_equation_name);

	const auto num_cells = grid.num_cells();
	this->average_pressures_.resize(num_cells);
	this->cell_index_to_is_shock_.resize(num_cells, false);
	this->cell_index_to_face_share_cell_indexes_table_ = grid.cell_index_to_face_share_cell_indexes_table_consider_pbdry();
}

void Shock_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution)
{
	if (this->criterion_solution_index_ < 0)
	{
		return; //it means that governing equation is linear advection
	}

	std::fill(this->cell_index_to_is_shock_.begin(), this->cell_index_to_is_shock_.end(), false);

	const auto P0_solutions = discrete_solution.calculate_P0_solutions();

	const auto num_cells = this->cell_index_to_is_shock_.size();
	for (uint i = 0; i < num_cells; ++i)
	{
		const auto& face_share_cell_indexes = this->cell_index_to_face_share_cell_indexes_table_[i];

		const auto my_value = P0_solutions[i][this->criterion_solution_index_];
		for (const auto face_share_cell_index : face_share_cell_indexes)
		{
			const auto other_value = P0_solutions[face_share_cell_index][this->criterion_solution_index_];

			if (std::abs(my_value - other_value) >= 0.1 * (std::min)(my_value, other_value))
			{
				this->cell_index_to_is_shock_[i] = true;
				break;
			}
		}
	}
}

bool Shock_Indicator::is_shock(const uint cell_index) const
{
	return this->cell_index_to_is_shock_[cell_index];
}

void Shock_Indicator::set_criterion_solution_index(const ushort space_dimension, const std::string& governing_equation_name)
{
	if (ms::compare_icase(governing_equation_name, "Linear_Advection"))
	{
		this->criterion_solution_index_ = -1;
	}
	else if (ms::compare_icase(governing_equation_name, "Burgers"))
	{
		this->criterion_solution_index_ = 0;
	}
	else if (ms::compare_icase(governing_equation_name, "Euler"))
	{
		if (space_dimension == 2)
		{
			this->criterion_solution_index_ = Euler_2D::pressure_index();
		}
		else if (space_dimension == 3)
		{
			this->criterion_solution_index_ = Euler_3D::pressure_index();
		}
	}
}

Discontinuity_Indicator::Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
	: criterion_solution_index_(criterion_solution_index)
{
	this->cell_index_to_face_share_cell_indexes_table_ = grid.cell_index_to_face_share_cell_indexes_table_consider_pbdry();
	this->cell_index_to_volume_table_ = grid.cell_index_to_volume_table();

	this->num_cells_ = grid.num_cells();
	this->cell_index_to_has_discontinuity_table_.resize(this->num_cells_, false);
	this->cell_index_to_QW_v_table_.resize(this->num_cells_);


	const auto cell_index_to_face_share_cell_indexes_table_ignore_pbdry = grid.cell_index_to_face_share_cell_indexes_table_ignore_pbdry();

	std::vector<std::vector<Euclidean_Vector>> cell_index_to_order_n_QPs_table(this->num_cells_);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto solution_degree = discrete_solution.solution_degree(cell_index);
		const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, solution_degree);

		cell_index_to_order_n_QPs_table[cell_index] = quadrature_rule.points;
		this->cell_index_to_QW_v_table_[cell_index] = quadrature_rule.weights;
	}


	const auto pbdry_pair_index_to_oc_nc_index_pair_table = grid.pbdry_pair_index_to_oc_nc_index_pair_table();
	const auto pbdry_pair_index_to_ocs_to_ncs_v_table = grid.pbdry_pair_index_to_ocs_to_ncs_v_table();
	const auto num_pbdry_pair = pbdry_pair_index_to_oc_nc_index_pair_table.size();

	std::vector<std::pair<std::vector<Euclidean_Vector>, std::vector<Euclidean_Vector>>> pbdry_pair_index_to_translated_oc_nc_order_n_QPs_pair_table(num_pbdry_pair);

	for (uint pbdry_pair_index = 0; pbdry_pair_index < num_pbdry_pair; ++pbdry_pair_index)
	{
		const auto [oc_index, nc_index] = pbdry_pair_index_to_oc_nc_index_pair_table[pbdry_pair_index];
		const auto& ocs_to_ncs_v = pbdry_pair_index_to_ocs_to_ncs_v_table.at(pbdry_pair_index);

		const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
		const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);

		const auto& oc_quadrature_rule = grid.get_cell_quadrature_rule(oc_index, oc_solution_degree);
		const auto& nc_quadrature_rule = grid.get_cell_quadrature_rule(nc_index, nc_solution_degree);

		auto translated_oc_QPs = oc_quadrature_rule.points;
		for (auto& QP : translated_oc_QPs)
		{
			QP += ocs_to_ncs_v;
		}

		auto translated_nc_QPs = nc_quadrature_rule.points;
		for (auto& QP : translated_nc_QPs)
		{
			QP -= ocs_to_ncs_v;
		}

		pbdry_pair_index_to_translated_oc_nc_order_n_QPs_pair_table[pbdry_pair_index] = { std::move(translated_oc_QPs),std::move(translated_nc_QPs) };
	}


	//for precalculation
	discrete_solution.precalculate_set_of_cell_index_to_target_cell_basis_QPs_m_(
		cell_index_to_order_n_QPs_table,
		cell_index_to_face_share_cell_indexes_table_ignore_pbdry,
		pbdry_pair_index_to_oc_nc_index_pair_table,
		pbdry_pair_index_to_translated_oc_nc_order_n_QPs_pair_table);
	//

	//test
	this->discontinuity_factor_.resize(this->num_cells_);
	//
}

void Discontinuity_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution)
{
	static constexpr ushort max_num_QPs = 100;
	std::array<double, max_num_QPs> nth_solution_at_target_cell_QPs = { 0 };

	//test
	std::fill(this->discontinuity_factor_.begin(), this->discontinuity_factor_.end(), 0.0);
	//

	for (uint i = 0; i < this->num_cells_; ++i)
	{
		double discontinuity_indicator_value = 0.0;

		const auto cell_volume = this->cell_index_to_volume_table_[i];
		const auto P0_value = discrete_solution.calculate_P0_nth_solution(i, this->criterion_solution_index_);

		const auto& QWs_v = this->cell_index_to_QW_v_table_[i];
		const auto num_QPs = static_cast<int>(QWs_v.size());

		const auto& face_share_cell_indexes = this->cell_index_to_face_share_cell_indexes_table_[i];
		const auto num_face_share_cells = face_share_cell_indexes.size();


		for (ushort j = 0; j < num_face_share_cells; ++j)
		{
			const auto my_cell_index = face_share_cell_indexes[j];

			discrete_solution.calculate_nth_solution_at_target_cell_QPs(nth_solution_at_target_cell_QPs.data(), i, my_cell_index, this->criterion_solution_index_);

			const auto P0_value_by_extrapolate = ms::BLAS::x_dot_y(num_QPs, nth_solution_at_target_cell_QPs.data(), QWs_v.data());

			//this->discontinuity_factor_[i] += std::abs(P0_value_by_extrapolate - P0_value * cell_volume);
			//this->discontinuity_factor_[i] += std::abs(P0_value_by_extrapolate / cell_volume - P0_value);
			discontinuity_indicator_value += std::abs(P0_value_by_extrapolate / cell_volume - P0_value);
		}

		//this->discontinuity_factor_[i] /= num_face_share_cells;
		discontinuity_indicator_value /= num_face_share_cells;
		const auto characteristic_length = std::sqrt(this->cell_index_to_volume_table_[i]);

		this->discontinuity_factor_[i] = std::log(discontinuity_indicator_value) / std::log(characteristic_length);
	}
}