#include "../INC/Indicator.h"

MLP_Criterion::MLP_Criterion(const Grid& grid, Discrete_Solution_DG& discrete_solution)
	: vnode_index_to_share_cell_index_set_(grid.get_vnode_index_to_share_cell_index_set_consider_pbdry())
{
	this->num_cells_ = static_cast<uint>(grid.num_cells());
	this->set_of_vnode_indexes_ = grid.cell_set_of_vertex_indexes();

	this->P0_values_.resize(this->num_cells_);

	this->set_of_allowable_min_max_criterion_values_.resize(this->num_cells_);
	for (uint i = 0; i < this->num_cells_; ++i)
	{
		const auto num_vertices = this->set_of_vnode_indexes_[i].size();
		this->set_of_allowable_min_max_criterion_values_[i].resize(num_vertices);
	}

	const auto num_vnode = this->vnode_index_to_share_cell_index_set_.size();
	this->vnode_index_to_allowable_min_max_criterion_value_.reserve(num_vnode);
	for (const auto& [vnode_index, share_cell_index_set] : this->vnode_index_to_share_cell_index_set_)
	{
		this->vnode_index_to_allowable_min_max_criterion_value_.emplace(vnode_index, std::pair<double, double>());
	}
}

void MLP_Criterion::caclulate(const Discrete_Solution_DG& discrete_solution)
{
	for (uint i = 0; i < this->num_cells_; ++i)
	{
		this->P0_values_[i] = discrete_solution.calculate_P0_nth_solution(i, this->criterion_equation_index_);
	}

	const auto start_iter = this->criterion_values_.begin();

	for (const auto& [vnode_index, share_cell_index_set] : this->vnode_index_to_share_cell_index_set_)
	{
		const auto num_share_cell = share_cell_index_set.size();
		std::fill(start_iter, start_iter + num_share_cell, 0.0);

		auto iter = start_iter;
		for (const auto cell_index : share_cell_index_set)
		{
			*iter = this->P0_values_[cell_index];
			++iter;

			//criterion_variables.push_back(this->P0_values_[cell_index]);
		}

		const auto min_criterion_value = *std::min_element(start_iter, iter);
		const auto max_criterion_value = *std::max_element(start_iter, iter);

		this->vnode_index_to_allowable_min_max_criterion_value_[vnode_index] = std::make_pair(min_criterion_value, max_criterion_value);
	}

	for (uint i = 0; i < this->num_cells_; ++i)
	{
		const auto& vnode_indexes = this->set_of_vnode_indexes_[i];
		const auto num_vnode_indexes = vnode_indexes.size();

		for (ushort j = 0; j < num_vnode_indexes; ++j)
		{
			const auto vnode_index = vnode_indexes[j];
			this->set_of_allowable_min_max_criterion_values_[i][j] = this->vnode_index_to_allowable_min_max_criterion_value_[vnode_index];
		}
	}
}

const std::vector<std::pair<double, double>>& MLP_Criterion::get_criterion_values(const uint cell_index) const
{
	return this->set_of_allowable_min_max_criterion_values_[cell_index];
}

double MLP_Criterion::get_P0_value(const uint cell_index) const
{
	return this->P0_values_[cell_index];
}

ushort MLP_Criterion::get_criterion_equation_index(void) const
{
	return this->criterion_equation_index_;
}

MLP_Indicator::MLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
{
	this->volumes_ = grid.cell_volumes();
	this->set_of_num_vertices_ = grid.cell_set_of_num_vertices();

	const auto set_of_verticies = grid.cell_set_of_verticies();
	discrete_solution.precalculate_cell_vertices_basis_values(set_of_verticies);
}


void MLP_Indicator::set_precalculated_result(const std::vector<double>* set_of_P1_projected_value_at_vertices_ptr)
{
	this->set_of_P1_projected_value_at_vertices_ptr_ = set_of_P1_projected_value_at_vertices_ptr;
}

cell_type MLP_Indicator::indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion& criterion) const
{
	const auto criterion_equation_index = criterion.get_criterion_equation_index();

	discrete_solution.calculate_P1_projected_nth_solution_at_vertices(this->P1_projected_value_at_vertices_.data(), cell_index, criterion_equation_index);
	discrete_solution.calculate_nth_solution_at_vertices(this->value_at_vertices_.data(), cell_index, criterion_equation_index);

	return this->check_cell_type(cell_index, this->P1_projected_value_at_vertices_.data(), this->value_at_vertices_.data(), criterion);


	//const auto P1_projected_value_at_vertices = discrete_solution.calculate_P1_projected_nth_solution_at_vertices(cell_index, criterion_equation_index);
	//const auto value_at_vertices = discrete_solution.calculate_nth_solution_at_vertices(cell_index, criterion_equation_index);

	//return this->check_cell_type(cell_index, P1_projected_value_at_vertices.data(), value_at_vertices.data(), criterion);
}

cell_type MLP_Indicator::indicate_opt(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion& criterion) const
{	
	const auto criterion_equation_index = criterion.get_criterion_equation_index();

	discrete_solution.calculate_nth_solution_at_vertices(this->value_at_vertices_.data(), cell_index, criterion_equation_index);

	return this->check_cell_type(cell_index, this->set_of_P1_projected_value_at_vertices_ptr_->data(), this->value_at_vertices_.data(), criterion);


	//const auto value_at_vertices = discrete_solution.calculate_nth_solution_at_vertices(cell_index, criterion_equation_index);

	//return this->check_cell_type(cell_index, this->set_of_P1_projected_value_at_vertices_ptr_->data(), value_at_vertices.data(), criterion);
}

cell_type MLP_Indicator::check_cell_type(const uint cell_index, const double* P1_projected_value_at_vertices, const double* value_at_vertices, const MLP_Criterion& criterion) const
{
	const auto P0_value = criterion.get_P0_value(cell_index);
	const auto& allowable_min_maxs = criterion.get_criterion_values(cell_index);

	bool has_smooth_extrma = false;

	for (ushort j = 0; j < this->set_of_num_vertices_[cell_index]; ++j)
	{
		const auto [allowable_min, allowable_max] = allowable_min_maxs[j];

		const auto criterion_value = value_at_vertices[j];
		const auto P1_projected_criterion_value = P1_projected_value_at_vertices[j];

		const auto higher_mode_criterion_value = criterion_value - P1_projected_criterion_value;
		const auto P1_mode_criterion_value = P1_projected_criterion_value - P0_value;

		if (!this->is_constant(criterion_value, P0_value, cell_index) && !this->is_satisfy_MLP_condition(P1_projected_criterion_value, allowable_min, allowable_max))
		{
			if (this->is_smooth_extrema(criterion_value, higher_mode_criterion_value, P1_mode_criterion_value, allowable_min, allowable_max))
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


bool MLP_Indicator::is_constant(const double solution, const double P0_solution, const uint cell_index) const
{
	const auto constant_criterion = (std::max)(1.0E-3 * std::abs(P0_solution), this->volumes_[cell_index]);
	return std::abs(solution - P0_solution) <= constant_criterion;
}

bool MLP_Indicator::is_satisfy_MLP_condition(const double P1_projected_value, const double allowable_min, const double allowable_max) const
{
	return allowable_min <= P1_projected_value && P1_projected_value <= allowable_max;
}

bool MLP_Indicator::is_smooth_extrema(const double solution, const double higher_mode_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) const
{
	if (P1_mode_solution > 0 && higher_mode_solution < 0 && solution > allowable_min)
	{
		return true;
	}
	else if (P1_mode_solution < 0 && higher_mode_solution > 0 && solution < allowable_max)
	{
		return true;
	}
	else
	{
		return false;
	}
}