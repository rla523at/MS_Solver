#include "../INC/Stability_Criterion.h"

MLP_Criterion_Base::MLP_Criterion_Base(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
	: vnode_index_to_share_cell_index_set_(grid.get_vnode_index_to_share_cell_index_set_consider_pbdry())
	, criterion_equation_index_(criterion_equation_index)
{
	this->num_cells_ = grid.num_cells();
	this->set_of_vertex_indexes_ = grid.cell_set_of_vertex_indexes();

	this->cell_index_to_P0_value_at_vertices_table_.resize(this->num_cells_);
	this->cell_index_to_allowable_min_max_value_at_vertices_table_.resize(this->num_cells_);
	for (uint i = 0; i < this->num_cells_; ++i)
	{
		const auto num_vertices = this->set_of_vertex_indexes_[i].size();
		this->cell_index_to_P0_value_at_vertices_table_[i].resize(num_vertices);
		this->cell_index_to_allowable_min_max_value_at_vertices_table_[i].resize(num_vertices);
	}

	const auto num_vnode = this->vnode_index_to_share_cell_index_set_.size();
	this->vnode_index_to_allowable_min_max_criterion_value_.reserve(num_vnode);

	for (const auto& [vnode_index, share_cell_index_set] : this->vnode_index_to_share_cell_index_set_)
	{
		this->vnode_index_to_allowable_min_max_criterion_value_.emplace(vnode_index, std::pair<double, double>());
	}
}

ushort MLP_Criterion_Base::get_criterion_equation_index(void) const
{
	return this->criterion_equation_index_;
}

const std::vector<double>& MLP_Criterion_Base::get_P0_value_at_vertices(const uint cell_index) const
{
	return this->cell_index_to_P0_value_at_vertices_table_[cell_index];
}

const std::vector<std::pair<double, double>>& MLP_Criterion_Base::get_allowable_min_max_value_at_vertices(const uint cell_index) const
{
	return this->cell_index_to_allowable_min_max_value_at_vertices_table_[cell_index];
}

void MLP_Criterion_Base::make_set_of_allowable_min_max_criterion_values(void)
{
	for (uint i = 0; i < this->num_cells_; ++i)
	{
		const auto& vertex_indexes = this->set_of_vertex_indexes_[i];
		const auto num_vertex_indexes = vertex_indexes.size();

		for (ushort j = 0; j < num_vertex_indexes; ++j)
		{
			const auto vnode_index = vertex_indexes[j];
			this->cell_index_to_allowable_min_max_value_at_vertices_table_[i][j] = this->vnode_index_to_allowable_min_max_criterion_value_[vnode_index];
		}
	}
}

