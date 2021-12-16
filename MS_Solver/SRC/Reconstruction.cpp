#include "../INC/Reconstruction.h"

hMLP_Reconstruction_DG::hMLP_Reconstruction_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution)
	:vnode_index_to_share_cell_index_set_(grid.get_vnode_index_to_share_cell_index_set_consider_pbdry())
{
	this->num_cells_ = static_cast<uint>(grid.num_cells());
	this->set_of_vnode_indexes_ = grid.cell_set_of_vertex_indexes();
	this->volumes_ = grid.cell_volumes();

	this->P0_criterion_values_.resize(this->num_cells_);
	this->set_of_P1_projected_criterion_values_at_verticies_.resize(this->num_cells_);//construction optimization할때 수정해야 됨
	
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

	//precalculate
	const auto set_of_verticies = grid.cell_set_of_verticies();
	discrete_solution.precalculate_cell_vertices_basis_values(set_of_verticies);

}

void hMLP_Reconstruction_DG::reconstruct(Discrete_Solution_DG& discrete_solution)
{
	this->precalculate(discrete_solution);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		auto solution_degree = discrete_solution.solution_degree(cell_index);

		const auto P0_criterion_value = this->P0_criterion_values_[cell_index];
		const auto& allowable_min_max_criterion_values = this->set_of_allowable_min_max_criterion_values_[cell_index];

		while (true)
		{
			if (this->is_trouble(discrete_solution, cell_index))
			{
				if (solution_degree <= 2)
				{
					const auto limiting_value = MLP_u1_Limiter::limiter_function(P0_criterion_value, this->set_of_P1_projected_criterion_values_at_verticies_[cell_index], allowable_min_max_criterion_values);

					solution_degree = 1;
					discrete_solution.project_to_Pn_space(cell_index, solution_degree);
					discrete_solution.limit_slope(cell_index, limiting_value);
					break;
				}
				else
				{
					discrete_solution.project_to_Pn_space(cell_index, --solution_degree);
				}
			}
			else
			{
				break;
			}
		}
	}
};

bool hMLP_Reconstruction_DG::is_trouble(const Discrete_Solution_DG& discrete_solution, const uint cell_index) const
{
	const auto  P0_criterion_value = this->P0_criterion_values_[cell_index];
	const auto& P1_projected_criterion_values_at_verticies = this->set_of_P1_projected_criterion_values_at_verticies_[cell_index];
	const auto  criterion_solution_at_vertices = discrete_solution.calculate_nth_solution_at_vertices(cell_index, this->criterion_variable_index_);

	const auto& allowable_min_max_criterion_values = this->set_of_allowable_min_max_criterion_values_[cell_index];
	const auto  volume = this->volumes_[cell_index];

	const auto num_vertices = allowable_min_max_criterion_values.size();

	for (ushort j = 0; j < num_vertices; ++j)
	{
		const auto [allowable_min, allowable_max] = allowable_min_max_criterion_values[j];

		const auto criterion_value = criterion_solution_at_vertices[j];
		const auto P1_projected_criterion_value = P1_projected_criterion_values_at_verticies[j];

		const auto higher_mode_criterion_value = criterion_value - P1_projected_criterion_value;
		const auto P1_mode_criterion_value = P1_projected_criterion_value - P0_criterion_value;

		if (!Constant_Region_Vertex_Indicator::is_constant(criterion_value, P0_criterion_value, volume) &&
			!P1_Projected_MLP_Trouble_Vertex_Indicator::is_satisfy(P1_projected_criterion_value, allowable_min, allowable_max) &&
			!MLP_Smooth_Extrema_Vertex_Indicator::is_smooth_extrema(criterion_value, higher_mode_criterion_value, P1_mode_criterion_value, allowable_min, allowable_max))
		{
			return true;
		}
	}

	return false;
}

void hMLP_Reconstruction_DG::precalculate(const Discrete_Solution_DG& discrete_solution)
{
	for (uint i = 0; i < this->num_cells_; ++i)
	{
		this->P0_criterion_values_[i] = discrete_solution.calculate_P0_nth_solution(i, this->criterion_variable_index_);
		this->set_of_P1_projected_criterion_values_at_verticies_[i] = discrete_solution.calculate_P1_projected_nth_solution_at_vertices(i, this->criterion_variable_index_);
	}

	for (const auto& [vnode_index, share_cell_index_set] : this->vnode_index_to_share_cell_index_set_)
	{
		const auto num_share_cell = share_cell_index_set.size();
		std::vector<double> criterion_variables;
		criterion_variables.reserve(num_share_cell);

		for (const auto cell_index : share_cell_index_set)
		{
			criterion_variables.push_back(this->P0_criterion_values_[cell_index]);
		}

		const auto min_criterion_value = *std::min_element(criterion_variables.begin(), criterion_variables.end());
		const auto max_criterion_value = *std::max_element(criterion_variables.begin(), criterion_variables.end());

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
