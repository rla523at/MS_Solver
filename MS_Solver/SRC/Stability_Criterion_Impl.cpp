#include "../INC/Stability_Criterion_Impl.h"

MLP_Criterion::MLP_Criterion(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
	: MLP_Criterion_Base(grid, discrete_solution, criterion_equation_index)
{
	this->P0_values_.resize(this->num_cells_);

	//precalculate
	discrete_solution.precalculate_cell_P0_basis_values();
}

void MLP_Criterion::precaclulate(const Discrete_Solution_DG& discrete_solution)
{
	for (uint i = 0; i < this->num_cells_; ++i)
	{
		this->P0_values_[i] = discrete_solution.calculate_P0_nth_solution(i, this->criterion_equation_index_);
		std::fill(this->cell_index_to_P0_value_at_vertices_table_[i].begin(), this->cell_index_to_P0_value_at_vertices_table_[i].end(), this->P0_values_[i]);
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
		}

		const auto min_criterion_value = *std::min_element(start_iter, iter);
		const auto max_criterion_value = *std::max_element(start_iter, iter);

		this->vnode_index_to_allowable_min_max_criterion_value_[vnode_index] = std::make_pair(min_criterion_value, max_criterion_value);
	}

	this->make_set_of_allowable_min_max_criterion_values();
}

Simplex_Decomposed_MLP_Criterion::Simplex_Decomposed_MLP_Criterion(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
	: MLP_Criterion_Base(grid, discrete_solution, criterion_equation_index)
{
	Profiler::set_time_point();

	this->vertex_index_to_matched_vertex_index_set_ = grid.vertex_index_to_peridoic_matched_vertex_index_set();
	this->set_of_vertex_index_to_simplex_P0_value_.resize(this->num_cells_);

	for (uint i = 0; i < this->num_cells_; ++i)
	{
		std::map<uint, double> vertex_index_to_simplex_P0_value;

		const auto& vertex_indexes = this->set_of_vertex_indexes_[i];
		const auto num_vertices = vertex_indexes.size();

		for (ushort j = 0; j < num_vertices; ++j)
		{
			vertex_index_to_simplex_P0_value.emplace(vertex_indexes[j], 0.0);
		}

		this->set_of_vertex_index_to_simplex_P0_value_[i] = std::move(vertex_index_to_simplex_P0_value);
	}

	//precalcuate
	discrete_solution.precalculate_set_of_simplex_P0_P1_projection_basis_vertices_m(grid);

	LOG << std::left << std::setw(50) << "@ Simplex Decomposed MLP Criterion precalculation" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << LOG.print_;
}

void Simplex_Decomposed_MLP_Criterion::precaclulate(const Discrete_Solution_DG& discrete_solution)
{
	for (uint i = 0; i < this->num_cells_; ++i)
	{
		const auto simplex_P0_value_at_vertices = discrete_solution.calculate_simplex_P0_projected_nth_solution_at_vertices(i, this->criterion_equation_index_);

		auto& vertex_index_to_simplex_P0_value = this->set_of_vertex_index_to_simplex_P0_value_[i];
		const auto& vertex_indexes = this->set_of_vertex_indexes_[i];

		for (ushort j = 0; j < vertex_indexes.size(); ++j)
		{
			vertex_index_to_simplex_P0_value.at(vertex_indexes[j]) = simplex_P0_value_at_vertices[j];
			this->cell_index_to_P0_value_at_vertices_table_[i][j] = simplex_P0_value_at_vertices[j];
		}
	}

	const auto start_iter = this->criterion_values_.begin();

	for (const auto& [vertex_index, share_cell_index_set] : this->vnode_index_to_share_cell_index_set_)
	{
		const auto num_share_cell = share_cell_index_set.size();
		std::fill(start_iter, start_iter + num_share_cell, 0.0);

		auto iter = start_iter;
		for (const auto cell_index : share_cell_index_set)
		{
			auto& vertex_index_to_simplex_P0_value = this->set_of_vertex_index_to_simplex_P0_value_[cell_index];

			if (vertex_index_to_simplex_P0_value.contains(vertex_index))
			{
				*iter = vertex_index_to_simplex_P0_value.at(vertex_index);
			}
			else
			{
				const auto& matched_vertex_index_set = this->vertex_index_to_matched_vertex_index_set_.at(vertex_index);

				for (const auto& matched_vertex_index : matched_vertex_index_set)
				{
					if (vertex_index_to_simplex_P0_value.contains(matched_vertex_index))
					{
						*iter = vertex_index_to_simplex_P0_value.at(matched_vertex_index);
						break;
					}

					if (matched_vertex_index == *matched_vertex_index_set.rbegin())
					{
						EXCEPTION("can not find matched vertex");
					}
				}
			}

			++iter;
		}

		const auto min_criterion_value = *std::min_element(start_iter, iter);
		const auto max_criterion_value = *std::max_element(start_iter, iter);

		this->vnode_index_to_allowable_min_max_criterion_value_[vertex_index] = std::make_pair(min_criterion_value, max_criterion_value);
	}

	this->make_set_of_allowable_min_max_criterion_values();
}

