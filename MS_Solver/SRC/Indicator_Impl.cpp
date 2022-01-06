#include "../INC/Cell_Indicator_Impl.h"

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

Cell_Type MLP_Indicator::indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& criterion) const
{
	const auto criterion_equation_index = criterion.get_criterion_equation_index();

	discrete_solution.calculate_nth_solution_at_vertices(this->value_at_vertices_.data(), cell_index, criterion_equation_index);

	return this->check_cell_type(cell_index, this->value_at_vertices_.data(), criterion);
}

Cell_Type MLP_Indicator::check_cell_type(const uint cell_index, const double* value_at_vertices, const MLP_Criterion_Base& criterion) const
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
				return Cell_Type::trouble;
			}
		}
	}

	if (has_smooth_extrma)
	{
		return Cell_Type::smooth_extrema;
	}
	else
	{
		return Cell_Type::normal;
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

Cell_Type hMLP_BD_Indicator::indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const
{
    auto cell_type = this->MLP_indicator_.indicate(discrete_solution, cell_index, stability_criterion);

    switch (cell_type)
    {
    case Cell_Type::normal:
    {
        if (this->subcell_oscillation_indicator_.is_typeI_cell(cell_index) && this->shock_indicator_->near_shock(cell_index))
        {
            cell_type = Cell_Type::typeI;
        }
        break;
    }
    case Cell_Type::smooth_extrema:
    {
        if (this->subcell_oscillation_indicator_.is_typeII_cell(cell_index) && this->discontinuity_indicator_->near_discontinuity(cell_index))
        {
            cell_type = Cell_Type::typeII;
        }
        break;
    }
    default:
        break;
    }

    return cell_type;
}