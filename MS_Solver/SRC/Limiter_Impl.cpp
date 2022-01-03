#include "../INC/Limiter_Impl.h"

void hMLP_Type_Limiter::apply_MLP_u1(const ushort cell_index, Discrete_Solution_DG& discrete_solution, const MLP_Criterion_Base& stability_criterion) const
{
	const auto projection_degree = 1;
	discrete_solution.project_to_Pn_space(cell_index, projection_degree);

	const auto limiting_value = this->limiting_function_.limiter_function(discrete_solution, cell_index, stability_criterion);
	discrete_solution.limit_slope(cell_index, limiting_value);
}

void hMLP_Limiter::limit(const ushort cell_index, const Cell_Type Cell_Type, Discrete_Solution_DG& discrete_solution, ushort& projection_degree, const MLP_Criterion_Base& stability_criterion) const
{
	switch (Cell_Type)
	{
	case Cell_Type::normal:
	case Cell_Type::smooth_extrema:
	{
		this->is_end_ = true;
		break;
	}
	case Cell_Type::trouble:
	{
		if (projection_degree <= 2)
		{
			this->apply_MLP_u1(cell_index, discrete_solution, stability_criterion);
			this->is_end_ = true;
		}
		else
		{
			discrete_solution.project_to_Pn_space(cell_index, --projection_degree);
			this->is_end_ = false;
		}

		break;
	}
	default:
		EXCEPTION("not supported cell type");
	}
}


void hMLP_BD_Limiter::limit(const ushort cell_index, const Cell_Type Cell_Type, Discrete_Solution_DG& discrete_solution, ushort& projection_degree, const MLP_Criterion_Base& stability_criterion) const
{
	switch (Cell_Type)
	{
	case Cell_Type::normal:
	case Cell_Type::smooth_extrema:
	{
		this->is_end_ = true;
		break;
	}
	case Cell_Type::trouble:
	case Cell_Type::typeII:
	{
		if (projection_degree <= 2)
		{
			this->apply_MLP_u1(cell_index, discrete_solution, stability_criterion);
			this->is_end_ = true;
		}
		else
		{
			discrete_solution.project_to_Pn_space(cell_index, --projection_degree);
			this->is_end_ = false;
		}

		break;
	}
	case Cell_Type::typeI:
	{
		this->apply_MLP_u1(cell_index, discrete_solution, stability_criterion);
		this->is_end_ = true;

		break;
	}
	default:
		EXCEPTION("not supported cell type");
	}
}