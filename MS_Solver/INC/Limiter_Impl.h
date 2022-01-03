#pragma once
#include "Limiter.h"
#include "Limiting_Function.h"

class hMLP_Type_Limiter : public Limiter
{
public:
	hMLP_Type_Limiter(const Grid& grid)
		:limiting_function_(grid) {};

public://Command
	void precalculate(const Discrete_Solution_DG& discrete_solution) override
	{
		this->limiting_function_.precalculate(discrete_solution);
	}

protected:
	void apply_MLP_u1(const ushort cell_index, Discrete_Solution_DG& discrete_solution, const MLP_Criterion_Base& stability_criterion) const
	{
		const auto projection_degree = 1;
		discrete_solution.project_to_Pn_space(cell_index, projection_degree);

		const auto limiting_value = this->limiting_function_.limiter_function(discrete_solution, cell_index, stability_criterion);
		discrete_solution.limit_slope(cell_index, limiting_value);
	}
protected:
	MLP_u1 limiting_function_;
};

class hMLP_Limiter : public hMLP_Type_Limiter
{
public:
	hMLP_Limiter(const Grid& grid)
		:hMLP_Type_Limiter(grid) {};

public://Query
    void limit(const ushort cell_index, const cell_type cell_type, Discrete_Solution_DG& discrete_solution, ushort& projection_degree, const MLP_Criterion_Base& stability_criterion) const override
    {
		switch (cell_type)
		{
		case cell_type::normal:
		case cell_type::smooth_extrema:
		{
			this->is_end_ = true;
			break;
		}
		case cell_type::trouble:
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
};

class hMLP_BD_Limiter : public hMLP_Type_Limiter
{
public:
	hMLP_BD_Limiter(const Grid& grid)
		:hMLP_Type_Limiter(grid) {};

public://Query
	void limit(const ushort cell_index, const cell_type cell_type, Discrete_Solution_DG& discrete_solution, ushort& projection_degree, const MLP_Criterion_Base& stability_criterion) const override
	{
		switch (cell_type)
		{
		case cell_type::normal:
		case cell_type::smooth_extrema:
		{
			this->is_end_ = true;
			break;
		}
		case cell_type::trouble:
		case cell_type::typeII:
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
		case cell_type::typeI:
		{
			this->apply_MLP_u1(cell_index, discrete_solution, stability_criterion);
			this->is_end_ = true;

			break;
		}
		default:
			EXCEPTION("not supported cell type");
		}
	}
};