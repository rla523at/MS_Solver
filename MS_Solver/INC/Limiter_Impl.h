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
	void apply_MLP_u1(const ushort cell_index, Discrete_Solution_DG& discrete_solution, const MLP_Criterion_Base& stability_criterion) const;

protected:
	MLP_u1 limiting_function_;
};

class hMLP_Limiter : public hMLP_Type_Limiter
{
public:
	hMLP_Limiter(const Grid& grid)
		:hMLP_Type_Limiter(grid) {};

public://Query
	void limit(const ushort cell_index, const Cell_Type Cell_Type, Discrete_Solution_DG& discrete_solution, ushort& projection_degree, const MLP_Criterion_Base& stability_criterion) const override;
};

class hMLP_BD_Limiter : public hMLP_Type_Limiter
{
public:
	hMLP_BD_Limiter(const Grid& grid)
		:hMLP_Type_Limiter(grid) {};

public://Query
	void limit(const ushort cell_index, const Cell_Type Cell_Type, Discrete_Solution_DG& discrete_solution, ushort& projection_degree, const MLP_Criterion_Base& stability_criterion) const override;
};