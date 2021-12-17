#pragma once
#include "Limiter_Function.h"

class Reconstruction_DG
{
public://Query
	virtual void reconstruct(Discrete_Solution_DG& discrete_solution) abstract;
};

class No_Reconstruction_DG : public Reconstruction_DG
{
public://Query
	void reconstruct(Discrete_Solution_DG& discrete_solution) override {};
};

class hMLP_Reconstruction_DG : public Reconstruction_DG
{
public:
    hMLP_Reconstruction_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
    void reconstruct(Discrete_Solution_DG& discrete_solution) override;

private:
    void precalculate(const Discrete_Solution_DG& discrete_solution);

private:
    uint num_cells_ = 0;
	
    MLP_Criterion stability_criterion_;
    MLP_Indicator indicator_;
    MLP_u1_Limiter limiter_;

    //precalculate
    std::vector<std::vector<double>> set_of_P1_projected_criterion_values_at_verticies_;
};

class Reconstruction_DG_Factory//static class
{
public://Query
    static std::unique_ptr<Reconstruction_DG> make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution);

private:
    Reconstruction_DG_Factory(void) = delete;
};