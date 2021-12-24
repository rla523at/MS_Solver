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
    hMLP_Reconstruction_DG(MLP_Criterion&& stability_criterion, MLP_Indicator&& MLP_indicator, MLP_u1_Limiter&& MLP_u1_limiter, const Grid& grid);

public:
    void reconstruct(Discrete_Solution_DG& discrete_solution) override;

private:
    void precalculate(const Discrete_Solution_DG& discrete_solution);

private:
    uint num_cells_ = 0;
	
    MLP_Criterion stability_criterion_;
    MLP_Indicator MLP_indicator_;
    MLP_u1_Limiter MLP_u1_limiter_;

    //precalculate
    std::vector<std::vector<double>> set_of_P1_projected_criterion_values_at_verticies_;
};

class hMLP_BD_Reconstruction_DG : public Reconstruction_DG
{
public:
    hMLP_BD_Reconstruction_DG(Simplex_Decomposed_MLP_Criterion&& stability_criterion, MLP_Indicator&& MLP_indicator, MLP_u1_Limiter&& MLP_u1_limiter,
        Subcell_Oscillation_Indicator&& boundary_indicator, Shock_Indicator&& shock_indicator,  const Grid& grid);

public:
    void reconstruct(Discrete_Solution_DG& discrete_solution) override;

private:
    void precalculate(const Discrete_Solution_DG& discrete_solution);
    void limit_solution(Discrete_Solution_DG& discrete_solution, const uint cell_index, ushort& solution_degree) const;

private:
    uint num_cells_ = 0;
    mutable bool is_end_limting_ = false;

    Simplex_Decomposed_MLP_Criterion stability_criterion_;
    MLP_Indicator MLP_indicator_;
    MLP_u1_Limiter MLP_u1_limiter_;
    Subcell_Oscillation_Indicator subcell_oscillation_indicator;
    Shock_Indicator shock_indicator_;

    //precalculate
    std::vector<std::vector<double>> set_of_P1_projected_criterion_values_at_verticies_;
};


#include "../INC/Post_Processor.h"
class Test_Reconstuction_DG : public Reconstruction_DG
{
public:
    Test_Reconstuction_DG(Discontinuity_Indicator&& discontinuity_indicator)
        : discontinuity_indicator_(discontinuity_indicator) {};

public:
    void reconstruct(Discrete_Solution_DG& discrete_solution) override;

private:
    Discontinuity_Indicator discontinuity_indicator_;
};

class Reconstruction_DG_Factory//static class
{
public://Query
    static std::unique_ptr<Reconstruction_DG> make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution);

private:
    Reconstruction_DG_Factory(void) = delete;
};