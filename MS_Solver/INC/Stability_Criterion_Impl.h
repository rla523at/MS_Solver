#pragma once
#include "Stability_Criterion.h"

class MLP_Criterion : public MLP_Criterion_Base
{
public:
    MLP_Criterion(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);

public://Command
    void precaclulate(const Discrete_Solution_DG& discrete_solution) override;

private:
    //precalculate
    std::vector<double> P0_values_;
};

class Simplex_Decomposed_MLP_Criterion : public MLP_Criterion_Base
{
public:
    Simplex_Decomposed_MLP_Criterion(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);

public://Command
    void precaclulate(const Discrete_Solution_DG& discrete_solution) override;

private:
    std::unordered_map<uint, std::set<uint>> vertex_index_to_matched_vertex_index_set_;

    //precalculate
    std::vector<std::map<uint, double>> set_of_vertex_index_to_simplex_P0_value_;
};