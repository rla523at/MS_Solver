#pragma once
#include "Boundary_Flux_Function.h"
#include "Discrete_Solution.h"
#include "Residual.h"

class Boundaries_DG
{
public:
    Boundaries_DG(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux);

public://Query
    void calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_soltuion) const;

private:
    uint num_boundaries_;
    ushort num_equations_;
    std::vector<uint> oc_indexes_;
    std::vector<std::unique_ptr<Boundary_Flux_Function>> boundary_flux_functions_;
    std::vector<Matrix> set_of_oc_side_basis_QPs_m_;
    std::vector<Matrix> set_of_oc_side_QWs_basis_m_;
    std::vector<std::vector<Euclidean_Vector>> set_of_normals_;
};