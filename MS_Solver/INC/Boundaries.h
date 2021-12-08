#pragma once
#include "Boundary_Flux_Function.h"
#include "Discrete_Solution.h"
#include "Residual.h"

class Boundaries_DG
{
public:
    Boundaries_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux);

public://Query
    void calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_soltuion) const;

private:
    void reset_residual_values(const uint oc_index) const;

private:
    uint num_boundaries_;
    std::vector<uint> oc_indexes_;
    std::vector<std::unique_ptr<Boundary_Flux_Function>> boundary_flux_functions_;
    std::vector<Matrix> set_of_oc_side_QWs_basis_m_;
    std::vector<std::vector<Euclidean_Vector>> set_of_normals_;

    //for construction optimization
    static constexpr ushort max_num_equation = 5;
    static constexpr ushort max_num_basis = 200;
    static constexpr ushort max_num_QPs = 200;

    mutable std::vector<Euclidean_Vector> solution_v_at_QPs_;
    mutable std::vector<Matrix> set_of_boundary_flux_QPs_m_;
    mutable std::array<double, max_num_equation * max_num_basis> residual_values_ = { 0 };
};