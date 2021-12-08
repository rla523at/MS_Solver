#pragma once
#include "Time_Step_Calculator.h"
#include "Discrete_Solution.h"
#include "Residual.h"

class Cells
{
public:
    Cells(const std::shared_ptr<Governing_Equation>& governing_equation, std::unique_ptr<Time_Step_Calculator>&& time_step_calculator, const Grid& grid);

protected:
    size_t num_cells_;
    std::shared_ptr<Governing_Equation> governing_equation_;
    std::unique_ptr<Time_Step_Calculator> time_step_calculator_;
};

class Cells_DG : public Cells
{
public:
    Cells_DG(const std::shared_ptr<Governing_Equation>& governing_equation, std::unique_ptr<Time_Step_Calculator>&& time_step_calculator, const Grid& grid, Discrete_Solution_DG& discrete_solution);

public://Query
    double calculate_time_step(const Discrete_Solution_DG& discrete_solution) const;
    void calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_solution) const;

private:
    void reset_residual_values(const uint oc_index) const;

private:
    std::vector<Matrix> set_of_QWs_gradient_basis_m_;

    //for construction optimization
    static constexpr ushort max_num_equation = 5;
    static constexpr ushort max_num_basis = 200;
    static constexpr ushort max_num_QPs = 200;

    std::vector<ushort> set_of_num_QPs_;
    mutable Matrix physical_flux_;
    //mutable std::vector<Euclidean_Vector> solution_v_at_QPs_;
    mutable std::vector<Matrix> set_of_flux_QPs_m_;
    mutable std::array<double, max_num_equation* max_num_basis> residual_values_ = { 0 };
    mutable std::array<Euclidean_Vector, max_num_QPs> solution_v_at_QPs_ = {};
};