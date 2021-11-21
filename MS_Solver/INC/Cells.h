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
    std::unique_ptr<Time_Step_Calculator> time_step_calculator_;
    std::shared_ptr<Governing_Equation> governing_equation_;
};

class Cells_DG : public Cells
{
public:
    Cells_DG(const std::shared_ptr<Governing_Equation>& governing_equation, std::unique_ptr<Time_Step_Calculator>&& time_step_calculator, const Grid& grid, Discrete_Solution_DG& discrete_solution);

public://Query
    double calculate_time_step(const Discrete_Solution_DG& discrete_solution) const;
    void calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_solution) const;

private:
    std::vector<Matrix> set_of_QWs_gradient_basis_m_;
};