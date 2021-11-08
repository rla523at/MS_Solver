#pragma once
#include "Time_Step_Calculator.h"
#include "Discrete_Solution.h"

class Cells
{
public:
    Cells(const Configuration& configuration, const Grid& grid);

public://Command
    virtual void update_solution(Euclidean_Vector&& updated_solution_v) abstract;

public://Query
    virtual double calculate_time_step(void) const abstract;
    virtual void calculate_RHS(double* RHS) const abstract;
    virtual const Euclidean_Vector& get_solution_vector(void) const abstract;
    virtual size_t num_solution_values(void) const abstract;

protected:
    size_t num_cells_;
    std::unique_ptr<Time_Step_Calculator> time_step_calculator_;
    std::unique_ptr<Governing_Equation> governing_equation_;
};

class Cells_HOM : public Cells
{
public:
    Cells_HOM(const Configuration& configuration, const Grid& grid);

public://Command
    void update_solution(Euclidean_Vector&& updated_solution_v) override;

public://Query
    double calculate_time_step(void) const override;
    void calculate_RHS(double* rhs) const override;
    const Euclidean_Vector& get_solution_vector(void) const override;
    size_t num_solution_values(void) const override;

private:
    void update_rhs(const uint cell_index, double* RHS, const Matrix& delta_rhs) const;

private:
    ushort num_equations_;
    ushort space_dimension_;

    Discrete_Solution_HOM discrete_solution_;
    std::vector<Matrix> set_of_QWs_gradient_basis_m_;

    std::vector<double> P0_basis_values_;
    std::vector<Matrix> set_of_basis_QPs_m_;
};


//#include "Cells_FVM.h"
//#include "Cells_HOM.h"
//#include "Spatial_Discrete_Method.h"
//
//
//template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
//class Cells;
//
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//class Cells<Governing_Equation, FVM, Reconstruction_Method> : public Cells_FVM<Governing_Equation>
//{   
//private:
//    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
//
//public:
//    Cells(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method) : Cells_FVM<Governing_Equation>(grid) {};
//};
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//class Cells<Governing_Equation, HOM, Reconstruction_Method> : public Cells_HOM<Governing_Equation, Reconstruction_Method>
//{
//private:
//    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
//
//public:
//    Cells(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method) : Cells_HOM<Governing_Equation, Reconstruction_Method>(grid, reconstruction_method) {};
//};