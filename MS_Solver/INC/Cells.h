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
    Cells_DG(const std::shared_ptr<Governing_Equation>& governing_equation, std::unique_ptr<Time_Step_Calculator>&& time_step_calculator, const Grid& grid, const Discrete_Solution_DG& discrete_solution);

public://Query
    double calculate_time_step(const Discrete_Solution_DG& discrete_solution) const;
    void calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_solution) const;

private:
    std::vector<double> P0_basis_values_;
    std::vector<Matrix> set_of_basis_QPs_m_;
    std::vector<Matrix> set_of_QWs_gradient_basis_m_;
};








































//class Cells
//{
//public:
//    Cells(const Configuration& configuration, const Grid& grid);
//
//public://Command
//    virtual void update_solution(Euclidean_Vector&& updated_solution_v) abstract;
//
//public://Query
//    virtual double calculate_time_step(void) const abstract;
//    virtual void calculate_RHS(double* RHS) const abstract;
//    virtual const Euclidean_Vector& get_solution_vector(void) const abstract;
//    virtual size_t num_solution_values(void) const abstract;
//
//protected:
//    size_t num_cells_;
//    std::unique_ptr<Time_Step_Calculator> time_step_calculator_;
//    std::unique_ptr<Governing_Equation> governing_equation_;
//};
//
//class Cells_DG : public Cells
//{
//public:
//    Cells_DG(const Configuration& configuration, const Grid& grid);
//
//public://Command
//    void update_solution(Euclidean_Vector&& updated_solution_v) override;
//
//public://Query
//    double calculate_time_step(void) const override;
//    void calculate_RHS(double* rhs) const override;
//    const Euclidean_Vector& get_solution_vector(void) const override;
//    size_t num_solution_values(void) const override;
//
//private:
//    void update_rhs(const uint cell_index, double* RHS, const Matrix& delta_rhs) const;
//    //std::vector<ushort> calculate_integrand_degrees(const std::vector<ushort>& solution_degrees) const;
//
//private:
//    Discrete_Solution_DG discrete_solution_;
//    std::vector<double> P0_basis_values_;
//    std::vector<Matrix> set_of_basis_QPs_m_;
//    std::vector<Matrix> set_of_QWs_gradient_basis_m_;
//};



























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