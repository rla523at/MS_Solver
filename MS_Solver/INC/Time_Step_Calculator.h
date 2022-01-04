#pragma once
#include "Governing_Equation.h"
#include "Grid.h"

class Time_Step_Calculator
{
public:
	virtual double calculate(const std::vector<Euclidean_Vector>& P0_solutions, const Governing_Equation& governing_equation) const abstract;
};

class CFL : public Time_Step_Calculator
{
public:
    CFL(const double cfl, const Grid& grid);

public://Query
    double calculate(const std::vector<Euclidean_Vector>& P0_solutions, const Governing_Equation& governing_equation) const override;

private:
    double cfl_;
    std::vector<double> cell_index_to_volume_reciprocal_table_;
    std::vector<std::vector<double>> cell_index_to_projected_volumes_table_;
};

class Constant_Dt : public Time_Step_Calculator
{
public:
    Constant_Dt(const double constant_dt) : constant_dt_(constant_dt) {};

public:
    double calculate(const std::vector<Euclidean_Vector>& P0_solutions, const Governing_Equation& governing_equation) const override;

private:
    double constant_dt_;
};

class Time_Step_Calculator_Factory
{
public:
    static std::unique_ptr<Time_Step_Calculator> make_unique(const Configuration& configuration, const Grid& grid);
};