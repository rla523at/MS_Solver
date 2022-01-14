#pragma once
#include "Governing_Equation.h"


class Scalar_Equation : public Governing_Equation
{
public:
    Scalar_Equation(void);

public:
    void extend_to_solution(Euclidean_Vector& governing_equation_solution) const override {};
    void extend_to_solution(const double* GE_solution_values, double* solution_values) const override { solution_values[0] = GE_solution_values[0]; };
    short pressure_index(void) const override { EXCEPTION("Scalar Equation doesn't have a notion of pressure"); return -1; };
};

class Linear_Advection : public Scalar_Equation
{
public:
    Linear_Advection(void) { this->type_ = Governing_Equation_Type::Linear_Advection; };

public://Query
    std::vector<std::vector<double>> calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
    //const Euclidean_Vector& get_advection_speed_vector(void) const { return this->advection_speeds_; };

protected:
    Euclidean_Vector advection_speeds_;
};

class Linear_Advection_2D : public Linear_Advection
{
public:
    Linear_Advection_2D(const double x_advection_speed, const double y_advection_speed);
    Linear_Advection_2D(const double* advection_speeds);

public:
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
    void calculate_physical_flux(Matrix& physical_flux, const Euclidean_Vector& solution) const override;
};

class Linear_Advection_3D : public Linear_Advection
{
public:
    Linear_Advection_3D(const double x_advection_speed, const double y_advection_speed, const double z_advection_speed);
    Linear_Advection_3D(const double* advection_speeds);

public:
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
    void calculate_physical_flux(Matrix& physical_flux, const Euclidean_Vector& solution) const override {};
};

class Burgers : public Scalar_Equation
{
public:
    Burgers(void) { this->type_ = Governing_Equation_Type::Burgers; };

public:
    std::vector<std::vector<double>> calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
    void calculate_physical_flux(Matrix& physical_flux, const Euclidean_Vector& solution) const override;
};

class Burgers_2D : public Burgers
{
public:
    Burgers_2D(void) { this->space_dimension_ = 2; };
};

class Burgers_3D : public Burgers
{
public:
    Burgers_3D(void) { this->space_dimension_ = 3; };

};

class Fluid_Governing_Equation : public Governing_Equation
{
protected:
    void check_non_physical_value(const double density, const double pressure) const
    {
        if (density < 0.0 || pressure < 0.0)
        {
            throw std::invalid_argument("density and pressure shold be positive");
        }
    }

protected:
    static constexpr double gamma_ = 1.4;
};

class Euler : public Fluid_Governing_Equation
{
public:
    Euler(void) { this->type_ = Governing_Equation_Type::Euler; };
};

class Euler_2D : public Euler
{
public:
    Euler_2D(void);

public://Query
    std::vector<std::vector<double>> calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
    void calculate_physical_flux(Matrix& physical_flux, const Euclidean_Vector& solution) const override;
    void extend_to_solution(Euclidean_Vector& governing_equation_solution) const override;
    void extend_to_solution(const double* GE_solution_values, double* solution_values) const override;
    short pressure_index(void) const override { return 6; };
};

class Euler_3D : public Euler
{
public:
    Euler_3D(void);

public://Query
    std::vector<std::vector<double>> calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
    void calculate_physical_flux(Matrix& physical_flux, const Euclidean_Vector& solution) const override;
    void extend_to_solution(Euclidean_Vector& governing_equation_solution) const override;
    void extend_to_solution(const double* GE_solution_values, double* solution_values) const override;
    short pressure_index(void) const override { return 8; };
};

class Governing_Equation_Factory//static class
{
public:
    static std::shared_ptr<Governing_Equation> make_shared(const Configuration& config);

private:
    Governing_Equation_Factory(void) = delete;
};