#pragma once
#include "Configuration.h"
#include "Euclidean_Vector.h"
#include "Matrix.h"

using ushort = unsigned short;

class Governing_Equation
{
public://Query    
    const std::vector<std::string>& get_solution_names(void) const;
    ushort num_equations(void) const;
    ushort num_solutions(void) const;
    ushort space_dimension(void) const;
    
    virtual std::vector<std::vector<double>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const abstract;
    virtual double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const abstract;
    virtual Matrix calculate_physical_flux(const Euclidean_Vector& solution) const abstract;
    virtual void extend_to_solution(Euclidean_Vector& governing_equation_solution) const abstract;
		
protected:
	ushort num_equations_;
    ushort num_solutions_;
    ushort space_dimension_;
    std::vector<std::string> solution_names_;
};

class Scalar_Equation : public Governing_Equation
{
public:
    Scalar_Equation(void);

public:
    void extend_to_solution(Euclidean_Vector& governing_equation_solution) const override {};
};

class Linear_Advection : public Scalar_Equation
{
public://Query
    std::vector<std::vector<double>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
    //Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
    const Euclidean_Vector& get_advection_speed_vector(void) const;

protected:
    Euclidean_Vector advection_speeds_;
};

class Linear_Advection_2D : public Linear_Advection
{
public:
    Linear_Advection_2D(const double x_advection_speed, const double y_advection_speed);
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
};

class Linear_Advection_3D : public Linear_Advection
{
public:
    Linear_Advection_3D(const double x_advection_speed, const double y_advection_speed, const double z_advection_speed);
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
};

class Burgers : public Scalar_Equation
{
public:
    std::vector<std::vector<double>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
};

class Burgers_2D : public Burgers
{
public:
    Burgers_2D(void);
};

class Burgers_3D : public Burgers
{
public:
    Burgers_3D(void);
};

class Euler_2D : public Governing_Equation
{
public:
	Euler_2D(void);

public://Query
    std::vector<std::vector<double>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
	Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
    void extend_to_solution(Euclidean_Vector& governing_equation_solution) const override;

private:
    static constexpr auto gamma_ = 1.4;
};

class Euler_3D : public Governing_Equation
{
public:
    Euler_3D(void);

public://Query
    std::vector<std::vector<double>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
    void extend_to_solution(Euclidean_Vector& governing_equation_solution) const override;

private:
    static constexpr auto gamma_ = 1.4;
};

class Governing_Equation_Factory//static class
{
public:
    static std::shared_ptr<Governing_Equation> make_shared(const Configuration& config);
    static std::unique_ptr<Governing_Equation> make_unique(const Configuration& config);

private:
	Governing_Equation_Factory(void) = delete;
};