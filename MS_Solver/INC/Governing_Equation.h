#pragma once
#include "Configuration.h"
#include "Euclidean_Vector.h"
#include "Matrix.h"

enum class Governing_Equation_Type
{
    Linear_Advection,
    Burgers,
    Euler,
    Not_Supported
};

class Governing_Equation
{
public://Query    
    const std::vector<std::string>& get_solution_names(void) const { return solution_names_; };
    ushort num_equations(void) const { return this->num_equations_; };
    ushort num_solutions(void) const { return this->num_solutions_; };
    ushort space_dimension(void) const { return this->space_dimension_; };
    Governing_Equation_Type type(void) const { return this->type_; };
    
    virtual std::vector<std::vector<double>> calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(const std::vector<Euclidean_Vector>& P0_solutions) const abstract;
    virtual double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const abstract;
    virtual Matrix calculate_physical_flux(const Euclidean_Vector& solution) const abstract;
    virtual void calculate_physical_flux(Matrix& physical_flux, const Euclidean_Vector& solution) const abstract;
    virtual void extend_to_solution(Euclidean_Vector& governing_equation_solution) const abstract;
    virtual void extend_to_solution(const double* GE_solution_values, double* solution_values) const abstract;		
    virtual short pressure_index(void) const abstract;

protected:
	ushort num_equations_ = 0;
    ushort num_solutions_ = 0;
    ushort space_dimension_ = 0;
    std::vector<std::string> solution_names_;
    Governing_Equation_Type type_ = Governing_Equation_Type::Not_Supported;
};