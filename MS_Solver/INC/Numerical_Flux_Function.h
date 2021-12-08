#pragma once
#include "Governing_Equation.h"


class Numerical_Flux_Function abstract
{
public:
	virtual Euclidean_Vector calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const abstract;
};

class LLF : public Numerical_Flux_Function
{
public:
    LLF(const std::shared_ptr<Governing_Equation>& governing_equation)
        : governing_equation_(governing_equation)
        , oc_physical_flux_(governing_equation->num_equations(), governing_equation->space_dimension())
        , nc_physical_flux_(governing_equation->num_equations(), governing_equation->space_dimension())
        , central_flux_(governing_equation->num_equations(), governing_equation->space_dimension())
        , central_normal_flux_(governing_equation->num_equations())
        , correction_flux_(governing_equation->num_equations()) {};

public:
    Euclidean_Vector calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const override;

private:
    std::shared_ptr<Governing_Equation> governing_equation_;

    //for construction optimization
    mutable Matrix oc_physical_flux_;
    mutable Matrix nc_physical_flux_;
    mutable Matrix central_flux_;
    mutable Euclidean_Vector central_normal_flux_;
    mutable Euclidean_Vector correction_flux_;
};

class Numerical_Flux_Function_Factory//static class
{
public:
    static std::shared_ptr<Numerical_Flux_Function> make_shared(const Configuration& configuration, const std::shared_ptr<Governing_Equation>& governing_equation);

private:
    Numerical_Flux_Function_Factory(void) = delete;
};