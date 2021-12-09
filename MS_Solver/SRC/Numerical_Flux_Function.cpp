#include "../INC/Numerical_Flux_Function.h"

#include <iostream>

Euclidean_Vector LLF::calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const 
{
    this->governing_equation_->calculate_physical_flux(this->oc_physical_flux_, oc_solution);
    this->governing_equation_->calculate_physical_flux(this->nc_physical_flux_, nc_solution);
    const auto inner_face_maximum_lambda = this->governing_equation_->calculate_inner_face_maximum_lambda(oc_solution, nc_solution, normal_vector);

    this->central_normal_flux_.initalize();
    ms::mpm(this->oc_physical_flux_, this->nc_physical_flux_, this->central_flux_.data());
    ms::mv(this->central_flux_, normal_vector, this->central_normal_flux_.data());
    this->central_normal_flux_ *= 0.5;

    ms::vmv(oc_solution, nc_solution, this->correction_flux_.data());
    this->correction_flux_ *= (0.5 * inner_face_maximum_lambda);

    return this->central_normal_flux_ + this->correction_flux_;
}

void LLF::calculate(double* numerical_flux_ptr, const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const
{
    this->governing_equation_->calculate_physical_flux(this->oc_physical_flux_, oc_solution);
    this->governing_equation_->calculate_physical_flux(this->nc_physical_flux_, nc_solution);
    const auto inner_face_maximum_lambda = this->governing_equation_->calculate_inner_face_maximum_lambda(oc_solution, nc_solution, normal_vector);

    this->central_normal_flux_.initalize();
    ms::mpm(this->oc_physical_flux_, this->nc_physical_flux_, this->central_flux_.data());
    ms::mv(this->central_flux_, normal_vector, this->central_normal_flux_.data());
    this->central_normal_flux_ *= 0.5;

    ms::vmv(oc_solution, nc_solution, this->correction_flux_.data());
    this->correction_flux_ *= (0.5 * inner_face_maximum_lambda);

    ms::vpv(this->central_normal_flux_, this->correction_flux_, numerical_flux_ptr);
}


std::shared_ptr<Numerical_Flux_Function> Numerical_Flux_Function_Factory::make_shared(const Configuration& configuration, const std::shared_ptr<Governing_Equation>& governing_equation)
{
    const auto name = configuration.get_numerical_flux();

    if (ms::contains_icase(name, "LLF"))
    {
        return std::make_shared<LLF>(governing_equation);
    }
    else
    {
        EXCEPTION("numerical flux in configuration file is not supported");
        return nullptr;
    }
}
