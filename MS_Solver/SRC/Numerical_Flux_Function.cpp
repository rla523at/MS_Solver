#include "../INC/Numerical_Flux_Function.h"

Euclidean_Vector LLF::calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const 
{
    const auto oc_physical_flux = this->governing_equation_->calculate_physical_flux(oc_solution);
    const auto nc_physical_flux = this->governing_equation_->calculate_physical_flux(nc_solution);
    const auto inner_face_maximum_lambda = this->governing_equation_->calculate_inner_face_maximum_lambda(oc_solution, nc_solution, normal_vector);

    return 0.5 * ((oc_physical_flux + nc_physical_flux) * normal_vector + inner_face_maximum_lambda * (oc_solution - nc_solution));
}

std::shared_ptr<Numerical_Flux_Function> Numerical_Flux_Function_Factory::make_shared(const Configuration& configuration, const std::shared_ptr<Governing_Equation>& governing_equation)
{
    const auto name = configuration.get("numerical_flux");

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
