#include "../INC/Numerical_Flux_Function.h"

#include <iostream>

Euclidean_Vector LLF::calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const 
{
    //const auto oc_physical_flux = this->governing_equation_->calculate_physical_flux(oc_solution);
    //const auto nc_physical_flux = this->governing_equation_->calculate_physical_flux(nc_solution);
    //const auto inner_face_maximum_lambda = this->governing_equation_->calculate_inner_face_maximum_lambda(oc_solution, nc_solution, normal_vector);

    //return 0.5 * ((oc_physical_flux + nc_physical_flux) * normal_vector + inner_face_maximum_lambda * (oc_solution - nc_solution));

    this->governing_equation_->calculate_physical_flux(this->oc_physical_flux_, oc_solution);
    this->governing_equation_->calculate_physical_flux(this->nc_physical_flux_, nc_solution);
    const auto inner_face_maximum_lambda = this->governing_equation_->calculate_inner_face_maximum_lambda(oc_solution, nc_solution, normal_vector);

    ms::mpm(this->oc_physical_flux_, this->nc_physical_flux_, this->central_flux_);
    ms::mv(this->central_flux_, normal_vector, this->central_normal_flux_);
    this->central_normal_flux_ *= 0.5;

    ms::vmv(oc_solution, nc_solution, this->correction_flux_);
    this->correction_flux_ *= (0.5 * inner_face_maximum_lambda);

    return this->central_normal_flux_ + this->correction_flux_;



    //const auto oc_physical_flux = this->governing_equation_->calculate_physical_flux(oc_solution);
    //const auto nc_physical_flux = this->governing_equation_->calculate_physical_flux(nc_solution);
    //const auto inner_face_maximum_lambda = this->governing_equation_->calculate_inner_face_maximum_lambda(oc_solution, nc_solution, normal_vector);

    //const auto central_normal_flux = 0.5 * (oc_physical_flux + nc_physical_flux) * normal_vector;
    //const auto correction_flux = 0.5 * inner_face_maximum_lambda * (oc_solution - nc_solution);


    //if (oc_physical_flux != this->oc_physical_flux_)
    //{
    //    std::cout << "oc physcial flux";
    //    std::exit(523);
    //}
    //else if (nc_physical_flux != this->nc_physical_flux_)
    //{
    //    std::cout << "nc physcial flux";
    //    std::exit(523);
    //}
    //else if (central_normal_flux != this->central_normal_flux_)
    //{
    //    std::cout << "central_flux * normal " << (oc_physical_flux + nc_physical_flux) * normal_vector;
    //    std::cout << "central_flux * normal " << this->central_flux_;


    //    //std::cout << "central_normal_flux";
    //    //std::cout << central_normal_flux;
    //    //std::cout << this->central_normal_flux_;
    //    std::exit(523);
    //}
    //else if (correction_flux != this->correction_flux_)
    //{
    //    std::cout << "correction_flux";
    //    std::exit(523);
    //}
    //else
    //    return central_normal_flux + correction_flux;

    //centeral normal flux

    //std::cout << "oc_physical_flux_ " << this->oc_physical_flux_ << "\n";
    //std::cout << "nc_physical_flux_ " << this->nc_physical_flux_ << "\n";
    //ms::mpm(this->oc_physical_flux_, this->nc_physical_flux_, this->central_flux_);    
    //std::cout << "central_flux_ " << this->central_flux_ << "\n";
    //std::cout << "normal_vector " << normal_vector << "\n";
    //ms::mv(this->central_flux_, normal_vector, this->central_normal_flux_);
    //std::cout << "centeral_normal_flux_ " << central_normal_flux_ << "\n";
    //this->central_normal_flux_ *= 0.5;
    //std::cout << "central_normal_flux_ " << this->central_normal_flux_ << "\n";

    //std::cout << "oc_solution " << oc_solution << "\n";
    //std::cout << "nc_solution " << nc_solution << "\n";
    //ms::vmv(oc_solution, nc_solution, this->correction_flux_);
    //std::cout << "correction_flux_ " << correction_flux_ << "\n";
    //std::cout << "constant " << (0.5 * inner_face_maximum_lambda) << "\n";
    //this->correction_flux_ *= (0.5 * inner_face_maximum_lambda);
    //std::cout << "correction_flux_ " << correction_flux_ << "\n";

    //return this->central_normal_flux_ + this->correction_flux_;
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
