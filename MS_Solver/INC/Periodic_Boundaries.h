#pragma once

#include "Periodic_Boundaries_FVM.h"
#include "Spatial_Discrete_Method.h"


template <typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
class Periodic_Boundaries;


template <typename Numerical_Flux_Function>
class Periodic_Boundaries<FVM, Constant_Reconstruction, Numerical_Flux_Function> : public Periodic_Boundaries_FVM_Constant<Numerical_Flux_Function>
{
private:
    Periodic_Boundaries(void) = delete;
};


template<typename Reconstruction_Method, typename Numerical_Flux_Function>
class Periodic_Boundaries<FVM, Reconstruction_Method, Numerical_Flux_Function> : public Periodic_Boundaries_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>
{
private:
    Periodic_Boundaries(void) = delete;
};