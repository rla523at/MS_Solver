#pragma once
#include "Boundaries_FVM.h"
#include "Boundaries_HOM.h"
#include "Spatial_Discrete_Method.h"

template <typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
class Boundaries;


template <typename Numerical_Flux_Function>
class Boundaries<FVM, Constant_Reconstruction, Numerical_Flux_Function> : public Boundaries_FVM_Constant<Numerical_Flux_Function>
{
private:
    static constexpr size_t space_dimension_ = Numerical_Flux_Function::space_dimension();

public:
    Boundaries(const Grid<space_dimension_>& grid, const Constant_Reconstruction& reconstruction_method)
        : Boundaries_FVM_Constant<Numerical_Flux_Function>(grid) {};
};


template <typename Reconstruction_Method, typename Numerical_Flux_Function>
class Boundaries<FVM, Reconstruction_Method, Numerical_Flux_Function> : public Boundaries_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>
{
private:
    static constexpr size_t space_dimension_ = Numerical_Flux_Function::space_dimension();

public:
    Boundaries(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
        : Boundaries_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>(grid, reconstruction_method) {};
};


template <typename Reconstruction_Method, typename Numerical_Flux_Function>
class Boundaries<HOM, Reconstruction_Method, Numerical_Flux_Function> : public Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>
{
private:
    static constexpr size_t space_dimension_ = Numerical_Flux_Function::space_dimension();

public:
    Boundaries(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
        : Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>(grid, reconstruction_method) {};
};