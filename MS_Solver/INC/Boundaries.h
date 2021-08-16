#pragma once
#include "Boundaries_FVM.h"
#include "Boundaries_HOM.h"
#include "Spatial_Discrete_Method.h"

template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Boundaries;


//template <typename Governing_Equation>
//class Boundaries<Governing_Equation, FVM, Constant_Reconstruction> : public Boundaries_FVM_Constant<Governing_Equation>
//{
//private:
//    Boundaries(void) = delete;
//};


template <typename Governing_Equation, typename Reconstruction_Method>
class Boundaries<Governing_Equation, FVM, Reconstruction_Method> : public Boundaries_FVM_Linear<Governing_Equation, Reconstruction_Method>
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

public:
    Boundaries(Grid<space_dimension_>&& grid, const Reconstruction_Method& reconstruction_method)
        : Boundaries_FVM_Linear<Governing_Equation, Reconstruction_Method>(std::move(grid), reconstruction_method) {};
};


template <typename Governing_Equation, typename Reconstruction_Method>
class Boundaries<Governing_Equation, HOM, Reconstruction_Method> : public Boundaries_HOM<Governing_Equation, Reconstruction_Method>
{
private:
    Boundaries(void) = delete;
};