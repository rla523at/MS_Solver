#pragma once
#include "Boundaries_FVM.h"
#include "Spatial_Discrete_Method.h"

template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Boundaries;


template <typename Governing_Equation>
class Boundaries<Governing_Equation, FVM, Constant_Reconstruction> : public Boundaries_FVM_Constant<Governing_Equation>
{
private:
    Boundaries(void) = delete;
};


template <typename Governing_Equation, typename Reconstruction_Method>
class Boundaries<Governing_Equation, FVM, Reconstruction_Method> : public Boundaries_FVM_Linear<Governing_Equation, Reconstruction_Method>
{
private:
    Boundaries(void) = delete;
};
