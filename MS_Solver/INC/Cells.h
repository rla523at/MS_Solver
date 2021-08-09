#pragma once
#include "Cells_FVM.h"
#include "Spatial_Discrete_Method.h"


template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Cells;


template <typename Governing_Equation, typename Reconstruction_Method>
class Cells<Governing_Equation, FVM, Reconstruction_Method> : public Cells_FVM<Governing_Equation::space_dimension()>
{   
private:
    Cells(void) = delete;
};