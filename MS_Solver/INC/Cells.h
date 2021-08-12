#pragma once
#include "Cells_FVM.h"
#include "Cells_HOM.h"
#include "Spatial_Discrete_Method.h"


template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Cells;


template <typename Governing_Equation, typename Reconstruction_Method>
class Cells<Governing_Equation, FVM, Reconstruction_Method> : public Cells_FVM<Governing_Equation>
{   
private:
    Cells(void) = delete;
};

template <typename Governing_Equation, typename Reconstruction_Method>
class Cells<Governing_Equation, HOM, Reconstruction_Method> : public Cells_HOM<Governing_Equation, Reconstruction_Method>
{
private:
    Cells(void) = delete;
};