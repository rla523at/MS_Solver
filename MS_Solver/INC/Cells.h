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
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

public:
    Cells(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method) : Cells_FVM<Governing_Equation>(grid) {};
};

template <typename Governing_Equation, typename Reconstruction_Method>
class Cells<Governing_Equation, HOM, Reconstruction_Method> : public Cells_HOM<Governing_Equation, Reconstruction_Method>
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

public:
    Cells(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method) : Cells_HOM<Governing_Equation>(grid, reconstruction_method) {};
};