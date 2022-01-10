#pragma once
#include "Measuring_Function_Impl.h"
#include "Trouble_Boundary_Indicator_Impl.h"

//TF
//T : Trouble Boundary Indicator Type
//F : Face Jump Measurer Type
class Trouble_Boundary_Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Trouble_Boundary_Indicator> make_11_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto face_jump_measurer = std::make_unique<Face_Jump_Measurer_Type1>(grid, discrete_solution, criterion_solution_index);
        return std::make_unique<Trouble_Boundary_Indicator_Type1>(grid, std::move(face_jump_measurer));
    }
    static std::unique_ptr<Trouble_Boundary_Indicator> make_12_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto face_jump_measurer = std::make_unique<Face_Jump_Measurer_Type2>(grid, discrete_solution, criterion_solution_index);
        return std::make_unique<Trouble_Boundary_Indicator_Type1>(grid, std::move(face_jump_measurer));
    }
    static std::unique_ptr<Trouble_Boundary_Indicator> make_21_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto face_jump_measurer = std::make_unique<Face_Jump_Measurer_Type1>(grid, discrete_solution, criterion_solution_index);
        return std::make_unique<Trouble_Boundary_Indicator_Type2>(grid, std::move(face_jump_measurer));
    }
    static std::unique_ptr<Trouble_Boundary_Indicator> make_22_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto face_jump_measurer = std::make_unique<Face_Jump_Measurer_Type2>(grid, discrete_solution, criterion_solution_index);
        return std::make_unique<Trouble_Boundary_Indicator_Type2>(grid, std::move(face_jump_measurer));
    }

private:
    Trouble_Boundary_Indicator_Factory(void) = delete;
};