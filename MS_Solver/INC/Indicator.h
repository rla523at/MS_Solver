#pragma once
#include "Stability_Criterion.h"
#include "Measuring_Function.h"

enum class Cell_Type
{
    normal = 0,
    smooth_extrema = 1,
    trouble = 2,
    typeI = 3,
    typeII = 4
};

class Cell_Indicator
{
public://Command
    virtual void precalculate(const Discrete_Solution_DG& discrete_solution) abstract;

public://Query
    virtual Cell_Type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const abstract;
};

class Subcell_Oscillation_Indicator
{
public:
    Subcell_Oscillation_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution);

public://Query
    bool is_typeI_cell(const uint cell_index) const { return typeI_threshold_number_ <= this->cell_index_to_num_troubled_faces_[cell_index]; };
    bool is_typeII_cell(const uint cell_index) const { return typeII_threshold_number_ <= this->cell_index_to_num_troubled_faces_[cell_index]; };

private:
    static constexpr auto typeI_threshold_number_ = 2;
    static constexpr auto typeII_threshold_number_ = 1;
    uint num_infcs_;

    std::vector<double> infc_index_to_characteristic_length_table_;
    std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;

    std::vector<uint> cell_index_to_num_troubled_faces_;
    Average_Solution_Jump_Measuring_Function measuring_function_;
};

class Discontinuity_Indicator
{
public://Command
    virtual void precalculate(const Discrete_Solution_DG& discrete_solution) abstract;

public://Query
    virtual bool near_discontinuity(const uint cell_index) const
    {
        return cell_index_to_near_discontinuity_table_[cell_index];
    }

protected:
    std::vector<bool> cell_index_to_near_discontinuity_table_;
};

class Shock_Indicator
{
public:
    virtual void precalculate(const Discrete_Solution_DG& discrete_solution) abstract;

public:
    virtual bool near_shock(const uint cell_index) const 
    {
        return cell_index_to_near_shock_table_[cell_index];
    }

protected:
    std::vector<bool> cell_index_to_near_shock_table_;
};

