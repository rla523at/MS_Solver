#pragma once
#include "Indicating_Function.h"
#include "Measuring_Function.h"

class Always_False_Discontinuity_Indicator : public Discontinuity_Indicating_Function
{
public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool has_discontinuity(const uint cell_index) const override { return false; };
};

class Always_True_Discontinuity_Indicator : public Discontinuity_Indicating_Function
{
public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool has_discontinuity(const uint cell_index) const override { return true; };
};


class Discontinuity_Indicating_Function_Base : public Discontinuity_Indicating_Function
{
public:
    Discontinuity_Indicating_Function_Base(const std::vector<double>& cell_index_to_threshold_value_table, std::unique_ptr<Discontinuity_Measuring_Function>&& measuring_function)
        : cell_index_to_threshold_value_table_(cell_index_to_threshold_value_table)
        , measuring_function_(std::move(measuring_function))
    {
        this->num_cells_ = cell_index_to_threshold_value_table_.size();
        this->cell_index_to_has_discontinuity_table_.resize(this->num_cells_, false);
    }

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override
    {
        std::fill(this->cell_index_to_has_discontinuity_table_.begin(), this->cell_index_to_has_discontinuity_table_.end(), false);

        const auto cell_index_to_measuring_value_table = measuring_function_->measure(discrete_solution);
        for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
        {
            const auto threshold_value = cell_index_to_threshold_value_table_[cell_index];
            const auto measuring_value = cell_index_to_measuring_value_table[cell_index];

            if (threshold_value < measuring_value)
            {
                this->cell_index_to_has_discontinuity_table_[cell_index] = true;
            }
        }
    }

public://Query
    bool has_discontinuity(const uint cell_index) const override
    {
        return this->cell_index_to_has_discontinuity_table_[cell_index];
    }

private:
    uint num_cells_;
    std::vector<double> cell_index_to_threshold_value_table_;
    std::unique_ptr<Discontinuity_Measuring_Function> measuring_function_;
    std::vector<bool> cell_index_to_has_discontinuity_table_;
};


//class Discontinuity_Indicator_Base : public Discontinuity_Indicator
//{
//public:
//    Discontinuity_Indicator_Base(const Grid& grid, std::unique_ptr<Discontinuity_Measure_Function>&& measure_function)
//        : num_cells_(grid.num_cells())
//        , cell_index_to_has_discontinuity_table_(num_cells_, false)
//        , measure_function_(std::move(measure_function)) {};
//
//public:
//    void precalculate(const Discrete_Solution_DG& discrete_solution) override
//    {
//        std::fill(this->cell_index_to_has_discontinuity_table_.begin(), this->cell_index_to_has_discontinuity_table_.end(), false);
//
//        const auto indicator_values = this->measure_function_->measure_indicator_values(discrete_solution);
//
//        for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
//        {
//
//        }
//    }
//
//public://Query
//    bool has_discontinuity(const uint cell_index) const override { return this->cell_index_to_has_discontinuity_table_[cell_index]; };
//
//protected:
//    uint num_cells_;
//    std::vector<bool> cell_index_to_has_discontinuity_table_;
//    std::unique_ptr<Discontinuity_Measure_Function> measure_function_;
//};


//class Discontinuity_Indicator_Base : public Discontinuity_Indicator
//{
//public:
//    Discontinuity_Indicator_Base(const Grid& grid, const ushort criterion_solution_index)
//        : criterion_solution_index_(criterion_solution_index)
//        , num_cells_(grid.num_cells())
//        , cell_index_to_face_share_cell_indexes_table_(grid.cell_index_to_face_share_cell_indexes_table_consider_pbdry())
//        , cell_index_to_has_discontinuity_table_(this->num_cells_, false) {};
//
//public://Query
//    bool has_discontinuity(const uint cell_index) const override { return this->cell_index_to_has_discontinuity_table_[cell_index]; };
//
//protected:
//    ushort criterion_solution_index_;
//    uint num_cells_;
//
//    std::vector<std::vector<uint>> cell_index_to_face_share_cell_indexes_table_;
//    std::vector<bool> cell_index_to_has_discontinuity_table_;
//};

//class Heuristic_Discontinuity_Indicator : public Discontinuity_Indicator_Base
//{
//public:
//    Heuristic_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
//        : Discontinuity_Indicator_Base(grid, criterion_solution_index)
//    {
//        //for precalculation
//        discrete_solution.precalculate_cell_P0_basis_values();
//    };
//
//public://Command
//    void precalculate(const Discrete_Solution_DG& discrete_solution) override;
//};

//class Extrapolation_Discontinuity_Indicator : public Discontinuity_Indicator_Base
//{
//public:
//    Extrapolation_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);
//
//public://Command
//    void precalculate(const Discrete_Solution_DG& discrete_solution) override;
//
//private:
//    std::vector<double> cell_index_to_volume_reciprocal_table_;
//    std::vector<double> cell_index_to_threshold_value_table_;
//    std::vector<Euclidean_Vector> cell_index_to_QW_v_table_;
//};