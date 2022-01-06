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
    Discontinuity_Indicating_Function_Base(const std::vector<double>& cell_index_to_threshold_value_table, std::unique_ptr<Discontinuity_Measuring_Function>&& measuring_function);

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override;

public://Query
    bool has_discontinuity(const uint cell_index) const override
    {
        return this->cell_index_to_has_discontinuity_table_[cell_index];
    }

protected:
    virtual bool compare(const double indicating_value, const double threshold_value) const abstract;

private:
    uint num_cells_ = 0;
    std::vector<double> cell_index_to_threshold_value_table_;
    std::unique_ptr<Discontinuity_Measuring_Function> measuring_function_;
    std::vector<bool> cell_index_to_has_discontinuity_table_;
};

class Less_is_Discontinuity : public Discontinuity_Indicating_Function_Base
{
public:
    Less_is_Discontinuity(const std::vector<double>& cell_index_to_threshold_value_table, std::unique_ptr<Discontinuity_Measuring_Function>&& measuring_function)
        :Discontinuity_Indicating_Function_Base(cell_index_to_threshold_value_table, std::move(measuring_function)) {};

private://Command
    bool compare(const double measuring_value, const double threshold_value) const
    {
        return measuring_value < threshold_value;
    };
};

class Greater_is_Discontinuity : public Discontinuity_Indicating_Function_Base
{
public:
    Greater_is_Discontinuity(const std::vector<double>& cell_index_to_threshold_value_table, std::unique_ptr<Discontinuity_Measuring_Function>&& measuring_function)
        :Discontinuity_Indicating_Function_Base(cell_index_to_threshold_value_table, std::move(measuring_function)) {};

private://Command
    bool compare(const double measuring_value, const double threshold_value) const
    {
        return threshold_value < measuring_value;
    };
};