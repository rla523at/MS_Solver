#pragma once
#include "Indicator.h"
#include "Measuring_Function.h"

class Always_False_Shock_Indicator : public Shock_Indicator
{
public:
    void check(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool is_near_shock(const uint cell_index) const override { return false; };
};

class Always_True_Shock_Indicator : public Shock_Indicator
{
public:
    void check(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool is_near_shock(const uint cell_index) const override { return true; };
};

// 0.1 <= max (|p - p_j| / min(|p|,|p_j|)) ==> Max scaled average difference가 0.1 이상이면 shock 으로 보겠다.
class Shock_Indicator_Type1 : public Shock_Indicator
{
public:
    Shock_Indicator_Type1(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort pressure_index);

public:
    void check(const Discrete_Solution_DG& discrete_solution) override;

private:
    static constexpr auto threshold_number_ = 0.1;

    uint num_inner_faces_;
    std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;
    Scaled_Average_Difference_Measurer measuring_function_;
};

class Shock_Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Shock_Indicator> make_always_false_indicator(void)
    {
        return std::make_unique<Always_False_Shock_Indicator>();
    }
    static std::unique_ptr<Shock_Indicator> make_type1_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Governing_Equation>& governing_equation)
    {
        const auto criterion_solution_index = find_criterion_index(governing_equation);
        return std::make_unique<Shock_Indicator_Type1>(grid, discrete_solution, criterion_solution_index);
    }
    static ushort find_criterion_index(const std::shared_ptr<Governing_Equation>& governing_equation)
    {
        const auto gov_eq_type = governing_equation->type();

        switch (gov_eq_type)
        {
        case Governing_Equation_Type::Burgers:
            return 0;
        case Governing_Equation_Type::Euler:
            return governing_equation->pressure_index();
        default:
            EXCEPTION("Wrong governing equation");
            return NULL;
        }
    }


private:
    Shock_Indicator_Factory(void) = delete;
};