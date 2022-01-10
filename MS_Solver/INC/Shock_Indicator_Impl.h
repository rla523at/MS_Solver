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

// 0.1 <= max (|p - p_j| / min(|p|,|p_j|)) ==> Max scaled average difference�� 0.1 �̻��̸� shock ���� ���ڴ�.
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
    static std::unique_ptr<Shock_Indicator> make_unique(const std::string& governing_equation_name, const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution);
    static std::unique_ptr<Shock_Indicator> make_always_false_indicator(void)
    {
        return std::make_unique<Always_False_Shock_Indicator>();
    }
    static std::unique_ptr<Shock_Indicator> make_type1_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::string& governing_equation_name)
    {
        const auto criterion_solution_index = find_pressure_index(governing_equation_name, grid.space_dimension());
        return std::make_unique<Shock_Indicator_Type1>(grid, discrete_solution, criterion_solution_index);
    }
    static ushort find_pressure_index(const std::string& governing_equation_name, const ushort space_dimension)
    {
        if (ms::compare_icase(governing_equation_name, "Burgers"))
        {
            return 0;
        }
        else if (ms::compare_icase(governing_equation_name, "Euler"))
        {
            if (space_dimension == 2)
            {
                return Euler_2D::pressure_index();
            }
            else if (space_dimension == 3)
            {
                return Euler_3D::pressure_index();
            }
            else
            {
                EXCEPTION("not supported dimension");
                return -1;
            }
        }
        else
        {
            EXCEPTION("not supported governing equation");
            return -1;
        }
    }

private:
    Shock_Indicator_Factory(void) = delete;
};