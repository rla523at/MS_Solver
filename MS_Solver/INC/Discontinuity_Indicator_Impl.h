#pragma once
#include "Indicator.h"
#include "Measuring_Function.h"

class Always_False_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool near_discontinuity(const uint cell_index) const override { return false; };
};

class Always_True_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool near_discontinuity(const uint cell_index) const override { return true; };
};

// 0.1 <= max (|rho - rho_j| / |rho|) ==> Max scaled average difference�� 0.1 �̻��̸� discontinuity�� ���ڴ�.
class Type1_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Type1_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override;

private:
    static constexpr auto rho_index = 0;
    static constexpr auto threshold_number_ = 0.1;

    uint num_inner_faces_;
    std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;
    Scaled_Average_Difference_Measuring_Function measuring_function_;    
};

// h <= max INT_{Omega} (q - q_j) / ||Omega|| ==> Max extrapolation differences�� h �̻��̸� discontinuity�� ���ڴ�.
// h = characteristic length of cell i
class Type2_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Type2_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override;

public:
    static constexpr auto rho_index = 0;    

    uint num_cells_;
    std::vector<double> cell_index_to_characteristic_length_;
    Extrapolation_Differences_Measuring_Function measuring_function_;
};

// Max(div(u)) <= 0 ==> Max div(u)�� 0���� ������ discontinuity�� ���ڴ�.
class Type3_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Type3_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override;

public:
    static constexpr auto threshold_number_ = 0.0;

    uint num_cells_;
    Divergence_Velocity_Measuring_Function measuring_function_;
};

class Discontinuity_Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Discontinuity_Indicator> make_unique(const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution);
    static std::unique_ptr<Discontinuity_Indicator> make_always_true_indicator(void)
    {
        return std::make_unique<Always_True_Discontinuity_Indicator>();
    }
    static std::unique_ptr<Discontinuity_Indicator> make_always_false_indicator(void)
    {
        return std::make_unique<Always_False_Discontinuity_Indicator>();
    }

private:
    Discontinuity_Indicator_Factory(void) = delete;
};