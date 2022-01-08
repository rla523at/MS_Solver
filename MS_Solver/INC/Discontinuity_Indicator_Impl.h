#pragma once
#include "Indicator.h"
#include "Measuring_Function.h"

class Always_False_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    void check(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool is_near_discontinuity(const uint cell_index) const override { return false; };
};

class Always_True_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    void check(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool is_near_discontinuity(const uint cell_index) const override { return true; };
};

// 0.1 <= max (|rho - rho_j| / |rho|) ==> Max scaled average difference가 0.1 이상이면 discontinuity로 보겠다.
class Type1_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Type1_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
    void check(const Discrete_Solution_DG& discrete_solution) override;

private:
    static constexpr auto rho_index = 0;
    static constexpr auto threshold_number_ = 0.1;

    uint num_inner_faces_;
    std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;
    Scaled_Average_Difference_Measurer measuring_function_;    
};

// h <= max INT_{Omega} (q - q_j) / ||Omega|| ==> Max extrapolation differences가 h 이상이면 discontinuity로 보겠다.
// h = characteristic length of cell i
class Type2_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Type2_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
    void check(const Discrete_Solution_DG& discrete_solution) override;

public:
    static constexpr auto rho_index = 0;    

    uint num_cells_;
    std::vector<double> cell_index_to_characteristic_length_;
    Extrapolation_Differences_Measurer measuring_function_;
};

// Max(div(u)) <= 0 ==> Max div(u)가 0보다 작으면 discontinuity로 보겠다.
class Type3_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Type3_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
    void check(const Discrete_Solution_DG& discrete_solution) override;

public:
    static constexpr auto threshold_number_ = 0.0;

    uint num_cells_;
    Divergence_Velocity_Measurer measuring_function_;
};

// h <= AVG INT_{w_j} (q - q_j) / ||w_j|| ==> average of average solution jump가 h 이상이면 discontinuity로 보겠다.
// h = characteristic length of cell i
class Type4_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Type4_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, std::unique_ptr<Average_Solution_Jump_Measurer>&& measurer);

public:
    void check(const Discrete_Solution_DG& discrete_solution) override;

public:
    static constexpr auto rho_index_ = 0;

    uint num_cells_ = 0;
    uint num_infcs_ = 0;
    std::vector<double> cell_index_to_characteristic_length_table_;
    std::vector<ushort> cell_index_to_num_infc_table_;
    std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;
    std::unique_ptr<Average_Solution_Jump_Measurer> measurer_;
};

