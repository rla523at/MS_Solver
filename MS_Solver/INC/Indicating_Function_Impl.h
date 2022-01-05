#pragma once
#include "Indicating_Function.h"

class Always_False_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool has_discontinuity(const uint cell_index) const override { return false; };
};

class Always_True_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override {};

public:
    bool has_discontinuity(const uint cell_index) const override { return true; };
};

//class Discontinuity_Measure_Function
//{
//public:
//    virtual std::vector<double> measure_indicator_values(const Discrete_Solution_DG& discrete_solution) const abstract;
//
//};
//
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


class Discontinuity_Indicator_Base : public Discontinuity_Indicator
{
public:
    Discontinuity_Indicator_Base(const Grid& grid, const ushort criterion_solution_index)
        : criterion_solution_index_(criterion_solution_index)
        , num_cells_(grid.num_cells())
        , cell_index_to_face_share_cell_indexes_table_(grid.cell_index_to_face_share_cell_indexes_table_consider_pbdry())
        , cell_index_to_has_discontinuity_table_(this->num_cells_, false) {};

public://Query
    bool has_discontinuity(const uint cell_index) const override { return this->cell_index_to_has_discontinuity_table_[cell_index]; };

protected:
    ushort criterion_solution_index_;
    uint num_cells_;

    std::vector<std::vector<uint>> cell_index_to_face_share_cell_indexes_table_;
    std::vector<bool> cell_index_to_has_discontinuity_table_;
};

class Heuristic_Discontinuity_Indicator : public Discontinuity_Indicator_Base
{
public:
    Heuristic_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
        : Discontinuity_Indicator_Base(grid, criterion_solution_index)
    {
        //for precalculation
        discrete_solution.precalculate_cell_P0_basis_values();
    };

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override;
};

class Extrapolation_Discontinuity_Indicator : public Discontinuity_Indicator_Base
{
public:
    Extrapolation_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override;

private:
    std::vector<double> cell_index_to_volume_reciprocal_table_;
    std::vector<double> cell_index_to_threshold_value_table_;
    std::vector<Euclidean_Vector> cell_index_to_QW_v_table_;
};

class Discontinuity_Indicator_Factory
{
public:
    static std::unique_ptr<Discontinuity_Indicator> make_unique(const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        if (ms::compare_icase(type_name, "Always_Flase"))
        {
            return std::make_unique<Always_False_Discontinuity_Indicator>();
        }
        else if (ms::compare_icase(type_name, "Always_True"))
        {
            return std::make_unique<Always_True_Discontinuity_Indicator>();
        }
        else if (ms::compare_icase(type_name, "Heuristic"))
        {
            return std::make_unique<Heuristic_Discontinuity_Indicator>(grid, discrete_solution, criterion_solution_index);
        }
        else if (ms::compare_icase(type_name, "Extrapolation"))
        {
            return std::make_unique<Extrapolation_Discontinuity_Indicator>(grid, discrete_solution, criterion_solution_index);
        }
        else
        {
            EXCEPTION(type_name + " is not supported discontinuity indicator type");
            return nullptr;
        }
    }

    static std::unique_ptr<Discontinuity_Indicator> make_always_false(void)
    {
        return std::make_unique<Always_False_Discontinuity_Indicator>();
    }

    static std::unique_ptr<Discontinuity_Indicator> make_always_true(void)
    {
        return std::make_unique<Always_True_Discontinuity_Indicator>();
    }

private:
    Discontinuity_Indicator_Factory(void) = delete;
};

class Shock_Indicator_Factory
{
public:
    static std::unique_ptr<Discontinuity_Indicator> make_unique(const std::string& governing_equation_name, const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        if (ms::compare_icase(governing_equation_name, "Linear_Advection"))
        {
            return Discontinuity_Indicator_Factory::make_always_false();
        }
        else if (ms::compare_icase(governing_equation_name, "Burgers"))
        {
            constexpr auto criterion_solution_index = 0;
            return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, criterion_solution_index);
        }
        else if (ms::compare_icase(governing_equation_name, "Euler"))
        {
            const auto space_dimension = grid.space_dimension();

            if (space_dimension == 2)
            {
                return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, Euler_2D::pressure_index());
            }
            else if (space_dimension == 3)
            {
                return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, Euler_3D::pressure_index());
            }
            else
            {
                EXCEPTION("not supported space_dimension");
                return nullptr;
            }
        }    
        else
        {
            EXCEPTION("not supported governing equation");
            return nullptr;
        }
    }

private:
    Shock_Indicator_Factory(void) = delete;
};

class Contact_Indicator_Factory
{
public:
    static std::unique_ptr<Discontinuity_Indicator> make_unique(const std::string& type_name, const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        constexpr auto solution_index = 0;
        return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, solution_index);
    }

private:
    Contact_Indicator_Factory(void) = delete;
};