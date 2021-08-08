#pragma once
#include "Grid_Builder.h"

//FVM이면 공통으로 사용하는 variable & method
template <ushort space_dimension>
class Cells_FVM
{
private:
    using This_         = Cells_FVM<space_dimension>;
    using SpaceVector_  = Euclidean_Vector<space_dimension>;

protected:
    inline static std::vector<SpaceVector_> centers_;
    inline static std::vector<double> volumes_;
    inline static std::vector<std::array<double, space_dimension>> coordinate_projected_volumes_;
    inline static std::vector<double> residual_scale_factors_;

private:
    Cells_FVM(void) = delete;

public:
    static void initialize(const Grid<space_dimension>& grid);
    static double calculate_time_step(const std::vector<std::array<double, space_dimension>>& coordinate_projected_maximum_lambdas, const double cfl);

    template <typename Residual, typename Solution>
    static void calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions);

    template <typename Initial_Condtion>
    static auto calculate_initial_solutions(void);

    template <typename Initial_Condition, typename Governing_Equation, typename Solution>
    static void estimate_error(const std::vector<Solution>& computed_solution, const double time);
};


//template definition
template <ushort space_dimension>
void Cells_FVM<space_dimension>::initialize(const Grid<space_dimension>& grid) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    This_::centers_.reserve(num_cell);
    This_::volumes_.reserve(num_cell);
    This_::coordinate_projected_volumes_.reserve(num_cell);
    This_::residual_scale_factors_.reserve(num_cell);

    for (const auto& cell_elemnt : cell_elements) {
        const auto& geometry = cell_elemnt.geometry_;
        const auto volume = geometry.volume();

        This_::centers_.push_back(geometry.center_node());
        This_::volumes_.push_back(volume);
        This_::coordinate_projected_volumes_.push_back(geometry.coordinate_projected_volume());
        This_::residual_scale_factors_.push_back(1.0 / volume);
    }

    Log::content_ << std::left << std::setw(50) << "@ Cells FVM precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
};

template <ushort space_dimension>
double Cells_FVM<space_dimension>::calculate_time_step(const std::vector<std::array<double, space_dimension>>& coordinate_projected_maximum_lambdas, const double cfl) {
    const auto num_cell = coordinate_projected_maximum_lambdas.size();

    std::vector<double> local_time_step(num_cell);
    for (size_t i = 0; i < num_cell; ++i) {
        const auto [x_projected_volume, y_projected_volume] = This_::coordinate_projected_volumes_[i];
        const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda] = coordinate_projected_maximum_lambdas[i];

        const auto x_radii = x_projected_volume * x_projeced_maximum_lambda;
        const auto y_radii = y_projected_volume * y_projeced_maximum_lambda;

        local_time_step[i] = cfl * This_::volumes_[i] / (x_radii + y_radii);
    }

    return *std::min_element(local_time_step.begin(), local_time_step.end());
}

template <ushort space_dimension>
template <typename Residual, typename Solution>
void Cells_FVM<space_dimension>::calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions){
    const auto num_cell = RHS.size();

    for (size_t i = 0; i < num_cell; ++i)
        RHS[i] *= This_::residual_scale_factors_[i];
}

template <ushort space_dimension>
template <typename Initial_Condtion>
auto Cells_FVM<space_dimension>::calculate_initial_solutions(void) {
    return Initial_Condtion::calculate_solutions(This_::centers_);
}

template <ushort space_dimension>
template <typename Initial_Condition, typename Governing_Equation, typename Solution>
void Cells_FVM<space_dimension>::estimate_error(const std::vector<Solution>& computed_solutions, const double time) {  

    Log::content_ << "================================================================================\n";
    Log::content_ << "\t\t\t\t Error Anlysis\n";
    Log::content_ << "================================================================================\n";

    if constexpr (std::is_same_v<Governing_Equation, Linear_Advection_2D>) {
        double global_L1_error = 0.0;
        double global_L2_error = 0.0;
        double global_Linf_error = 0.0;

        const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Governing_Equation>(This_::centers_, time);
        const auto num_solutions = computed_solutions.size();

        for (size_t i = 0; i < num_solutions; ++i) {
            const auto local_error = (exact_solutions[i] - computed_solutions[i]).L1_norm();
            global_L1_error += local_error;
            global_L2_error += local_error * local_error;
            global_Linf_error = max(global_Linf_error, local_error);
        }

        global_L1_error = global_L1_error / num_solutions;
        global_L2_error = global_L2_error / num_solutions;

        global_L2_error = std::sqrt(global_L2_error);

        Log::content_ << "L1 error \t\tL2 error \t\tLinf error \n";
        Log::content_ << ms::double_to_string(global_L1_error) << "\t" << ms::double_to_string(global_L2_error) << "\t" << ms::double_to_string(global_Linf_error) << "\n\n";

    }
    else
        Log::content_ << Governing_Equation::name() << " does not provide error analysis result.\n\n";

    Log::print();
}
