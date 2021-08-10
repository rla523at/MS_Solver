#pragma once
#include "Grid_Builder.h"
#include "Governing_Equation.h"

//FVM이면 공통으로 사용하는 variable & method
template <typename Governing_Equation>
class Cells_FVM
{
private:
    static constexpr ushort space_dimension_    = Governing_Equation::space_dimension();
    static constexpr ushort num_equation_       = Governing_Equation::num_equation();

    using This_         = Cells_FVM<Governing_Equation>;
    using Space_Vector_  = Euclidean_Vector<space_dimension_>;    

protected:
    inline static std::vector<Space_Vector_> centers_;
    inline static std::vector<double> volumes_;
    inline static std::vector<std::array<double, space_dimension_>> coordinate_projected_volumes_;
    inline static std::vector<double> residual_scale_factors_;

private:
    Cells_FVM(void) = delete;

public:
    static void initialize(const Grid<space_dimension_>& grid);
    static double calculate_time_step(const std::vector<Euclidean_Vector<num_equation_>>& solutions, const double cfl);
    static void calculate_RHS(std::vector<Euclidean_Vector<num_equation_>>& RHS, const std::vector<Euclidean_Vector<num_equation_>>& solutions);

    template <typename Initial_Condition>
    static auto calculate_initial_solutions(void);

    template <typename Initial_Condition>
    static void estimate_error(const std::vector<Euclidean_Vector<num_equation_>>& computed_solutions, const double time);
};


//template definition
template <typename Governing_Equation>
void Cells_FVM<Governing_Equation>::initialize(const Grid<space_dimension_>& grid) {
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

template <typename Governing_Equation>
double Cells_FVM<Governing_Equation>::calculate_time_step(const std::vector<Euclidean_Vector<num_equation_>>& solutions, const double cfl) {
    const auto projected_maximum_lambdas = Governing_Equation::calculate_coordinate_projected_maximum_lambdas(solutions);
    const auto num_cell = projected_maximum_lambdas.size();

    std::vector<double> local_time_step(num_cell);
    for (size_t i = 0; i < num_cell; ++i) {
        const auto [x_projected_volume, y_projected_volume] = This_::coordinate_projected_volumes_[i];
        const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda] = projected_maximum_lambdas[i];

        const auto x_radii = x_projected_volume * x_projeced_maximum_lambda;
        const auto y_radii = y_projected_volume * y_projeced_maximum_lambda;

        if (!std::isnormal(x_radii) || !std::isnormal(y_radii))
            std::cout << "not normal!";

        local_time_step[i] = cfl * This_::volumes_[i] / (x_radii + y_radii);
    }

    return *std::min_element(local_time_step.begin(), local_time_step.end());
}

template <typename Governing_Equation>
void Cells_FVM<Governing_Equation>::calculate_RHS(std::vector<Euclidean_Vector<num_equation_>>& RHS, const std::vector<Euclidean_Vector<num_equation_>>& solutions) {
    const auto num_cell = RHS.size();

    for (size_t i = 0; i < num_cell; ++i)
        RHS[i] *= This_::residual_scale_factors_[i];
}

template <typename Governing_Equation>
template <typename Initial_Condition>
auto Cells_FVM<Governing_Equation>::calculate_initial_solutions(void) {
    return Initial_Condition::calculate_solutions(This_::centers_);
}

template <typename Governing_Equation>
template <typename Initial_Condition>
void Cells_FVM<Governing_Equation>::estimate_error(const std::vector<Euclidean_Vector<num_equation_>>& computed_solutions, const double time) {

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
