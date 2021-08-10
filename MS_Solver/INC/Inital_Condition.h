#pragma once
#include "Governing_Equation.h"
#include "Setting.h"

#include <numbers>


class IC {}; //Initial Condition

class Sine_Wave_2D : public IC {
    static constexpr ushort num_eqation_ = 1;
    static constexpr ushort dimension_ = 2;
    static constexpr double x_wave_number_ = 2 * std::numbers::pi / static_cast<double>(X_WAVE_LENGTH);
    static constexpr double y_wave_number_ = 2 * std::numbers::pi / static_cast<double>(Y_WAVE_LENGTH);

    using This_             = Sine_Wave_2D;
    using Space_Vector_     = Euclidean_Vector<dimension_>;
    using Solution_         = Euclidean_Vector<num_eqation_>;

private:
    Sine_Wave_2D(void) = delete;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& space_vectors);
    static std::string name(void);


    template <typename Governing_Equation>
    static std::vector<Solution_> calculate_exact_solutions(const std::vector<Space_Vector_>& cell_centers, const double end_time);
};


template <>
std::vector<Sine_Wave_2D::Solution_> Sine_Wave_2D::calculate_exact_solutions<Linear_Advection_2D>(const std::vector<Space_Vector_>& cell_centers, const double end_time);


class Square_Wave_2D : public IC {
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution = Euclidean_Vector<num_eqation_>;
public:
    static std::vector<Solution> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
    static std::string name(void) { return "Square_Wave_2D"; };
    
    template <typename Governing_Equation>
    static std::vector<Solution> calculate_exact_solutions(const std::vector<Space_Vector_>& cell_centers, const double end_time);

private:
    Square_Wave_2D(void) = delete;
};


template <>
std::vector<Square_Wave_2D::Solution> Square_Wave_2D::calculate_exact_solutions<Linear_Advection_2D>(const std::vector<Space_Vector_>& cell_centers, const double end_time);


class Modified_SOD_2D : public IC {
    static constexpr size_t num_eqation_ = 4;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;
public:
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
    static std::string name(void) { return "Modifid_SOD"; };

private:
    Modified_SOD_2D(void) = delete;
};


namespace ms {
    template <typename T>
    inline constexpr bool is_initial_condition = std::is_base_of_v<IC, T>;
}