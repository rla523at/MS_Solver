#pragma once
#include "Governing_Equation.h"

#include <numbers>


class IC {}; //Initial Condition

template <double x_wave_length, double y_wave_length>
class Sine_Wave_2D : public IC {
private:
    static constexpr ushort num_eqation_ = 1;
    static constexpr ushort dimension_ = 2;
    static constexpr double x_wave_number_ = 2 * std::numbers::pi / x_wave_length;
    static constexpr double y_wave_number_ = 2 * std::numbers::pi / y_wave_length;

    using This_             = Sine_Wave_2D;
    using Space_Vector_     = Euclidean_Vector<dimension_>;
    using Solution_         = Euclidean_Vector<num_eqation_>;

private:
    Sine_Wave_2D(void) = delete;

public:
    static Euclidean_Vector<1> calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Euclidean_Vector<1>> calculate_solutions(const std::vector<Space_Vector_>& space_vectors);
    static std::string name(void);

    template <typename Governing_Equation>
    static std::vector<Euclidean_Vector<1>> calculate_exact_solutions(const std::vector<Space_Vector_>& cell_centers, const double end_time);
};


class Square_Wave_2D : public IC {
private:
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_ = Euclidean_Vector<num_eqation_>;

private:
    Square_Wave_2D(void) = delete;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
    static std::string name(void) { return "Square_Wave_2D"; };
    
    template <typename Governing_Equation>
    static std::vector<Solution_> calculate_exact_solutions(const std::vector<Space_Vector_>& cell_centers, const double end_time);
};


//template <>
//std::vector<Square_Wave_2D::Solution> Square_Wave_2D::calculate_exact_solutions<Linear_Advection_2D>(const std::vector<Space_Vector_>& cell_centers, const double end_time);


class Modified_SOD_2D : public IC {
    static constexpr size_t num_eqation_ = 4;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;
public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
    static std::string name(void) { return "Modifid_SOD"; };

private:
    Modified_SOD_2D(void) = delete;
};


namespace ms {
    template <typename T>
    inline constexpr bool is_initial_condition = std::is_base_of_v<IC, T>;
}


//template definition part
template <double x_wave_length, double y_wave_length>
Euclidean_Vector<1> Sine_Wave_2D<x_wave_length, y_wave_length>::calculate_solution(const Space_Vector_& space_vector) {
    const auto x_coord = space_vector.at(0);
    const auto y_coord = space_vector.at(1);

    return { std::sin(This_::x_wave_number_ * x_coord) * std::sin(This_::y_wave_number_ * y_coord) };
}

template <double x_wave_length, double y_wave_length>
std::vector<Euclidean_Vector<1>> Sine_Wave_2D<x_wave_length, y_wave_length>::calculate_solutions(const std::vector<Space_Vector_>& space_vectors) {
    const auto num_cell = space_vectors.size();

    std::vector<Solution_> solutions_(num_cell);
    for (size_t i = 0; i < num_cell; ++i) {
        const auto x_coord = space_vectors[i].at(0);
        const auto y_coord = space_vectors[i].at(1);

        solutions_[i] = std::sin(This_::x_wave_number_ * x_coord) * std::sin(This_::y_wave_number_ * y_coord);
    }

    return solutions_;
}

template <double x_wave_length, double y_wave_length>
std::string Sine_Wave_2D<x_wave_length, y_wave_length>::name(void) {
    return "Sine_Wave_2D";
};

template <double x_wave_length, double y_wave_length>
template <typename Governing_Equation>
std::vector<Euclidean_Vector<1>> Sine_Wave_2D<x_wave_length, y_wave_length>::calculate_exact_solutions(const std::vector<Space_Vector_>& cell_centers, const double end_time) {
    static_require(ms::is_Linear_Advection_2D<Governing_Equation>, "excat solution can be calculated when governing equation is linear advection");
    
    const auto num_cell = cell_centers.size();
    const auto [x_advection_speed, y_advection_speed] = Governing_Equation::advection_speed();

    std::vector<Solution_> exact_solutions_(num_cell);
    for (size_t i = 0; i < num_cell; ++i) {
        const auto x_coord = cell_centers[i].at(0);
        const auto y_coord = cell_centers[i].at(1);
        exact_solutions_[i] = std::sin(This_::x_wave_number_ * (x_coord - x_advection_speed * end_time)) * std::sin(This_::y_wave_number_ * (y_coord - y_advection_speed * end_time));
    }

    return exact_solutions_;
}