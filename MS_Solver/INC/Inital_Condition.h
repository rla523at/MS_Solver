#pragma once
#include "Governing_Equation.h"

#include <numbers>


class IC {}; //Initial Condition


class Constant1_2D : public IC {
private:
    Constant1_2D(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;
    static constexpr ushort dimension_ = 2;

    using This_ = Constant1_2D;
    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_ = Euclidean_Vector<num_eqation_>;

public:
    static std::string name(void);

public:
    static Euclidean_Vector<1> calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Euclidean_Vector<1>> calculate_solutions(const std::vector<Space_Vector_>& space_vectors);
};


class Sine_Wave_2D : public IC 
{
private:
    Sine_Wave_2D(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;
    static constexpr ushort dimension_ = 2;
    inline static double x_wave_number_;
    inline static double y_wave_number_;

    using This_             = Sine_Wave_2D;
    using Space_Vector_     = Euclidean_Vector<dimension_>;
    using Solution_         = Euclidean_Vector<num_eqation_>;

public:
    static void initialize(const double x_wave_length, const double y_wave_length);
    static std::string name(void);

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& space_vectors);
    static std::vector<Euclidean_Vector<1>> calculate_exact_solutions(const std::vector<Space_Vector_>& space_vectors, const double end_time);
};


class Square_Wave_2D : public IC 
{
private:
    Square_Wave_2D(void) = delete;

private:
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_ = Euclidean_Vector<num_eqation_>;

public:
    static std::string name(void) { return "Square_Wave_2D"; };

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
};

class Circle_Wave_2D : public IC {
private:
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_ = Euclidean_Vector<num_eqation_>;

private:
    Circle_Wave_2D(void) = delete;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
    static std::string name(void) { return "Circle_Wave_2D"; };
};

class Gaussian_Wave_2D : public IC {
private:
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t dimension_ = 2;
    static constexpr double beta_ = 20.0;

    using This_ = Gaussian_Wave_2D;
    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_ = Euclidean_Vector<num_eqation_>;

private:
    Gaussian_Wave_2D(void) = delete;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
    static std::string name(void) { return "Gaussian_Wave_2D"; };
};


class SOD_2D : public IC 
{
private:
    SOD_2D(void) = delete;

private:
    static constexpr size_t num_eqation_ = 4;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;
    using This_         = SOD_2D;

public:
    static std::string name(void) { return "SOD"; };
    
public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);//for FVM
};


class Modified_SOD_2D : public IC 
{
private:
    Modified_SOD_2D(void) = delete;

private:
    static constexpr size_t num_eqation_ = 4;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;
    using This_         = Modified_SOD_2D;

public:
    static std::string name(void) { return "Modifid_SOD"; };

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);//for FVM
};


class Shu_Osher_2D : public IC 
{
private:
    Shu_Osher_2D(void) = delete;

private:
    static constexpr size_t num_eqation_ = 4;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = Euclidean_Vector<dimension_>;
    using Solution_ = Euclidean_Vector<num_eqation_>;
    using This_ = Shu_Osher_2D;

public:
    static std::string name(void) { return "Shu_Osher"; };

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);//for FVM
};


namespace ms {
    template <typename T>
    inline constexpr bool is_initial_condition = std::is_base_of_v<IC, T>;
}