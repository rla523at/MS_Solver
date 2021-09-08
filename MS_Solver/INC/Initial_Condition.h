#pragma once
#include "Governing_Equation.h"

#include <numbers>


class IC {}; //Initial Condition


template <ushort space_dimension_>
class Constant1 : public IC 
{
private:
    Constant1(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_         = Constant1<space_dimension_>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Euclidean_Vector<1> calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Euclidean_Vector<1>> calculate_solutions(const std::vector<Space_Vector_>& space_vectors);

public:
    static std::string name(void);
    static constexpr ushort space_dimension(void);
};


template<ushort space_dimension_>
class Sine_Wave : public IC 
{
// sin(kx)sin(ky)sin(kz)
private:
    Sine_Wave(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_             = Sine_Wave<space_dimension_>;
    using Space_Vector_     = Euclidean_Vector<space_dimension_>;

private:
    inline static std::array<double, space_dimension_> wave_numbers_;

public:
    static void initialize(const std::array<double, space_dimension_>& wave_lengths);

public:
    static Euclidean_Vector<1> calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Euclidean_Vector<1>> calculate_solutions(const std::vector<Space_Vector_>& space_vectors);
    static std::vector<Euclidean_Vector<1>> calculate_exact_solutions(const std::vector<Space_Vector_>& space_vectors, const double end_time);

public:
    static std::string name(void);
    static constexpr ushort space_dimension(void);
};


template <ushort space_dimension_>
class Square_Wave : public IC 
{
private:
    Square_Wave(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_         = Square_Wave<space_dimension_>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Euclidean_Vector<1> calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Euclidean_Vector<1>> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);

public:
    static std::string name(void);
    static constexpr ushort space_dimension(void);
};


template <ushort space_dimension_>
class Circle_Wave : public IC 
{
private:
    Circle_Wave(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_         = Circle_Wave<space_dimension_>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Euclidean_Vector<1> calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Euclidean_Vector<1>> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);

public:
    static std::string name(void);
    static constexpr ushort space_dimension(void);
};


template <ushort space_dimension_>
class Gaussian_Wave : public IC 
{
private:
    Gaussian_Wave(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_         = Gaussian_Wave<space_dimension_>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Euclidean_Vector<1> calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Euclidean_Vector<1>> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);

public:
    static std::string name(void);
    static constexpr ushort space_dimension(void);
};


template <ushort space_dimension_>
class SOD : public IC 
{
private:
    SOD(void) = delete;

private:
    static constexpr ushort num_eqation_ = 2 + space_dimension_;

    using This_ = SOD;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;
    
public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);//for FVM

public:
    static std::string name(void);
    static constexpr ushort space_dimension(void);
};


template <ushort space_dimension_>
class Modified_SOD : public IC 
{
private:
    Modified_SOD(void) = delete;

private:
    static constexpr ushort num_eqation_ = 2 + space_dimension_;

    using This_         = Modified_SOD<space_dimension_>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);//for FVM

public:
    static std::string name(void);
    static constexpr ushort space_dimension(void);
};

template <ushort space_dimension_>
class Shu_Osher : public IC 
{
private:
    Shu_Osher(void) = delete;

private:
    static constexpr size_t num_eqation_ = 2 + space_dimension_;

    using This_         = Shu_Osher<space_dimension_>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);//for FVM

public:
    static std::string name(void);
    static constexpr ushort space_dimension(void);
};


namespace ms {
    template <typename T>
    inline constexpr bool is_initial_condition = std::is_base_of_v<IC, T>;

    template <typename T>
    inline constexpr bool is_sine_wave = std::is_same_v<T, Sine_Wave<T::space_dimension()>>;
}


//Template Definition Part
template <ushort space_dimension_>
Euclidean_Vector<1> Constant1<space_dimension_>::calculate_solution(const Space_Vector_& space_vector) {
    return { 1 };
};

template <ushort space_dimension_>
std::vector<Euclidean_Vector<1>> Constant1<space_dimension_>::calculate_solutions(const std::vector<Space_Vector_>& space_vectors) {
    std::vector<Euclidean_Vector<1>> result(space_vectors.size(), { 1 });
    return result;
};

template <ushort space_dimension_>
std::string Constant1<space_dimension_>::name(void) {
    return "Constant1";
}

template <ushort space_dimension_>
constexpr ushort Constant1<space_dimension_>::space_dimension(void) {
    return space_dimension_;
}

template <ushort space_dimension_>
void Sine_Wave<space_dimension_>::initialize(const std::array<double, space_dimension_>& wave_lengths) {
    for (ushort i = 0; i < space_dimension_; ++i)
        This_::wave_numbers_[i] = 2 * std::numbers::pi / wave_lengths[i];
}

template <ushort space_dimension_>
Euclidean_Vector<1> Sine_Wave<space_dimension_>::calculate_solution(const Space_Vector_& space_vector) {
    double solution = 1.0;

    for (ushort i = 0; i < space_dimension_; ++i) 
        solution *= std::sin(This_::wave_numbers_[i] * space_vector[i]);

    return { solution };
}

template<ushort space_dimension_>
std::vector<Euclidean_Vector<1>> Sine_Wave<space_dimension_>::calculate_solutions(const std::vector<Space_Vector_>& space_vectors) {
    const auto num_cell = space_vectors.size();

    std::vector<Euclidean_Vector<1>> solutions(num_cell);
    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(space_vectors[i]);

    return solutions;
}

template<ushort space_dimension_>
std::vector<Euclidean_Vector<1>> Sine_Wave<space_dimension_>::calculate_exact_solutions(const std::vector<Space_Vector_>& space_vectors, const double end_time) {
    const auto advection_speed = Linear_Advection<space_dimension_>::advection_speed();

    const auto num_cell = space_vectors.size();
    std::vector<Euclidean_Vector<1>> exact_solutions(num_cell);

    for (size_t i = 0; i < num_cell; ++i) {
        const auto& cell_center = space_vectors[i];
        
        double exact_solution = 1.0;

        for (ushort j = 0; j < space_dimension_; ++j)
            exact_solution *= std::sin(This_::wave_numbers_[j] * (cell_center[j] - advection_speed[j] * end_time));

        exact_solutions[i] = exact_solution;
    }

    return exact_solutions;
}

template <ushort space_dimension_>
std::string Sine_Wave<space_dimension_>::name(void) {
    return "Sine_Wave";
};

template <ushort space_dimension_>
constexpr ushort Sine_Wave<space_dimension_>::space_dimension(void) {
    return space_dimension_;
}

template <ushort space_dimension_>
Euclidean_Vector<1> Square_Wave<space_dimension_>::calculate_solution(const Space_Vector_& space_vector) {
    for (ushort i = 0; i < space_dimension_; ++i) {
        if (space_vector[i] < 0.25 || 0.75 < space_vector[i])
            return { 0 };
    }
    return { 1 };    
}

template <ushort space_dimension_>
std::vector<Euclidean_Vector<1>> Square_Wave<space_dimension_>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();
    std::vector<Solution_> solutions(num_cell);

    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template <ushort space_dimension_>
std::string Square_Wave<space_dimension_>::name(void) {
    return "Square_Wave";
};

template <ushort space_dimension_>
constexpr ushort Square_Wave<space_dimension_>::space_dimension(void) {
    return space_dimension_;
}

template <ushort space_dimension_>
Euclidean_Vector<1> Circle_Wave<space_dimension_>::calculate_solution(const Space_Vector_& space_vector) {
    double temp = 0.0;
    for (ushort i = 0; i < space_dimension_; ++i) 
        temp += (space_vector[i] - 0.5) * (space_vector[i] - 0.5);

    if (temp <= 0.25 * 0.25)
        return { 1 };
    else
        return { 0 };
}

template <ushort space_dimension_>
std::vector<Euclidean_Vector<1>> Circle_Wave<space_dimension_>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();
    std::vector<Solution_> solutions(num_cell);

    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template <ushort space_dimension_>
std::string Circle_Wave<space_dimension_>::name(void) {
    return "Circle_Wave";
};

template <ushort space_dimension_>
constexpr ushort Circle_Wave<space_dimension_>::space_dimension(void) {
    return space_dimension_;
}

template <ushort space_dimension_>
Euclidean_Vector<1> Gaussian_Wave<space_dimension_>::calculate_solution(const Space_Vector_& space_vector) {
    constexpr double beta = 20.0;

    double temp = 0.0;
    for (ushort i = 0; i < space_dimension_; ++i)
        temp += (space_vector[i] - 0.5) * (space_vector[i] - 0.5);

    return { std::exp(-beta * temp) };
}

template <ushort space_dimension_>
std::vector<Euclidean_Vector<1>> Gaussian_Wave<space_dimension_>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();
    std::vector<Solution_> solutions(num_cell);

    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template <ushort space_dimension_>
std::string Gaussian_Wave<space_dimension_>::name(void) {
    return "Gaussian_Wave";
};

template <ushort space_dimension_>
constexpr ushort Gaussian_Wave<space_dimension_>::space_dimension(void) {
    return space_dimension_;
}

template <ushort space_dimension_>
SOD<space_dimension_>::Solution_ SOD<space_dimension_>::calculate_solution(const Space_Vector_& space_vector) {
    constexpr auto gamma = 1.4;
    constexpr auto c = 1 / (gamma - 1);
    constexpr auto discontinuity_location = 0.5;

    const auto x_coordinate = space_vector.at(0);

    if constexpr (space_dimension_ == 2) {
        if (x_coordinate <= discontinuity_location) {
            constexpr auto rho = 1.0;
            constexpr auto u = 0.0;
            constexpr auto v = 0.0;
            constexpr auto p = 1.0;

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

            return { rho, rhou, rhov, rhoE };
        }
        else {
            constexpr auto rho = 0.125;
            constexpr auto u = 0.0;
            constexpr auto v = 0.0;
            constexpr auto p = 0.1;

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

            return { rho, rhou, rhov, rhoE };
        }
    }
    else if constexpr (space_dimension_ ==3) {
        if (x_coordinate <= discontinuity_location) {
            constexpr auto rho = 1.0;
            constexpr auto u = 0.0;
            constexpr auto v = 0.0;
            constexpr auto w = 0.0;
            constexpr auto p = 1.0;

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhow = rho * w;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v + rhow * w);

            return { rho, rhou, rhov, rhow, rhoE };
        }
        else {
            constexpr auto rho = 0.125;
            constexpr auto u = 0.0;
            constexpr auto v = 0.0;
            constexpr auto w = 0.0;
            constexpr auto p = 0.1;

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhow = rho * w;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v + rhow * w);

            return { rho, rhou, rhov, rhow, rhoE };
        }
    }
    else {
        throw std::runtime_error("not supported space dimension");
        return {};
    }
}

template <ushort space_dimension_>
std::vector<typename SOD<space_dimension_>::Solution_> SOD<space_dimension_>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();

    std::vector<Solution_> solutions(num_cell);
    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template <ushort space_dimension_>
std::string SOD<space_dimension_>::name(void) {
    return "SOD";
};

template <ushort space_dimension_>
constexpr ushort SOD<space_dimension_>::space_dimension(void) {
    return space_dimension_;
}

template <ushort space_dimension_>
Modified_SOD<space_dimension_>::Solution_ Modified_SOD<space_dimension_>::calculate_solution(const Space_Vector_& space_vector) {
    constexpr auto gamma = 1.4;
    constexpr auto c = 1.0 / (gamma - 1.0);
    constexpr auto discontinuity_location = 0.3;

    const auto x_coordinate = space_vector.at(0);

    if constexpr (space_dimension_ == 2) {
        if (x_coordinate <= discontinuity_location) {
            constexpr auto rho = 1.0;
            constexpr auto u = 0.75;
            constexpr auto v = 0.0;
            constexpr auto p = 1.0;

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

            return { rho, rhou, rhov, rhoE };
        }
        else {
            constexpr auto rho = 0.125;
            constexpr auto u = 0.0;
            constexpr auto v = 0.0;
            constexpr auto p = 0.1;

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

            return { rho, rhou, rhov, rhoE };
        }
    }
    else if constexpr (space_dimension_ == 3) {
        if (x_coordinate <= discontinuity_location) {
            constexpr auto rho = 1.0;
            constexpr auto u = 0.75;
            constexpr auto v = 0.0;
            constexpr auto w = 0.0;
            constexpr auto p = 1.0;

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhow = rho * w;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v + rhow * w);

            return { rho, rhou, rhov, rhow, rhoE };
        }
        else {
            constexpr auto rho = 0.125;
            constexpr auto u = 0.0;
            constexpr auto v = 0.0;
            constexpr auto w = 0.0;
            constexpr auto p = 0.1;

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhow = rho * w;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v + rhow * w);

            return { rho, rhou, rhov, rhow, rhoE };
        }
    }
    else {
        throw std::runtime_error("not supported space dimension");
        return {};
    }
}

template <ushort space_dimension_>
std::vector<typename Modified_SOD<space_dimension_>::Solution_> Modified_SOD<space_dimension_>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();

    std::vector<Solution_> solutions(num_cell);
    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template <ushort space_dimension_>
std::string Modified_SOD<space_dimension_>::name(void) {
    return "Modified_SOD";
};

template <ushort space_dimension_>
constexpr ushort Modified_SOD<space_dimension_>::space_dimension(void) {
    return space_dimension_;
}

template <ushort space_dimension_>
Shu_Osher<space_dimension_>::Solution_ Shu_Osher<space_dimension_>::calculate_solution(const Space_Vector_& space_vector) {
    constexpr auto gamma = 1.4;
    constexpr auto c = 1 / (gamma - 1);
    constexpr auto discontinuity_location = -4.0;

    const auto x_coordinate = space_vector.at(0);

    if constexpr (space_dimension_ == 2) {
        if (x_coordinate <= discontinuity_location) {
            constexpr auto rho = 3.857143;
            constexpr auto u = 2.629369;    //rhou = 10.1418522328      
            constexpr auto v = 0.0;
            constexpr auto p = 10.333333;   //rhoE = 39.1666684317   

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

            return { rho, rhou, rhov, rhoE };
        }
        else {
            const auto rho = 1 + 0.2 * std::sin(5 * x_coordinate);
            constexpr auto u = 0.0;
            constexpr auto v = 0.0;
            constexpr auto p = 1.0;

            const auto rhou = rho * u;
            const auto rhov = rho * v;
            const auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

            return { rho, rhou, rhov, rhoE };
        }
    }
    else if constexpr (space_dimension_ == 3) {
        if (x_coordinate <= discontinuity_location) {
            constexpr auto rho = 3.857143;
            constexpr auto u = 2.629369;    //rhou = 10.1418522328 
            constexpr auto v = 0.0;
            constexpr auto w = 0.0;
            constexpr auto p = 10.333333;   //rhoE = 39.1666684317

            constexpr auto rhou = rho * u;
            constexpr auto rhov = rho * v;
            constexpr auto rhow = rho * w;
            constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v + rhow * w);

            return { rho, rhou, rhov, rhow, rhoE };
        }
        else {
            const auto rho = 1 + 0.2 * std::sin(5 * x_coordinate);
            constexpr auto u = 0.0;
            constexpr auto v = 0.0;
            constexpr auto w = 0.0;
            constexpr auto p = 1.0;

            const auto rhou = rho * u;
            const auto rhov = rho * v;
            const auto rhow = rho * w;
            const auto rhoE = p * c + 0.5 * (rhou * u + rhov * v + rhow * w);

            return { rho, rhou, rhov, rhow, rhoE };
        }
    }
    else {
        throw std::runtime_error("not supported space dimension");
        return {};
    }
}

template <ushort space_dimension_>
std::vector<typename Shu_Osher<space_dimension_>::Solution_> Shu_Osher<space_dimension_>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();

    std::vector<Solution_> solutions(num_cell);
    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template <ushort space_dimension_>
std::string Shu_Osher<space_dimension_>::name(void) {
    return "Shu_Osher";
};

template <ushort space_dimension_>
constexpr ushort Shu_Osher<space_dimension_>::space_dimension(void) {
    return space_dimension_;
}