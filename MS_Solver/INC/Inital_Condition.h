#pragma once
#include "Governing_Equation.h"

#include <numbers>


class IC {}; //Initial Condition


template<ushort space_dimension>
class Constant1 : public IC 
{
private:
    Constant1(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_         = Constant1<space_dimension>;
    using Space_Vector_ = Euclidean_Vector<space_dimension>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Euclidean_Vector<1> calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Euclidean_Vector<1>> calculate_solutions(const std::vector<Space_Vector_>& space_vectors);

public:
    static std::string name(void);
};


template<ushort space_dimension>
class Sine_Wave : public IC 
{
// sin(kx)sin(ky)sin(kz)
private:
    Sine_Wave(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_             = Sine_Wave<space_dimension>;
    using Space_Vector_     = Euclidean_Vector<space_dimension>;
    using Solution_         = Euclidean_Vector<num_eqation_>;

private:
    inline static std::array<double, space_dimension> wave_numbers_;

public:
    static void initialize(const std::array<double, space_dimension>& wave_lengths);

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& space_vectors);
    static std::vector<Euclidean_Vector<1>> calculate_exact_solutions(const std::vector<Space_Vector_>& space_vectors, const double end_time);

public:
    static std::string name(void);
};


template <ushort space_dimension>
class Square_Wave : public IC 
{
private:
    Square_Wave(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_         = Square_Wave<space_dimension>;
    using Space_Vector_ = Euclidean_Vector<space_dimension>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);

public:
    static std::string name(void);
};


template <ushort space_dimension>
class Circle_Wave : public IC 
{
private:
    Circle_Wave(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_         = Circle_Wave<space_dimension>;
    using Space_Vector_ = Euclidean_Vector<space_dimension>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);

public:
    static std::string name(void);
};


template <ushort space_dimension>
class Gaussian_Wave : public IC 
{
private:
    Gaussian_Wave(void) = delete;

private:
    static constexpr ushort num_eqation_ = 1;

    using This_         = Gaussian_Wave<space_dimension>;
    using Space_Vector_ = Euclidean_Vector<space_dimension>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);

public:
    static std::string name(void);
};


template <ushort space_dimension>
class SOD : public IC 
{
private:
    SOD(void) = delete;

private:
    static constexpr ushort num_eqation_ = 2 + space_dimension;

    using This_ = SOD;
    using Space_Vector_ = Euclidean_Vector<space_dimension>;
    using Solution_     = Euclidean_Vector<num_eqation_>;
    
public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);//for FVM

public:
    static std::string name(void);
};


template <ushort space_dimension>
class Modified_SOD : public IC 
{
private:
    Modified_SOD(void) = delete;

private:
    static constexpr ushort num_eqation_ = 2 + space_dimension;

    using This_         = Modified_SOD<space_dimension>;
    using Space_Vector_ = Euclidean_Vector<space_dimension>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);//for FVM

public:
    static std::string name(void);
};

template <ushort space_dimension>
class Shu_Osher : public IC 
{
private:
    Shu_Osher(void) = delete;

private:
    static constexpr size_t num_eqation_ = 2 + space_dimension;

    using This_         = Shu_Osher<space_dimension>;
    using Space_Vector_ = Euclidean_Vector<space_dimension>;
    using Solution_     = Euclidean_Vector<num_eqation_>;

public:
    static Solution_ calculate_solution(const Space_Vector_& space_vector);
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);//for FVM

public:
    static std::string name(void);
};


namespace ms {
    template <typename T>
    inline constexpr bool is_initial_condition = std::is_base_of_v<IC, T>;
}


//Template Definition Part
template<ushort space_dimension>
Euclidean_Vector<1> Constant1<space_dimension>::calculate_solution(const Space_Vector_& space_vector) {
    return { 1 };
};

template<ushort space_dimension>
std::vector<Euclidean_Vector<1>> Constant1<space_dimension>::calculate_solutions(const std::vector<Space_Vector_>& space_vectors) {
    std::vector<Euclidean_Vector<1>> result(space_vectors.size(), { 1 });
    return result;
};

template<ushort space_dimension>
std::string Constant1<space_dimension>::name(void) {
    return "Constant" + std::to_string(space_dimension) + "D";
}

template<ushort space_dimension>
void Sine_Wave<space_dimension>::initialize(const std::array<double, space_dimension>& wave_lengths) {
    for (ushort i = 0; i < space_dimension; ++i)
        This_::wave_numbers_[i] = 2 * std::numbers::pi / wave_lengths[i];
}

template<ushort space_dimension>
Euclidean_Vector<1> Sine_Wave<space_dimension>::calculate_solution(const Space_Vector_& space_vector) {
    double solution = 1.0;

    for (ushort i = 0; i < space_dimension; ++i) 
        solution *= std::sin(This_::wave_numbers_[i] * space_vector[i]);

    return { solution };
}

template<ushort space_dimension>
std::vector<Euclidean_Vector<1>> Sine_Wave<space_dimension>::calculate_solutions(const std::vector<Space_Vector_>& space_vectors) {
    const auto num_cell = space_vectors.size();

    std::vector<Solution_> solutions(num_cell);
    for (size_t i = 0; i < num_cell; ++i) 
        solutions[i]= This_::calculate_solution(space_vectors[i]);

    return solutions;
}

template<ushort space_dimension>
std::vector<Euclidean_Vector<1>> Sine_Wave<space_dimension>::calculate_exact_solutions(const std::vector<Space_Vector_>& cell_centers, const double end_time) {
    const auto advection_speed = Linear_Advection<space_dimension>::advection_speed();

    const auto num_cell = cell_centers.size();
    std::vector<Solution_> exact_solutions(num_cell);

    for (size_t i = 0; i < num_cell; ++i) {
        const auto& cell_center = cell_centers[i];
        
        double exact_solution = 1.0;

        for (ushort j = 0; j < space_dimension; ++j)
            exact_solution *= std::sin(This_::wave_numbers_[j] * (cell_center[j] - advection_speed[j] * end_time));

        exact_solutions[i] = exact_solution;
    }

    return exact_solutions;
}

template<ushort space_dimension>
std::string Sine_Wave<space_dimension>::name(void) {
    return "Sine_Wave" + std::to_string(space_dimension) + "D";
};

template <ushort space_dimension>
Euclidean_Vector<1> Square_Wave<space_dimension>::calculate_solution(const Space_Vector_& space_vector) {
    for (ushort i = 0; i < space_dimension; ++i) {
        if (space_vector[i] < 0.25 || 0.75 < space_vector[i])
            return { 0 };
    }
    return { 1 };    
}

template <ushort space_dimension>
std::vector<Euclidean_Vector<1>> Square_Wave<space_dimension>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();
    std::vector<Solution_> solutions(num_cell);

    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template<ushort space_dimension>
std::string Square_Wave<space_dimension>::name(void) {
    return "Square_Wave" + std::to_string(space_dimension) + "D";
};

template <ushort space_dimension>
Euclidean_Vector<1> Circle_Wave<space_dimension>::calculate_solution(const Space_Vector_& space_vector) {
    double temp = 0.0;
    for (ushort i = 0; i < space_dimension; ++i) 
        temp += (space_vector[i] - 0.5) * (space_vector[i] - 0.5);

    if (temp <= 0.25 * 0.25)
        return { 1 };
    else
        return { 0 };
}

template <ushort space_dimension>
std::vector<Euclidean_Vector<1>> Circle_Wave<space_dimension>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();
    std::vector<Solution_> solutions(num_cell);

    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template<ushort space_dimension>
std::string Circle_Wave<space_dimension>::name(void) {
    return "Circle_Wave" + std::to_string(space_dimension) + "D";
};

template <ushort space_dimension>
Euclidean_Vector<1> Gaussian_Wave<space_dimension>::calculate_solution(const Space_Vector_& space_vector) {
    constexpr double beta = 20.0;

    double temp = 0.0;
    for (ushort i = 0; i < space_dimension; ++i)
        temp += (space_vector[i] - 0.5) * (space_vector[i] - 0.5);

    return { std::exp(-beta * temp) };
}

template <ushort space_dimension>
std::vector<Euclidean_Vector<1>> Gaussian_Wave<space_dimension>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();
    std::vector<Solution_> solutions(num_cell);

    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template <ushort space_dimension>
std::string Gaussian_Wave<space_dimension>::name(void) {
    return "Gaussian_Wave" + std::to_string(space_dimension) + "D";
};

template <ushort space_dimension>
SOD<space_dimension>::Solution_ SOD<space_dimension>::calculate_solution(const Space_Vector_& space_vector) {
    constexpr auto gamma = 1.4;
    constexpr auto c = 1 / (1.4 - 1);
    constexpr auto discontinuity_location = 0.5;

    const auto x_coordinate = space_vector.at(0);

    if constexpr (space_dimension == 2) {
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
    else if constexpr (space_dimension ==3) {
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
}

template <ushort space_dimension>
std::vector<typename SOD<space_dimension>::Solution_> SOD<space_dimension>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();

    std::vector<Solution_> solutions(num_cell);
    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template <ushort space_dimension>
std::string SOD<space_dimension>::name(void) {
    return "SOD" + std::to_string(space_dimension) + "D";
};

template <ushort space_dimension>
Modified_SOD<space_dimension>::Solution_ Modified_SOD<space_dimension>::calculate_solution(const Space_Vector_& space_vector) {
    constexpr auto gamma = 1.4;
    constexpr auto c = 1.0 / (gamma - 1.0);
    constexpr auto discontinuity_location = 0.3;

    const auto x_coordinate = space_vector.at(0);

    if constexpr (space_dimension == 2) {
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
    else if constexpr (space_dimension == 3) {
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
}

template <ushort space_dimension>
std::vector<typename Modified_SOD<space_dimension>::Solution_> Modified_SOD<space_dimension>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();

    std::vector<Solution_> solutions(num_cell);
    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}

template <ushort space_dimension>
std::string Modified_SOD<space_dimension>::name(void) {
    return "Modified_SOD" + std::to_string(space_dimension) + "D";
};


template <ushort space_dimension>
Shu_Osher<space_dimension>::Solution_ Shu_Osher<space_dimension>::calculate_solution(const Space_Vector_& space_vector) {
    constexpr auto gamma = 1.4;
    constexpr auto c = 1 / (gamma - 1);
    constexpr auto discontinuity_location = -4.0;

    const auto x_coordinate = space_vector.at(0);

    if (x_coordinate < discontinuity_location) {
        constexpr auto rho = 3.857143;
        constexpr auto u = 2.629369;	
        constexpr auto v = 0.0;			
        constexpr auto p = 10.333333;	

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

    if constexpr (space_dimension == 2) {
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
    else if constexpr (space_dimension == 3) {
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
}

template <ushort space_dimension>
std::vector<typename Shu_Osher<space_dimension>::Solution_> Shu_Osher<space_dimension>::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
    const auto num_cell = cell_centers.size();

    std::vector<Solution_> solutions(num_cell);
    for (size_t i = 0; i < num_cell; ++i)
        solutions[i] = This_::calculate_solution(cell_centers[i]);

    return solutions;
}
