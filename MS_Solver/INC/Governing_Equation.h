#pragma once
#include "Matrix.h"


using ushort = unsigned short;

class GE {}; // Governing Equation


class SCL_2D : public GE    // 2D Scalar Conservation Law 
{
protected:
    static constexpr ushort num_equation_ = 1;
    static constexpr ushort space_dimension_ = 2;

private:
    SCL_2D(void) = delete;

public:
    using Space_Vector_  = Euclidean_Vector<space_dimension_>;
    using Solution_      = Euclidean_Vector<num_equation_>;
    using Physical_Flux_ = Matrix<num_equation_, space_dimension_>;

    static constexpr ushort space_dimension(void) { return space_dimension_; };
    static constexpr ushort num_equation(void) { return num_equation_; };
};


template <double x_advection_speed, double y_advection_speed>
class Linear_Advection_2D : public SCL_2D
{
private:
    using This_ = Linear_Advection_2D;

private:
    static constexpr std::array<double, space_dimension_> advection_speeds_ = { x_advection_speed, y_advection_speed };

private:
    Linear_Advection_2D(void) = delete;

public:
    static constexpr std::array<double, space_dimension_> advection_speed(void) { return advection_speeds_; };
    static auto physical_flux(const Solution_& solution);
    static auto physical_fluxes(const std::vector<Solution_>& solutions);
    static auto calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions);
    static double inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector);
    static std::string name(void) { return "Linear_Advection_2D"; }; 
};


class Burgers_2D : public SCL_2D
{
private:
    using This_ = Burgers_2D;

private:
    Burgers_2D(void) = delete;

public:
    static Physical_Flux_ physical_flux(const Solution_& solution);
    static std::vector<Physical_Flux_> physical_fluxes(const std::vector<Solution_>& solutions);
    static std::vector<std::array<double, space_dimension_>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions);
    static double inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector);
    static std::string name(void) { return "Burgers_2D"; };    
};


class Euler_2D : public GE
{
private:
    static constexpr ushort num_equation_ = 4;
    static constexpr ushort space_dimension_ = 2;

    using This_                 = Euler_2D;
    using Space_Vector_         = Euclidean_Vector<space_dimension_>;
    using Solution_             = Euclidean_Vector<num_equation_>;
    using Physical_Flux_        = Matrix<num_equation_, space_dimension_>;

private:
    Euler_2D(void) = delete;

public:
    static Solution_ conservative_to_primitive(const Solution_& conservative_variable);
    static std::vector<std::array<double, space_dimension_>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& conservative_variables);
    static Physical_Flux_ physical_flux(const Solution_& conservative_variable);
    static Physical_Flux_ physical_flux(const Solution_& conservative_variable, const Solution_& primitivie_variable);
    static std::vector<Physical_Flux_> physical_fluxes(const std::vector<Solution_>& conservative_variables, const std::vector<Solution_>& primitive_variables);
    static double inner_face_maximum_lambda(const Solution_& oc_primitive_variable, const Solution_& nc_primitive_variable, const Space_Vector_& nomal_vector);
    
    static constexpr ushort space_dimension(void) { return space_dimension_; };
    static constexpr ushort num_equation(void) { return num_equation_; };
    static std::string name(void) { return "Euler_2D"; };
};

namespace ms {
    template <typename T>
    inline constexpr bool is_governing_equation = std::is_base_of_v<GE, T>;
    template <typename T>
    inline constexpr bool is_SCL_2D = std::is_base_of_v<SCL_2D, T>;
    template <typename T>
    inline constexpr bool is_Linear_Advection_2D = std::is_same_v<T, Linear_Advection_2D<T::advection_speed()[0], T::advection_speed()[1]>>;
}


//template definition part
template <double x_advection_speed, double y_advection_speed>
auto Linear_Advection_2D<x_advection_speed, y_advection_speed>::physical_flux(const Solution_& solution) {
    const auto sol = solution.at(0);	//scalar

    Physical_Flux_ physical_flux = { x_advection_speed * sol , y_advection_speed * sol };
    return physical_flux;
}

template <double x_advection_speed, double y_advection_speed>
auto Linear_Advection_2D<x_advection_speed, y_advection_speed>::physical_fluxes(const std::vector<Solution_>& solutions) {
    //static size_t num_solution = solutions.size();
    const size_t num_solution = solutions.size();

    std::vector<Physical_Flux_> physical_fluxes(num_solution);
    for (size_t i = 0; i < num_solution; ++i) {
        const auto sol = solutions[i].at(0);	//scalar
        physical_fluxes[i] = { x_advection_speed * sol , y_advection_speed * sol };
    }

    return physical_fluxes;
}

template <double x_advection_speed, double y_advection_speed>
auto Linear_Advection_2D<x_advection_speed, y_advection_speed>::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions) {
    //static size_t num_solution = solutions.size();
    const size_t num_solution = solutions.size();

    static double absolute_x_advection_speed = std::abs(x_advection_speed);
    static double absolute_y_advection_speed = std::abs(y_advection_speed);

    std::vector<std::array<double, This_::space_dimension_>> projected_maximum_lambdas(num_solution, { absolute_x_advection_speed,absolute_y_advection_speed });
    return projected_maximum_lambdas;
}

template <double x_advection_speed, double y_advection_speed>
double Linear_Advection_2D<x_advection_speed, y_advection_speed>::inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector) {
    return std::abs(nomal_vector.inner_product(This_::advection_speeds_));
}