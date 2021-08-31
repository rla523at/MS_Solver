#pragma once
#include "Matrix.h"


class GE {}; // Governing Equation


template <ushort space_dimension_>
class SCL : public GE    //Scalar Conservation Law 
{
private:
    SCL(void) = delete;

protected:
    static constexpr ushort num_equation_ = 1;

public:
    static constexpr ushort space_dimension(void) { return space_dimension_; };
    static constexpr ushort num_equation(void) { return num_equation_; };
};


template <ushort space_dimension>
class Linear_Advection : public SCL<space_dimension>
{
private:
    Linear_Advection(void) = delete;

private:
    using This_             = Linear_Advection<space_dimension>;
    using Space_Vector_     = Euclidean_Vector<space_dimension>;
    using Solution_         = Euclidean_Vector<This_::num_equation_>;
    using Physical_Flux_    = Matrix<This_::num_equation_, space_dimension>;

private:
    inline static Euclidean_Vector<space_dimension> advection_speeds_;

public:
    static void initialize(const Euclidean_Vector<space_dimension>& advection_speed) { advection_speeds_ = advection_speed; };

public:
    static Euclidean_Vector<space_dimension> advection_speed(void) { return advection_speeds_; };
    static auto physical_flux(const This_::Solution_& solution);
    static auto physical_fluxes(const std::vector<This_::Solution_>& solutions);
    static std::vector<std::array<double, space_dimension>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions);
    static double inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector);

public:
    static std::string name(void) { return "Linear_Advection_2D"; }; 
};


template <ushort space_dimension>
class Burgers : public SCL<space_dimension>
{
private:
    Burgers(void) = delete;

private:
    using This_             = Burgers<space_dimension>;
    using Space_Vector_     = Euclidean_Vector<space_dimension>;
    using Solution_         = Euclidean_Vector<This_::num_equation_>;
    using Physical_Flux_    = Matrix<This_::num_equation_, space_dimension>;

public:
    static auto physical_flux(const Solution_& solution);
    static auto physical_fluxes(const std::vector<Solution_>& solutions);
    static std::vector<std::array<double, space_dimension>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions);
    static double inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector);

public:
    static std::string name(void) { return "Burgers<space_dimension>"; };    
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
    inline constexpr bool is_SCL = std::is_base_of_v<SCL, T>;
}


//Template Definition Part

template <ushort space_dimension>
auto Linear_Advection<space_dimension>::physical_flux(const Solution_& solution) {
    const auto sol = solution.at(0);	//scalar

    std::array<double, space_dimension> flux_values;
    for (ushort i = 0; i < space_dimension; ++i) 
        flux_values[i] = sol * This_::advection_speeds_[i];
    
    return Matrix<1, space_dimension>(flux_values);
}

template <ushort space_dimension>
auto Linear_Advection<space_dimension>::physical_fluxes(const std::vector<Solution_>& solutions) {
    const auto num_solution = solutions.size();    
    std::vector<This_::Physical_Flux_> physical_fluxes(num_solution);

    for (size_t i = 0; i < num_solution; ++i)
        physical_fluxes[i] = This_::physical_flux(solutions[i]);

    return physical_fluxes;
}

template <ushort space_dimension>
std::vector<std::array<double, space_dimension>> Linear_Advection<space_dimension>::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions) {
    const auto num_solution = solutions.size();
    
    std::array<double, space_dimension> projected_maximum_lambda;
    for (ushort i = 0; i < space_dimension; ++i)
        projected_maximum_lambda[i] = std::abs(This_::advection_speeds_[i]);
    
    return std::vector<std::array<double, space_dimension>>(num_solution, projected_maximum_lambda);
}

template <ushort space_dimension>
double Linear_Advection<space_dimension>::inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector) {
    return std::abs(This_::advection_speeds_.inner_product(nomal_vector));
}

template <ushort space_dimension>
auto Burgers<space_dimension>::physical_flux(const Solution_& solution) {
    const auto sol = solution.at(0); //scalar
    const auto flux_value = 0.5 * sol * sol;

    std::array<double, space_dimension> flux_values;
    flux_values.fill(flux_value);

    return Physical_Flux_(flux_values);
}

template <ushort space_dimension>
auto Burgers<space_dimension>::physical_fluxes(const std::vector<Solution_>& solutions) {
    const auto num_solution = solutions.size();
    std::vector<This_::Physical_Flux_> physical_fluxes(num_solution);

    for (size_t i = 0; i < num_solution; ++i)
        physical_fluxes[i] = This_::physical_flux(solutions[i]);

    return physical_fluxes;
}

template <ushort space_dimension>
std::vector<std::array<double, space_dimension>> Burgers<space_dimension>::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions) {
    const auto num_solution = solutions.size();

    std::vector<std::array<double, space_dimension>> projected_maximum_lambdas(num_solution);
    for (size_t i = 0; i < num_solution; ++i) {
        const auto maximum_lambdas = std::abs(solutions[i].at(0));	//scalar
        projected_maximum_lambdas[i].fill(maximum_lambdas);
    }

    return projected_maximum_lambdas;
}

template <ushort space_dimension>
double Burgers<space_dimension>::inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& normal_vector) {
    double normal_component_sum = 0.0;
    for (ushort i = 0; i < space_dimension; ++i)
        normal_component_sum += normal_vector[i];
            
    return std::max(std::abs(solution_o.at(0) * normal_component_sum), std::abs(solution_n.at(0) * normal_component_sum));
}

