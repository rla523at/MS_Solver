#pragma once
#include "Governing_Equation.h"


using uint = unsigned int;


class NFF {};    // Numerical Flux Function


template <typename Governing_Equation>
class LLF : public NFF  // Local Lax Fridrich method
{
private:
    static_require(ms::is_governing_equation<Governing_Equation>, "It should be Governing Equation");

    static constexpr ushort space_dimension_    = Governing_Equation::space_dimension();
    static constexpr ushort num_equation_       = Governing_Equation::num_equation();

    using This_             = LLF<Governing_Equation>;
    using Space_Vector_     = Euclidean_Vector<space_dimension_>;
    using Solution_         = Euclidean_Vector<num_equation_>;
    using Numerical_Flux_   = Euclidean_Vector<num_equation_>;

public:
    static std::string name(void) { return "LLF"; };
    static constexpr ushort space_dimension(void) { return This_::space_dimension_; };
    static constexpr ushort num_equation(void) { return This_::num_equation_; };

private:
    LLF(void) = delete;

public:
    static auto calculate(const std::vector<Solution_>& solutions, const std::vector<Space_Vector_>& normals, const std::vector<std::pair<uint, uint>>& oc_nc_index_pairs);
    static auto calculate(const Solution_& oc_side_solution, const Solution_& nc_side_solution, const Space_Vector_& normal);
};


template<>
class LLF<Euler_2D> : public NFF
{
private:
    static constexpr ushort space_dimension_    = Euler_2D::space_dimension();
    static constexpr ushort num_equation_       = Euler_2D::num_equation();

    using This_             = LLF<Euler_2D>;
    using Space_Vector_     = Euclidean_Vector<space_dimension_>;
    using Solution_         = Euclidean_Vector<num_equation_>;
    using Numerical_Flux_   = Euclidean_Vector<num_equation_>;

public:
    static std::string name(void) { return "LLF"; };
    static constexpr ushort space_dimension(void) { return This_::space_dimension_; };
    static constexpr ushort num_equation(void) { return This_::num_equation_; };

private:
    LLF(void) = delete;

public:
    static Numerical_Flux_ calculate(const Solution_& oc_side_cvariable, const Solution_& nc_side_cvariable, const Space_Vector_& normal);
    static std::vector<Numerical_Flux_> calculate(const std::vector<Solution_>& conservative_variables, const std::vector<Space_Vector_>& normals, const std::vector<std::pair<uint, uint>>& oc_nc_index_pairs);
};


namespace ms {
    template <typename T>
    inline constexpr bool is_numeirical_flux_function = std::is_base_of_v<NFF, T>;
}


//template definition part
template <typename Governing_Equation>
auto LLF<Governing_Equation>::calculate(const std::vector<Solution_>& solutions, const std::vector<Space_Vector_>& normals, const std::vector<std::pair<uint, uint>>& oc_nc_index_pairs) {
    const auto num_inner_face = normals.size();
    const auto physical_fluxes = Governing_Equation::physical_fluxes(solutions);

    std::vector<Numerical_Flux_> inner_face_numerical_fluxes(num_inner_face);
    for (size_t i = 0; i < num_inner_face; ++i) {
        const auto [oc_index, nc_index] = oc_nc_index_pairs[i];
        const auto oc_physical_flux = physical_fluxes[oc_index];
        const auto nc_physical_flux = physical_fluxes[nc_index];

        const auto oc_solution = solutions[oc_index];
        const auto nc_solution = solutions[nc_index];
        const auto normal = normals[i];
        const auto inner_face_maximum_lambda = Governing_Equation::inner_face_maximum_lambda(oc_solution, nc_solution, normal);

        inner_face_numerical_fluxes[i] = 0.5 * ((oc_physical_flux + nc_physical_flux) * normal + inner_face_maximum_lambda * (oc_solution - nc_solution));
    }
    return inner_face_numerical_fluxes;
};

template <typename Governing_Equation>
auto LLF<Governing_Equation>::calculate(const Solution_& oc_side_solution, const Solution_& nc_side_solution, const Space_Vector_& normal) {
    const auto oc_physical_flux = Governing_Equation::physical_flux(oc_side_solution);
    const auto nc_physical_flux = Governing_Equation::physical_flux(nc_side_solution);
    const auto inner_face_maximum_lambda = Governing_Equation::inner_face_maximum_lambda(oc_side_solution, nc_side_solution, normal);

    Numerical_Flux_ LLF_flux = 0.5 * ((oc_physical_flux + nc_physical_flux) * normal + inner_face_maximum_lambda * (oc_side_solution - nc_side_solution));
    return LLF_flux;
}