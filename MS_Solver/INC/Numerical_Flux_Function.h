#pragma once
#include "Governing_Equation.h"


using uint = unsigned int;


class NFF {};    // Numerical Flux Function


template <typename Governing_Equation>
class LLF : public NFF  // Local Lax Fridrich method
{
private:
    LLF(void) = delete;

public:
    using Governing_Equation_ = Governing_Equation;

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

public:
    static Numerical_Flux_ calculate(const Solution_& oc_side_solution, const Solution_& nc_side_solution, const Space_Vector_& normal);
    static std::vector<Numerical_Flux_> calculate(const std::vector<Solution_>& solutions, const std::vector<Space_Vector_>& normals, const std::vector<std::pair<uint, uint>>& oc_nc_index_pairs);
};


namespace ms {
    template <typename T>
    inline constexpr bool is_numeirical_flux_function = std::is_base_of_v<NFF, T>;
}


//template definition part
template <typename Governing_Equation>
LLF<Governing_Equation>::Numerical_Flux_ LLF<Governing_Equation>::calculate(const Solution_& oc_side_solution, const Solution_& nc_side_solution, const Space_Vector_& normal) {
    if constexpr (ms::is_SCL<Governing_Equation>) {
        const auto oc_physical_flux = Governing_Equation::physical_flux(oc_side_solution);
        const auto nc_physical_flux = Governing_Equation::physical_flux(nc_side_solution);
        const auto inner_face_maximum_lambda = Governing_Equation::inner_face_maximum_lambda(oc_side_solution, nc_side_solution, normal);

        return 0.5 * ((oc_physical_flux + nc_physical_flux) * normal + inner_face_maximum_lambda * (oc_side_solution - nc_side_solution));
    }
    else if constexpr (ms::is_Euler<Governing_Equation>) {
        const auto oc_side_pvariable = Governing_Equation::conservative_to_primitive(oc_side_solution);
        const auto nc_side_pvariable = Governing_Equation::conservative_to_primitive(nc_side_solution);

        const auto oc_physical_flux = Governing_Equation::physical_flux(oc_side_solution, oc_side_pvariable);
        const auto nc_physical_flux = Governing_Equation::physical_flux(nc_side_solution, nc_side_pvariable);
        const auto inner_face_maximum_lambda = Governing_Equation::inner_face_maximum_lambda(oc_side_pvariable, nc_side_pvariable, normal);

        return 0.5 * ((oc_physical_flux + nc_physical_flux) * normal + inner_face_maximum_lambda * (oc_side_solution - nc_side_solution));
    }
    else {
        throw std::invalid_argument("not supported governing equation");
        return {};
    }
}

template <typename Governing_Equation>
std::vector<typename LLF<Governing_Equation>::Numerical_Flux_> LLF<Governing_Equation>::calculate(const std::vector<Solution_>& solutions, const std::vector<Space_Vector_>& normals, const std::vector<std::pair<uint, uint>>& oc_nc_index_pairs) {
      const auto num_inner_face = normals.size();
    std::vector<Numerical_Flux_> inner_face_numerical_fluxes(num_inner_face);

    for (uint i = 0; i < num_inner_face; ++i) {
        const auto [oc_index, nc_index] = oc_nc_index_pairs[i];
        const auto oc_solution = solutions[oc_index];
        const auto nc_solution = solutions[nc_index];
        const auto normal = normals[i];

        inner_face_numerical_fluxes[i] = This_::calculate(oc_solution, nc_solution, normal);
    }

    return inner_face_numerical_fluxes;
};