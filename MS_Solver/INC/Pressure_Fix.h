#pragma once
#include "Governing_Equation.h"
#include "Reconstruction_Method_HOM.h"
#include "Setting.h"

namespace ms {
	template <typename Governing_Equation, typename Spatial_Discrete_Method>
	inline constexpr bool can_use_pressure_fix = std::is_same_v<Governing_Equation, Euler_2D> && std::is_same_v<Spatial_Discrete_Method, HOM>;
}

class Pressure_Fix
{
private:
	using This_ = Pressure_Fix;

private:
	inline static double fix_rate = 0.5;
	inline static std::vector<Dynamic_Matrix> set_of_basis_qnodes_;

private:
	Pressure_Fix(void) = delete;	

public:
	template <typename Reconstruction_Method>
	static void initialize(const Grid<Reconstruction_Method::space_dimension()>& grid) {
#ifdef PRESSURE_FIX_MODE

		SET_TIME_POINT;

		const auto& cell_elements = grid.elements.cell_elements;

		const auto num_cell = cell_elements.size();
		This_::set_of_basis_qnodes_.reserve(num_cell); // 2 * solution_order

		constexpr auto solution_order = Reconstruction_Method::solution_order();
		constexpr auto integrand_order = 2 * solution_order + 1;

		for (uint i = 0; i < num_cell; ++i) {
			const auto& geometry = cell_elements[i].geometry_;

			const auto& quadrature_rule = geometry.get_quadrature_rule(integrand_order);
			This_::set_of_basis_qnodes_.push_back(Reconstruction_Method::calculate_basis_nodes(i, quadrature_rule.points));
		}

		Log::content_ << std::left << std::setw(50) << "@ Pressure Fix precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
		Log::print();

#endif

	}

	template <ushort num_equation, ushort num_basis>
	static void reconstruct(std::vector<Matrix<num_equation, num_basis>>& solution_coefficients) {
#ifdef PRESSURE_FIX_MODE

		const auto num_solution = solution_coefficients.size();

		for (uint i = 0; i < num_solution; ++i) {
			ushort fix_count = 0;

			auto solution_qnodes = solution_coefficients[i] * This_::set_of_basis_qnodes_[i];
			const auto [temp, num_qnode] = solution_qnodes.size();

			for (ushort j = 0; j < num_qnode; ++j) {
				while (true) {
					dynamic_require(fix_count < 10, "More then 10 attemps to fix pressure is meaningless");

					const auto cvariable = solution_qnodes.column<num_equation>(j);
					const auto pvariable = Euler_2D::conservative_to_primitive(cvariable);
					const auto pressure = pvariable[2];

					if (pressure < 0.0) {
						fix_count++;

						std::array<double, num_basis> limiting_values;
						limiting_values.fill(fix_rate);
						limiting_values[0] = 1.0; //preserve P0 values
						const auto limiting_matrix = Matrix<num_basis, num_basis>::diagonal_matrix(limiting_values);

						solution_coefficients[i] *= limiting_matrix;
						solution_qnodes = solution_coefficients[i] * This_::set_of_basis_qnodes_[i];
						continue;
					}

					break;
				}
			}			
		}

#endif
	}
};