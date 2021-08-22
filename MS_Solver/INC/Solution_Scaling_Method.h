#pragma once
#include "Governing_Equation.h"
#include "Reconstruction_Method_HOM.h"

namespace ms {
	template <typename Governing_Equation, typename Spatial_Discrete_Method>
	inline constexpr bool can_use_scaliling_method = std::is_same_v<Governing_Equation, Euler_2D> && std::is_same_v<Spatial_Discrete_Method, HOM>;
}

class Solution_Scaler
{
private:
	using This_ = Solution_Scaler;

private:
	inline static double fix_rate = 0.5;
	inline static std::vector<Dynamic_Matrix> set_of_cell_basis_qnodes_;
	inline static std::vector<std::pair<uint, Dynamic_Matrix>> cell_index_basis_qnodes_pairs;

private:
	Solution_Scaler(void) = delete;	

public:
	static void record_cell_basis_qnodes(const std::vector<Dynamic_Matrix>& set_of_cell_basis_qnodes) {
		This_::set_of_cell_basis_qnodes_ = set_of_cell_basis_qnodes;
	}

	static void record_face_basis_qnodes(const uint cell_index, const Dynamic_Matrix& basis_qnodes) {
		This_::cell_index_basis_qnodes_pairs.push_back({ cell_index, basis_qnodes });
	}

	template <ushort num_equation, ushort num_basis>
	static void inspect_and_scale(std::vector<Matrix<num_equation, num_basis>>& solution_coefficients) {
		//check cell
		const auto num_cell = This_::set_of_cell_basis_qnodes_.size();

		for (uint i = 0; i < num_cell; ++i) {
			ushort fix_count = 0;

			auto solution_qnodes = solution_coefficients[i] * This_::set_of_cell_basis_qnodes_[i];
			const auto [temp, num_qnode] = solution_qnodes.size();

			for (ushort j = 0; j < num_qnode; ++j) {
				while (true) {
					dynamic_require(fix_count < 10, "More then 10 attemps to fix pressure is meaningless");

					const auto cvariable = solution_qnodes.column<num_equation>(j);
					const auto pvariable = Euler_2D::conservative_to_primitive(cvariable);
					const auto rho = cvariable[0];
					const auto pressure = pvariable[2];

					if (rho <= 0.0 || pressure <= 0.0) {
						fix_count++;
						const auto fix_matrix = This_::calculate_fix_matrix<num_basis>();
						solution_coefficients[i] *= fix_matrix;

						solution_qnodes = solution_coefficients[i] * This_::set_of_cell_basis_qnodes_[i];
						continue;
					}

					break;
				}
			}			
		}


		//check face
		for (const auto& [cell_index, basis_qnodes] : This_::cell_index_basis_qnodes_pairs) {
			ushort fix_count = 0;

			auto solution_qnodes = solution_coefficients[cell_index] * basis_qnodes;
			const auto [temp, num_qnode] = solution_qnodes.size();

			for (ushort j = 0; j < num_qnode; ++j) {
				while (true) {
					dynamic_require(fix_count < 10, "More then 10 attemps to fix pressure is meaningless");

					const auto cvariable = solution_qnodes.column<num_equation>(j);
					const auto pvariable = Euler_2D::conservative_to_primitive(cvariable);
					const auto rho = cvariable[0];
					const auto pressure = pvariable[2];

					if (rho <= 0.0 || pressure <= 0.0) {
						fix_count++;
						const auto fix_matrix = This_::calculate_fix_matrix<num_basis>();
						solution_coefficients[cell_index] *= fix_matrix;

						solution_qnodes = solution_coefficients[cell_index] * basis_qnodes;
						continue;
					}

					break;
				}
			}
		}
	}

private:
	template <ushort num_basis>
	static Matrix<num_basis, num_basis> calculate_fix_matrix(void) {
		std::array<double, num_basis> limiting_values;
		limiting_values.fill(fix_rate);
		limiting_values[0] = 1.0; //preserve P0 values
		return limiting_values;
	}
};