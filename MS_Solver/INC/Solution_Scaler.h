#pragma once
#include "Governing_Equation.h"
#include "Reconstruction_Method_HOM.h"

namespace ms {
	template <typename Governing_Equation, typename Spatial_Discrete_Method>
	inline constexpr bool can_use_scaliling_method = ms::is_Euler<Governing_Equation> && std::is_same_v<Spatial_Discrete_Method, HOM>;
}

template <ushort space_dimension>
class Solution_Scaler
{
private:
	Solution_Scaler(void) = delete;

private:
	using This_ = Solution_Scaler;

private:
	inline static double fix_rate = 0.5;
	inline static std::vector<Dynamic_Matrix> set_of_cell_basis_qnodes_;
	inline static std::vector<std::pair<uint, Dynamic_Matrix>> cell_index_basis_qnodes_pairs;


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
					const auto cvariable = solution_qnodes.column<num_equation>(j);
					const auto pvariable = Euler<space_dimension>::conservative_to_primitive(cvariable);
					const auto rho = cvariable[0];
					const auto p = pvariable[num_equation - 2];

					if (rho <=0 || p <=0) {
						dynamic_require(fix_count < 10, "More then 10 attemps to fix is meaningless");
						fix_count++;						
						solution_coefficients[i] *= This_::calculate_fix_matrix<num_basis>();

						solution_qnodes = solution_coefficients[i] * This_::set_of_cell_basis_qnodes_[i];

						Log::content_ << "\n" << i << " cell has negative density or pressure, fix " << fix_count << " time\n";
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
					const auto cvariable = solution_qnodes.column<num_equation>(j);
					const auto pvariable = Euler<space_dimension>::conservative_to_primitive(cvariable);
					const auto rho = cvariable[0];
					const auto p = pvariable[num_equation - 2];


					if (rho <= 0 || p <= 0) {
						dynamic_require(fix_count < 10, "More then 10 attemps to fix is meaningless");
						fix_count++;
						solution_coefficients[cell_index] *= This_::calculate_fix_matrix<num_basis>();

						solution_qnodes = solution_coefficients[cell_index] * basis_qnodes;

						Log::content_ << "\n" << cell_index << " cell has negative density or pressure, fix " << fix_count << " time\n";
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