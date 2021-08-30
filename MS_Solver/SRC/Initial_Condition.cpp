#include "../INC/Inital_Condition.h"

Square_Wave_2D::Solution_ Square_Wave_2D::calculate_solution(const Space_Vector_& space_vector) {
	const auto x_coord = space_vector.at(0);
	const auto y_coord = space_vector.at(1);

	if (0.25 <= x_coord && x_coord <= 0.75 && 0.25 <= y_coord && y_coord <= 0.75)
		return { 1 };
	else
		return { 0 };
}

std::vector<Square_Wave_2D::Solution_> Square_Wave_2D::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
	const auto num_cell = cell_centers.size();

	std::vector<Solution_> solutions_(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coord = cell_centers[i].at(0);
		const auto y_coord = cell_centers[i].at(1);

		if (0.25 <= x_coord && x_coord <= 0.75 && 0.25 <= y_coord && y_coord <= 0.75)
			solutions_[i] = 1;
		else
			solutions_[i] = 0;
	}

	return solutions_;
}

Circle_Wave_2D::Solution_ Circle_Wave_2D::calculate_solution(const Space_Vector_& space_vector) {
	const auto x_coord = space_vector.at(0);
	const auto y_coord = space_vector.at(1);

	if ( (x_coord - 0.5) * (x_coord - 0.5) + (y_coord - 0.5) * (y_coord - 0.5) <= 0.25 * 0.25 )
		return { 1 };
	else
		return { 0 };
}

std::vector<Circle_Wave_2D::Solution_> Circle_Wave_2D::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
	const auto num_cell = cell_centers.size();

	std::vector<Solution_> solutions_(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coord = cell_centers[i].at(0);
		const auto y_coord = cell_centers[i].at(1);

		if ((x_coord - 0.5) * (x_coord - 0.5) + (y_coord - 0.5) * (y_coord - 0.5) <= 0.25 * 0.25 )
			solutions_[i] = 1;
		else
			solutions_[i] = 0;
	}

	return solutions_;
}

Euclidean_Vector<1> Gaussian_Wave_2D::calculate_solution(const Space_Vector_& space_vector) {
	const auto x_coord = space_vector.at(0);
	const auto y_coord = space_vector.at(1);

	return { exp(- This_::beta_*( (x_coord - 0.5) * (x_coord - 0.5) + (y_coord - 0.5) * (y_coord - 0.5))) };
}

std::vector<Euclidean_Vector<1>> Gaussian_Wave_2D::calculate_solutions(const std::vector<Space_Vector_>& space_vectors) {
	const auto num_cell = space_vectors.size();

	std::vector<Solution_> solutions_(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coord = space_vectors[i].at(0);
		const auto y_coord = space_vectors[i].at(1);

		solutions_[i] = exp(-This_::beta_ * ((x_coord - 0.5) * (x_coord - 0.5) + (y_coord - 0.5) * (y_coord - 0.5)));
	}

	return solutions_;
}
//template <>
//std::vector<Square_Wave_2D::Solution> Square_Wave_2D::calculate_exact_solutions<Linear_Advection_2D>(const std::vector<Space_Vector_>& space_vector, const double end_time) {
//	const auto num_cell = space_vector.size();
//	const auto [x_advection_speed, y_advection_speed] = Linear_Advection_2D::advection_speed();
//
//	std::vector<Solution> exact_solutions_(num_cell);
//	for (size_t i = 0; i < num_cell; ++i) {
//		const auto x_coord = space_vector[i].at(0);
//		const auto y_coord = space_vector[i].at(1);
//
//		//Assume that domian [0,1] x [0,1]
//		const auto exact_x_start	= 0.25 + x_advection_speed * end_time - static_cast<int>(0.25 + x_advection_speed * end_time);
//		const auto exact_x_end		= 0.75 + x_advection_speed * end_time - static_cast<int>(0.75 + x_advection_speed * end_time);
//		const auto exact_y_start	= 0.25 + y_advection_speed * end_time - static_cast<int>(0.25 + y_advection_speed * end_time);
//		const auto exact_y_end		= 0.75 + y_advection_speed * end_time - static_cast<int>(0.75 + y_advection_speed * end_time);
//
//		if (exact_x_start <= x_coord && x_coord <= exact_x_end &&
//			exact_y_start <= y_coord && y_coord <= exact_y_end)
//			exact_solutions_[i] = 1;
//		else
//			exact_solutions_[i] = 0;
//	}
//
//	return exact_solutions_;
//}
SOD_2D::Solution_ SOD_2D::calculate_solution(const Space_Vector_& space_vector) {
	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (1.4 - 1);
	constexpr auto discontinuity_location = 0.5;

	const auto x_coordinate = space_vector.at(0);

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


Modified_SOD_2D::Solution_ Modified_SOD_2D::calculate_solution(const Space_Vector_& space_vector) {
	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (1.4 - 1);

	const auto x_coordinate = space_vector.at(0);

	if (x_coordinate <= 0.3) {
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


std::vector<Modified_SOD_2D::Solution_> Modified_SOD_2D::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
	const auto num_cell = cell_centers.size();
	
	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (1.4 - 1);
	std::vector<Solution_> solutions(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coordinate = cell_centers[i].at(0);

		if (x_coordinate <= 0.3) {
			const auto rho = 1.0;
			const auto u = 0.75;
			const auto v = 0.0;
			const auto p = 1.0;

			const auto rhou = rho * u;
			const auto rhov = rho * v;
			const auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

			solutions[i] = { rho, rhou, rhov, rhoE };
		}
		else {
			const auto rho = 0.125;
			const auto u = 0.0;
			const auto v = 0.0;
			const auto p = 0.1;

			const auto rhou = rho * u;
			const auto rhov = rho * v;
			const auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

			solutions[i] = { rho, rhou, rhov, rhoE };
		}
	}

	return solutions;
}
