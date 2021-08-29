#include "../INC/Inital_Condition.h"

Euclidean_Vector<1> Constant1_2D::calculate_solution(const Space_Vector_& space_vector) {
	return { 1 };
};

std::vector<Euclidean_Vector<1>> Constant1_2D::calculate_solutions(const std::vector<Space_Vector_>& space_vectors) {
	std::vector<Euclidean_Vector<1>> result(space_vectors.size(), { 1 });
	return result;
};

std::string Constant1_2D::name(void) {
	return "Constant_2D";
}

void Sine_Wave_2D::initialize(const double x_wave_length, const double y_wave_length) {
	This_::x_wave_number_ = 2 * std::numbers::pi / x_wave_length;
	This_::y_wave_number_ = 2 * std::numbers::pi / y_wave_length;
}

Euclidean_Vector<1> Sine_Wave_2D::calculate_solution(const Space_Vector_& space_vector) {
	const auto x_coord = space_vector.at(0);
	const auto y_coord = space_vector.at(1);

	return { std::sin(This_::x_wave_number_ * x_coord) * std::sin(This_::y_wave_number_ * y_coord) };
}

std::vector<Euclidean_Vector<1>> Sine_Wave_2D::calculate_solutions(const std::vector<Space_Vector_>& space_vectors) {
	const auto num_cell = space_vectors.size();

	std::vector<Solution_> solutions_(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coord = space_vectors[i].at(0);
		const auto y_coord = space_vectors[i].at(1);

		solutions_[i] = std::sin(This_::x_wave_number_ * x_coord) * std::sin(This_::y_wave_number_ * y_coord);
	}

	return solutions_;
}

std::string Sine_Wave_2D::name(void) {
	return "Sine_Wave_2D";
};

template <typename Governing_Equation>
std::vector<Euclidean_Vector<1>> Sine_Wave_2D::calculate_exact_solutions(const std::vector<Space_Vector_>& cell_centers, const double end_time) {
	static_require(ms::is_Linear_Advection_2D<Governing_Equation>, "excat solution can be calculated when governing equation is linear advection");

	const auto num_cell = cell_centers.size();
	const auto [x_advection_speed, y_advection_speed] = Governing_Equation::advection_speed();

	std::vector<Solution_> exact_solutions_(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coord = cell_centers[i].at(0);
		const auto y_coord = cell_centers[i].at(1);
		exact_solutions_[i] = std::sin(This_::x_wave_number_ * (x_coord - x_advection_speed * end_time)) * std::sin(This_::y_wave_number_ * (y_coord - y_advection_speed * end_time));
	}

	return exact_solutions_;
}


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

std::vector<SOD_2D::Solution_> SOD_2D::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
	const auto num_cell = cell_centers.size();

	std::vector<Solution_> solutions(num_cell);
	for (size_t i = 0; i < num_cell; ++i)
		solutions[i] = This_::calculate_solution(cell_centers[i]);

	return solutions;
}

Modified_SOD_2D::Solution_ Modified_SOD_2D::calculate_solution(const Space_Vector_& space_vector) {
	constexpr auto gamma = 1.4;
	constexpr auto c = 1.0 / (gamma - 1.0);

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

	std::vector<Solution_> solutions(num_cell);
	for (size_t i = 0; i < num_cell; ++i)
		solutions[i] = This_::calculate_solution(cell_centers[i]);

	return solutions;
}

Shu_Osher_2D::Solution_ Shu_Osher_2D::calculate_solution(const Space_Vector_& space_vector) {
	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);
	constexpr auto discontinuity_location = -4.0;

	const auto x_coordinate = space_vector.at(0);

	if (x_coordinate < discontinuity_location) {
		constexpr auto rho = 3.857143;
		constexpr auto u = 2.629369;	//rhou = 10.1418522328
		constexpr auto v = 0.0;			//rhov = 0.0
		constexpr auto p = 10.333333;	//rhoE = 39.1666684317

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

std::vector<Shu_Osher_2D::Solution_> Shu_Osher_2D::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
	const auto num_cell = cell_centers.size();

	std::vector<Solution_> solutions(num_cell);
	for (size_t i = 0; i < num_cell; ++i)
		solutions[i] = This_::calculate_solution(cell_centers[i]);	

	return solutions;
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