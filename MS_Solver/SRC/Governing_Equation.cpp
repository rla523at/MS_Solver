#include "../INC/Governing_Equation.h"

Linear_Advection_2D::Physical_Flux_ Linear_Advection_2D::physical_flux(const Solution_& solution) {
	const auto [x_advection_speed, y_advection_speed] = Linear_Advection_2D::advection_speeds_;
	const auto sol = solution.at(0);	//scalar

	Physical_Flux_ physical_flux = { x_advection_speed * sol , y_advection_speed * sol };
	return physical_flux;
}

std::vector<Linear_Advection_2D::Physical_Flux_> Linear_Advection_2D::physical_fluxes(const std::vector<Solution_>& solutions) {
	//static size_t num_solution = solutions.size();
	const size_t num_solution = solutions.size();

	const auto [x_advection_speed, y_advection_speed] = Linear_Advection_2D::advection_speeds_;
	std::vector<Physical_Flux_> physical_fluxes(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto sol = solutions[i].at(0);	//scalar
		physical_fluxes[i] = { x_advection_speed * sol , y_advection_speed * sol };
	}

	return physical_fluxes;
}

//Dynamic_Matrix_ Linear_Advection_2D::flux_nodes(const Dynamic_Matrix_& solution_nodes) {
//	const auto [num_eq, num_node] = solution_nodes.size();
//
//	Dynamic_Matrix_ flux_nodes(num_eq, This_::space_dimension_ * num_node);	
//	
//	for (size_t i = 0; i < num_node; ++i) {
//		const auto physical_flux = This_::physical_flux(solution_nodes.column<This_::num_equation_>(i));
//		flux_nodes.change_columns(i * This_::space_dimension_, physical_flux);
//	}
//
//	return flux_nodes;
//}

std::vector<std::array<double, Linear_Advection_2D::space_dimension_>> Linear_Advection_2D::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions) {
	//static size_t num_solution = solutions.size();
	const size_t num_solution = solutions.size();

	static double absolute_x_advection_speed = std::abs(advection_speeds_[0]);
	static double absolute_y_advection_speed = std::abs(advection_speeds_[1]);

	std::vector<std::array<double, Linear_Advection_2D::space_dimension_>> projected_maximum_lambdas(num_solution, { absolute_x_advection_speed,absolute_y_advection_speed });
	return projected_maximum_lambdas;
}

double Linear_Advection_2D::inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector) {
	return std::abs(nomal_vector.inner_product(advection_speeds_));
}


Burgers_2D::Physical_Flux_ Burgers_2D::physical_flux(const Solution_& solution) {
	const auto sol = solution.at(0); //scalar

	const auto temp_val = 0.5 * sol * sol;
	return { temp_val, temp_val };
}

std::vector<Burgers_2D::Physical_Flux_> Burgers_2D::physical_fluxes(const std::vector<Solution_>& solutions) {
	static size_t num_solution = solutions.size();


	std::vector<Physical_Flux_> physical_fluxes(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto sol = solutions[i].at(0);	//scalar
		const auto temp_val = 0.5 * sol * sol;
		physical_fluxes[i] = { temp_val, temp_val };
	}

	return physical_fluxes;
}

//Dynamic_Matrix_ Burgers_2D::flux_nodes(const Dynamic_Matrix_& solution_nodes) {
//	const auto [num_eq, num_node] = solution_nodes.size();
//
//	Dynamic_Matrix_ flux_nodes(num_eq, This_::space_dimension_ * num_node);
//
//	for (size_t i = 0; i < num_node; ++i) {
//		const auto physical_flux = This_::physical_flux(solution_nodes.column<This_::num_equation_>(i));
//		flux_nodes.change_columns(i * This_::space_dimension_, physical_flux);
//	}
//
//	return flux_nodes;
//}

std::vector<std::array<double, Burgers_2D::space_dimension_>> Burgers_2D::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions) {
	static size_t num_solution = solutions.size();

	std::vector<std::array<double, Burgers_2D::space_dimension_>> projected_maximum_lambdas(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto maximum_lambdas = std::abs(solutions[i].at(0));	//scalar
		projected_maximum_lambdas[i] = { maximum_lambdas, maximum_lambdas };
	}

	return projected_maximum_lambdas;
}

double Burgers_2D::inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector) {
	const auto normal_component_sum = nomal_vector.at(0) + nomal_vector.at(1);
	return std::max(std::abs(solution_o.at(0) * normal_component_sum), std::abs(solution_n.at(0) * normal_component_sum));
}

Euler_2D::Solution_ Euler_2D::conservative_to_primitive(const Solution_& conservative_variable) {
	constexpr auto gamma = 1.4;
	
	const auto rho = conservative_variable.at(0);
	const auto rhou = conservative_variable.at(1);
	const auto rhov = conservative_variable.at(2);
	const auto rhoE = conservative_variable.at(3);

	const auto one_over_rho = 1.0 / rho;	

	const auto u = rhou * one_over_rho;
	const auto v = rhov * one_over_rho;
	const auto p = (rhoE - 0.5 * (rhou * u + rhov * v)) * (gamma - 1);
	const auto a = std::sqrt(gamma * p * one_over_rho);

	return { u,v,p,a };
}

std::vector<std::array<double, Euler_2D::space_dimension_>> Euler_2D::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& conservative_variables) {
	static auto num_solution = conservative_variables.size();

	std::vector<Solution_> primitive_variables;
	primitive_variables.reserve(num_solution);

	for (size_t i = 0; i < num_solution; ++i)
		primitive_variables.push_back(Euler_2D::conservative_to_primitive(conservative_variables[i]));

	std::vector<std::array<double,space_dimension_>> coordinate_projected_maximum_lambdas(num_solution);

	for (size_t i = 0; i < num_solution; ++i) {
		const auto u = primitive_variables[i].at(0);
		const auto v = primitive_variables[i].at(1);
		const auto a = primitive_variables[i].at(3);

		const auto x_projected_maximum_lambda = std::abs(u) + a;
		const auto y_projected_maximum_lambda = std::abs(v) + a;

		coordinate_projected_maximum_lambdas[i] = { x_projected_maximum_lambda, y_projected_maximum_lambda };
	}

	return coordinate_projected_maximum_lambdas;
}

Euler_2D::Physical_Flux_ Euler_2D::physical_flux(const Solution_& cvariable) {
	const auto pvariable = conservative_to_primitive(cvariable);
	return physical_flux(cvariable, pvariable);
}

Euler_2D::Physical_Flux_ Euler_2D::physical_flux(const Solution_& conservative_variable, const Solution_& primitivie_variable) {
	const auto rho = conservative_variable.at(0);
	const auto rhou = conservative_variable.at(1);
	const auto rhov = conservative_variable.at(2);
	const auto rhoE = conservative_variable.at(3);
	const auto u = primitivie_variable.at(0);
	const auto v = primitivie_variable.at(1);
	const auto p = primitivie_variable.at(2);
	const auto a = primitivie_variable.at(3);
	const auto rhouv = rhou * v;

	return 
	{
		rhou,				rhov,
		rhou * u + p,		rhouv,
		rhouv,				rhov * v + p,
		(rhoE + p) * u,		(rhoE + p) * v
	};
}

std::vector<Euler_2D::Physical_Flux_> Euler_2D::physical_fluxes(const std::vector<Solution_>& conservative_variables, const std::vector<Solution_>& primitive_variables) {
	static const size_t num_solution = conservative_variables.size();
	
	std::vector<Physical_Flux_> physical_fluxes;
	physical_fluxes.reserve(num_solution);

	for (size_t i = 0; i < num_solution; ++i) 
		physical_fluxes.push_back(physical_flux(conservative_variables[i], primitive_variables[i]));
	
	return physical_fluxes;
}

//Dynamic_Matrix_ Euler_2D::flux_nodes(const Dynamic_Matrix_& solution_nodes) {
//	const auto [num_eq, num_node] = solution_nodes.size();
//
//	Dynamic_Matrix_ flux_nodes(num_eq, This_::space_dimension_ * num_node);
//
//	for (size_t i = 0; i < num_node; ++i) {
//		const auto physical_flux = This_::physical_flux(solution_nodes.column<This_::num_equation_>(i));
//		flux_nodes.change_columns(i * This_::space_dimension_, physical_flux);
//	}
//
//	return flux_nodes;
//}

double Euler_2D::inner_face_maximum_lambda(const Solution_& oc_primitive_variable, const Solution_& nc_primitive_variable, const Space_Vector_& nomal_vector) {
	const auto oc_u = oc_primitive_variable.at(0);
	const auto oc_v = oc_primitive_variable.at(1);
	const auto oc_a = oc_primitive_variable.at(3);	
	const auto oc_side_face_maximum_lambda = std::abs(oc_u * nomal_vector.at(0) + oc_v * nomal_vector.at(1)) + oc_a;

	const auto nc_u = nc_primitive_variable.at(0);
	const auto nc_v = nc_primitive_variable.at(1);
	const auto nc_a = nc_primitive_variable.at(3);
	const auto nc_side_face_maximum_lambda = std::abs(nc_u * nomal_vector.at(0) + nc_v * nomal_vector.at(1)) + nc_a;

	return std::max(oc_side_face_maximum_lambda, nc_side_face_maximum_lambda);
}