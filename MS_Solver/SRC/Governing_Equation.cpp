#include "../INC/Governing_Equation.h"

std::vector<std::array<double, Euler_2D::space_dimension_>> Euler_2D::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& conservative_variables) {
	auto num_solution = conservative_variables.size();

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
	const size_t num_solution = conservative_variables.size();
	
	std::vector<Physical_Flux_> physical_fluxes;
	physical_fluxes.reserve(num_solution);

	for (size_t i = 0; i < num_solution; ++i) 
		physical_fluxes.push_back(physical_flux(conservative_variables[i], primitive_variables[i]));
	
	return physical_fluxes;
}

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