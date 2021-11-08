#include "../INC/Governing_Equation.h"

const std::vector<std::string>& Governing_Equation::get_variable_names(void) const 
{
	return variable_names_;
}

ushort Governing_Equation::num_equations(void) const 
{
	return this->num_equations_;
}

ushort Governing_Equation::space_dimension(void) const 
{
	return this->space_dimension_;
}

Euler2D::Euler2D(void) 
{
	this->space_dimension_ = 2;
	this->num_equations_ = 4;
	this->variable_names_ = { "rho", "rhou", "rhov", "rhoE", "u", "v", "p", "e" };
}

Matrix Euler2D::calculate_physical_flux(const Euclidean_Vector& solution) const
{
	const auto rho = solution.at(0);
	const auto rhou = solution.at(1);
	const auto rhov = solution.at(2);
	const auto rhoE = solution.at(3);

	const auto velocity = this->calculate_velocity(solution);
	const auto u = velocity.at(0);
	const auto v = velocity.at(1);

	const auto p = this->calculate_pressure(solution, velocity);
	const auto rhouv = rhou * v;

	REQUIRE(rho >= 0 && p >= 0, "density and pressure shold be positive");
		
	return { this->num_equations_, this->space_dimension_,
		{
		rhou,				rhov,
		rhou * u + p,		rhouv,
		rhouv,				rhov * v + p,
		(rhoE + p) * u,		(rhoE + p) * v
		} };
}

Euclidean_Vector Euler2D::calculate_velocity(const Euclidean_Vector& solution) const
{
	const auto rho = solution[0];
	const auto one_over_rho = 1.0 / rho;

	const auto rhou = solution[1];
	const auto rhov = solution[2];

	return { rhou * one_over_rho, rhov * one_over_rho };
}

double Euler2D::calculate_pressure(const Euclidean_Vector& solution, const Euclidean_Vector& velocity) const
{
	constexpr auto gamma = 1.4;
	
	const auto rho = solution[0];
	const auto rhou = solution[1];
	const auto rhov = solution[2];
	const auto rhoE = solution[3];
	const auto u = velocity[0];
	const auto v = velocity[1];

	return (rhoE - 0.5 * (rhou * u + rhov * v)) * (gamma - 1);
}
