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

Euclidean_Vector Euler2D::calculate_solution(const Euclidean_Vector& GE_solution) const 
{
	const auto rho = GE_solution.at(0);
	const auto rhou = GE_solution.at(1);
	const auto rhov = GE_solution.at(2);
	const auto rhoE = GE_solution.at(3);

	const auto one_over_rho = 1.0 / rho;

	const auto u = rhou * one_over_rho;
	const auto v = rhov * one_over_rho;
	const auto p = (rhoE - 0.5 * (rhou * u + rhov * v)) * (this->gamma_ - 1);
	const auto a = std::sqrt(this->gamma_ * p * one_over_rho);

	return { rho, rhou, rhov, rhoE, u, v, p, a };
}

std::vector<std::vector<double>> Euler2D::calculate_coordinate_projected_maximum_lambda(const std::vector<Euclidean_Vector>& P0_solutions) const
{
	auto num_solution = P0_solutions.size();

	std::vector<std::vector<double>> coordinate_projected_maximum_lambdas(num_solution);

	for (size_t i = 0; i < num_solution; ++i)
	{
		const auto& solution = P0_solutions[i];
		const auto u = solution[4];
		const auto v = solution[5];
		const auto a = solution[7];

		const auto x_projected_maximum_lambda = std::abs(u) + a;
		const auto y_projected_maximum_lambda = std::abs(v) + a;

		coordinate_projected_maximum_lambdas[i] = { x_projected_maximum_lambda, y_projected_maximum_lambda };
	}

	return coordinate_projected_maximum_lambdas;
}

double Euler2D::calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const
{
	const auto oc_u = oc_solution[4];
	const auto oc_v = oc_solution[5];
	const auto oc_a = oc_solution[7];
	const auto oc_side_face_maximum_lambda = std::abs(oc_u * normal_vector[0] + oc_v * normal_vector[1]) + oc_a;

	const auto nc_u = nc_solution[4];
	const auto nc_v = nc_solution[5];
	const auto nc_a = nc_solution[7];
	const auto nc_side_face_maximum_lambda = std::abs(nc_u * normal_vector[0] + nc_v * normal_vector[1]) + nc_a;

	return (std::max)(oc_side_face_maximum_lambda, nc_side_face_maximum_lambda);	
}


Matrix Euler2D::calculate_physical_flux(const Euclidean_Vector& solution) const
{
	const auto rho = solution[0];
	const auto rhou = solution[1];
	const auto rhov = solution[2];
	const auto rhoE = solution[3];
	const auto u = solution[4];
	const auto v = solution[5];
	const auto p = solution[6];

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

std::shared_ptr<Governing_Equation> Governing_Equation_Factory::make_shared(const Configuration& config)
{
	const auto governing_equation = config.get("Governing_Equation");
	const auto space_dimension = config.get<ushort>("space_dimension");

	if (ms::contains_icase(governing_equation, "Euler"))
	{
		if (space_dimension == 2)
		{
			return std::make_shared<Euler2D>();
		}
	}
}

std::unique_ptr<Governing_Equation> Governing_Equation_Factory::make_unique(const Configuration& config)
{
	const auto governing_equation = config.get("Governing_Equation");
	const auto space_dimension = config.get<ushort>("space_dimension");

	if (ms::contains_icase(governing_equation, "Euler")) 
	{
		if (space_dimension == 2)
		{
			return std::make_unique<Euler2D>();
		}
	}
}