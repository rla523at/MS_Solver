#include "../INC/Governing_Equation.h"

const std::vector<std::string>& Governing_Equation::get_solution_names(void) const 
{
	return solution_names_;
}

ushort Governing_Equation::num_equations(void) const 
{
	return this->num_equations_;
}

ushort Governing_Equation::num_solutions(void) const
{
	return this->num_solutions_;
}

ushort Governing_Equation::space_dimension(void) const 
{
	return this->space_dimension_;
}

Scalar_Equation::Scalar_Equation(void)
{
	this->num_solutions_ = 1;
	this->num_equations_ = 1;
	this->solution_names_ = { "q" };
}

std::vector<std::vector<double>> Linear_Advection::calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const
{
	const auto num_solution = P0_solutions.size();

	std::vector<double> coordinate_projected_maximum_lambda(this->space_dimension_);

	for (int i = 0; i < this->space_dimension_; ++i)
	{
		coordinate_projected_maximum_lambda[i] = std::abs(this->advection_speeds_[i]);
	}		

	return std::vector<std::vector<double>>(num_solution, coordinate_projected_maximum_lambda);
}

double Linear_Advection::calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const
{
	return std::abs(this->advection_speeds_.inner_product(normal_vector));
}

Matrix Linear_Advection::calculate_physical_flux(const Euclidean_Vector& solution) const
{
	std::vector<double> values(this->space_dimension_);

	for (int i = 0; i < this->space_dimension_; ++i)
	{
		values[i] = this->advection_speeds_[i] * solution[0];
	}

	return { this->num_equations_, this->space_dimension_, std::move(values) };
}

Linear_Advection_2D::Linear_Advection_2D(const double x_advection_speed, const double y_advection_speed)
{
	this->space_dimension_ = 2;
	this->advection_speeds_ = { x_advection_speed, y_advection_speed };
}

Linear_Advection_3D::Linear_Advection_3D(const double x_advection_speed, const double y_advection_speed, const double z_advection_speed)
{
	this->space_dimension_ = 3;
	this->advection_speeds_ = { x_advection_speed, y_advection_speed, z_advection_speed };
}

std::vector<std::vector<double>> Burgers::calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const
{
	const auto num_solution = P0_solutions.size();

	std::vector<std::vector<double>> projected_maximum_lambdas(num_solution);
	for (int i = 0; i < num_solution; ++i) 
	{
		const auto maximum_lambda = std::abs(P0_solutions[i][0]);		
		projected_maximum_lambdas[i] = std::vector<double>(this->space_dimension_, maximum_lambda);
	}

	return projected_maximum_lambdas;
}
double Burgers::calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const
{
	double normal_component_sum = 0.0;
	for (ushort i = 0; i < this->space_dimension_; ++i)
		normal_component_sum += normal_vector[i];

	return (std::max)(std::abs(oc_solution[0] * normal_component_sum), std::abs(nc_solution[0] * normal_component_sum));
}
Matrix Burgers::calculate_physical_flux(const Euclidean_Vector& solution) const
{
	const double flux_value = 0.5 * std::pow(solution[0], 2.0);

	std::vector<double> flux_values(this->space_dimension_, flux_value);

	return { this->num_equations_, this->space_dimension_, std::move(flux_values) };
}

Burgers_2D::Burgers_2D(void)
{
	this->space_dimension_ = 2;
}

Burgers_3D::Burgers_3D(void)
{
	this->space_dimension_ = 3;
}

Euler_2D::Euler_2D(void) 
{
	this->space_dimension_ = 2;
	this->num_solutions_ = 8;
	this->num_equations_ = 4;
	this->solution_names_ = { "rho", "rhou", "rhov", "rhoE", "u", "v", "p", "a" };
}

std::vector<std::vector<double>> Euler_2D::calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const
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

double Euler_2D::calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& normal_vector) const
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

Matrix Euler_2D::calculate_physical_flux(const Euclidean_Vector& solution) const
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

void Euler_2D::extend_to_solution(Euclidean_Vector& GE_solution) const
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

	GE_solution = { rho, rhou, rhov, rhoE, u, v, p, a };
}

std::shared_ptr<Governing_Equation> Governing_Equation_Factory::make_shared(const Configuration& config)
{
	const auto governing_equation = config.get("Governing_Equation");
	const auto space_dimension = config.get<ushort>("space_dimension");

	if (ms::contains_icase(governing_equation, "Linear_Advection"))
	{
		if (space_dimension == 2)
		{
			const auto x_advection_speed = config.get<double>("x_advection_speed");
			const auto y_advection_speed = config.get<double>("y_advection_speed");
			return std::make_shared<Linear_Advection_2D>(x_advection_speed, y_advection_speed);
		}
		else
		{
			EXCEPTION("Linear advection does not support space dimension in configuration file");
			return nullptr;
		}
	}
	else if (ms::contains_icase(governing_equation, "Euler"))
	{
		if (space_dimension == 2)
		{
			return std::make_shared<Euler_2D>();
		}
		else
		{
			EXCEPTION("Euler does not support space dimension in configuration file");
			return nullptr;
		}
	}
	else
	{
		EXCEPTION("governing equation in configuration file is not supported");
		return nullptr;
	}
}

std::unique_ptr<Governing_Equation> Governing_Equation_Factory::make_unique(const Configuration& config)
{
	const auto governing_equation = config.get("Governing_Equation");
	const auto space_dimension = config.get<ushort>("space_dimension");

	if (ms::contains_icase(governing_equation, "Linear_Advection"))
	{
		if (space_dimension == 2)
		{
			const auto x_advection_speed = config.get<double>("x_advection_speed");
			const auto y_advection_speed = config.get<double>("y_advection_speed");
			return std::make_unique<Linear_Advection_2D>(x_advection_speed,y_advection_speed);
		}
		else
		{
			EXCEPTION("Linear advection does not support space dimension in configuration file");
			return nullptr;
		}
	}
	else if (ms::contains_icase(governing_equation, "Euler")) 
	{
		if (space_dimension == 2)
		{
			return std::make_unique<Euler_2D>();
		}
		else
		{
			EXCEPTION("Euler does not support space dimension in configuration file");
			return nullptr;
		}
	}
	else
	{
		EXCEPTION("governing equation in configuration file is not supported");
		return nullptr;
	}
}