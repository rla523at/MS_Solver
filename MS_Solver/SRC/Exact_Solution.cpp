#include "../INC/Exact_Solution.h"

Euclidean_Vector Linear_Advection_Exact_Soltuion::calculate_exact_solution_vector(const Euclidean_Vector& point, const double time) const
{
	std::vector<double> exact_values(num_equation_);

	const auto initial_point = point - time * this->advection_speed_vector_;	
	const auto initial_point_consider_pbdry = this->consider_pbdry(initial_point);

	return this->initial_condition_->calculate_solution(initial_point_consider_pbdry);
}

Euclidean_Vector Linear_Advection_Exact_Soltuion::consider_pbdry(const Euclidean_Vector& point) const
{
	// Requirement : domain should be [0, this->periodic_length[0]] x [0, this->periodic_length[1]] x [0, this->periodic_length[2]]

	const auto space_dimension = point.size();
	std::vector<double> coordinate_values_consider_pbdry(space_dimension);

	for (ushort i = 0; i < space_dimension; ++i)
	{
		coordinate_values_consider_pbdry[i] = ms::remainder(point[i], this->periodic_lengths_[i]);
	}

	return coordinate_values_consider_pbdry;
}

std::unique_ptr<Exact_Solution> Exact_Solution_Factory::make_unique(const Configuration& configuration)
{
	const auto governing_equation_name = configuration.get("governing_equation");
	const auto initial_condition_name = configuration.get("initial_condition");

	if (ms::contains_icase(governing_equation_name, "Linear", "Advection") &&
		ms::contains_icase(initial_condition_name, "sine", "wave"))
	{
		const auto space_dimension = configuration.get<ushort>("space_dimension");
		const auto initial_condition = Initial_Condition_Factory::make_unique(configuration);

		std::vector<double> periodic_lengths(space_dimension);
		std::vector<double> advection_speeds(space_dimension);

		if (space_dimension == 2)
		{
			periodic_lengths[0] = configuration.get<double>("x_periodic_lenght");
			periodic_lengths[1] = configuration.get<double>("y_periodic_lenght");

			advection_speeds[0] = configuration.get<double>("x_advection_speed");
			advection_speeds[1] = configuration.get<double>("y_advection_speed");
		}
		else if (space_dimension == 3)
		{
			periodic_lengths[0] = configuration.get<double>("x_periodic_lenght");
			periodic_lengths[1] = configuration.get<double>("y_periodic_lenght");
			periodic_lengths[2] = configuration.get<double>("z_periodic_lenght");

			advection_speeds[0] = configuration.get<double>("x_advection_speed");
			advection_speeds[1] = configuration.get<double>("y_advection_speed");
			advection_speeds[2] = configuration.get<double>("z_advection_speed");
		}
		else
		{
			EXCEPTION("space dimension in configuration file is not supported");
		}

		return std::make_unique<Linear_Advection_Exact_Soltuion>(std::move(advection_speeds), std::move(initial_condition), std::move(periodic_lengths));
	}
	else
	{
		EXCEPTION("Exact solution does not provided in this configuration");
		return nullptr;
	}
}


namespace ms
{
	double remainder(const double dividend, const double divisor)
	{
		REQUIRE(divisor >= 0, "divisor should be positive");

		const auto quotient = std::floor(dividend / divisor);
		return dividend - divisor * quotient;
	}

	int quotient(const double dividend, const double divisor)
	{
		REQUIRE(divisor >= 0, "divisor should be positive");

		return static_cast<int>(std::floor(dividend / divisor));
	}
}