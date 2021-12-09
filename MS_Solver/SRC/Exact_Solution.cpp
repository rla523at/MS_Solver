#include "../INC/Exact_Solution.h"

Euclidean_Vector Linear_Advection_Exact_Soltuion::calculate_exact_solution_vector(const Euclidean_Vector& point, const double time) const
{
	const auto initial_point = point - time * this->advection_speed_vector_;	
	const auto initial_point_consider_pbdry = this->consider_pbdry(initial_point);

	return this->initial_condition_->calculate_solution(initial_point_consider_pbdry);
}

std::vector<Euclidean_Vector> Linear_Advection_Exact_Soltuion::calculate_exact_solution_vectors(const std::vector<Euclidean_Vector>& points, const double time) const
{
	const auto num_points = points.size();
	std::vector<Euclidean_Vector> exact_solution_vs(num_points);

	for (ushort i = 0; i < num_points; ++i)
	{
		exact_solution_vs[i] = this->calculate_exact_solution_vector(points[i], time);
	}

	return exact_solution_vs;
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
	const auto& governing_equation_name = configuration.get_governing_equation();

	if (ms::contains_icase(governing_equation_name, "Linear", "Advection"))
	{
		const auto space_dimension = configuration.space_dimension();
		const auto periodic_lengths_ptr = configuration.periodic_lengths_ptr();
		const auto advection_speeds_ptr = configuration.advection_speeds_ptr();

		std::vector<double> advection_speeds(advection_speeds_ptr, advection_speeds_ptr + space_dimension);
		std::vector<double> periodic_lengths(periodic_lengths_ptr, periodic_lengths_ptr + space_dimension);
		auto initial_condition = Initial_Condition_Factory::make_unique(configuration);

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