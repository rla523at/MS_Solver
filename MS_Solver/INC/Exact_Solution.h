#pragma once
#include "Configuration.h"
#include "Initial_Condition.h"
#include "Initial_Condition_Impl.h"

class Exact_Solution
{
public:
	virtual Euclidean_Vector calculate_exact_solution_vector(const Euclidean_Vector& point, const double time) const abstract;
	virtual std::vector<Euclidean_Vector> calculate_exact_solution_vectors(const std::vector<Euclidean_Vector>& points, const double time) const abstract;
};

class Linear_Advection_Exact_Soltuion : public Exact_Solution
{
public:
	Linear_Advection_Exact_Soltuion(std::vector<double>&& advection_speed_vector, std::unique_ptr<Initial_Condition>&& initial_condition, std::vector<double>&& periodic_length)
		: advection_speed_vector_(std::move(advection_speed_vector))
		, initial_condition_(std::move(initial_condition))
		, periodic_lengths_(std::move(periodic_length)) {};

public:
	Euclidean_Vector calculate_exact_solution_vector(const Euclidean_Vector& point, const double time) const override;
	std::vector<Euclidean_Vector> calculate_exact_solution_vectors(const std::vector<Euclidean_Vector>& points, const double time) const override;

private:
	Euclidean_Vector consider_pbdry(const Euclidean_Vector& point) const;

private:
	static constexpr ushort num_equation_ = 1;
	Euclidean_Vector advection_speed_vector_;
	std::unique_ptr<Initial_Condition> initial_condition_;
	std::vector<double> periodic_lengths_;
};

class Exact_Solution_Factory//static class
{
public:
	static std::unique_ptr<Exact_Solution> make_unique(const Configuration& configuration);

private:
	Exact_Solution_Factory(void) = delete;
};

namespace ms
{
	double remainder(const double dividend, const double divisor);
	int quotient(const double dividend, const double divisor);
}

