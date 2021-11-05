#pragma once
#include "Configuration.h"
#include "Grid.h"
#include "Governing_Equation.h"
#include "Initial_Condition.h"

using ushort = unsigned short;
using uint = unsigned int;

class Discretized_Solution
{
public:
	Discretized_Solution(const Governing_Equation& governing_equation);

public://Command
	virtual void set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) abstract;

public://Query
	virtual std::vector<std::vector<Euclidean_Vector>> calculate_set_of_post_point_solutions(void) const abstract;
	const std::vector<std::string>& get_variable_names(void) const;

protected:
	ushort num_equations_;
	std::vector<std::string> solution_variable_names_;
	size_t num_cells_;
	std::vector<double> discretized_solutions_;
};


class Discretized_Solution_FVM : public Discretized_Solution
{
public:
	Discretized_Solution_FVM(const Governing_Equation& governing_equation) 
		: Discretized_Solution(governing_equation) {};

public://Command
	void set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) override;

public: //Query
	std::vector<std::vector<Euclidean_Vector>> calculate_set_of_post_point_solutions(void) const override;

private:
	Euclidean_Vector calculate_solution_at_center(const uint icell) const;
};


class Discretized_Solution_HOM : public Discretized_Solution
{
public:
	Discretized_Solution_HOM(const Governing_Equation& governing_equation, const ushort solution_order)
		:Discretized_Solution(governing_equation),
		solution_order_(solution_order) {};

public://Command
	void set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) override;

public:
	std::vector<std::vector<Euclidean_Vector>> calculate_set_of_post_point_solutions(void) const override;

private:
	Matrix_Wrapper get_coefficient_matrix_wrapper(const uint icell) const;
	std::vector<Euclidean_Vector> calculate_post_point_solutions(const uint icell) const;

private:
	ushort solution_order_;
	std::vector<Vector_Function<Polynomial>> basis_vector_functions_;


	size_t num_post_points_;	
	std::vector<size_t> set_of_data_start_index_;

	std::vector<Matrix> set_of_basis_post_points_;
};

class Discretized_Solution_Factory
{
public:
	static std::unique_ptr<Discretized_Solution> make(const Configuration& configuration) {
		const auto governing_eqaution = Governing_Equation_Factory::make(configuration);

		const auto spatial_discretization_method = configuration.get("spatial_discretization_method");
		
		if (ms::contains_icase(spatial_discretization_method, "FVM")) 
		{
			return std::make_unique<Discretized_Solution_FVM>(*governing_eqaution);
		}
		else if (ms::contains_icase(spatial_discretization_method, "HOM")) 
		{
			const auto solution_order = configuration.get("solution_order");
			return std::make_unique<Discretized_Solution_HOM>(*governing_eqaution, solution_order);
		}
		else
			EXCEPTION("discretization method in configuration file was not supported");
	}
};


namespace ms {
	template <typename T>
	void insert(std::vector<T>& vec, const std::vector<T>& other_vec) {
		vec.insert(vec.end(), other_vec.begin(), other_vec.end());
	}
}