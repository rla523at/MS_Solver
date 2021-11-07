#pragma once
#include "Configuration.h"
#include "Grid.h"
#include "Governing_Equation.h"
#include "Initial_Condition.h"

using ushort = unsigned short;
using uint = unsigned int;

class Discrete_Solution
{
public:
	Discrete_Solution(const Governing_Equation& governing_equation, const Grid& grid);

public://Command
	virtual void set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) abstract;

public://Query
	virtual std::vector<Euclidean_Vector> calculate_P0_solutions(void) const abstract;
	virtual std::vector<std::vector<Euclidean_Vector>> calculate_set_of_post_point_solutions(void) const abstract;
	const std::vector<std::string>& get_variable_names(void) const;

protected:
	ushort num_equations_;
	std::vector<std::string> solution_variable_names_;
	size_t num_cells_;
	std::vector<double> discretized_solutions_;
};


class Discretized_Solution_FVM : public Discrete_Solution
{
public:
	Discretized_Solution_FVM(const Governing_Equation& governing_equation, const Grid& grid, const Initial_Condition& initial_condition);

public://Command
	void set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) override;

public: //Query
	std::vector<std::vector<Euclidean_Vector>> calculate_set_of_post_point_solutions(void) const override;

private:
	Euclidean_Vector calculate_solution_at_center(const uint icell) const;
};


class Discrete_Solution_HOM : public Discrete_Solution
{
public:
	Discrete_Solution_HOM(const Governing_Equation& governing_equation, const Grid& grid, const Initial_Condition& initial_condition, const ushort solution_order);

public://Command
	void set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) override;

public://Query
	std::vector<std::vector<Euclidean_Vector>> calculate_set_of_post_point_solutions(void) const override;

private:
	double* coefficient_pointer(const uint icell);

	Matrix_Wrapper get_coefficient_matrix(const uint icell) const;
	std::vector<Euclidean_Vector> calculate_post_point_solutions(const uint icell) const;
	Euclidean_Vector calculate_basis_vector_value(const uint icell, const Euclidean_Vector& node) const;
	const double* coefficient_pointer(const uint icell) const;	
	size_t num_total_basis(void) const;

private:
	ushort solution_order_;
	std::vector<Vector_Function<Polynomial>> basis_vector_functions_;
	std::vector<size_t> set_of_num_basis_;
	std::vector<size_t> coeffcieint_start_indexes_;

	size_t num_post_points_;	
	std::vector<Matrix> set_of_basis_post_points_m_;
};

//class Discretized_Solution_Factory
//{
//public:
//	static std::unique_ptr<Discretized_Solution> make(const Configuration& configuration, const Grid& grid, const Initial_Condition& initial_condition) {
//		const auto governing_eqaution = Governing_Equation_Factory::make(configuration);
//
//		const auto spatial_discretization_method = configuration.get("spatial_discretization_method");
//		
//		if (ms::contains_icase(spatial_discretization_method, "FVM")) 
//		{
//			return std::make_unique<Discretized_Solution_FVM>(*governing_eqaution);
//		}
//		else if (ms::contains_icase(spatial_discretization_method, "HOM")) 
//		{
//			const auto solution_order = configuration.get<ushort>("solution_order");
//			return std::make_unique<Discretized_Solution_HOM>(*governing_eqaution, solution_order);
//		}
//		else
//			EXCEPTION("discretization method in configuration file was not supported");
//	}
//};


namespace ms {
	template <typename T>
	void insert(std::vector<T>& vec, const std::vector<T>& other_vec) {
		vec.insert(vec.end(), other_vec.begin(), other_vec.end());
	}
}