#pragma once
#include "Grid.h"
#include "Governing_Equation.h"
#include "Initial_Condition.h"

using ushort = unsigned short;
using uint = unsigned int;

class Discrete_Solution
{
public:
	Discrete_Solution(const Grid& grid, const std::shared_ptr<Governing_Equation>& governing_equation);

public://Command
	void update_solution(Euclidean_Vector&& updated_solution);

public://Query
	const Euclidean_Vector& get_solution_vector(void) const;
	size_t num_values(void) const;

protected:
	ushort space_dimension_;
	ushort num_equations_;
	size_t num_cells_;
	Euclidean_Vector value_v_;
	std::shared_ptr<Governing_Equation> governing_equation_;
};

class Discrete_Solution_DG : public Discrete_Solution
{
public:	
	Discrete_Solution_DG(const Grid& grid, const std::shared_ptr<Governing_Equation>& governing_equation, const Initial_Condition& initial_condition, const ushort solution_degree);

public://Query
	double calculate_P0_basis_value(const uint cell_index) const;
	Matrix calculate_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const;
	Matrix_Function<Polynomial> calculate_tranposed_gradient_basis(const uint cell_index) const;
	std::vector<Euclidean_Vector> calculate_P0_solutions(const std::vector<double>& P0_basis_values) const;
	std::vector<Euclidean_Vector> calculate_solution_at_points(const uint cell_index, const Matrix& basis_points_m) const;

	size_t coefficient_start_index(const uint cell_index) const;
	ushort num_basis(const uint cell_index) const;
	ushort solution_degree(const uint cell_index) const;

	const std::vector<size_t>& get_coefficient_start_indexes(void) const;
	const std::vector<ushort>& get_solution_degrees(void) const;
	const std::vector<ushort>& get_set_of_num_basis(void) const;

private:
	Euclidean_Vector calculate_basis_vector_value(const uint cell_index, const Euclidean_Vector& node) const;

	const double* coefficient_pointer(const uint cell_index) const;
	Matrix_Wrapper coefficient_matrix(const uint cell_index) const;
	size_t num_total_basis(void) const;
	Euclidean_Vector P0_coefficient_vector(const uint cell_index) const;
	Euclidean_Vector calculate_initial_values(const Grid& grid, const Initial_Condition& initial_condition) const;
	

private:
	std::vector<ushort> solution_degrees_;
	std::vector<Vector_Function<Polynomial>> basis_vector_functions_;
	std::vector<ushort> set_of_num_basis_;
	std::vector<size_t> coefficieint_start_indexes_;
};

//class Discretized_Solution_FVM : public Discrete_Solution
//{
//public:
//	Discretized_Solution_FVM(const Governing_Equation& governing_equation, const Grid& grid, const Initial_Condition& initial_condition);
//
//public://Command
//	void set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition);
//
//public: //Query
//	std::vector<std::vector<Euclidean_Vector>> calculate_set_of_post_point_solutions(void) const;
//
//private:
//	Euclidean_Vector calculate_solution_at_center(const uint cell_index) const;
//};


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

//class Discrete_Solution
//{
//public:
//	Discrete_Solution(const Governing_Equation& governing_equation, const Grid& grid);
//
//	//public://Command
//		//virtual void set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) abstract;
//
//	//public://Query
//		//virtual std::vector<Euclidean_Vector> calculate_P0_solutions(void) const abstract;
//		//virtual std::vector<std::vector<Euclidean_Vector>> calculate_set_of_post_point_solutions(void) const abstract;
//
//public://Command
//	void update_solution(Euclidean_Vector&& updated_solution);
//
//public://Query
//	const std::vector<std::string>& get_variable_names(void) const;
//	const Euclidean_Vector& get_solution_vector(void) const;
//	size_t num_values(void) const;
//
//protected:
//	ushort num_equations_;
//	std::vector<std::string> solution_variable_names_;
//	size_t num_cells_;
//	Euclidean_Vector value_v_;
//};


namespace ms {
	template <typename T>
	void insert(std::vector<T>& vec, const std::vector<T>& other_vec) {
		vec.insert(vec.end(), other_vec.begin(), other_vec.end());
	}
}