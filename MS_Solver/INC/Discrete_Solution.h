#pragma once
#include "Grid.h"
#include "Governing_Equation.h"
#include "Initial_Condition.h"

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
	Euclidean_Vector calculate_basis_point_v(const uint cell_index, const Euclidean_Vector& point) const;
	Matrix_Function<Polynomial> calculate_tranposed_gradient_basis(const uint cell_index) const;

	std::vector<Euclidean_Vector> calculate_P0_solutions(const std::vector<double>& P0_basis_values) const;
	std::vector<Euclidean_Vector> calculate_solution_at_points(const uint cell_index, const Matrix& basis_points_m) const;

	size_t coefficient_start_index(const uint cell_index) const;
	ushort num_basis(const uint cell_index) const;
	ushort num_equations(void) const;
	ushort solution_degree(const uint cell_index) const;

	const std::vector<size_t>& get_coefficient_start_indexes(void) const;
	const std::vector<ushort>& get_solution_degrees(void) const;
	const std::vector<ushort>& get_set_of_num_basis(void) const;

private:
	Euclidean_Vector P0_coefficient_v(const uint cell_index) const;
	const double* coefficient_pointer(const uint cell_index) const;
	Matrix_Wrapper coefficient_m(const uint cell_index) const;
	size_t num_total_basis(void) const;
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


//namespace ms 
//{
//	template <typename T>
//	void insert(std::vector<T>& vec, const std::vector<T>& other_vec) 
//	{
//		vec.insert(vec.end(), other_vec.begin(), other_vec.end());
//	}
//}