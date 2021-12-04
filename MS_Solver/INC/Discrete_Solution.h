#pragma once
#include "Grid.h"
#include "Governing_Equation.h"
#include "Initial_Condition.h"

class Discrete_Solution
{
public:
	Discrete_Solution(const std::shared_ptr<Governing_Equation>& governing_equation, const Grid& grid);

public://Command
	virtual void precalculate_post_elements(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_element_center_points) abstract;
	virtual void precalculate_post_points(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_points) abstract;
	void update_solution(const Euclidean_Vector& updated_solution_v);
	void update_solution(Euclidean_Vector&& updated_solution_v);

public://Query
	virtual std::vector<Euclidean_Vector> calculate_solution_at_post_element_centers(const uint cell_index) const abstract;
	virtual std::vector<Euclidean_Vector> calculate_solution_at_post_points(const uint cell_index) const abstract;
	const std::vector<std::string>& get_solution_names(void) const;
	Euclidean_Vector solution_vector(void) const;
	Euclidean_Vector_Constant_Wrapper solution_vector_constant_wrapper(void) const;
	Euclidean_Vector_Wrapper solution_vector_wrapper(void);
	ushort num_equations(void) const;
	ushort num_solutions(void) const;
	size_t num_values(void) const;

protected:
	ushort space_dimension_;
	ushort num_equations_;	
	size_t num_cells_;
	std::vector<double> values_;
	std::shared_ptr<Governing_Equation> governing_equation_;
};

class Discrete_Solution_DG : public Discrete_Solution
{
public:	
	Discrete_Solution_DG(const std::shared_ptr<Governing_Equation>& governing_equation, const Grid& grid, const Initial_Condition& initial_condition, const ushort solution_degree);

public://Command
	void precalculate_post_elements(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_element_center_points) override;
	void precalculate_post_points(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_points) override;

	void precalculate_basis_bdry_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<Quadrature_Rule>& quadrature_rules);
	void precalcualte_cell_QPs_basis_values(const std::vector<Quadrature_Rule>& quadrature_rules);
	void precalculate_cell_P0_basis_values(void);
	void precalculate_infs_ocs_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ocs_QPs);
	void precalculate_infs_ncs_QPs_basis_values(const std::vector<uint>& nc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ncs_QPs);

public://Query	
	std::vector<Euclidean_Vector> calculate_solution_at_post_element_centers(const uint cell_index) const override;
	std::vector<Euclidean_Vector> calculate_solution_at_post_points(const uint cell_index) const override;

	std::vector<Euclidean_Vector> calculate_P0_solutions(void) const;
	std::vector<Euclidean_Vector> calculate_solution_at_bdry_QPs(const uint bdry_index, const uint oc_index) const;
	std::vector<Euclidean_Vector> calculate_solution_at_cell_QPs(const uint cell_index) const;
	std::vector<Euclidean_Vector> calculate_solution_at_infc_ocs_QPs(const uint infs_index, const uint oc_index) const;
	std::vector<Euclidean_Vector> calculate_solution_at_infc_ncs_QPs(const uint infs_index, const uint nc_index) const;

	Euclidean_Vector calculate_basis_point_v(const uint cell_index, const Euclidean_Vector& point) const;
	Matrix_Function<Polynomial> calculate_tranposed_gradient_basis(const uint cell_index) const;

	const std::vector<size_t>& get_coefficient_start_indexes(void) const;
	ushort solution_degree(const uint cell_index) const;
	ushort num_basis(const uint cell_index) const;
	ushort maximum_solution_degree(void) const;


private:
	std::vector<double> calculate_initial_values(const Grid& grid, const Initial_Condition& initial_condition) const;

	double calculate_P0_basis_value(const uint cell_index) const;
	Matrix calculate_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const;

	Euclidean_Vector calculate_P0_solution_precalculated(const uint cell_index) const;
	std::vector<Euclidean_Vector> calculate_solution_at_precalulated_points(const uint cell_index, const Matrix& basis_points_m) const;

	size_t coefficient_start_index(const uint cell_index) const;
	size_t num_total_basis(void) const;

	const std::vector<ushort>& get_solution_degrees(void) const;
	const std::vector<ushort>& get_set_of_num_basis(void) const;

	const double* coefficient_pointer(const uint cell_index) const;
	Euclidean_Vector P0_coefficient_v(const uint cell_index) const;	
	Matrix_Constant_Wrapper coefficient_matrix_contant_wrapper(const uint cell_index) const;
	
private:
	std::vector<ushort> solution_degrees_;
	std::vector<Vector_Function<Polynomial>> basis_vector_functions_;
	std::vector<ushort> set_of_num_basis_;
	std::vector<size_t> coefficieint_start_indexes_;

	//precalculated
	std::vector<Matrix> set_of_basis_post_element_center_points_m_;
	std::vector<Matrix> set_of_basis_post_points_m_;

	std::vector<Matrix> set_of_bdry_basis_QPs_m_;		//boudnary quadrature point basis value matrix
	std::vector<double> cell_P0_basis_values_;
	std::vector<Matrix> set_of_cell_basis_QPs_m_;		//cell quadrature point basis value matrix
	std::vector<Matrix> set_of_infc_basis_ocs_QPs_m_;	//inner face owner cell side quadratue point basis value matrix
	std::vector<Matrix> set_of_infc_basis_ncs_QPs_m_;	//inner face neighbor cell side quadratue point basis value matrix
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