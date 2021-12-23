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

	Euclidean_Vector_Wrapper discrete_solution_vector_wrapper(void);

public://Query
	virtual std::vector<Euclidean_Vector> calculate_solution_at_post_element_centers(const uint cell_index) const abstract;
	virtual std::vector<Euclidean_Vector> calculate_solution_at_post_points(const uint cell_index) const abstract;

	Euclidean_Vector discrete_solution_vector(void) const;
	Constant_Euclidean_Vector_Wrapper discrete_solution_constant_vector_wrapper(void) const;
	const std::vector<std::string>& get_solution_names(void) const;
	ushort num_equations(void) const;
	ushort num_solutions(void) const;
	size_t num_total_values(void) const;

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

	void limit_slope(const uint cell_index, const double limiting_value);
	void project_to_Pn_space(const uint cell_index, const ushort Pn);

	//for time step calculation
	void precalculate_cell_P0_basis_values(void);

	//for DG RHS calculation
	void precalculate_basis_bdry_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<Quadrature_Rule>& quadrature_rules);
	void precalcualte_cell_QPs_basis_values(const std::vector<Quadrature_Rule>& quadrature_rules);
	void precalculate_infs_ocs_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ocs_QPs);
	void precalculate_infs_ncs_QPs_basis_values(const std::vector<uint>& nc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ncs_QPs);
	
	//for simplex decomposed MLP criterion
	void precalculate_set_of_simplex_P0_P1_projection_basis_vertices_m(const Grid& grid);

	//for MLP Indicator
	void precalculate_cell_vertices_basis_values(const std::vector<std::vector<Euclidean_Vector>>& set_of_verticies);

	//for troubled boundary indicator
	void precalculate_infs_ocs_jump_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ocs_jump_QPs);
	void precalculate_infs_ncs_jump_QPs_basis_values(const std::vector<uint>& nc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ncs_jump_QPs);	

	//for discontinuity indicator
	void precalculate_set_of_cell_index_to_target_cell_basis_QPs_m_(const std::vector<Quadrature_Rule>& quadrature_rules, const std::vector<std::vector<uint>>& set_of_face_neighbor_cell_indexes);

public://Query	
	std::vector<Euclidean_Vector> calculate_solution_at_post_element_centers(const uint cell_index) const override;
	std::vector<Euclidean_Vector> calculate_solution_at_post_points(const uint cell_index) const override;

	Euclidean_Vector calculate_basis_point_v(const uint cell_index, const Euclidean_Vector& point) const;
	Matrix_Function<Polynomial> calculate_tranposed_gradient_basis(const uint cell_index) const;
	const std::vector<size_t>& get_coefficient_start_indexes(void) const;
	ushort solution_degree(const uint cell_index) const;
	ushort num_basis(const uint cell_index) const;
	ushort num_values(const uint cell_index) const;
	ushort maximum_solution_degree(void) const;

	//for time step calculation
	std::vector<Euclidean_Vector> calculate_P0_solutions(void) const;	

	//for DG RHS calculation
	std::vector<Euclidean_Vector> calculate_solution_at_bdry_QPs(const uint bdry_index, const uint oc_index) const;
	std::vector<Euclidean_Vector> calculate_solution_at_cell_QPs(const uint cell_index) const;
	std::vector<Euclidean_Vector> calculate_solution_at_infc_ocs_QPs(const uint infs_index, const uint oc_index) const;
	std::vector<Euclidean_Vector> calculate_solution_at_infc_ncs_QPs(const uint infs_index, const uint nc_index) const;
	void calculate_solution_at_bdry_QPs(Euclidean_Vector* solution_at_QPs, const uint bdry_index, const uint oc_index) const;
	void calculate_solution_at_cell_QPs(Euclidean_Vector* solution_at_QPs, const uint cell_index) const;
	void calculate_solution_at_infc_ocs_QPs(Euclidean_Vector* solution_at_infc_ocs_QPs, const uint infs_index, const uint oc_index) const;
	void calculate_solution_at_infc_ncs_QPs(Euclidean_Vector* solution_at_infc_ncs_QPs, const uint infs_index, const uint nc_index) const;

	//for MLP criterion
	double calculate_P0_nth_solution(const uint cell_index, const ushort equation_index) const;
	std::vector<double> calculate_P1_projected_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const;
	std::vector<double> calculate_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const;
	void calculate_P1_projected_nth_solution_at_vertices(double* P1_projected_nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const;
	void calculate_nth_solution_at_vertices(double* nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const;

	//for simplex decomposed MLP criterion
	std::vector<double> calculate_simplex_P0_projected_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const;
	std::vector<double> calculate_simplex_P1_projected_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const;
	void calculate_simplex_P0_projected_nth_solution_at_vertices(double* simplex_P0_projected_nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const;
	void calculate_simplex_P1_projected_nth_solution_at_vertices(double* simplex_P1_projected_nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const;

	//for subcell oscillation indicator
	void calculate_nth_solution_at_infc_ocs_jump_QPs(double* nth_solution_at_infc_ocs_jump_QPs, const uint infc_index, const uint oc_index, const ushort equation_index) const;
	void calculate_nth_solution_at_infc_ncs_jump_QPs(double* nth_solution_at_infc_ncs_jump_QPs, const uint infc_index, const uint nc_index, const ushort equation_index) const;

	//for discontinuity indicator	
	void calculate_nth_solution_at_target_cell_QPs(double* nth_solution_at_target_cell_QPs, const uint target_cell_index, const uint my_cell_index, const ushort equation_index) const;

private:	
	std::vector<double> calculate_initial_values(const Grid& grid, const Initial_Condition& initial_condition) const;

	Matrix calculate_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const;
	Vector_Function<Polynomial> calculate_simplex_Pn_projection_basis_vector_function(const uint cell_index, const ushort Pn, const Geometry& sub_simplex_geometry) const;
	Constant_Matrix_Wrapper Pn_projected_basis_constant_matrix_wrapper(const Constant_Matrix_Wrapper& basis_points_m, const ushort Pn) const;
	size_t num_total_basis(void) const;

	double* coefficient_pointer(const uint cell_index);
	const double* coefficient_pointer(const uint cell_index) const;
	Matrix_Wrapper coefficient_matrix_wrapper(const uint cell_index);
	double P0_nth_coefficient(const uint cell_index, const ushort equation_index) const;
	Euclidean_Vector P0_coefficient_v(const uint cell_index) const;	
	Constant_Matrix_Wrapper coefficient_matrix_contant_wrapper(const uint cell_index) const;
	Constant_Matrix_Wrapper nth_coefficient_matrix_contant_wrapper(const uint cell_index, const ushort equation_index) const;
	Constant_Matrix_Wrapper Pn_projected_mth_coefficient_contant_matrix_wrapper(const uint cell_index, const ushort Pn, const ushort equation_index) const;

	std::vector<Euclidean_Vector> calculate_solution_at_precalulated_points(const uint cell_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	std::vector<double> calculate_nth_solution_at_precalulated_points(const uint cell_index, const ushort equation_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	std::vector<double> calculate_Pn_projected_mth_solution_at_precalulated_points(const uint cell_index, const ushort Pn, const ushort equation_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	void calculate_solution_at_precalulated_points(Euclidean_Vector* solution_v_at_points, const uint cell_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	void calculate_nth_solution_at_precalulated_points(double* nth_solution_at_points, const uint cell_index, const ushort equation_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	void calculate_Pn_projected_mth_solution_at_precalulated_points(double* Pn_projected_mth_solution_at_points, const uint cell_index, const ushort Pn, const ushort equation_index, const Constant_Matrix_Wrapper& basis_points_m) const;

private:
	std::vector<ushort> set_of_num_values_;
	std::vector<ushort> solution_degrees_;
	std::vector<Vector_Function<Polynomial>> basis_vector_functions_;
	std::vector<ushort> set_of_num_basis_;
	std::vector<size_t> coefficieint_start_indexes_;

	//precalculated
	static constexpr ushort max_solution_degree = 20;
	std::array<ushort, max_solution_degree> degree_to_num_basis_table = { 0 };

	std::vector<Matrix> set_of_basis_post_element_center_points_m_;
	std::vector<Matrix> set_of_basis_post_points_m_;

	std::vector<Matrix> set_of_bdry_basis_QPs_m_;		//boudnary quadrature point basis value matrix
	std::vector<double> cell_P0_basis_values_;
	std::vector<Matrix> set_of_cell_basis_QPs_m_;		//cell quadrature point basis value matrix
	std::vector<Matrix> set_of_cell_basis_vertices_m_;	
	std::vector<Matrix> set_of_infc_basis_ocs_flux_QPs_m_;	//inner face owner cell side flux quadratue point basis value matrix
	std::vector<Matrix> set_of_infc_basis_ncs_flux_QPs_m_;	//inner face neighbor cell side flux quadratue point basis value matrix

	std::vector<Matrix> set_of_simplex_P0_projected_basis_vertices_m_;
	std::vector<Matrix> set_of_simplex_P1_projected_basis_vertices_m_;
	std::vector<Matrix> set_of_infc_basis_ocs_jump_QPs_m_;	//inner face owner cell side jump quadratue point basis value matrix
	std::vector<Matrix> set_of_infc_basis_ncs_jump_QPs_m_;	//inner face neighbor cell side jump quadratue point basis value matrix

	std::vector<std::map<uint, Matrix>> set_of_cell_index_to_target_cell_basis_QPs_m_;

	//optimize construction
	static constexpr ushort max_num_equation = 5;
	static constexpr ushort max_num_precalculated_points = 200;
	mutable Euclidean_Vector GE_soluion_;
	mutable Euclidean_Vector solution_;
	mutable std::array<double, max_num_equation * max_num_precalculated_points> solution_at_points_values_ = { 0 };
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