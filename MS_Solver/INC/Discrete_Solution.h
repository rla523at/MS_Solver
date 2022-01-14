#pragma once
#include "Grid.h"
#include "Governing_Equation.h"
#include "Initial_Condition.h"

class Discrete_Solution
{
public:
	Discrete_Solution(const std::shared_ptr<Governing_Equation>& governing_equation, const Grid& grid);

public://Command
	virtual void precalculate_cell_post_element_center_points_basis_values(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_element_center_points) abstract;
	virtual void precalculate_cell_post_point_basis_values(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_points) abstract;	

	Euclidean_Vector_Wrapper discrete_solution_vector_wrapper(void);

public://Query
	virtual std::vector<Euclidean_Vector> calculate_solution_at_post_element_centers(const uint cell_index) const abstract;
	virtual std::vector<Euclidean_Vector> calculate_solution_at_post_points(const uint cell_index) const abstract;
	virtual void calculate_solution_at_post_points(Euclidean_Vector* solution_at_post_points, const uint cell_index) const abstract;

	Euclidean_Vector discrete_solution_vector(void) const;
	Constant_Euclidean_Vector_Wrapper discrete_solution_constant_vector_wrapper(void) const;
	const std::vector<std::string>& get_solution_names(void) const;
	uint num_cells(void) const { return this->num_cells_; };
	ushort num_equations(void) const;
	ushort num_solutions(void) const;
	size_t num_total_values(void) const;

protected:
	ushort space_dimension_;
	ushort num_equations_;	
	uint num_cells_;
	std::vector<double> values_;
	std::shared_ptr<Governing_Equation> governing_equation_;
};

class Discrete_Solution_DG : public Discrete_Solution
{
public:	
	Discrete_Solution_DG(void) = default;
	Discrete_Solution_DG(const std::shared_ptr<Governing_Equation>& governing_equation, const Grid& grid, const Initial_Condition& initial_condition, const ushort solution_degree);

public://Command

	void limit_slope(const uint cell_index, const double limiting_value);
	void project_to_Pn_space(const uint cell_index, const ushort Pn);

	//for boundary
	void precalculate_bdry_RHS_QPs_basis_values(const Grid& grid);

	//for cell
	void precalculate_cell_post_element_center_points_basis_values(const std::vector<std::vector<Euclidean_Vector>>& cell_index_to_post_element_center_points) override;
	void precalculate_cell_post_point_basis_values(const std::vector<std::vector<Euclidean_Vector>>& cell_index_to_post_points) override;
	void precalculate_cell_P0_basis_values(void);
	void precalculate_cell_RHS_QPs_basis_values(const std::vector<Quadrature_Rule>& quadrature_rules);
	void precalculate_cell_RHS_QPs_ddx_basis_values(const std::vector<std::vector<Euclidean_Vector>>& cell_index_to_QPs);
	void precalculate_cell_RHS_QPs_ddy_basis_values(const std::vector<std::vector<Euclidean_Vector>>& cell_index_to_QPs);
	void precalculate_cell_jump_QPS_basis_values(const Grid& grid);
	void precalculate_cell_jump_QPs_face_neighbor_cell_basis_values(const Grid& grid);
	void precalculate_cell_vertices_basis_values(const std::vector<std::vector<Euclidean_Vector>>& set_of_verticies);
	void precalculate_set_of_simplex_P0_P1_projection_basis_vertices_m(const Grid& grid);

	//for inner face
	void precalculate_infs_ocs_RHS_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ocs_QPs);
	void precalculate_infs_ncs_RHS_QPs_basis_values(const std::vector<uint>& nc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ncs_QPs);
	void precalculate_infs_ocs_jump_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ocs_jump_QPs);
	void precalculate_infs_ncs_jump_QPs_basis_values(const std::vector<uint>& nc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ncs_jump_QPs);

public://Query	
	Euclidean_Vector calculate_basis_point_v(const uint cell_index, const Euclidean_Vector& point) const;
	Matrix_Function<Polynomial> calculate_tranposed_gradient_basis(const uint cell_index) const;
	const std::vector<size_t>& get_coefficient_start_indexes(void) const;
	ushort solution_degree(const uint cell_index) const;
	ushort num_basis(const uint cell_index) const;
	ushort num_values(const uint cell_index) const;
	ushort maximum_solution_degree(void) const;


	double calculate_P0_nth_solution(const uint cell_index, const ushort solution_index) const;
	std::vector<Euclidean_Vector> calculate_P0_solutions(void) const;

	//at post
	std::vector<Euclidean_Vector> calculate_solution_at_post_element_centers(const uint cell_index) const override;
	std::vector<Euclidean_Vector> calculate_solution_at_post_points(const uint cell_index) const override;
	void calculate_solution_at_post_points(Euclidean_Vector* solution_at_post_points, const uint cell_index) const override;

	//at boundary
	std::vector<Euclidean_Vector> calculate_solution_at_bdry_RHS_QPs(const uint bdry_index, const uint oc_index) const;
	void calculate_solution_at_bdry_RHS_QPs(Euclidean_Vector* solution_at_QPs, const uint bdry_index, const uint oc_index);

	//at cell
	std::vector<Euclidean_Vector> calculate_solution_at_cell_RHS_QPs(const uint cell_index) const;
	void calculate_solution_at_cell_RHS_QPs(Euclidean_Vector* solution_at_QPs, const uint cell_index);
	void calculate_solution_at_cell_RHS_QPs(Euclidean_Vector* solution_at_QPs, const uint cell_index) const;
	void calculate_nth_solution_at_cell_jump_QPs(double* solution_at_cell_jump_QPs, const uint cell_index, const ushort solution_index) const;
	void calculate_nth_extrapolated_solution_at_cell_jump_QPs(double* nth_solution_at_target_cell_QPs, const uint face_share_cell_index, const uint target_cell_index, const ushort solution_index) const;
	void calculate_ddx_GE_solution_at_cell_RHS_QPs(Euclidean_Vector* solution_at_QPs, const uint cell_index) const;
	void calculate_ddy_GE_solution_at_cell_RHS_QPs(Euclidean_Vector* solution_at_QPs, const uint cell_index) const;

	//at inner face
	std::vector<Euclidean_Vector> calculate_solution_at_infc_ocs_RHS_QPs(const uint infs_index, const uint oc_index) const;
	std::vector<Euclidean_Vector> calculate_solution_at_infc_ncs_RHS_QPs(const uint infs_index, const uint nc_index) const;
	void calculate_solution_at_infc_ocs_RHS_QPs(Euclidean_Vector* solution_at_infc_ocs_QPs, const uint infs_index, const uint oc_index);
	void calculate_solution_at_infc_ncs_RHS_QPs(Euclidean_Vector* solution_at_infc_ncs_QPs, const uint infs_index, const uint nc_index);
	void calculate_nth_solution_at_infc_ocs_jump_QPs(double* nth_solution_at_infc_ocs_jump_QPs, const uint infc_index, const uint oc_index, const ushort equation_index) const;
	void calculate_nth_solution_at_infc_ncs_jump_QPs(double* nth_solution_at_infc_ncs_jump_QPs, const uint infc_index, const uint nc_index, const ushort equation_index) const;

	//at vertices
	std::vector<double> calculate_P1_projected_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const;
	std::vector<double> calculate_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const;
	std::vector<double> calculate_simplex_P0_projected_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const;
	std::vector<double> calculate_simplex_P1_projected_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const;
	void calculate_P1_projected_nth_solution_at_vertices(double* P1_projected_nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const;
	void calculate_nth_solution_at_vertices(double* nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const;
	void calculate_simplex_P0_projected_nth_solution_at_vertices(double* simplex_P0_projected_nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const;
	void calculate_simplex_P1_projected_nth_solution_at_vertices(double* simplex_P1_projected_nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const;

private:	
	std::vector<double> calculate_initial_values(const Grid& grid, const Initial_Condition& initial_condition) const;

	Matrix calculate_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const;
	Matrix calculate_ddx_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const;
	Matrix calculate_ddy_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const;

	Vector_Function<Polynomial> calculate_simplex_Pn_projection_basis_vector_function(const uint cell_index, const ushort Pn, const Geometry& sub_simplex_geometry) const;
	Constant_Matrix_Wrapper Pn_projected_basis_constant_matrix_wrapper(const Constant_Matrix_Wrapper& basis_points_m, const ushort Pn) const;
	size_t num_total_basis(void) const;

	double* coefficient_pointer(const uint cell_index);
	const double* coefficient_pointer(const uint cell_index) const;
	Matrix_Wrapper coefficient_matrix_wrapper(const uint cell_index);
	double P0_nth_coefficient(const uint cell_index, const ushort equation_index) const;
	//Euclidean_Vector P0_coefficient_v(const uint cell_index) const;	
	void calculate_P0_coefficient_v(double* ptr, const uint cell_index) const;

	Constant_Matrix_Wrapper coefficient_constant_matrix_wrapper(const uint cell_index) const;
	Constant_Matrix_Wrapper nth_coefficient_matrix_contant_wrapper(const uint cell_index, const ushort equation_index) const;
	Constant_Matrix_Wrapper Pn_projected_mth_coefficient_contant_matrix_wrapper(const uint cell_index, const ushort Pn, const ushort equation_index) const;

	Constant_Matrix_Wrapper calculate_GE_solution_points_constant_matrix_wrapper(const uint cell_index, const Constant_Matrix_Wrapper & basis_points_m) const;

	std::vector<Euclidean_Vector> calculate_solution_at_precalulated_points(const uint cell_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	std::vector<double> calculate_nth_solution_at_precalulated_points(const uint cell_index, const ushort solution_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	std::vector<double> calculate_Pn_projected_mth_solution_at_precalulated_points(const uint cell_index, const ushort Pn, const ushort solution_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	void calculate_GE_solution_at_precalulated_points(Euclidean_Vector* solution_v_at_points, const uint cell_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	void calculate_solution_at_precalulated_points(Euclidean_Vector* solution_v_at_points, const uint cell_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	void calculate_nth_solution_at_precalulated_points(double* nth_solution_at_points, const uint cell_index, const ushort solution_index, const Constant_Matrix_Wrapper& basis_points_m) const;
	void calculate_Pn_projected_mth_solution_at_precalulated_points(double* Pn_projected_mth_solution_at_points, const uint cell_index, const ushort Pn, const ushort equation_index, const Constant_Matrix_Wrapper& basis_points_m) const;

	void scailing(const uint cell_index, const double scail_factor);


private:
	std::vector<ushort> set_of_num_values_;
	std::vector<ushort> cell_index_to_solution_degree_table_;
	std::vector<Vector_Function<Polynomial>> basis_vector_functions_;
	std::vector<ushort> set_of_num_basis_;
	std::vector<size_t> coefficieint_start_indexes_;

	double scaling_factor_ = 0.5;

	//precalculated
	static constexpr ushort max_solution_degree = 20;
	std::array<ushort, max_solution_degree> degree_to_num_basis_table = { 0 };

	std::vector<Matrix> set_of_basis_post_element_center_points_m_;

	std::vector<Matrix> boundary_index_to_basis_RHS_QPs_m_;		//boudnary quadrature point basis value matrix
	
	std::vector<Matrix> inner_face_index_to_basis_ocs_RHS_QPs_m_;	//inner face owner cell side flux quadratue point basis value matrix
	std::vector<Matrix> inner_face_index_to_basis_ncs_RHS_QPs_m_;	//inner face neighbor cell side flux quadratue point basis value matrix
	std::vector<Matrix> inner_face_index_to_basis_ocs_jump_QPs_m_;	//inner face owner cell side jump quadratue point basis value matrix
	std::vector<Matrix> inner_face_index_to_basis_ncs_jump_QPs_m_;	//inner face neighbor cell side jump quadratue point basis value matrix

	std::vector<double> cell_index_to_P0_basis_value_;
	std::vector<Matrix> cell_index_to_basis_post_points_m_table_;
	std::vector<Matrix> cell_index_to_basis_RHS_QPs_m_table_;		
	std::vector<Matrix> cell_index_to_basis_jump_QPs_m_table_;
	std::vector<Matrix> cell_index_to_basis_vertices_m_table_;
	std::vector<Matrix> cell_index_to_ddx_basis_QPs_m_table_;
	std::vector<Matrix> cell_index_to_ddy_basis_QPs_m_table_;	
	std::vector<std::map<uint, Matrix>> cell_index_to__face_neighbor_cell_index_to_basis_my_QPs_m__table_;
	

	std::vector<Matrix> set_of_simplex_P0_projected_basis_vertices_m_;
	std::vector<Matrix> set_of_simplex_P1_projected_basis_vertices_m_;

	//optimize construction
	static constexpr ushort max_num_equation = 5;
	static constexpr ushort max_num_precalculated_points = 200;
	mutable Euclidean_Vector P0_coefficient_v_;
	mutable Euclidean_Vector GE_soluion_v_;
	mutable Euclidean_Vector solution_v_;
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