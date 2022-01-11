#include "../INC/Discrete_Solution.h"

Discrete_Solution::Discrete_Solution(const std::shared_ptr<Governing_Equation>& governing_equation, const Grid& grid)
	:governing_equation_(governing_equation)
{
	this->num_cells_ = grid.num_cells();
	this->space_dimension_ = this->governing_equation_->space_dimension();
	this->num_equations_ = this->governing_equation_->num_equations();
}

Euclidean_Vector Discrete_Solution::discrete_solution_vector(void) const
{
	return this->values_;
}

Constant_Euclidean_Vector_Wrapper Discrete_Solution::discrete_solution_constant_vector_wrapper(void) const
{
	return { this->values_.size(), this->values_.data() };
}

Euclidean_Vector_Wrapper Discrete_Solution::discrete_solution_vector_wrapper(void)
{
	return { this->values_.size(), this->values_.data() };
}


const std::vector<std::string>& Discrete_Solution::get_solution_names(void) const
{
	return this->governing_equation_->get_solution_names();
}

ushort Discrete_Solution::num_equations(void) const
{
	return this->num_equations_;
}

ushort Discrete_Solution::num_solutions(void) const
{
	return this->governing_equation_->num_solutions();
}

size_t Discrete_Solution::num_total_values(void) const
{
	return this->values_.size();
}


//Discretized_Solution_FVM::Discretized_Solution_FVM(const Governing_Equation& governing_equation, const Grid& grid, const Initial_Condition& initial_condition)
//	:Discrete_Solution(grid, governing_equation)
//{
//	this->set_initial_condition(grid, initial_condition);
//}
//
//void Discretized_Solution_FVM::set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) 
//{
//	//Celss�� �־�� �ɰŰ�������
//	std::vector<double> initial_values(this->num_cells_ * this->num_equations_);
//	const auto cell_center_nodes = grid.cell_center_nodes();
//
//	for (size_t i = 0; i < this->num_cells_; ++i) {
//		const auto initial_solution = initial_condition.calculate_solution(cell_center_nodes[i]);
//
//		for (ushort j = 0; j < this->num_equations_; ++j) {
//			const auto index = i * this->num_equations_ + j;
//			initial_values[index] = initial_solution[j];
//		}
//	}
//
//	this->values_ = std::move(initial_values);
//}
//
//std::vector<std::vector<Euclidean_Vector>> Discretized_Solution_FVM::calculate_set_of_post_point_solutions(void) const 
//{
//	std::vector<std::vector<Euclidean_Vector>> set_of_post_point_solutions(this->num_cells_);
//
//	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index) {
//		auto solution = this->calculate_solution_at_center(cell_index);
//		set_of_post_point_solutions[cell_index].push_back(std::move(solution));
//	}
//
//	return set_of_post_point_solutions;
//}
//
//Euclidean_Vector Discretized_Solution_FVM::calculate_solution_at_center(const uint cell_index) const 
//{
//	const auto start_index = cell_index * this->num_equations_;
//	const auto num_variable = this->num_equations_;
//
//	const auto start_iter = this->values_.begin() + start_index;
//
//	std::vector<double> values = { start_iter, start_iter + num_variable };
//	return values;
//}
//

Discrete_Solution_DG::Discrete_Solution_DG(const std::shared_ptr<Governing_Equation>& governing_equation, const Grid& grid, const Initial_Condition& initial_condition, const ushort solution_degree)
	: Discrete_Solution(governing_equation, grid)
	, P0_coefficient_v_(this->governing_equation_->num_equations())
	, GE_soluion_v_(this->governing_equation_->num_equations())
	, solution_v_(this->governing_equation_->num_solutions())
{
	Profiler::set_time_point();

	for (ushort i = 0; i < this->max_solution_degree; ++i)
	{
		this->degree_to_num_basis_table[i] = ms::combination_with_repetition(1 + this->space_dimension_, i);
	}

	this->cell_index_to_solution_degree_table_.resize(this->num_cells_, solution_degree);
	
	this->basis_vector_functions_.resize(this->num_cells_);
	this->set_of_num_basis_.resize(this->num_cells_);
	this->set_of_num_values_.resize(this->num_cells_);
	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		this->basis_vector_functions_[cell_index] = grid.cell_basis_vector_function(cell_index, this->cell_index_to_solution_degree_table_[cell_index]);
		this->set_of_num_basis_[cell_index] = static_cast<ushort>(this->basis_vector_functions_[cell_index].size());		
		this->set_of_num_values_[cell_index] = this->num_equations_ * this->set_of_num_basis_[cell_index];
	}

	this->coefficieint_start_indexes_.resize(this->num_cells_);
	for (size_t i = 0; i < this->num_cells_ - 1; ++i)
	{
		this->coefficieint_start_indexes_[i + 1] = this->coefficieint_start_indexes_[i] + this->set_of_num_values_[i];
	}

	this->values_ = this->calculate_initial_values(grid, initial_condition);	

	LOG << std::left << std::setw(50) << "@ Discrete Solution DG precalculation" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << LOG.print_;
}

void Discrete_Solution_DG::precalculate_cell_post_element_center_points_basis_values(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_element_center_points)
{
	REQUIRE(set_of_post_element_center_points.size() == this->num_cells_, "size of set should be same with number of cells");

	this->set_of_basis_post_element_center_points_m_.reserve(this->num_cells_);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto& points = set_of_post_element_center_points[cell_index];
		this->set_of_basis_post_element_center_points_m_.push_back(this->calculate_basis_points_m(cell_index, points));
	}
}

void Discrete_Solution_DG::precalculate_cell_post_point_basis_values(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_points)
{
	REQUIRE(set_of_post_points.size() == this->num_cells_, "size of set should be same with number of cells");

	this->cell_index_to_basis_post_points_m_table_.reserve(this->num_cells_);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto& points = set_of_post_points[cell_index];
		this->cell_index_to_basis_post_points_m_table_.push_back(this->calculate_basis_points_m(cell_index, points));
	}
}

void Discrete_Solution_DG::precalculate_bdry_RHS_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<Quadrature_Rule>& quadrature_rules)
{
	const auto num_bdrys = oc_indexes.size();
	this->boundary_index_to_basis_RHS_QPs_m_.reserve(num_bdrys);

	for (uint bdry_index = 0; bdry_index < num_bdrys; ++bdry_index)
	{
		const auto& QPs = quadrature_rules[bdry_index].points;
		this->boundary_index_to_basis_RHS_QPs_m_.push_back(this->calculate_basis_points_m(oc_indexes[bdry_index], QPs));
	}
}

void Discrete_Solution_DG::precalculate_cell_RHS_QPs_basis_values(const std::vector<Quadrature_Rule>& quadrature_rules)
{
	this->cell_index_to_basis_RHS_QPs_m_table_.resize(this->num_cells_);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto& QPs = quadrature_rules[cell_index].points;
		this->cell_index_to_basis_RHS_QPs_m_table_[cell_index] = this->calculate_basis_points_m(cell_index, QPs);
	}
}

void Discrete_Solution_DG::precalculate_cell_vertices_basis_values(const std::vector<std::vector<Euclidean_Vector>>& set_of_verticies)
{
	this->cell_index_to_basis_vertices_m_table_.resize(this->num_cells_);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto& vertices = set_of_verticies[cell_index];
		this->cell_index_to_basis_vertices_m_table_[cell_index] = this->calculate_basis_points_m(cell_index, vertices);
	}
}

void Discrete_Solution_DG::precalculate_cell_P0_basis_values(void)
{
	if (this->cell_index_to_P0_basis_value_.empty())
	{
		this->cell_index_to_P0_basis_value_.resize(this->num_cells_);

		for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
		{
			const auto& basis_vector_function = this->basis_vector_functions_[cell_index];
			const auto& P0_basis_function = basis_vector_function[0];

			this->cell_index_to_P0_basis_value_[cell_index] = P0_basis_function.to_constant();
		}
	}
}

void Discrete_Solution_DG::precalculate_infs_ocs_RHS_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ocs_QPs)
{
	const auto num_infcs = oc_indexes.size();
	this->inner_face_index_to_basis_ocs_RHS_QPs_m_.reserve(num_infcs);

	for (uint infc_index = 0; infc_index < num_infcs; ++infc_index)
	{
		const auto oc_index = oc_indexes[infc_index];
		const auto& QPs = set_of_ocs_QPs[infc_index];
		this->inner_face_index_to_basis_ocs_RHS_QPs_m_.push_back(this->calculate_basis_points_m(oc_index, QPs));
	}
}

void Discrete_Solution_DG::precalculate_infs_ncs_RHS_QPs_basis_values(const std::vector<uint>& nc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ncs_QPs)
{
	const auto num_infcs = nc_indexes.size();
	this->inner_face_index_to_basis_ncs_RHS_QPs_m_.reserve(num_infcs);

	for (uint infc_index = 0; infc_index < num_infcs; ++infc_index)
	{
		const auto nc_index = nc_indexes[infc_index];
		const auto& QPs = set_of_ncs_QPs[infc_index];
		this->inner_face_index_to_basis_ncs_RHS_QPs_m_.push_back(this->calculate_basis_points_m(nc_index, QPs));
	}
}

void Discrete_Solution_DG::precalculate_set_of_simplex_P0_P1_projection_basis_vertices_m(const Grid& grid)
{
	this->set_of_simplex_P0_projected_basis_vertices_m_.resize(this->num_cells_);
	this->set_of_simplex_P1_projected_basis_vertices_m_.resize(this->num_cells_);

	constexpr ushort P0 = 0;
	constexpr ushort P1 = 1;

	const auto num_P0_basis = this->degree_to_num_basis_table[P0];
	const auto num_P1_basis = this->degree_to_num_basis_table[P1];

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto num_basis = this->set_of_num_basis_[cell_index];
		const auto vertices = grid.cell_vertices(cell_index);

		if (grid.cell_is_simplex(cell_index))
		{
			auto basis_vertices_m_ = this->calculate_basis_points_m(cell_index, vertices);

			basis_vertices_m_.change_rows(num_P1_basis, num_basis, 0.0);
			this->set_of_simplex_P1_projected_basis_vertices_m_[cell_index] = basis_vertices_m_;

			basis_vertices_m_.change_rows(num_P0_basis, num_P1_basis, 0.0);
			this->set_of_simplex_P0_projected_basis_vertices_m_[cell_index] = std::move(basis_vertices_m_);
		}
		else
		{
			const auto sub_simplex_geometries = grid.cell_sub_simplex_geometries(cell_index);

			const auto num_vertices = vertices.size();
			Matrix simplex_P0_projected_basis_vnodes(num_basis, num_vertices);
			Matrix simplex_P1_projected_basis_vnodes(num_basis, num_vertices);

			for (ushort j = 0; j < num_vertices; ++j)
			{
				const auto simplex_P0_projection_basis_vector_function = this->calculate_simplex_Pn_projection_basis_vector_function(cell_index, P0, sub_simplex_geometries[j]);
				const auto simplex_P1_projection_basis_vector_function = this->calculate_simplex_Pn_projection_basis_vector_function(cell_index, P1, sub_simplex_geometries[j]);

				const auto& vertex = vertices[j];

				simplex_P0_projected_basis_vnodes.change_column(j, simplex_P0_projection_basis_vector_function(vertex));
				simplex_P1_projected_basis_vnodes.change_column(j, simplex_P1_projection_basis_vector_function(vertex));
			}

			this->set_of_simplex_P0_projected_basis_vertices_m_[cell_index] = std::move(simplex_P0_projected_basis_vnodes);
			this->set_of_simplex_P1_projected_basis_vertices_m_[cell_index] = std::move(simplex_P1_projected_basis_vnodes);
		}
	}
}

void Discrete_Solution_DG::precalculate_infs_ocs_jump_QPs_basis_values(const std::vector<uint>& oc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ocs_jump_QPs)
{
	const auto num_infcs = oc_indexes.size();
	this->inner_face_index_to_basis_ocs_jump_QPs_m_.reserve(num_infcs);

	for (int infc_index = 0; infc_index < num_infcs; ++infc_index)
	{
		const auto oc_index = oc_indexes[infc_index];
		const auto& QPs = set_of_ocs_jump_QPs[infc_index];
		this->inner_face_index_to_basis_ocs_jump_QPs_m_.push_back(this->calculate_basis_points_m(oc_index, QPs));
	}
}

void Discrete_Solution_DG::precalculate_infs_ncs_jump_QPs_basis_values(const std::vector<uint>& nc_indexes, const std::vector<std::vector<Euclidean_Vector>>& set_of_ncs_jump_QPs)
{
	const auto num_infcs = nc_indexes.size();
	this->inner_face_index_to_basis_ncs_jump_QPs_m_.reserve(num_infcs);

	for (int infc_index = 0; infc_index < num_infcs; ++infc_index)
	{
		const auto nc_index = nc_indexes[infc_index];
		const auto& QPs = set_of_ncs_jump_QPs[infc_index];
		this->inner_face_index_to_basis_ncs_jump_QPs_m_.push_back(this->calculate_basis_points_m(nc_index, QPs));
	}
}

void Discrete_Solution_DG::precalculate_cell_jump_QPS_basis_values(const Grid& grid)
{
	this->cell_index_to_basis_jump_QPs_m_table_.resize(this->num_cells_);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto solution_degree = this->cell_index_to_solution_degree_table_[cell_index];
		
		const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, solution_degree);

		this->cell_index_to_basis_jump_QPs_m_table_[cell_index] = this->calculate_basis_points_m(cell_index, quadrature_rule.points);
	}
}

void Discrete_Solution_DG::precalculate_cell_jump_QPs_face_neighbor_cell_basis_values(const Grid& grid)
{
	this->cell_index_to__face_neighbor_cell_index_to_basis_my_QPs_m__table_.resize(this->num_cells_);

	//inter cell face
	const auto num_inter_cell_faces = grid.num_inter_cell_faces();
	const auto inter_cell_face_index_to_oc_nc_index_pair_table = grid.inter_cell_face_index_to_oc_nc_index_pair_table();

	for (size_t i = 0; i < num_inter_cell_faces; ++i)
	{
		const auto [oc_index, nc_index] = inter_cell_face_index_to_oc_nc_index_pair_table[i];

		const auto oc_solution_degree = this->cell_index_to_solution_degree_table_[oc_index];
		const auto nc_solution_degree = this->cell_index_to_solution_degree_table_[nc_index];
		const auto max_degree = (std::max)(oc_solution_degree, nc_solution_degree);

		const auto& oc_quadrature_rule = grid.get_cell_quadrature_rule(oc_index, max_degree);
		const auto& nc_quadrature_rule = grid.get_cell_quadrature_rule(nc_index, max_degree);

		this->cell_index_to__face_neighbor_cell_index_to_basis_my_QPs_m__table_[oc_index].emplace(nc_index, this->calculate_basis_points_m(nc_index, oc_quadrature_rule.points));
		this->cell_index_to__face_neighbor_cell_index_to_basis_my_QPs_m__table_[nc_index].emplace(oc_index, this->calculate_basis_points_m(oc_index, nc_quadrature_rule.points));
	}

	//periodic boundary face
	const auto pbdry_pair_index_to_oc_nc_index_pair_table = grid.pbdry_pair_index_to_oc_nc_index_pair_table();
	const auto pbdry_pair_index_to_ocs_to_ncs_v_table = grid.pbdry_pair_index_to_ocs_to_ncs_v_table();
	const auto num_pbdry_pair = pbdry_pair_index_to_oc_nc_index_pair_table.size();

	for (uint pbdry_pair_index = 0; pbdry_pair_index < num_pbdry_pair; ++pbdry_pair_index)
	{
		const auto [oc_index, nc_index] = pbdry_pair_index_to_oc_nc_index_pair_table[pbdry_pair_index];
		const auto& ocs_to_ncs_v = pbdry_pair_index_to_ocs_to_ncs_v_table.at(pbdry_pair_index);

		const auto oc_solution_degree = this->cell_index_to_solution_degree_table_[oc_index];
		const auto nc_solution_degree = this->cell_index_to_solution_degree_table_[nc_index];

		const auto& oc_quadrature_rule = grid.get_cell_quadrature_rule(oc_index, oc_solution_degree);
		const auto& nc_quadrature_rule = grid.get_cell_quadrature_rule(nc_index, nc_solution_degree);

		auto translated_oc_QPs = oc_quadrature_rule.points;
		for (auto& QP : translated_oc_QPs)
		{
			QP += ocs_to_ncs_v;
		}

		auto translated_nc_QPs = nc_quadrature_rule.points;
		for (auto& QP : translated_nc_QPs)
		{
			QP -= ocs_to_ncs_v;
		}

		this->cell_index_to__face_neighbor_cell_index_to_basis_my_QPs_m__table_[oc_index].emplace(nc_index, this->calculate_basis_points_m(nc_index, translated_oc_QPs));
		this->cell_index_to__face_neighbor_cell_index_to_basis_my_QPs_m__table_[nc_index].emplace(oc_index, this->calculate_basis_points_m(oc_index, translated_nc_QPs));
	}
}

void Discrete_Solution_DG::precalculate_cell_RHS_QPs_ddx_basis_values(const std::vector<std::vector<Euclidean_Vector>>& cell_index_to_QPs)
{
	this->cell_index_to_ddx_basis_QPs_m_table_.resize(this->num_cells_);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		this->cell_index_to_ddx_basis_QPs_m_table_[cell_index] = this->calculate_ddx_basis_points_m(cell_index, cell_index_to_QPs[cell_index]);
	}
}

void Discrete_Solution_DG::precalculate_cell_RHS_QPs_ddy_basis_values(const std::vector<std::vector<Euclidean_Vector>>& cell_index_to_QPs)
{
	this->cell_index_to_ddy_basis_QPs_m_table_.resize(this->num_cells_);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		this->cell_index_to_ddy_basis_QPs_m_table_[cell_index] = this->calculate_ddy_basis_points_m(cell_index, cell_index_to_QPs[cell_index]);
	}
}

void Discrete_Solution_DG::project_to_Pn_space(const uint cell_index, const ushort Pn)
{	
	REQUIRE(Pn < this->solution_degree(cell_index), "projection degree can not exceed given range");

	const auto num_projected_basis = this->degree_to_num_basis_table[Pn];
	const auto num_basis = this->num_basis(cell_index);
	constexpr auto zero = 0.0;

	auto coefficient_mw = this->coefficient_matrix_wrapper(cell_index);
	coefficient_mw.change_columns(num_projected_basis, num_basis, zero);
}

void Discrete_Solution_DG::scailing(const uint cell_index, const double scaling_factor)
{
	//LOG << "Scailing " << cell_index << " cell by factor " << scaling_factor << "\n" << Log::print_;
	LOG << "Scailing " << cell_index << " cell by factor " << scaling_factor << "\n";


	constexpr auto P0 = 0;
	const auto num_P0_basis = this->degree_to_num_basis_table[P0];

	const auto solution_degree = this->cell_index_to_solution_degree_table_[cell_index];
	const auto num_basis = this->degree_to_num_basis_table[solution_degree];

	auto coefficient_mw = this->coefficient_matrix_wrapper(cell_index);	
	coefficient_mw.scalar_multiplcation_at_columns(num_P0_basis, num_basis, scaling_factor);
}

void Discrete_Solution_DG::limit_slope(const uint cell_index, const double limiting_value)
{	
	constexpr auto P0 = 0;
	constexpr auto P1 = 1;
	const auto num_P0_basis = this->degree_to_num_basis_table[P0];
	const auto num_P1_basis = this->degree_to_num_basis_table[P1];

	auto coefficient_mw = this->coefficient_matrix_wrapper(cell_index);
	coefficient_mw.scalar_multiplcation_at_columns(num_P0_basis, num_P1_basis, limiting_value);
}

std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_solution_at_post_element_centers(const uint cell_index) const 
{
	REQUIRE(!this->set_of_basis_post_element_center_points_m_.empty(), "basis value should be precalculated");
	return this->calculate_solution_at_precalulated_points(cell_index, this->set_of_basis_post_element_center_points_m_[cell_index]);
}

std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_solution_at_post_points(const uint cell_index) const 
{
	REQUIRE(!this->cell_index_to_basis_post_points_m_table_.empty(), "basis value should be precalculated");
	//return this->calculate_solution_at_precalulated_points(cell_index, this->set_of_basis_post_points_m_[cell_index]);

	const auto& basis_points_m = this->cell_index_to_basis_post_points_m_table_[cell_index];

	const auto num_post_points = basis_points_m.num_column();
	const auto GE_solution_points_m = this->coefficient_constant_matrix_wrapper(cell_index) * basis_points_m;

	std::vector<Euclidean_Vector> solution_at_points(num_post_points);

	for (int i = 0; i < num_post_points; ++i)
	{
		solution_at_points[i] = GE_solution_points_m.column(i);
		try
		{
			this->governing_equation_->extend_to_solution(solution_at_points[i]);
		}
		catch (...) 
		{
			//ignore non-physical value exception
		}
	}

	return solution_at_points;
}

void Discrete_Solution_DG::calculate_solution_at_post_points(Euclidean_Vector* solution_at_post_points, const uint cell_index) const
{
	REQUIRE(!this->cell_index_to_basis_post_points_m_table_.empty(), "basis value should be precalculated");

	const auto& basis_points_m = this->cell_index_to_basis_post_points_m_table_[cell_index];
	const auto GE_solution_points_mcw = this->calculate_GE_solution_points_constant_matrix_wrapper(cell_index, basis_points_m);
	const auto num_points = GE_solution_points_mcw.num_column();

	for (size_t i = 0; i < num_points; ++i)
	{
		GE_solution_points_mcw.column(i, this->GE_soluion_v_.data());
		
		try
		{
			this->governing_equation_->extend_to_solution(this->GE_soluion_v_.data(), solution_at_post_points[i].data());
		}
		catch (...)
		{
			//ignore non-physical value exception
		}
	}
}



std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_P0_solutions(void) const
{
	std::vector<Euclidean_Vector> P0_solutions(this->num_cells_);

	const auto n = static_cast<int>(this->num_equations_);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		this->calculate_P0_coefficient_v(this->GE_soluion_v_.data(), cell_index);
		const auto P0_basis = this->cell_index_to_P0_basis_value_[cell_index];				

		ms::BLAS::cx(P0_basis, n, this->GE_soluion_v_.data());

		this->governing_equation_->extend_to_solution(this->GE_soluion_v_.data(), this->solution_v_.data());
		
		P0_solutions[cell_index] = this->solution_v_;
	}

	return P0_solutions;
}

std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_solution_at_bdry_RHS_QPs(const uint bdry_index, const uint oc_index) const
{
	REQUIRE(!this->boundary_index_to_basis_RHS_QPs_m_.empty(), "boundary basis QPs should be precalculated");
	return this->calculate_solution_at_precalulated_points(oc_index, this->boundary_index_to_basis_RHS_QPs_m_[bdry_index]);
}

std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_solution_at_cell_RHS_QPs(const uint cell_index) const
{
	REQUIRE(!this->cell_index_to_basis_RHS_QPs_m_table_.empty(), "cell basis QPs should be precalculated");
	return this->calculate_solution_at_precalulated_points(cell_index, this->cell_index_to_basis_RHS_QPs_m_table_[cell_index]);
}

std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_solution_at_infc_ocs_RHS_QPs(const uint infs_index, const uint oc_index) const
{
	REQUIRE(!this->inner_face_index_to_basis_ocs_RHS_QPs_m_.empty(), "inner face basis owner cell side QPs should be precalculated");
	return this->calculate_solution_at_precalulated_points(oc_index, this->inner_face_index_to_basis_ocs_RHS_QPs_m_[infs_index]);
}

std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_solution_at_infc_ncs_RHS_QPs(const uint infs_index, const uint nc_index) const
{
	REQUIRE(!this->inner_face_index_to_basis_ncs_RHS_QPs_m_.empty(), "inner face basis neighbor cell side QPs should be precalculated");
	return this->calculate_solution_at_precalulated_points(nc_index, this->inner_face_index_to_basis_ncs_RHS_QPs_m_[infs_index]);
}

std::vector<double> Discrete_Solution_DG::calculate_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const
{
	REQUIRE(!this->cell_index_to_basis_vertices_m_table_.empty(), "cell basis vertices should be precalculated");
	return this->calculate_nth_solution_at_precalulated_points(cell_index, equation_index, this->cell_index_to_basis_vertices_m_table_[cell_index]);
}

double Discrete_Solution_DG::calculate_P0_nth_solution(const uint cell_index, const ushort solution_index) const
{
	REQUIRE(!this->cell_index_to_P0_basis_value_.empty(), "cell P0 basis values should be precalculated");	

	if (solution_index <= this->governing_equation_->num_equations())
	{
		return this->P0_nth_coefficient(cell_index, solution_index) * this->cell_index_to_P0_basis_value_[cell_index];
	}
	else if (solution_index <= this->governing_equation_->num_solutions())
	{
		const auto n = static_cast<int>(this->num_equations_);

		this->calculate_P0_coefficient_v(this->GE_soluion_v_.data(), cell_index);
		const auto P0_basis = this->cell_index_to_P0_basis_value_[cell_index];

		ms::BLAS::cx(P0_basis, n, this->GE_soluion_v_.data());

		this->governing_equation_->extend_to_solution(this->GE_soluion_v_.data(), this->solution_v_.data());

		return this->solution_v_[solution_index];
	}
	else
	{
		EXCEPTION("equation index exceed given range");
		return NULL;
	}
}

std::vector<double> Discrete_Solution_DG::calculate_simplex_P0_projected_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const
{
	REQUIRE(!this->set_of_simplex_P0_projected_basis_vertices_m_.empty(), "simplex P0 projected basis vertices should be precalculated");
	return this->calculate_nth_solution_at_precalulated_points(cell_index, equation_index, this->set_of_simplex_P0_projected_basis_vertices_m_[cell_index]);
}

std::vector<double> Discrete_Solution_DG::calculate_simplex_P1_projected_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const
{
	REQUIRE(!this->set_of_simplex_P1_projected_basis_vertices_m_.empty(), "simplex P1 projected basis vertices should be precalculated");
	return this->calculate_nth_solution_at_precalulated_points(cell_index, equation_index, this->set_of_simplex_P1_projected_basis_vertices_m_[cell_index]);
}

std::vector<double> Discrete_Solution_DG::calculate_P1_projected_nth_solution_at_vertices(const uint cell_index, const ushort equation_index) const
{
	REQUIRE(!this->cell_index_to_basis_vertices_m_table_.empty(), "cell basis vertices should be precalculated");
	constexpr auto P1 = 1;
	return this->calculate_Pn_projected_mth_solution_at_precalulated_points(cell_index, P1, equation_index, this->cell_index_to_basis_vertices_m_table_[cell_index]);
}

void Discrete_Solution_DG::calculate_solution_at_bdry_RHS_QPs(Euclidean_Vector* solution_at_QPs, const uint bdry_index, const uint oc_index)
{
	REQUIRE(!this->boundary_index_to_basis_RHS_QPs_m_.empty(), "basis value should be precalculated");
	try
	{
		this->calculate_solution_at_precalulated_points(solution_at_QPs, oc_index, this->boundary_index_to_basis_RHS_QPs_m_[bdry_index]);
	}
	catch (...)
	{
		this->scailing(oc_index, this->scaling_factor_);
		this->calculate_solution_at_bdry_RHS_QPs(solution_at_QPs, bdry_index, oc_index);
	}
}

void Discrete_Solution_DG::calculate_solution_at_cell_RHS_QPs(Euclidean_Vector* solution_at_QPs_ptr, const uint cell_index)
{
	REQUIRE(!this->cell_index_to_basis_RHS_QPs_m_table_.empty(), "basis value should be precalculated");
	try
	{
		this->calculate_solution_at_precalulated_points(solution_at_QPs_ptr, cell_index, this->cell_index_to_basis_RHS_QPs_m_table_[cell_index]);
	}
	catch (...)
	{
		this->scailing(cell_index, this->scaling_factor_);
		this->calculate_solution_at_cell_RHS_QPs(solution_at_QPs_ptr, cell_index);
	}
}

void Discrete_Solution_DG::calculate_solution_at_cell_RHS_QPs(Euclidean_Vector* solution_at_QPs_ptr, const uint cell_index) const
{
	REQUIRE(!this->cell_index_to_basis_RHS_QPs_m_table_.empty(), "basis value should be precalculated");
	try
	{
		this->calculate_solution_at_precalulated_points(solution_at_QPs_ptr, cell_index, this->cell_index_to_basis_RHS_QPs_m_table_[cell_index]);
	}
	catch (...)
	{
		//ignore negative density pressure error
	}
}

void Discrete_Solution_DG::calculate_solution_at_infc_ocs_RHS_QPs(Euclidean_Vector* solution_at_infc_ocs_QPs, const uint infs_index, const uint oc_index)
{
	REQUIRE(!this->inner_face_index_to_basis_ocs_RHS_QPs_m_.empty(), "basis value should be precalculated");
	try
	{
		this->calculate_solution_at_precalulated_points(solution_at_infc_ocs_QPs, oc_index, this->inner_face_index_to_basis_ocs_RHS_QPs_m_[infs_index]);
	}
	catch (...)
	{
		this->scailing(oc_index, this->scaling_factor_);
		this->calculate_solution_at_infc_ocs_RHS_QPs(solution_at_infc_ocs_QPs, infs_index, oc_index);
	}
}

void Discrete_Solution_DG::calculate_solution_at_infc_ncs_RHS_QPs(Euclidean_Vector* solution_at_infc_ncs_QPs, const uint infs_index, const uint nc_index)
{
	REQUIRE(!this->inner_face_index_to_basis_ncs_RHS_QPs_m_.empty(), "basis value should be precalculated");
	try
	{
		this->calculate_solution_at_precalulated_points(solution_at_infc_ncs_QPs, nc_index, this->inner_face_index_to_basis_ncs_RHS_QPs_m_[infs_index]);
	}
	catch (...)
	{
		this->scailing(nc_index, this->scaling_factor_);
		this->calculate_solution_at_infc_ncs_RHS_QPs(solution_at_infc_ncs_QPs, infs_index, nc_index);
	}
}

void Discrete_Solution_DG::calculate_simplex_P0_projected_nth_solution_at_vertices(double* simplex_P0_projected_nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const
{
	REQUIRE(!this->set_of_simplex_P0_projected_basis_vertices_m_.empty(), "simplex P0 projected basis vertices should be precalculated");
	this->calculate_nth_solution_at_precalulated_points(simplex_P0_projected_nth_solution_at_vertices, cell_index, equation_index, this->set_of_simplex_P0_projected_basis_vertices_m_[cell_index]);
}

void Discrete_Solution_DG::calculate_simplex_P1_projected_nth_solution_at_vertices(double* simplex_P1_projected_nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const
{
	REQUIRE(!this->set_of_simplex_P1_projected_basis_vertices_m_.empty(), "simplex P1 projected basis vertices should be precalculated");
	this->calculate_nth_solution_at_precalulated_points(simplex_P1_projected_nth_solution_at_vertices, cell_index, equation_index, this->set_of_simplex_P1_projected_basis_vertices_m_[cell_index]);
}

void Discrete_Solution_DG::calculate_nth_solution_at_vertices(double* nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const
{
	REQUIRE(!this->cell_index_to_basis_vertices_m_table_.empty(), "cell basis QPs should be precalculated");
	this->calculate_nth_solution_at_precalulated_points(nth_solution_at_vertices, cell_index, equation_index, this->cell_index_to_basis_vertices_m_table_[cell_index]);
}

void Discrete_Solution_DG::calculate_P1_projected_nth_solution_at_vertices(double* P1_projected_nth_solution_at_vertices, const uint cell_index, const ushort equation_index) const
{
	REQUIRE(!this->cell_index_to_basis_vertices_m_table_.empty(), "cell basis vertices should be precalculated");
	constexpr auto P1 = 1;
	this->calculate_Pn_projected_mth_solution_at_precalulated_points(P1_projected_nth_solution_at_vertices, cell_index, P1, equation_index, this->cell_index_to_basis_vertices_m_table_[cell_index]);
}

void Discrete_Solution_DG::calculate_nth_solution_at_infc_ocs_jump_QPs(double* nth_solution_at_infc_ocs_jump_QPs, const uint infc_index, const uint oc_index, const ushort equation_index) const
{
	REQUIRE(!this->inner_face_index_to_basis_ocs_jump_QPs_m_.empty(), "basis value should be precalculated");
	this->calculate_nth_solution_at_precalulated_points(nth_solution_at_infc_ocs_jump_QPs, oc_index, equation_index, this->inner_face_index_to_basis_ocs_jump_QPs_m_[infc_index]);
}

void Discrete_Solution_DG::calculate_nth_solution_at_infc_ncs_jump_QPs(double* nth_solution_at_infc_ncs_jump_QPs, const uint infc_index, const uint nc_index, const ushort equation_index) const
{
	REQUIRE(!this->inner_face_index_to_basis_ncs_jump_QPs_m_.empty(), "basis value should be precalculated");
	this->calculate_nth_solution_at_precalulated_points(nth_solution_at_infc_ncs_jump_QPs, nc_index, equation_index, this->inner_face_index_to_basis_ncs_jump_QPs_m_[infc_index]);
}

void Discrete_Solution_DG::calculate_nth_solution_at_cell_jump_QPs(double* nth_solution_at_cell_jump_QPs, const uint cell_index, const ushort solution_index) const
{
	REQUIRE(!this->cell_index_to_basis_jump_QPs_m_table_.empty(), "basis value should be precalculated");
	const auto basis_jump_QPs_m = this->cell_index_to_basis_jump_QPs_m_table_[cell_index];
	this->calculate_nth_solution_at_precalulated_points(nth_solution_at_cell_jump_QPs, cell_index, solution_index, basis_jump_QPs_m);
}

void Discrete_Solution_DG::calculate_nth_extrapolated_solution_at_cell_jump_QPs(double* nth_solution_at_target_cell_QPs, const uint face_share_cell_index, const uint target_cell_index, const ushort solution_index) const
{
	REQUIRE(!this->cell_index_to__face_neighbor_cell_index_to_basis_my_QPs_m__table_.empty(), "basis value should be precalculated");
	const auto& basis_points_m = this->cell_index_to__face_neighbor_cell_index_to_basis_my_QPs_m__table_[target_cell_index].at(face_share_cell_index);
	this->calculate_nth_solution_at_precalulated_points(nth_solution_at_target_cell_QPs, face_share_cell_index, solution_index, basis_points_m);
}

void Discrete_Solution_DG::calculate_ddx_GE_solution_at_cell_RHS_QPs(Euclidean_Vector* solution_at_QPs, const uint cell_index) const
{
	REQUIRE(!this->cell_index_to_ddx_basis_QPs_m_table_.empty(), "basis value shold be precalculated");
	const auto& basis_points_m = this->cell_index_to_ddx_basis_QPs_m_table_[cell_index];
	this->calculate_GE_solution_at_precalulated_points(solution_at_QPs, cell_index, basis_points_m);
}

void Discrete_Solution_DG::calculate_ddy_GE_solution_at_cell_RHS_QPs(Euclidean_Vector* solution_at_QPs, const uint cell_index) const
{
	REQUIRE(!this->cell_index_to_ddy_basis_QPs_m_table_.empty(), "basis value shold be precalculated");
	const auto& basis_points_m = this->cell_index_to_ddy_basis_QPs_m_table_[cell_index];
	this->calculate_GE_solution_at_precalulated_points(solution_at_QPs, cell_index, basis_points_m);
}

Matrix Discrete_Solution_DG::calculate_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const
{
	const auto& basis_functions = this->basis_vector_functions_[cell_index];

	const auto num_basis = this->set_of_num_basis_[cell_index];
	const auto num_points = points.size();
	Matrix basis_points_m(num_basis, num_points);

	for (ushort j = 0; j < num_points; ++j)
	{
		basis_points_m.change_column(j, basis_functions(points[j]));
	}

	return basis_points_m;
}

Matrix Discrete_Solution_DG::calculate_ddx_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const
{
	const auto& basis_vf = this->basis_vector_functions_[cell_index];
	
	constexpr auto x = 0;
	const auto dx_basis_vf = basis_vf.get_differentiate(x);

	const auto num_basis = this->set_of_num_basis_[cell_index];
	const auto num_points = points.size();
	Matrix ddx_basis_points_m(num_basis, num_points);

	for (ushort j = 0; j < num_points; ++j)
	{
		ddx_basis_points_m.change_column(j, dx_basis_vf(points[j]));
	}

	return ddx_basis_points_m;
}

Matrix Discrete_Solution_DG::calculate_ddy_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const
{
	const auto& basis_vf = this->basis_vector_functions_[cell_index];

	constexpr auto y = 1;
	const auto dy_basis_vf = basis_vf.get_differentiate(y);

	const auto num_basis = this->set_of_num_basis_[cell_index];
	const auto num_points = points.size();
	Matrix ddy_basis_points_m(num_basis, num_points);

	for (ushort j = 0; j < num_points; ++j)
	{
		ddy_basis_points_m.change_column(j, dy_basis_vf(points[j]));
	}

	return ddy_basis_points_m;
}

Vector_Function<Polynomial> Discrete_Solution_DG::calculate_simplex_Pn_projection_basis_vector_function(const uint cell_index, const ushort Pn, const Geometry& sub_simplex_geometry) const
{
	const auto& basis_vector_function = this->basis_vector_functions_[cell_index];
	const auto Pn_simplex_basis_vector_function = sub_simplex_geometry.orthonormal_basis_vector_function(Pn);

	const auto num_basis = this->set_of_num_basis_[cell_index];

	std::vector<Polynomial> simplex_Pn_projection_basis_functions(num_basis);

	for (ushort i = 0; i < num_basis; ++i)
	{
		for (ushort j = 0; j < this->degree_to_num_basis_table[Pn]; ++j)
		{
			simplex_Pn_projection_basis_functions[i] += ms::inner_product(basis_vector_function[i], Pn_simplex_basis_vector_function[j], sub_simplex_geometry) * Pn_simplex_basis_vector_function[j];
		}
	}

	return simplex_Pn_projection_basis_functions;
}

Euclidean_Vector Discrete_Solution_DG::calculate_basis_point_v(const uint cell_index, const Euclidean_Vector& node) const
{
	return this->basis_vector_functions_[cell_index](node);
}

Matrix_Function<Polynomial> Discrete_Solution_DG::calculate_tranposed_gradient_basis(const uint cell_index) const
{
	const auto num_basis = this->set_of_num_basis_[cell_index];
	const auto& basis_function = this->basis_vector_functions_[cell_index];

	Matrix_Function<Polynomial> transposed_gradient_basis(this->space_dimension_, num_basis);
	for (ushort i = 0; i < num_basis; ++i)
	{
		transposed_gradient_basis.change_column(i, basis_function[i].gradient(this->space_dimension_));
	}

	return transposed_gradient_basis;
}

Constant_Matrix_Wrapper Discrete_Solution_DG::calculate_GE_solution_points_constant_matrix_wrapper(const uint cell_index, const Constant_Matrix_Wrapper& basis_points_m) const
{
	const auto num_points = basis_points_m.num_column();
	std::fill(this->solution_at_points_values_.begin(), this->solution_at_points_values_.begin() + this->num_equations_ * num_points, 0.0);

	ms::gemm(this->coefficient_constant_matrix_wrapper(cell_index), basis_points_m, this->solution_at_points_values_.data());
	return { this->num_equations_, num_points, this->solution_at_points_values_.data() };
}

std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_solution_at_precalulated_points(const uint cell_index, const Constant_Matrix_Wrapper& basis_points_m) const
{
	const auto num_points = basis_points_m.num_column();
	const auto GE_solution_points_m = this->coefficient_constant_matrix_wrapper(cell_index) * basis_points_m;

	std::vector<Euclidean_Vector> solution_at_points(num_points);

	for (int i = 0; i < num_points; ++i)
	{
		solution_at_points[i] = GE_solution_points_m.column(i);
		this->governing_equation_->extend_to_solution(solution_at_points[i]);
	}

	return solution_at_points;
}

std::vector<double> Discrete_Solution_DG::calculate_nth_solution_at_precalulated_points(const uint cell_index, const ushort equation_index, const Constant_Matrix_Wrapper& basis_points_m) const
{
	const auto num_points = basis_points_m.num_column();
	
	std::vector<double> nth_solution_at_points(num_points);
	ms::gemm(this->nth_coefficient_matrix_contant_wrapper(cell_index, equation_index), basis_points_m, nth_solution_at_points.data());
	
	return nth_solution_at_points;
}

std::vector<double> Discrete_Solution_DG::calculate_Pn_projected_mth_solution_at_precalulated_points(const uint cell_index, const ushort Pn, const ushort equation_index, const Constant_Matrix_Wrapper& basis_points_m) const
{
	const auto num_points = basis_points_m.num_column();
	std::vector<double> nth_Pn_projected_solution_at_points(num_points);

	const auto Pn_projected_nth_coefficient_mcw = this->Pn_projected_mth_coefficient_contant_matrix_wrapper(cell_index, Pn, equation_index);
	const auto Pn_projected_basis_mcw = this->Pn_projected_basis_constant_matrix_wrapper(basis_points_m, Pn);
	ms::gemm(Pn_projected_nth_coefficient_mcw, Pn_projected_basis_mcw, nth_Pn_projected_solution_at_points.data());

	return nth_Pn_projected_solution_at_points;	
}

void Discrete_Solution_DG::calculate_GE_solution_at_precalulated_points(Euclidean_Vector* solution_v_at_points_ptr, const uint cell_index, const Constant_Matrix_Wrapper& basis_points_m) const
{
	const auto GE_solution_points_mcw = this->calculate_GE_solution_points_constant_matrix_wrapper(cell_index, basis_points_m);
	const auto num_points = GE_solution_points_mcw.num_column();

	for (size_t i = 0; i < num_points; ++i)
	{
		GE_solution_points_mcw.column(i, solution_v_at_points_ptr[i].data());
	}
}

void Discrete_Solution_DG::calculate_solution_at_precalulated_points(Euclidean_Vector* solution_v_at_points_ptr, const uint cell_index, const Constant_Matrix_Wrapper& basis_points_m) const
{
	const auto GE_solution_points_mcw = this->calculate_GE_solution_points_constant_matrix_wrapper(cell_index, basis_points_m);
	const auto num_points = GE_solution_points_mcw.num_column();

	for (size_t i = 0; i < num_points; ++i)
	{
		GE_solution_points_mcw.column(i, this->GE_soluion_v_.data());
		this->governing_equation_->extend_to_solution(GE_soluion_v_.data(), solution_v_at_points_ptr[i].data());
	}
}

void Discrete_Solution_DG::calculate_nth_solution_at_precalulated_points(double* nth_solution_at_points, const uint cell_index, const ushort solution_index, const Constant_Matrix_Wrapper& basis_points_m) const
{
	if (solution_index < this->governing_equation_->num_equations())
	{
		const auto num_points = basis_points_m.num_column();
		std::fill(nth_solution_at_points, nth_solution_at_points + num_points, 0.0);

		ms::gemm(this->nth_coefficient_matrix_contant_wrapper(cell_index, solution_index), basis_points_m, nth_solution_at_points);
	}
	else if (solution_index <= this->governing_equation_->num_solutions())
	{
		const auto num_points = basis_points_m.num_column();
		std::fill(this->solution_at_points_values_.begin(), this->solution_at_points_values_.begin() + this->num_equations_ * num_points, 0.0);

		ms::gemm(this->coefficient_constant_matrix_wrapper(cell_index), basis_points_m, this->solution_at_points_values_.data());
		Constant_Matrix_Wrapper GE_solution_points_cmw(this->num_equations_, num_points, this->solution_at_points_values_.data());

		for (int i = 0; i < num_points; ++i)
		{
			GE_solution_points_cmw.column(i, this->GE_soluion_v_.data());			
			this->governing_equation_->extend_to_solution(this->GE_soluion_v_.data(), this->solution_v_.data());
			nth_solution_at_points[i] = this->solution_v_[solution_index];
		}
	}
	else
	{
		EXCEPTION("equation index exceed given range");
	}
}

void Discrete_Solution_DG::calculate_Pn_projected_mth_solution_at_precalulated_points(double* Pn_projected_mth_solution_at_points, const uint cell_index, const ushort Pn, const ushort equation_index, const Constant_Matrix_Wrapper& basis_points_m) const
{
	const auto num_points = basis_points_m.num_column();
	std::fill(Pn_projected_mth_solution_at_points, Pn_projected_mth_solution_at_points + num_points, 0.0);

	const auto Pn_projected_nth_coefficient_cmw = this->Pn_projected_mth_coefficient_contant_matrix_wrapper(cell_index, Pn, equation_index);
	const auto Pn_projected_basis_cmw = this->Pn_projected_basis_constant_matrix_wrapper(basis_points_m, Pn);
	ms::gemm(Pn_projected_nth_coefficient_cmw, Pn_projected_basis_cmw, Pn_projected_mth_solution_at_points);
}

ushort Discrete_Solution_DG::num_basis(const uint cell_index) const
{
	return this->set_of_num_basis_[cell_index];
}

ushort Discrete_Solution_DG::num_values(const uint cell_index) const
{
	return this->set_of_num_values_[cell_index];
}


ushort Discrete_Solution_DG::maximum_solution_degree(void) const
{
	return *std::max_element(this->cell_index_to_solution_degree_table_.begin(), this->cell_index_to_solution_degree_table_.end());
}

ushort Discrete_Solution_DG::solution_degree(const uint cell_index) const
{
	return this->cell_index_to_solution_degree_table_[cell_index];
}

const std::vector<size_t>& Discrete_Solution_DG::get_coefficient_start_indexes(void) const
{
	return this->coefficieint_start_indexes_;
}

double* Discrete_Solution_DG::coefficient_pointer(const uint cell_index) 
{
	return this->values_.data() + this->coefficieint_start_indexes_[cell_index];
}

Matrix_Wrapper Discrete_Solution_DG::coefficient_matrix_wrapper(const uint cell_index)
{
	return { this->num_equations_, this->set_of_num_basis_[cell_index], this->coefficient_pointer(cell_index) };
}

const double* Discrete_Solution_DG::coefficient_pointer(const uint cell_index) const
{
	return this->values_.data() + this->coefficieint_start_indexes_[cell_index];
}

double Discrete_Solution_DG::P0_nth_coefficient(const uint cell_index, const ushort equation_index) const
{
	const auto num_basis = this->set_of_num_basis_[cell_index];
	const auto jump_index = equation_index * num_basis;

	return *(this->coefficient_pointer(cell_index) + jump_index);
}

//Euclidean_Vector Discrete_Solution_DG::P0_coefficient_v(const uint cell_index) const
//{
//	const auto coefficient_mcw = this->coefficient_constant_matrix_wrapper(cell_index);
//
//	constexpr auto P0_column_index = 0;
//	return coefficient_mcw.column(P0_column_index);
//}

void Discrete_Solution_DG::calculate_P0_coefficient_v(double* ptr, const uint cell_index) const
{
	const auto coefficient_mcw = this->coefficient_constant_matrix_wrapper(cell_index);

	constexpr auto P0_column_index = 0;
	coefficient_mcw.column(P0_column_index, ptr);
}

Constant_Matrix_Wrapper Discrete_Solution_DG::coefficient_constant_matrix_wrapper(const uint cell_index) const
{
	return { this->num_equations_, this->set_of_num_basis_[cell_index], this->coefficient_pointer(cell_index) };
}

Constant_Matrix_Wrapper Discrete_Solution_DG::nth_coefficient_matrix_contant_wrapper(const uint cell_index, const ushort equation_index) const
{
	const auto num_basis = this->set_of_num_basis_[cell_index];
	const auto jump_index = equation_index * num_basis;

	return { 1, num_basis , this->coefficient_pointer(cell_index) + jump_index };
}

Constant_Matrix_Wrapper Discrete_Solution_DG::Pn_projected_mth_coefficient_contant_matrix_wrapper(const uint cell_index, const ushort Pn, const ushort equation_index) const
{
	const auto num_basis = this->set_of_num_basis_[cell_index];
	const auto jump_index = equation_index * num_basis;
	const auto num_projected_basis = this->degree_to_num_basis_table[Pn];

	return { 1, num_projected_basis , this->coefficient_pointer(cell_index) + jump_index };
}


Constant_Matrix_Wrapper Discrete_Solution_DG::Pn_projected_basis_constant_matrix_wrapper(const Constant_Matrix_Wrapper& basis_points_m, const ushort Pn) const
{
	const auto num_Pn_basis = this->degree_to_num_basis_table[Pn];
	const auto num_points = basis_points_m.num_column();
	return { num_Pn_basis, num_points, basis_points_m.data() };
}

size_t Discrete_Solution_DG::num_total_basis(void) const
{
	size_t num_total_basis = 0;
	for (const auto num_basis : this->set_of_num_basis_)
	{
		num_total_basis += num_basis;
	}
	
	return num_total_basis;
}

std::vector<double> Discrete_Solution_DG::calculate_initial_values(const Grid& grid, const Initial_Condition& initial_condition) const
{
	std::vector<double> initial_values(this->num_total_basis() * this->num_equations_);

	for (uint i = 0; i < this->num_cells_; ++i)
	{
		const auto integrand_degree = this->cell_index_to_solution_degree_table_[i] * 2;

		const auto quadrature_rule = grid.get_cell_quadrature_rule(i, integrand_degree);
		const auto& qnodes = quadrature_rule.points;
		const auto& qweights = quadrature_rule.weights;

		const auto num_qnode = qnodes.size();

		Matrix initial_solution_qnodes_m(this->num_equations_, num_qnode);
		Matrix basis_weight_m(num_qnode, this->set_of_num_basis_[i]);

		for (ushort q = 0; q < num_qnode; ++q)
		{
			initial_solution_qnodes_m.change_column(q, initial_condition.calculate_solution(qnodes[q]));
			basis_weight_m.change_row(q, this->calculate_basis_point_v(i, qnodes[q]) * qweights[q]);
		}

		auto ptr = initial_values.data() + this->coefficieint_start_indexes_[i];
		ms::gemm(initial_solution_qnodes_m, basis_weight_m, ptr);
	}

	return initial_values;
}
