#include "../INC/Discrete_Solution.h"

Discrete_Solution::Discrete_Solution(const Grid& grid, const std::shared_ptr<Governing_Equation>& governing_equation)
	:governing_equation_(governing_equation)
{
	this->num_cells_ = grid.num_cells();
	this->space_dimension_ = this->governing_equation_->space_dimension();
	this->num_equations_ = this->governing_equation_->num_equations();
}

void Discrete_Solution::update_solution(Euclidean_Vector&& updated_solution)
{
	this->value_v_ = std::move(updated_solution);
}

const Euclidean_Vector& Discrete_Solution::get_solution_vector(void) const
{
	return this->value_v_;
}

size_t Discrete_Solution::num_values(void) const
{
	return this->value_v_.size();
}


//Discretized_Solution_FVM::Discretized_Solution_FVM(const Governing_Equation& governing_equation, const Grid& grid, const Initial_Condition& initial_condition)
//	:Discrete_Solution(grid, governing_equation)
//{
//	this->set_initial_condition(grid, initial_condition);
//}
//
//void Discretized_Solution_FVM::set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) 
//{
//	//Celss¿¡ ÀÖ¾î¾ß µÉ°Å°°Àºµª¼õ
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
//	this->value_v_ = std::move(initial_values);
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
//	const auto start_iter = this->value_v_.begin() + start_index;
//
//	std::vector<double> values = { start_iter, start_iter + num_variable };
//	return values;
//}
//

Discrete_Solution_DG::Discrete_Solution_DG(const Grid& grid, const std::shared_ptr<Governing_Equation>& governing_equation, const Initial_Condition& initial_condition, const ushort solution_degree)
	: Discrete_Solution(grid, governing_equation)
{
	this->solution_degrees_.resize(this->num_cells_, solution_degree);

	this->basis_vector_functions_ = grid.cell_basis_vector_functions(this->solution_degrees_);

	this->set_of_num_basis_.resize(this->num_cells_);
	for (size_t i = 0; i < this->num_cells_; ++i)
		this->set_of_num_basis_[i] = this->basis_vector_functions_[i].size();

	this->coefficieint_start_indexes_.resize(this->num_cells_);
	for (size_t i = 0; i < this->num_cells_ - 1; ++i)
		this->coefficieint_start_indexes_[i + 1] = this->num_equations_ * this->set_of_num_basis_[i];

	this->value_v_ = this->calculate_initial_values(grid, initial_condition);
}

double Discrete_Solution_DG::calculate_P0_basis_value(const uint cell_index) const
{
	const auto& basis_vector_function = this->basis_vector_functions_[cell_index];
	const auto& P0_basis_function = basis_vector_function[0];

	return P0_basis_function.to_constant();
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
		transposed_gradient_basis.change_column(i, basis_function[i].gradient());
	}

	return transposed_gradient_basis;
}

std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_P0_solutions(const std::vector<double>& P0_basis_values) const
{
	std::vector<Euclidean_Vector> P0_solutions(this->num_cells_);

	for (size_t i = 0; i < this->num_cells_; ++i)
	{
		const auto P0_coefficient_v = this->P0_coefficient_v(i);
		const auto P0_basis = P0_basis_values[i];
		const auto GE_solution = P0_coefficient_v * P0_basis;
		P0_solutions[i] = this->governing_equation_->calculate_solution(GE_solution);
	}

	return P0_solutions;
}

std::vector<Euclidean_Vector> Discrete_Solution_DG::calculate_solution_at_points(const uint cell_index, const Matrix& basis_points_m) const
{
	const auto num_points = basis_points_m.num_column();
	const auto GE_solution_points_m = this->coefficient_m(cell_index) * basis_points_m;

	std::vector<Euclidean_Vector> solution_at_points(num_points);
	
	for (int i = 0; i < num_points; ++i)
	{
		const auto GE_solution = GE_solution_points_m.column(i);
		solution_at_points[i] = this->governing_equation_->calculate_solution(GE_solution);
	}

	return solution_at_points;
}


size_t Discrete_Solution_DG::coefficient_start_index(const uint cell_index) const
{
	return this->coefficieint_start_indexes_[cell_index];
}

ushort Discrete_Solution_DG::num_basis(const uint cell_index) const
{
	return this->set_of_num_basis_[cell_index];
}

ushort Discrete_Solution_DG::num_equations(void) const
{
	return this->num_equations_;
}

ushort Discrete_Solution_DG::solution_degree(const uint cell_index) const
{
	return this->solution_degrees_[cell_index];
}

const std::vector<size_t>& Discrete_Solution_DG::get_coefficient_start_indexes(void) const
{
	return this->coefficieint_start_indexes_;
}

const std::vector<ushort>& Discrete_Solution_DG::get_solution_degrees(void) const
{
	return this->solution_degrees_;
}

const std::vector<ushort>& Discrete_Solution_DG::get_set_of_num_basis(void) const
{
	return this->set_of_num_basis_;
}

Euclidean_Vector Discrete_Solution_DG::P0_coefficient_v(const uint cell_index) const
{
	const auto coefficient_matrix = this->coefficient_m(cell_index);

	constexpr auto P0_column_index = 0;
	return coefficient_matrix.column(P0_column_index);
}

const double* Discrete_Solution_DG::coefficient_pointer(const uint cell_index) const
{
	return this->value_v_.data() + this->coefficieint_start_indexes_[cell_index];
}

Matrix_Wrapper Discrete_Solution_DG::coefficient_m(const uint cell_index) const
{
	return { this->num_equations_, this->set_of_num_basis_[cell_index], this->coefficient_pointer(cell_index) };
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

Euclidean_Vector Discrete_Solution_DG::calculate_initial_values(const Grid& grid, const Initial_Condition& initial_condition) const
{
	std::vector<double> initial_values(this->num_total_basis() * this->num_equations_);

	for (uint i = 0; i < this->num_cells_; ++i)
	{
		const auto quadrature_rule = grid.cell_quadrature_rule(i, this->solution_degrees_[i] * 2);
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
