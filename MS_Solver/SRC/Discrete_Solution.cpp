#include "../INC/Discrete_Solution.h"

Discrete_Solution::Discrete_Solution(const Governing_Equation& governing_equation, const Grid& grid) 
{
	this->num_equations_ = governing_equation.num_equations();
	this->solution_variable_names_ = governing_equation.get_variable_names();	
	this->num_cells_ = grid.num_cells();
}

void Discrete_Solution::update_solution(Euclidean_Vector&& updated_solution)
{
	this->value_v_ = std::move(updated_solution);
}

const std::vector<std::string>& Discrete_Solution::get_variable_names(void) const 
{
	return this->solution_variable_names_;
}

const Euclidean_Vector& Discrete_Solution::get_solution_vector(void) const
{
	return this->value_v_;
}

size_t Discrete_Solution::num_values(void) const
{
	return this->value_v_.size();
}

Discretized_Solution_FVM::Discretized_Solution_FVM(const Governing_Equation& governing_equation, const Grid& grid, const Initial_Condition& initial_condition)
	:Discrete_Solution(governing_equation, grid)
{
	this->set_initial_condition(grid, initial_condition);
}

void Discretized_Solution_FVM::set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) 
{
	//Celss¿¡ ÀÖ¾î¾ß µÉ°Å°°Àºµª¼õ
	std::vector<double> initial_values(this->num_cells_ * this->num_equations_);
	const auto cell_center_nodes = grid.cell_center_nodes();

	for (size_t i = 0; i < this->num_cells_; ++i) {
		const auto initial_solution = initial_condition.calculate_solution(cell_center_nodes[i]);

		for (ushort j = 0; j < this->num_equations_; ++j) {
			const auto index = i * this->num_equations_ + j;
			initial_values[index] = initial_solution[j];
		}
	}

	this->value_v_ = std::move(initial_values);
}

std::vector<std::vector<Euclidean_Vector>> Discretized_Solution_FVM::calculate_set_of_post_point_solutions(void) const 
{
	std::vector<std::vector<Euclidean_Vector>> set_of_post_point_solutions(this->num_cells_);

	for (uint icell = 0; icell < this->num_cells_; ++icell) {
		auto solution = this->calculate_solution_at_center(icell);
		set_of_post_point_solutions[icell].push_back(std::move(solution));
	}

	return set_of_post_point_solutions;
}

Euclidean_Vector Discretized_Solution_FVM::calculate_solution_at_center(const uint icell) const 
{
	const auto start_index = icell * this->num_equations_;
	const auto num_variable = this->num_equations_;

	const auto start_iter = this->value_v_.begin() + start_index;

	std::vector<double> values = { start_iter, start_iter + num_variable };
	return values;
}

Discrete_Solution_HOM::Discrete_Solution_HOM(const Configuration& configuration, const Governing_Equation& governing_equation, const Grid& grid)
	:Discrete_Solution(governing_equation, grid)
{
	const auto solution_degree = configuration.get<ushort>("solution_degree");
	this->solution_degrees_.resize(this->num_cells_, solution_degree);

	const auto initial_condition = Initial_Condition_Factory::make(configuration);
	this->set_initial_condition(grid, *initial_condition);
}

void Discrete_Solution_HOM::set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) 
{	
	this->basis_vector_functions_ = grid.cell_basis_vector_functions(this->solution_degrees_);

	this->set_of_num_basis_.resize(this->num_cells_);
	for (size_t i = 0; i < this->num_cells_; ++i)
		this->set_of_num_basis_[i] = this->basis_vector_functions_[i].range_dimension();

	this->coefficieint_start_indexes_.resize(this->num_cells_);
	for (size_t i = 0; i < this->num_cells_ - 1; ++i) 		
		this->coefficieint_start_indexes_[i + 1] = this->num_equations_ * this->set_of_num_basis_[i];

	this->value_v_ = this->calculate_initial_values(grid, initial_condition);
}

std::vector<double> Discrete_Solution_HOM::calculate_P0_basis_values(void) const
{
	std::vector<double> P0_basis_values(this->num_cells_);

	for (int i = 0; i < this->num_cells_; ++i)
		P0_basis_values[i] = this->calculate_P0_basis_value(i);

	return P0_basis_values;
}

std::vector<Euclidean_Vector> Discrete_Solution_HOM::calculate_P0_solutions(const std::vector<double>& P0_basis_values) const
{
	std::vector<Euclidean_Vector> P0_solutions(this->num_cells_);

	for (size_t i = 0; i < this->num_cells_; ++i)
	{
		const auto P0_coefficient_v = this->P0_coefficient_vector(i);
		const auto P0_basis = P0_basis_values[i];
		P0_solutions[i] = P0_coefficient_v * P0_basis;
	}

	return P0_solutions;
}

Matrix Discrete_Solution_HOM::caclulate_basis_points_m(const uint cell_index, const std::vector<Euclidean_Vector>& points) const
{
	const auto& basis_functions = this->basis_vector_functions_[cell_index];

	const auto num_basis = this->set_of_num_basis_[cell_index];
	const auto num_points = points.size();
	Matrix basis_points_m(num_basis, num_points);

	for (ushort j = 0; j < num_points; ++j)
		basis_points_m.change_column(j, basis_functions(points[j]));

	return basis_points_m;
}

std::vector<Euclidean_Vector> Discrete_Solution_HOM::calculate_solution_at_points(const uint icell, const Matrix& basis_points_m) const
{
	const auto num_points = basis_points_m.num_column();
	const auto solution_points_m = this->coefficient_matrix(icell) * basis_points_m;

	std::vector<Euclidean_Vector> solution_at_points(num_points);
	
	for (int i = 0; i < num_points; ++i)
	{
		solution_at_points[i] = solution_points_m.column(i);
	}

	return solution_at_points;
}


size_t Discrete_Solution_HOM::coefficient_start_index(const uint icell) const
{
	return this->coefficieint_start_indexes_[icell];
}

const std::vector<ushort>& Discrete_Solution_HOM::get_solution_degrees(void) const
{
	return this->solution_degrees_;
}


Matrix_Wrapper Discrete_Solution_HOM::coefficient_matrix(const uint icell) const 
{
	return { this->num_equations_, this->set_of_num_basis_[icell], this->coefficient_pointer(icell) };
}

Euclidean_Vector Discrete_Solution_HOM::calculate_basis_vector_value(const uint icell, const Euclidean_Vector& node) const
{
	return this->basis_vector_functions_[icell](node);
}

const double* Discrete_Solution_HOM::coefficient_pointer(const uint icell) const
{
	return this->value_v_.data() + this->coefficieint_start_indexes_[icell];
}

size_t Discrete_Solution_HOM::num_total_basis(void) const
{
	size_t num_total_basis = 0;
	for (const auto num_basis : this->set_of_num_basis_)
		num_total_basis += num_basis;
	
	return num_total_basis;
}

double Discrete_Solution_HOM::calculate_P0_basis_value(const uint icell) const
{
	const auto& basis_vector_function = this->basis_vector_functions_[icell];
	const auto& P0_basis_function = basis_vector_function[0];

	return P0_basis_function.to_constant();
}

Euclidean_Vector Discrete_Solution_HOM::P0_coefficient_vector(const uint icell) const
{
	const auto coefficient_matrix = this->coefficient_matrix(icell);

	constexpr auto P0_column_index = 0;
	return coefficient_matrix.column(P0_column_index);
}

Euclidean_Vector Discrete_Solution_HOM::calculate_initial_values(const Grid& grid, const Initial_Condition& initial_condition) const
{
	std::vector<double> initial_values(this->num_total_basis() * this->num_equations_);

	std::vector<ushort> cell_integral_degree(this->num_cells_);
	for (int i = 0; i < this->num_cells_; ++i)
		cell_integral_degree[i] = this->solution_degrees_[i] * 2;

	const auto quadrature_rules = grid.cell_quadrature_rules(cell_integral_degree);

	for (uint i = 0; i < this->num_cells_; ++i)
	{
		const auto& qnodes = quadrature_rules[i].nodes;
		const auto& qweights = quadrature_rules[i].weights;

		const auto num_qnode = qnodes.size();

		Matrix initial_solution_qnodes_m(this->num_equations_, num_qnode);
		Matrix basis_weight_m(num_qnode, this->set_of_num_basis_[i]);

		for (ushort q = 0; q < num_qnode; ++q)
		{
			initial_solution_qnodes_m.change_column(q, initial_condition.calculate_solution(qnodes[q]));
			basis_weight_m.change_row(q, this->calculate_basis_vector_value(i, qnodes[q]) * qweights[q]);
		}

		auto ptr = initial_values.data() + this->coefficieint_start_indexes_[i];
		ms::gemm(initial_solution_qnodes_m, basis_weight_m, ptr);
	}

	return initial_values;
}
