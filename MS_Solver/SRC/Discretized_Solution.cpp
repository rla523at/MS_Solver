#include "../INC/Discretized_Solution.h"

Discretized_Solution::Discretized_Solution(const Governing_Equation& governing_equation, const Grid& grid) 
{
	this->num_equations_ = governing_equation.num_equations();
	this->solution_variable_names_ = governing_equation.get_variable_names();	
	this->num_cells_ = grid.num_cells();
}

const std::vector<std::string>& Discretized_Solution::get_variable_names(void) const 
{
	return this->solution_variable_names_;
}

Discretized_Solution_FVM::Discretized_Solution_FVM(const Governing_Equation& governing_equation, const Grid& grid, const Initial_Condition& initial_condition)
	:Discretized_Solution(governing_equation, grid)
{
	this->set_initial_condition(grid, initial_condition);
}

void Discretized_Solution_FVM::set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) 
{
	this->discretized_solutions_.resize(this->num_cells_ * this->num_equations_);

	const auto cell_center_nodes = grid.cell_center_nodes();

	for (size_t i = 0; i < this->num_cells_; ++i) {
		const auto initial_solution = initial_condition.calculate_solution(cell_center_nodes[i]);

		for (ushort j = 0; j < this->num_equations_; ++j) {
			const auto index = i * this->num_equations_ + j;
			this->discretized_solutions_[index] = initial_solution[j];
		}
	}
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

	const auto start_iter = this->discretized_solutions_.begin() + start_index;

	std::vector<double> values = { start_iter, start_iter + num_variable };
	return values;
}

Discretized_Solution_HOM::Discretized_Solution_HOM(const Governing_Equation& governing_equation, const Grid& grid, const Initial_Condition& initial_condition, const ushort solution_order)
	:Discretized_Solution(governing_equation, grid),
	solution_order_(solution_order) 
{
	this->set_initial_condition(grid, initial_condition);
}

void Discretized_Solution_HOM::set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) 
{	
	this->basis_vector_functions_ = grid.cell_basis_vector_functions(this->solution_order_);

	this->set_of_num_basis_.resize(this->num_cells_);
	for (size_t i = 0; i < this->num_cells_; ++i)
		this->set_of_num_basis_[i] = this->basis_vector_functions_[i].range_dimension();

	this->coeffcieint_start_indexes_.resize(this->num_cells_);
	for (size_t i = 0; i < this->num_cells_ - 1; ++i) 		
		this->coeffcieint_start_indexes_[i + 1] = this->set_of_num_basis_[i];

	this->discretized_solutions_.resize(this->num_total_basis() * this->num_equations_);

	const auto quadrature_rules = grid.cell_quadrature_rules(this->solution_order_ * 2);

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
			basis_weight_m.change_row(q, this->calculate_basis_vector_value(i,qnodes[q]) * qweights[q]);
		}

		ms::gemm(initial_solution_qnodes_m, basis_weight_m, this->coefficient_pointer(i));
	}
}

double* Discretized_Solution_HOM::coefficient_pointer(const uint icell)
{
	return this->discretized_solutions_.data() + this->coeffcieint_start_indexes_[icell];
}

std::vector<std::vector<Euclidean_Vector>> Discretized_Solution_HOM::calculate_set_of_post_point_solutions(void) const 
{
	std::vector<std::vector<Euclidean_Vector>> set_of_post_point_solutions(this->num_cells_);

	for (uint icell = 0; icell < this->num_cells_; ++icell) 
		set_of_post_point_solutions[icell] = this->calculate_post_point_solutions(icell);

	return set_of_post_point_solutions;
}

Matrix_Wrapper Discretized_Solution_HOM::get_coefficient_matrix(const uint icell) const 
{
	return { this->num_equations_, this->set_of_num_basis_[icell], this->coefficient_pointer(icell) };
}

std::vector<Euclidean_Vector> Discretized_Solution_HOM::calculate_post_point_solutions(const uint icell) const 
{
	const auto solution_post_points_m = this->get_coefficient_matrix(icell) * this->set_of_basis_post_points_m_[icell];	//E x PP
	const auto [num_equations, num_post_points] = solution_post_points_m.size();

	std::vector<Euclidean_Vector> solution_at_post_points;
	solution_at_post_points.reserve(num_post_points);

	for (ushort i = 0; i < num_post_points; ++i)
		solution_at_post_points.push_back(solution_post_points_m.column(i));

	return solution_at_post_points;
}

Euclidean_Vector Discretized_Solution_HOM::calculate_basis_vector_value(const uint icell, const Euclidean_Vector& node) const
{
	return this->basis_vector_functions_[icell](node);
}

const double* Discretized_Solution_HOM::coefficient_pointer(const uint icell) const
{
	return this->discretized_solutions_.data() + this->coeffcieint_start_indexes_[icell];
}

size_t Discretized_Solution_HOM::num_total_basis(void) const
{
	size_t num_total_basis = 0;
	for (const auto num_basis : this->set_of_num_basis_)
		num_total_basis += num_basis;
	
	return num_total_basis;
}