#include "../INC/Discretized_Solution.h"

Discretized_Solution::Discretized_Solution(const Governing_Equation& governing_equation) {
	this->num_equations_ = governing_equation.num_equations();
	this->solution_variable_names_ = governing_equation.get_variable_names();	
}

const std::vector<std::string>& Discretized_Solution::get_variable_names(void) const {
	return this->solution_variable_names_;
}

void Discretized_Solution_FVM::set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) {
	this->num_cells_ = grid.num_cells();
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

std::vector<std::vector<Euclidean_Vector>> Discretized_Solution_FVM::calculate_set_of_post_point_solutions(void) const {
	std::vector<std::vector<Euclidean_Vector>> set_of_post_point_solutions(this->num_cells_);

	for (uint icell = 0; icell < this->num_cells_; ++icell) {
		auto solution = this->calculate_solution_at_center(icell);
		set_of_post_point_solutions[icell].push_back(std::move(solution));
	}

	return set_of_post_point_solutions;
}

Euclidean_Vector Discretized_Solution_FVM::calculate_solution_at_center(const uint icell) const {
	const auto start_index = icell * this->num_equations_;
	const auto num_variable = this->num_equations_;

	const auto start_iter = this->discretized_solutions_.begin() + start_index;

	std::vector<double> values = { start_iter, start_iter + num_variable };
	return values;
}

void Discretized_Solution_HOM::set_initial_condition(const Grid& grid, const Initial_Condition& initial_condition) {
	
	//this->basis_vector_functions_ = grid.cell_basi


	this->num_cells_ = grid.num_cells();

}
std::vector<std::vector<Euclidean_Vector>> Discretized_Solution_HOM::calculate_set_of_post_point_solutions(void) const {
	std::vector<std::vector<Euclidean_Vector>> set_of_post_point_solutions(this->num_cells_);

	for (uint icell = 0; icell < this->num_cells_; ++icell) 
		set_of_post_point_solutions[icell] = this->calculate_post_point_solutions(icell);

	return set_of_post_point_solutions;
}


Matrix_Wrapper Discretized_Solution_HOM::get_coefficient_matrix_wrapper(const uint icell) const {
	return { this->num_equations_, this->set_of_num_basis_[icell], this->discretized_solutions_.data() + this->set_of_data_start_index_[icell] };
}

std::vector<Euclidean_Vector> Discretized_Solution_HOM::calculate_post_point_solutions(const uint icell) const {
	const auto solution_post_points_m = this->get_coefficient_matrix_wrapper(icell) * this->set_of_basis_post_points_[icell];	//E x PP
	const auto [num_equations, num_post_points] = solution_post_points_m.size();

	std::vector<Euclidean_Vector> solution_at_post_points;
	solution_at_post_points.reserve(num_post_points);

	for (ushort i = 0; i < num_post_points; ++i)
		solution_at_post_points.push_back(solution_post_points_m.column(i));

	return solution_at_post_points;
}