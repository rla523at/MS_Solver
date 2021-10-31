#pragma once
#include "Governing_Equation.h"
#include "Euclidean_Vector.h"
#include "Grid.h"
#include "Initial_Condition.h"

using ushort = unsigned short;
using uint = unsigned int;

namespace ms {
	template <typename T>
	void insert(std::vector<T>& vec, const std::vector<T>& other_vec) {
		vec.insert(vec.end(), other_vec.begin(), other_vec.end());
	}
}


class Discretized_Solution
{
protected:
	Discretized_Solution(const Grid& grid, const Initial_Condition& initial_condition) {
		this->num_equations_ = initial_condition.num_equations();
		this->num_cells_ = grid.num_cells();
		this->solution_variable_names_ = initial_condition.solution_variable_names();
	}

public:	//Query
	virtual std::vector<std::vector<double>> calculate_post_point_solutions_by_variable(void) const abstract;
	const std::vector<std::string>& get_variable_names(void) const {
		return this->solution_variable_names_;
	}

protected:
	ushort num_equations_;
	size_t num_cells_;
	std::vector<std::string> solution_variable_names_;
	std::vector<double> discretized_solutions_;
};


class Discretized_Solution_FVM : public Discretized_Solution
{
public:
	Discretized_Solution_FVM(const Grid& grid, const Initial_Condition& initial_condition)
		: Discretized_Solution(grid, initial_condition) {
		this->discretized_solutions_.resize(num_cells_ * num_equations_);


	}

public: //Query
	std::vector<std::vector<double>> calculate_post_point_solutions_by_variable(void) const override {
		std::vector<std::vector<double>> post_point_solutions_by_variable(num_equations_);
		for (auto& vec : post_point_solutions_by_variable)
			vec.reserve(this->num_cells_);
				
		for (uint icell = 0; icell < this->num_cells_; ++icell) {
			const auto solution = this->calculate_solution_at_post_point(icell);
		
			for (ushort j = 0; j < this->num_equations_; ++j) 
				post_point_solutions_by_variable[j].push_back(solution[j]);
		}

		return post_point_solutions_by_variable;
	}

private:
	Euclidean_Vector calculate_solution_at_post_point(const uint icell) const {
		const auto start_index = icell * this->num_equations_;
		const auto num_variable = this->num_equations_;

		const auto start_iter = this->discretized_solutions_.begin() + start_index;

		std::vector<double> values = { start_iter, start_iter + num_variable };

		return std::move(values);
	}
};


class Discretized_Solution_HOM : public Discretized_Solution
{
public:
	std::vector<std::vector<double>> calculate_post_point_solutions_by_variable(void) const override {
		std::vector<std::vector<double>> post_point_solutions_by_variable(num_equations_);
		for (auto& vec : post_point_solutions_by_variable)
			vec.reserve(this->num_post_points_);

		const auto num_cell = this->set_of_num_basis_.size();

		for (uint i = 0; i < num_cell; ++i) {
			const auto solution_at_post_points = this->calculate_solution_at_post_points(i);
			const auto num_post_points = solution_at_post_points.size();

			for (ushort q = 0; q < num_post_points; ++q) {
				const auto& solution = solution_at_post_points[q];

				for (ushort j = 0; j < this->num_equations_; ++j)
					post_point_solutions_by_variable[j].push_back(solution[j]);
			}
		}

		return post_point_solutions_by_variable;
	}

private:
	Matrix_Wrapper get_coefficient_matrix_wrapper(const uint icell) const {
		return { this->num_equations_, this->set_of_num_basis_[icell], this->discretized_solutions_.data() + this->set_of_data_start_index_[icell] };
	}
	std::vector<Euclidean_Vector> calculate_solution_at_post_points(const uint icell) const {
		const auto solution_post_points_m = this->get_coefficient_matrix_wrapper(icell) * this->set_of_basis_post_points_[icell];	//E x PP
		const auto [num_equations, num_post_points] = solution_post_points_m.size();

		std::vector<Euclidean_Vector> solution_at_post_points;
		solution_at_post_points.reserve(num_post_points);

		for (ushort i = 0; i < num_post_points; ++i)
			solution_at_post_points.push_back(solution_post_points_m.column(i));

		return solution_at_post_points;
	}

private:
	size_t num_post_points_;

	std::vector<ushort> set_of_num_basis_;
	std::vector<size_t> set_of_data_start_index_;

	std::vector<Matrix> set_of_basis_post_points_;
};