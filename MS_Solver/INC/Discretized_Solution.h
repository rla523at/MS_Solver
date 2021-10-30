#pragma once
#include <vector>
#include <string_view>

#include "Governing_Equation.h"


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
public:	
	virtual std::vector<std::vector<double>> calculate_post_point_solutions_by_variable(void) const abstract;

	const std::vector<std::string>& get_variable_names(void) const {
		return this->governing_equation_->get_variable_names();
	}

protected:
	ushort num_equation_;
	std::vector<double> discretized_solutions_;
	std::unique_ptr<Governing_Equation> governing_equation_;
};


class Discretized_Solution_FVM : public Discretized_Solution
{
public:
	std::vector<std::vector<double>> calculate_post_point_solutions_by_variable(void) const override {
		const auto num_cell = this->discretized_solutions_.size() / num_equation_;
		
		std::vector<std::vector<double>> post_point_solutions_by_variable(num_equation_);

		for (ushort j = 0; j < this->num_equation_; ++j) {
			post_point_solutions_by_variable[j].resize(num_cell);
			for (uint i = 0; i < num_cell; ++i) {
				post_point_solutions_by_variable[j][i] = this->discretized_solutions_[i * this->num_equation_ + j];
			}
		}

		return post_point_solutions_by_variable;
	}
};


class Discretized_Solution_HOM : public Discretized_Solution
{
public:
	std::vector<std::vector<double>> calculate_post_point_solutions_by_variable(void) const override {
		const auto num_cell = this->set_of_num_basis_.size();
		
		std::vector<std::vector<double>> post_point_solutions_by_variable(num_equation_);

		for (uint i = 0; i < num_cell; ++i) {
			const auto solution_at_post_points = this->get_cell_coefficient(i) * this->set_of_basis_post_points_[i];	//E x PP
			for (ushort j = 0; j < this->num_equation_; ++j)
				ms::insert(post_point_solutions_by_variable[j], solution_at_post_points.row(j));
		}

		return post_point_solutions_by_variable;
	}

private:
	Matrix_Wrapper get_cell_coefficient(const uint icell) const {
		return { this->num_equation_, this->set_of_num_basis_[icell], this->discretized_solutions_.data() + this->set_of_data_start_index_[icell] };
	}

private:
	std::vector<ushort> set_of_num_basis_;
	std::vector<size_t> set_of_data_start_index_;

	std::vector<Matrix> set_of_basis_post_points_;
};