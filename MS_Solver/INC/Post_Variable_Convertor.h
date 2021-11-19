#pragma once
#include "Discrete_Solution.h"
#include "Grid.h"

class Post_Variable_Convertor
{
public:
	Post_Variable_Convertor(const Grid& grid, const ushort post_order);

public://Query
	std::vector<double> convert_cell_center_values(const std::vector<double>& values) const;
	std::vector<double> convert_values(const std::vector<double>& values) const;

	virtual std::vector<std::vector<double>> calculate_set_of_post_point_solution_values(void) const abstract;
	virtual std::vector<std::vector<double>> calculate_set_of_cell_center_solution_values(void) const abstract;
	virtual const std::vector<std::string>& get_solution_names(void) const abstract;

protected:
	size_t num_cells_ = 0;
	size_t num_post_points_ = 0;
	size_t num_post_elements_ = 0;
	std::vector<ushort> set_of_num_post_elements_;
	std::vector<ushort> set_of_num_post_points_;
};

class Post_Variable_Convertor_DG : public Post_Variable_Convertor
{
public:
	Post_Variable_Convertor_DG(const Grid& grid, const ushort post_order, const Discrete_Solution_DG& discrete_solution);

public:
	std::vector<std::vector<double>> calculate_set_of_post_point_solution_values(void) const override;
	std::vector<std::vector<double>> calculate_set_of_cell_center_solution_values(void) const override;
	const std::vector<std::string>& get_solution_names(void) const override;

private:
	std::vector<double> P0_basis_values_;
	std::vector<Matrix> set_of_basis_post_points_m_;
	const Discrete_Solution_DG& discrete_solution_;
};


class Post_Variable_Converter_Factory//static class
{
public:
	static std::unique_ptr<Post_Variable_Convertor> make_unique(const Grid& grid, const ushort post_order, const Discrete_Solution_DG& discrete_solution);

private:
	Post_Variable_Converter_Factory(void) = delete;
};
