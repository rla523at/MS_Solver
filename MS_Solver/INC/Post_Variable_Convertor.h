#pragma once
#include "Discrete_Solution.h"
#include "Grid.h"

class Post_Variable_Convertor
{
public:
	Post_Variable_Convertor(const Grid& grid, const ushort post_order, const Discrete_Solution& discrete_solution);

public://Query
	virtual std::vector<double> convert_values(const std::vector<double>& values) const abstract;
	virtual std::vector<std::vector<double>> calculate_set_of_post_point_solution_values(void) const abstract;
	virtual std::string variable_location_str(void) const abstract;
	const std::vector<std::string>& get_solution_names(void) const;

protected:
	size_t num_cells_ = 0;
	const Discrete_Solution& discrete_solution_;
};

class Node_Base_Convertor : public Post_Variable_Convertor
{
public:
	Node_Base_Convertor(const Grid& grid, const ushort post_order, Discrete_Solution& discrete_solution);

public:
	std::vector<double> convert_values(const std::vector<double>& values) const override;
	std::vector<std::vector<double>> calculate_set_of_post_point_solution_values(void) const override;
	std::string variable_location_str(void) const override;

private:
	size_t num_post_points_ = 0;
	std::vector<ushort> set_of_num_post_points_;
};

class Center_Base_Convertor : public Post_Variable_Convertor
{
public:
	Center_Base_Convertor(const Grid& grid, const ushort post_order, Discrete_Solution& discrete_solution);

public:
	std::vector<double> convert_values(const std::vector<double>& values) const override;
	std::vector<std::vector<double>> calculate_set_of_post_point_solution_values(void) const override;
	std::string variable_location_str(void) const override;

private:
	size_t num_post_elements_ = 0;
	std::vector<ushort> set_of_num_post_elements_;
};

class Post_Variable_Converter_Factory//static class
{
public:
	static std::unique_ptr<Post_Variable_Convertor> make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution& discrete_solution);

private:
	Post_Variable_Converter_Factory(void) = delete;
};
