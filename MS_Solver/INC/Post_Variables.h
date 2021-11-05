#pragma once
#include "Grid.h"
#include "Configuration.h"

class Post_Variables
{
public://Command
	void record_grid_data(const Grid& grid, const ushort post_order);
	void record_variable(const std::string_view name, const std::vector<double>& values);
	void set_post_variable_format(const std::string& post_variable_location_str);
	void syncronize_solution_time(const double& get_solution_time);

public://Query
	Zone_Type zone_type(void) const;
	size_t num_post_node(void) const;
	size_t num_post_element(void) const;
	ushort num_grid_variable(void) const;
	ushort num_solution_variable(void) const;
	double solution_time(void) const;
	std::string grid_variable_str(void) const;

	const std::vector<std::vector<double>>& get_post_nodes_by_axis(void) const;
	const std::vector<std::vector<int>>& get_connectivities(void) const;

	std::string solution_variable_str(void) const;
	std::string solution_variable_location_str(void) const {
		return value_format_convertor_->solution_variable_location_str(this->solution_variable_name_to_value_.size());
	}
	std::vector<std::vector<double>> calculate_set_of_solution_datas(void) const;

private:
	std::vector<std::vector<double>> make_post_coordinate_blocks(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_nodes) const;

private:
	std::unique_ptr<Post_Value_Format_Converter> value_format_convertor_;

	const double* solution_time_ptr_ = nullptr;

	Zone_Type zone_type_;
	size_t num_post_nodes_ = 0;
	size_t num_post_elements_ = 0;
	ushort num_grid_variables_ = 0;

	std::vector<std::vector<int>> set_of_connectivities_;
	std::vector<std::vector<double>> post_nodes_by_axis_;
	std::map<std::string, std::vector<double>> solution_variable_name_to_value_;
};

enum class Post_Variable_Location {
	node,
	cell_center
};

enum class Zone_Type {
	FETriangle = 2,
	FETetrahedron = 4
};

class Post_Value_Format_Converter
{
public://Command
	virtual void set_grid_data(const Grid& grid, const ushort post_order) abstract;

public://Query
	virtual std::vector<double> convert_to_post_variable_format(const std::vector<double>& values) const abstract;
	virtual std::string solution_variable_location_str(const size_t get_num_solution_variable) const abstract;

protected:
	size_t num_cells_;
	size_t num_post_nodes_;
};

class Cell_Center_Format_Convertor : public Post_Value_Format_Converter
{
public://Command
	void set_grid_data(const Grid& grid, const ushort post_order) override;

public://Query
	std::vector<double> convert_to_post_variable_format(const std::vector<double>& values) const override;
	std::string solution_variable_location_str(const size_t get_num_solution_variable) const override;
};

class Node_Format_Convertor : public Post_Value_Format_Converter
{
public://Command
	void set_grid_data(const Grid& grid, const ushort post_order) override;

public://Query
	std::vector<double> convert_to_post_variable_format(const std::vector<double>& values) const override;
	std::string solution_variable_location_str(const size_t get_num_solution_variable) const override;

private:
	std::vector<ushort> cell_set_of_num_post_nodes_;
};

//static class
class Post_Variable_Format_Converter_Factory
{
public:
	static std::unique_ptr<Post_Value_Format_Converter> make(const std::string_view type) {
		if (ms::contains_icase(type.data(), "cell", "center"))
			return std::make_unique<Cell_Center_Format_Convertor>();
		else if (ms::contains_icase(type.data(), "node"))
			return std::make_unique<Node_Format_Convertor>();
		else
			EXCEPTION("wrong type");
	}

private:
	Post_Variable_Format_Converter_Factory(void) = delete;
};

namespace ms {
	template <typename T>
	size_t size_of_vvec(const std::vector<std::vector<T>>& vvec) {
		size_t size = 0;
		for (const auto& vec : vvec)
			size += vec.size();
		return size;
	};
}