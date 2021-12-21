#pragma once
#include "Post_Variable_Convertor.h"

enum class Zone_Type
{
	FETriangle = 2,
	FETetrahedron = 4,
	Noting = -1
};

class Post_Variables
{
public:
	Post_Variables(const Grid& grid, std::unique_ptr<Post_Variable_Convertor>&& post_variable_convertor, const ushort post_order);

public://Command
	void clear_variables(void);
	void record_solution(void);
	void record_variable(const std::string& name, const std::vector<double>& values);
	void syncronize_solution_time(const double& compute);	

public://Query
	const std::vector<std::vector<double>>& get_post_coordinate_blocks(void) const;
	const std::vector<std::vector<int>>& get_connectivities(void) const;
	const std::vector<std::vector<double>>& get_set_of_post_variable_values(void) const;
	const std::vector<std::string>& get_post_variable_names(void) const;
	std::string grid_variable_str(void) const;
	size_t num_post_points(void) const;
	size_t num_post_elements(void) const;
	ushort num_grid_variable(void) const;
	ushort num_post_variables(void) const;
	double solution_time(void) const;
	std::string variable_location_str(void) const;
	Zone_Type zone_type(void) const;

	bool is_emptry(void) const;
private:
	std::vector<std::vector<double>> make_post_coordinate_blocks(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_nodes) const;

private:	
	std::unique_ptr<Post_Variable_Convertor> post_variable_convertor_;

	ushort space_dimension_ = 0;
	Zone_Type zone_type_;
	size_t num_post_points_ = 0;
	size_t num_post_elements_ = 0;
	const double* solution_time_ptr_ = nullptr;

	std::vector<std::vector<int>> set_of_connectivities_;
	std::vector<std::vector<double>> post_coordinate_blocks;

	std::vector<std::string> post_variable_names_;
	std::vector<std::vector<double>> set_of_post_variable_values_;
};

namespace ms 
{
	template <typename T>	size_t size_of_vvec(const std::vector<std::vector<T>>& vvec) 
	{
		size_t size = 0;
		for (const auto& vec : vvec)
		{
			size += vec.size();
		}

		return size;
	};
	template <typename T>	bool contains(const std::vector<T>& vec, const T& val)
	{
		return ms::find(vec, val) != vec.end();
	}

}