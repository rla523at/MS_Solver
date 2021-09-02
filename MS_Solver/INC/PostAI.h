#pragma once
#include "Grid_Builder.h"

class Post_AI_Data
{
private:
	using This_ = Post_AI_Data;

public:
	static inline bool post_condition_ = false;

public:
	inline static std::string path_;
	inline static size_t num_data_;
	inline static Text comment_;
	inline static std::vector<Text> ai_data_text_set_;
	inline static std::vector<std::vector<size_t>> vertex_share_cell_indexes_set_;
	inline static std::vector<size_t> target_cell_indexes_;

private:
	Post_AI_Data(void) = delete;
	
public:
	static void set_path(const std::string& path) { This_::path_ = path; }

	template <size_t space_dimension>
	static void intialize(const Grid<space_dimension>& grid);

	template <size_t num_equation, size_t space_dimension>
	static void record_solution_datas(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::vector<Matrix<num_equation, space_dimension>>& solution_gradients);

	template <size_t num_equation, size_t space_dimension>
	static void conditionally_record_solution_datas(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::vector<Matrix<num_equation, space_dimension>>& solution_gradients);

	template <ushort num_equation>
	static void record_limiting_value(const size_t index, const std::array<double,num_equation>& limiting_value);

	template <ushort num_equation>
	static void conditionally_record_limiting_value(const size_t index, const std::array<double, num_equation>& limiting_value);

	static void post(void);
	
	static void conditionally_post(void);

	static void post_scatter_data(const std::vector<double>& limiting_values);



//private: //for test
	template <size_t space_dimension>
	static auto calculate_face_share_cell_indexes_set(const Grid<space_dimension>& grid);

	template <size_t space_dimension>
	static auto calculate_vertex_nodes_coordinate_string_set(const Grid<space_dimension>& grid);

	template <size_t num_equation>
	static auto convert_to_solution_strings(const std::vector<Euclidean_Vector<num_equation>>& solutions);

	template <size_t num_equation, size_t space_dimension>
	static std::vector<std::string> convert_to_solution_gradient_strings(const std::vector<Matrix<num_equation, space_dimension>>& solution_gradients);
};


//template definition part
template <size_t space_dimension>
void Post_AI_Data::intialize(const Grid<space_dimension>& grid) {

	const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_index_set;
	const auto& cell_elements = grid.elements.cell_elements;
	This_::num_data_ = cell_elements.size();

	This_::vertex_share_cell_indexes_set_.reserve(num_data_);
	This_::ai_data_text_set_.resize(num_data_);
	This_::target_cell_indexes_.reserve(num_data_);

	auto file_path = path_ + "AI_Solver_Data_1.txt";	

	const auto parsed_path = ms::parse(path_, '/');
	const auto GE_comments = parsed_path[4];
	const auto IC_comments = parsed_path[5];
	const auto grid_comments = ms::parse(parsed_path[7], '_');

	comment_ << "***********************************\n"
		"***********************************\n"
		"G.E : " + GE_comments + "\n" +
		"I.C : " + IC_comments + "\n" +
		"Grid : " + grid_comments[0] + "\n" +
		"***********************************\n"
		"***********************************\n\n";	

	comment_.add_write(file_path);	

	const auto face_share_cell_indexes_set = calculate_face_share_cell_indexes_set(grid);
	const auto vnodes_coordinate_string_set = calculate_vertex_nodes_coordinate_string_set(grid);

	for (size_t i = 0; i < num_data_; ++i) {
		const auto& cell_element = cell_elements[i];
		const auto& cell_geometry = cell_element.geometry_;

		//vertex share cell indexes temp
		std::set<size_t> vertex_share_cell_indexes_temp;

		const auto vnode_indexes = cell_element.vertex_node_indexes();
		for (const auto& vnode_index : vnode_indexes) {
			const auto& vnode_share_cell_indexes = vnode_index_to_share_cell_indexes.at(vnode_index);
			vertex_share_cell_indexes_temp.insert(vnode_share_cell_indexes.begin(), vnode_share_cell_indexes.end());
		}

		//chunk edge connectivities //quad3������ ����� �۵����� �ʴ� algorithm
		std::set<std::set<size_t>> face_share_cell_index_pairs;

		for (const auto chunk_cell_index : vertex_share_cell_indexes_temp) {
			const auto& face_share_cell_indexes = face_share_cell_indexes_set.at(chunk_cell_index);

			std::vector<size_t> face_share_cell_indexes_in_chunk;
			std::set_intersection(vertex_share_cell_indexes_temp.begin(), vertex_share_cell_indexes_temp.end(), face_share_cell_indexes.begin(), face_share_cell_indexes.end(), std::back_inserter(face_share_cell_indexes_in_chunk));

			for (const auto face_share_cell_index_in_chunk : face_share_cell_indexes_in_chunk)
				face_share_cell_index_pairs.insert({ chunk_cell_index, face_share_cell_index_in_chunk });
		}

		//vertex_share_cell_indexes
		vertex_share_cell_indexes_temp.erase(i);

		std::vector<size_t> vertex_share_cell_indexes;
		vertex_share_cell_indexes.push_back(i);
		vertex_share_cell_indexes.insert(vertex_share_cell_indexes.end(), vertex_share_cell_indexes_temp.begin(), vertex_share_cell_indexes_temp.end());

		// header string
		ai_data_text_set_[i] << "temporary header text";

		// node number string
		const auto num_node = vertex_share_cell_indexes.size();
		ai_data_text_set_[i] << "@nodeNumber\n" + std::to_string(num_node);

		// node index order string
		std::string node_index_order_string = "@nodeIndexOrder\n";
		for (const auto& node_index : vertex_share_cell_indexes)
			node_index_order_string += std::to_string(node_index) + "\t";
		ai_data_text_set_[i] << std::move(node_index_order_string);

		//edge number string
		const auto num_edge = face_share_cell_index_pairs.size();
		ai_data_text_set_[i] << "@edgeNumber\n" + std::to_string(num_edge);

		//connectivity string
		std::string node_connectivity_string = "@connectivity\n";
		for (const auto& chunk_edge_connectivity : face_share_cell_index_pairs) {
			for (const auto& node_index : chunk_edge_connectivity)
				node_connectivity_string += std::to_string(node_index) + "\t";
			node_connectivity_string += "\n";
		}
		node_connectivity_string.pop_back();
		ai_data_text_set_[i] << std::move(node_connectivity_string);

		//cell coords string
		std::string cell_coords_string = "@cellCoords\n";
		for (const auto vertex_share_cell_index : vertex_share_cell_indexes)
			cell_coords_string += vnodes_coordinate_string_set[vertex_share_cell_index];
		cell_coords_string.pop_back();

		ai_data_text_set_[i] << std::move(cell_coords_string);

		//vertex_share_cell_indexes
		vertex_share_cell_indexes_set_.push_back(std::move(vertex_share_cell_indexes));
	}

}


template <size_t num_equation, size_t space_dimension>
void Post_AI_Data::record_solution_datas(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::vector<Matrix<num_equation, space_dimension>>& solution_gradients) {

	dynamic_require(num_data_ == solutions.size(),			"number of solution should be same with number of data");
	dynamic_require(num_data_ == solution_gradients.size(), "number of solution gradient should be same with number of data");

	const auto solution_strings = convert_to_solution_strings(solutions);	//?
	const auto solution_gradient_strings = convert_to_solution_gradient_strings(solution_gradients);

	std::string cell_average_string;
	std::string cell_gradient_string;
	for (size_t i = 0; i < num_data_; ++i) {
		const auto& vertex_share_cell_indexes = vertex_share_cell_indexes_set_[i];
		const auto num_vertex_share_cell = vertex_share_cell_indexes.size();


		std::vector<double> vertex_share_cell_solutions(num_vertex_share_cell);
		for (size_t i = 0; i < num_vertex_share_cell; ++i)
			vertex_share_cell_solutions[i] = solutions[vertex_share_cell_indexes[i]].at(0);

		const auto min_solution = *std::min_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
		const auto max_solution = *std::max_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
		const auto solution_diff = max_solution - min_solution;

		if (solution_diff < 0.01)
			continue;
			
		target_cell_indexes_.push_back(i);

		cell_average_string = "@cellAverage\n";
		cell_gradient_string = "@cellGradient\n";
		for (const auto vertex_share_cell_index : vertex_share_cell_indexes) {
			cell_average_string += solution_strings[vertex_share_cell_index] + "\n";
			cell_gradient_string += solution_gradient_strings[vertex_share_cell_index] + "\n";
		}
		cell_average_string.pop_back();
		cell_gradient_string.pop_back();

		ai_data_text_set_[i] << std::move(cell_average_string) << std::move(cell_gradient_string);
	}

}

template <size_t num_equation, size_t space_dimension>
void Post_AI_Data::conditionally_record_solution_datas(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::vector<Matrix<num_equation, space_dimension>>& solution_gradients) {
	if (This_::post_condition_)
		This_::record_solution_datas(solutions, solution_gradients);
}



template <ushort num_equation>
void Post_AI_Data::record_limiting_value(const size_t index, const std::array<double, num_equation>& limiting_value) {

	if (std::find(target_cell_indexes_.begin(), target_cell_indexes_.end(), index) == target_cell_indexes_.end())
		return;

	std::string limiting_value_string = "@limiterFunction\n";

	for (size_t i = 0; i < num_equation; ++i)
		limiting_value_string += ms::double_to_str_sp(limiting_value[i]) + "\t";
	limiting_value_string += "\n";

	ai_data_text_set_[index] << std::move(limiting_value_string);
}

template <ushort num_equation>
void Post_AI_Data::conditionally_record_limiting_value(const size_t index, const std::array<double, num_equation>& limiting_value) {
	if (This_::post_condition_)
		This_::record_limiting_value(index, limiting_value);
}




template <size_t space_dimension>
auto Post_AI_Data::calculate_face_share_cell_indexes_set(const Grid<space_dimension>& grid) {
	const auto& vnode_index_to_share_cell_index_set = grid.connectivity.vnode_index_to_share_cell_index_set;
	const auto& cell_elements = grid.elements.cell_elements;
	const auto num_cell = cell_elements.size();

	//face share cell indexes set
	std::vector<std::set<size_t>> face_share_cell_indexes_set;
	face_share_cell_indexes_set.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = cell_elements[i].geometry_;

		const auto face_vnode_indexes_set = element.face_vertex_node_indexes_set();
		const auto num_face = face_vnode_indexes_set.size();

		std::set<size_t> face_share_cell_indexes;

		for (const auto& face_vnode_indexes : face_vnode_indexes_set) {
			std::vector<size_t> this_face_share_cell_indexes;

			const auto num_face_vnode = face_vnode_indexes.size();

			const auto& set_0 = vnode_index_to_share_cell_index_set.at(face_vnode_indexes[0]);
			const auto& set_1 = vnode_index_to_share_cell_index_set.at(face_vnode_indexes[1]);
			std::set_intersection(set_0.begin(), set_0.end(), set_1.begin(), set_1.end(), std::back_inserter(this_face_share_cell_indexes));

			if (2 < num_face_vnode) {
				std::vector<size_t> buffer;
				for (size_t i = 2; i < num_face_vnode; ++i) {
					const auto& set_i = vnode_index_to_share_cell_index_set.at(face_vnode_indexes[i]);

					buffer.clear();
					std::set_intersection(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), set_i.begin(), set_i.end(), std::back_inserter(buffer));
					std::swap(this_face_share_cell_indexes, buffer);
				}
			}

			const auto my_index_pos_iter = std::find(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), i);
			dynamic_require(my_index_pos_iter != this_face_share_cell_indexes.end(), "my index should be included in this face share cell indexes");

			this_face_share_cell_indexes.erase(my_index_pos_iter);
			dynamic_require(this_face_share_cell_indexes.size() == 1, "face share cell should be unique");

			face_share_cell_indexes.insert(this_face_share_cell_indexes.front());
		}

		face_share_cell_indexes_set.push_back(std::move(face_share_cell_indexes));
	}

	return face_share_cell_indexes_set;
}


template <size_t space_dimension>
auto Post_AI_Data::calculate_vertex_nodes_coordinate_string_set(const Grid<space_dimension>& grid) {
	const auto& cell_elements = grid.elements.cell_elements;
	const auto num_cell = cell_elements.size();

	std::vector<std::string> vnodes_coordinate_string_set;
	vnodes_coordinate_string_set.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = element.geometry_;

		const auto vnodes = geometry.vertex_nodes();
		const auto num_vnode = vnodes.size();
		std::string vnodes_coordinate_string = std::to_string(num_vnode) + "\n";

		for (const auto& vnode : vnodes) {
			for (size_t i = 0; i < space_dimension; ++i)
				vnodes_coordinate_string += ms::double_to_str_sp(vnode.at(i)) + "\t";

			vnodes_coordinate_string += "\n";
		}

		vnodes_coordinate_string_set.push_back(std::move(vnodes_coordinate_string));
	}

	return vnodes_coordinate_string_set;
}


template <size_t num_equation>
auto Post_AI_Data::convert_to_solution_strings(const std::vector<Euclidean_Vector<num_equation>>& solutions) {
	const auto num_solution = solutions.size();
	
	std::vector<std::string> solution_strings;
	solution_strings.reserve(num_solution);

	std::string solution_string;
	for (size_t i = 0; i < num_solution; ++i) {

		const auto& solution = solutions[i];
		for (size_t j = 0; j < num_equation; ++j)
			solution_string += ms::double_to_str_sp(solution.at(j)) + "\t";

		solution_strings.push_back(std::move(solution_string));
	}

	return solution_strings;
}


template <size_t num_equation, size_t space_dimension>
std::vector<std::string> Post_AI_Data::convert_to_solution_gradient_strings(const std::vector<Matrix<num_equation, space_dimension>>& solution_gradients) {
	const auto num_solution = solution_gradients.size();
	//const auto [num_equation, space_dimension] = solution_gradients.front().size();

	std::vector<std::string> solution_gradient_strings;
	solution_gradient_strings.reserve(num_solution);

	std::string solution_gradient_string;
	for (size_t i = 0; i < num_solution; ++i) {

		const auto& solution_gradient = solution_gradients[i];
		for (size_t j = 0; j < num_equation; ++j)
			for (size_t k = 0; k < space_dimension; ++k)
				solution_gradient_string += ms::double_to_str_sp(solution_gradient.at(j, k)) + "\t";

		solution_gradient_strings.push_back(std::move(solution_gradient_string));
	}

	return solution_gradient_strings;
}


