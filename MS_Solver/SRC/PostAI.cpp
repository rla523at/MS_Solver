#include "../INC/PostAI.h"

void Post_AI_Data::set_path(const std::string& path) {
#ifdef POST_AI_DATA

	path_ = path;

#endif
}

void Post_AI_Data::post(void) {
#ifdef POST_AI_DATA

	static size_t num_post = 1;
	static size_t num_post_data = 1;

	auto file_name = "AI_Solver_Data_" + std::to_string(num_post) + ".txt";
	auto file_path = path_ + file_name;

	constexpr size_t num_solution_str = 3;
	for (const auto target_cell_index : target_cell_indexes_) {
		auto& data_text = ai_data_text_set_.at(target_cell_index);

		data_text.front() = "#########################" + std::to_string(num_post_data++);

		data_text.add_write(file_path);
		data_text.erase(data_text.end() - num_solution_str, data_text.end());

		if (10000 < num_post_data) {
			num_post_data = 1;

			file_name = "AI_Solver_Data_" + std::to_string(++num_post) + ".txt";
			file_path = path_ + file_name;
		}
	}

	target_cell_indexes_.clear();

#endif
}

std::vector<std::string> Post_AI_Data::convert_to_solution_gradient_strings(const std::vector<Dynamic_Matrix_>& solution_gradients) {
	const auto num_solution = solution_gradients.size();
	const auto [num_equation, space_dimension] = solution_gradients.front().size();

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

void Post_AI_Data::record_limiting_value(const size_t index, const std::vector<double>& limiting_value) {
#ifdef POST_AI_DATA
	size_t num_equation = limiting_value.size();

	if (std::find(target_cell_indexes_.begin(), target_cell_indexes_.end(), index) == target_cell_indexes_.end())
		return;

	std::string limiting_value_string = "@limiterFunction\n";

	for (size_t i = 0; i < num_equation; ++i)
		limiting_value_string += ms::double_to_str_sp(limiting_value[i]) + "\t";
	limiting_value_string += "\n";

	ai_data_text_set_[index] << std::move(limiting_value_string);

#endif
}