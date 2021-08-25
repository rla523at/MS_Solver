#include "../INC/PostAI.h"

void Post_AI_Data::set_path(const std::string& path) {
#ifdef POST_AI_DATA_MODE

	path_ = path;

#endif
}

void Post_AI_Data::post(void) {
#ifdef POST_AI_DATA_MODE

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
			comment_.add_write(file_path);
		}
	}

	target_cell_indexes_.clear();
	//target_cell_indexes_.reset();		

#endif
}


