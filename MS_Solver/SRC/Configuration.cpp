#include "../INC/Configuration.h"

Configuration::Configuration(const std::string_view file_path)
{
	this->date_time_str_ = ms::date_time_string();
	const auto config_text = this->read_file(file_path);
	this->set_value(config_text);
};

std::string Configuration::get(const std::string_view config_name) const
{
	REQUIRE(this->name_to_value_str_.contains(ms::get_upper_case(config_name.data())), std::string(config_name.data()) + " should be exist in configuration file");
	return name_to_value_str_.at(ms::get_upper_case(config_name.data()));
};

std::string Configuration::post_folder_path_str(void) const
{
	const auto post_base_path = this->get("Post_base_path");
	REQUIRE(post_base_path.back() == '/', "post base path should be end with /");

	const auto governing_equation_name = this->get("governing_equation");
	const auto initial_condition_name = this->get("initial_condition");
	const auto grid_file_name = this->grid_file_name();

	return post_base_path + governing_equation_name + "/" + initial_condition_name + "/" + grid_file_name + "_" + this->date_time_str_ + "/";
}

std::string Configuration::configuration_str(void) const
{
	const auto space_dimension = this->get("space_dimension");
	const auto governing_equation_str = this->governing_equation_str();
	const auto initial_condition_str = this->initial_condition_str();
	const auto grid_file_name = this->grid_file_name();
	const auto numerical_flux_name = this->get("numerical_flux");
	const auto time_discrete_scheme = this->get("time_discrete_scheme");	
	const auto time_step_str = this->time_step_str();
	const auto solve_end_condition_str = this->solve_end_condition_str();

	std::ostringstream os;

	os << "current date : " << this->date_time_str_ << "\n\n";
	os << "================================================================================\n";
	os << "\t\t\t\t SETTING \n";
	os << "================================================================================\n";
	os << std::left << std::setw(40) << "Space dimension" << space_dimension << "\n";
	os << std::left << std::setw(40) << "Governing Equation" << governing_equation_str << "\n";
	os << std::left << std::setw(40) << "Initial Condtion" << initial_condition_str << "\n";
	os << std::left << std::setw(40) << "Grid" << grid_file_name << "\n";
	//os << "Spatial Discrete Method" << SPATIAL_DISCRETE_METHOD::name() << "\n";
	//if constexpr (SCAILING_METHOD_FLAG)
	//	os << "Reconstruction Method" << RECONSTRUCTION_METHOD::name() << " with scailing method\n";
	//else
	//	os << "Reconstruction Method" << RECONSTRUCTION_METHOD::name() << "\n";

	os << std::left << std::setw(40) << "Numeraical Flux Function" << numerical_flux_name << "\n";
	os << std::left << std::setw(40) << "Time Discrete Scheme" << time_discrete_scheme << "\n";
	os << std::left << std::setw(40) << "Time Step Method" << time_step_str << "\n";
	os << std::left << std::setw(40) << "Solve End Condtion" << solve_end_condition_str << "\n";
	os << "================================================================================\n";
	os << "================================================================================\n\n";

	return os.str();
}

std::string Configuration::governing_equation_str(void) const
{
	const auto governing_equation_name = this->get("governing_equation");

	if (ms::contains_icase(governing_equation_name, "Linear", "Advection"))
	{
		const auto space_dimension = this->get<short>("space_dimension");

		if (space_dimension == 2)
		{
			const auto x_advection_speed = this->get("x_advection_speed");
			const auto y_advection_speed = this->get("y_advection_speed");

			return governing_equation_name + "(" + x_advection_speed + "," + y_advection_speed + ")";
		}
		else if (space_dimension == 3)
		{
			const auto x_advection_speed = this->get("x_advection_speed");
			const auto y_advection_speed = this->get("y_advection_speed");
			const auto z_advection_speed = this->get("z_advection_speed");

			return governing_equation_name + "(" + x_advection_speed + "," + y_advection_speed + "," + z_advection_speed + ")";
		}
		else
		{
			EXCEPTION("space dimension in configuration file is not supported");
			return "";
		}
	}
	else
	{
		return governing_equation_name;
	}
}


std::string Configuration::grid_file_name(void) const
{
	const auto grid_file_path = this->get("grid_file_path");
	const auto parsed_strs = ms::parse(grid_file_path, '/');
	return ms::parse(parsed_strs.back(), '.').front();
}

std::string Configuration::initial_condition_str(void) const
{
	const auto initial_condition_name = this->get("initial_condition");

	if (ms::contains_icase(initial_condition_name, "Sine", "Wave"))
	{
		const auto space_dimension = this->get<short>("space_dimension");

		if (space_dimension == 2)
		{
			const auto x_wave_length = this->get("x_wave_length");
			const auto y_wave_length = this->get("y_wave_length");

			return initial_condition_name + "(" + x_wave_length + "," + y_wave_length + ")";
		}
		else if (space_dimension == 3)
		{
			const auto x_wave_length = this->get("x_wave_length");
			const auto y_wave_length = this->get("y_wave_length");
			const auto z_wave_length = this->get("z_wave_length");

			return initial_condition_name + "(" + x_wave_length + "," + y_wave_length + "," + z_wave_length + ")";
		}
		else
		{
			EXCEPTION("space dimension in configuration file is not supported");
			return "";
		}
	}
	else
	{
		return initial_condition_name;
	}
}

Text Configuration::read_file(const std::string_view file_path) const
{
	Text config_text;
	config_text.read(file_path);

	const std::vector<char> trim_character_set = { ' ', '\t' };
	const std::string comment_tag = "//";

	for (auto& sentence : config_text)
	{
		sentence.remove(trim_character_set);

		const auto comment_pos = sentence.find_position(comment_tag);

		if (comment_pos != std::string::npos)
		{
			sentence.remove_from_here(comment_pos);
		}
	}

	config_text.remove_empty_line();
	return config_text;
};

void Configuration::set_value(const Text& config_text)
{
	const auto delimiter = '=';

	for (const auto& sentence : config_text)
	{
		auto parsed_sentences = sentence.parse(delimiter);

		auto name_str = parsed_sentences.front().get_upper_case().get_string();
		auto value_str = parsed_sentences.back().get_string();

		this->name_to_value_str_.emplace(std::move(name_str), std::move(value_str));
	}
}

std::string Configuration::solve_end_condition_str(void) const
{
	const auto solve_end_controller_type = this->get("solve_end_controller_type");

	if (ms::contains_icase(solve_end_controller_type, "Time"))
	{
		return solve_end_controller_type + "(" + this->get("solve_end_time") + ")";
	}
	else if (ms::contains_icase(solve_end_controller_type, "Iter"))
	{
		return solve_end_controller_type + "(" + this->get("solve_end_iter") + ")";
	}
	else
	{
		EXCEPTION("time step calculator type in configuration file is not supported");
		return "";
	}
}

std::string Configuration::time_step_str(void) const
{
	const auto time_step_calculator_name = this->get("time_step_calculator_type");

	if (ms::contains_icase(time_step_calculator_name, "CFL"))
	{
		return time_step_calculator_name + "(" + this->get("CFL_number") + ")";
	}
	else if (ms::contains_icase(time_step_calculator_name, "Constant_dt"))
	{
		return time_step_calculator_name + "(" + this->get("Constant_dt") + ")";
	}
	else
	{
		EXCEPTION("time step calculator type in configuration file is not supported");
		return "";
	}
}


namespace ms
{
	std::string date_time_string(void)
	{
		time_t date_calculator = time(nullptr);
		tm date;
		localtime_s(&date, &date_calculator);

		return std::to_string(date.tm_year + 1900) + "." + std::to_string(date.tm_mon + 1) + "." + std::to_string(date.tm_mday)
			+ "__(" + std::to_string(date.tm_hour) + "£º" + std::to_string(date.tm_min) + ")";
	}
}