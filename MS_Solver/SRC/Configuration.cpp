#include "../INC/Configuration.h"

Configuration::Configuration(const std::string_view file_path)
{
	this->date_time_str_ = ms::date_time_string();
	const auto config_text = this->read_file(file_path);
	this->set_value(config_text);
};

ushort Configuration::space_dimension(void) const
{
	return this->space_dimension_;
}

const std::string& Configuration::get_governing_equation(void) const
{
	return this->governing_equation_;
}

const std::string& Configuration::get_initial_condition(void) const
{
	return this->initial_condition_;
}

const std::string& Configuration::get_grid_file_path(void) const
{
	return this->grid_file_path_;
}

const std::string& Configuration::get_grid_file_type(void) const
{
	return this->grid_file_type_;
}

const std::string& Configuration::get_spatial_discrete_scheme(void) const
{
	return this->spatial_discrete_scheme_;
}

ushort Configuration::solution_degree(void) const
{
	return this->solution_degree_;
}

const std::string& Configuration::get_numerical_flux(void) const
{
	return this->numerical_flux_;
}

const std::string& Configuration::get_time_discrete_scheme(void) const
{
	return this->time_discrete_scheme_;
}

const std::string& Configuration::get_time_step_method(void) const
{
	return this->time_step_method_;
}

double Configuration::CFL_number(void) const
{
	return this->CFL_number_;
}

double Configuration::constant_dt(void) const
{
	return this->constant_dt_;
}

const std::string& Configuration::get_solve_end_controller_type(void) const
{
	return this->solve_end_controller_type_;
}

double Configuration::solve_end_time(void) const
{
	return this->solve_end_time_;
}

size_t Configuration::solve_end_iter(void) const
{
	return this->solve_end_iter_;
}

const std::string& Configuration::get_solve_post_controller_type(void) const
{
	return this->solve_post_controller_type_;
}

double Configuration::solve_post_time_step(void) const
{
	return this->solve_post_time_step_;
}

size_t Configuration::solve_post_iter_unit(void) const
{
	return this->solve_post_iter_unit_;
}

const std::string& Configuration::get_post_base_path(void) const
{
	return this->post_base_path_;
}

const std::string& Configuration::get_post_file_format(void) const
{
	return this->post_file_format_;
}

const std::string& Configuration::get_post_point_location(void) const
{
	return this->post_point_location_;
}

ushort Configuration::post_order(void) const
{
	return this->post_order_;
}

const double* Configuration::advection_speeds_ptr(void) const
{
	return this->advection_speeds_.data();
}

const double* Configuration::wave_lengths_ptr(void) const
{
	return this->wave_lengths_.data();
}

const double* Configuration::periodic_lengths_ptr(void) const
{
	return this->periodic_lengths_.data();
}

std::string Configuration::post_folder_path_str(void) const
{
	REQUIRE(this->post_base_path_.back() == '/', "post base path should be end with /");
	return this->post_base_path_ + this->governing_equation_str() + "/" + this->initial_condition_str() + "/" + this->grid_file_name() + "_" + this->date_time_str_ + "/";
}

std::string Configuration::configuration_str(void) const
{
	const auto governing_equation_str = this->governing_equation_str();
	const auto initial_condition_str = this->initial_condition_str();
	const auto grid_file_name = this->grid_file_name();
	const auto time_step_str = this->time_step_str();
	const auto solve_end_condition_str = this->solve_end_condition_str();

	std::ostringstream os;

	os << "current date : " << this->date_time_str_ << "\n\n";
	os << "================================================================================\n";
	os << "\t\t\t\t Configuration \n";
	os << "================================================================================\n";
	os << std::left << std::setw(40) << "Space dimension" << this->space_dimension_ << "\n";
	os << std::left << std::setw(40) << "Governing Equation" << governing_equation_str << "\n";
	os << std::left << std::setw(40) << "Initial Condtion" << initial_condition_str << "\n";
	os << std::left << std::setw(40) << "Grid" << grid_file_name << "\n";
	//os << "Spatial Discrete Method" << SPATIAL_DISCRETE_METHOD::name() << "\n";
	//if constexpr (SCAILING_METHOD_FLAG)
	//	os << "Reconstruction Method" << RECONSTRUCTION_METHOD::name() << " with scailing method\n";
	//else
	//	os << "Reconstruction Method" << RECONSTRUCTION_METHOD::name() << "\n";

	os << std::left << std::setw(40) << "Numeraical Flux Function" << this->numerical_flux_ << "\n";
	os << std::left << std::setw(40) << "Time Discrete Scheme" << this->time_discrete_scheme_ << "\n";
	os << std::left << std::setw(40) << "Time Step Method" << time_step_str << "\n";
	os << std::left << std::setw(40) << "Solve End Condtion" << solve_end_condition_str << "\n";
	os << "================================================================================\n";
	os << "================================================================================\n\n";

	return os.str();
}

std::string Configuration::governing_equation_str(void) const
{
	auto result = this->governing_equation_;

	if (ms::contains_icase(this->governing_equation_, "Linear", "Advection"))
	{
		result += "(";

		for (ushort i = 0; i < this->space_dimension_; ++i)
		{
			result += std::to_string(this->advection_speeds_[i]) + ", ";
		}

		result.pop_back();
		result.pop_back();

		result += ")";
	}

	return result;
}


std::string Configuration::grid_file_name(void) const
{
	const auto parsed_strs = ms::parse(this->grid_file_path_, '/');
	return ms::parse(parsed_strs.back(), '.').front();
}

std::string Configuration::initial_condition_str(void) const
{
	std::string result = this->initial_condition_;

	if (ms::contains_icase(this->initial_condition_, "Sine", "Wave"))
	{
		result += "(";

		for (ushort i = 0; i < this->space_dimension_; ++i)
		{
			result += std::to_string(this->wave_lengths_[i]) + ", ";
		}

		result.pop_back();
		result.pop_back();

		result += ")";
	}
		
	return result;
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
	std::map<Sentence, Sentence> name_to_value;

	for (const auto& sentence : config_text)
	{
		auto parsed_sentences = sentence.parse('=');

		REQUIRE(parsed_sentences.size() == 2, "configuration file has wrong format");

		parsed_sentences.front().upper_case(); //name part

		name_to_value.emplace(std::move(parsed_sentences.front()), std::move(parsed_sentences.back()));
	}

	this->space_dimension_ = this->get<ushort>(name_to_value, "space_dimension");

	this->governing_equation_ = this->get<std::string>(name_to_value, "governing_equation");
	this->initial_condition_ = this->get<std::string>(name_to_value, "initial_condition");

	this->grid_file_path_ = this->get<std::string>(name_to_value, "grid_file_path");
	this->grid_file_type_ = this->get<std::string>(name_to_value, "grid_file_type");

	this->spatial_discrete_scheme_ = this->get<std::string>(name_to_value, "spatial_discrete_scheme");
	this->solution_degree_ = this->get<ushort>(name_to_value, "solution_degree");

	this->numerical_flux_ = this->get<std::string>(name_to_value, "numerical_flux");

	this->time_discrete_scheme_ = this->get<std::string>(name_to_value, "time_discrete_scheme");

	this->time_step_method_ = this->get<std::string>(name_to_value, "time_step_method");
	this->CFL_number_ = this->get<double>(name_to_value, "CFL_number");
	this->constant_dt_ = this->get<double>(name_to_value, "constant_dt");

	this->solve_end_controller_type_ = this->get<std::string>(name_to_value, "solve_end_controller_type");
	this->solve_end_time_ = this->get<double>(name_to_value, "solve_end_time");
	this->solve_end_iter_ = this->get<size_t>(name_to_value, "solve_end_iter");

	this->solve_post_controller_type_ = this->get<std::string>(name_to_value, "solve_post_controller_type");
	this->solve_post_time_step_ = this->get<double>(name_to_value, "solve_post_time_step");
	this->solve_post_iter_unit_ = this->get<size_t>(name_to_value, "solve_post_iter_unit");

	this->post_base_path_ = this->get<std::string>(name_to_value, "post_base_path");
	this->post_file_format_ = this->get<std::string>(name_to_value, "post_file_format");
	this->post_point_location_ = this->get<std::string>(name_to_value, "post_point_location");
	this->post_order_ = this->get<ushort>(name_to_value, "post_order");


	const auto x_advection_speed = this->get<double>(name_to_value, "x_advection_speed");
	const auto y_advection_speed = this->get<double>(name_to_value, "y_advection_speed");
	const auto z_advection_speed = this->get<double>(name_to_value, "z_advection_speed");
	this->advection_speeds_ = { x_advection_speed, y_advection_speed, z_advection_speed };

	const auto x_wave_length = this->get<double>(name_to_value, "x_wave_length");
	const auto y_wave_length = this->get<double>(name_to_value, "y_wave_length");
	const auto z_wave_length = this->get<double>(name_to_value, "z_wave_length");
	this->wave_lengths_ = { x_wave_length, y_wave_length, z_wave_length };

	const auto x_periodic_length = this->get<double>(name_to_value, "x_periodic_length");
	const auto y_periodic_length = this->get<double>(name_to_value, "y_periodic_length");
	const auto z_periodic_length = this->get<double>(name_to_value, "z_periodic_length");
	this->periodic_lengths_ = { x_periodic_length, y_periodic_length, z_periodic_length };
}

std::string Configuration::solve_end_condition_str(void) const
{
	if (ms::contains_icase(this->solve_end_controller_type_, "Time"))
	{
		return this->solve_end_controller_type_ + "(" + std::to_string(this->solve_end_time_) + ")";
	}
	else if (ms::contains_icase(this->solve_end_controller_type_, "Iter"))
	{
		return this->solve_end_controller_type_ + "(" + std::to_string(this->solve_end_iter_) + ")";
	}
	else
	{
		EXCEPTION("time step calculator type in configuration file is not supported");
		return "";
	}
}

std::string Configuration::time_step_str(void) const
{
	if (ms::contains_icase(this->time_step_method_, "CFL"))
	{
		return this->time_step_method_ + "(" + std::to_string(this->CFL_number_) + ")";
	}
	else if (ms::contains_icase(this->time_step_method_, "Constant_dt"))
	{
		return this->time_step_method_ + "(" + std::to_string(this->constant_dt_) + ")";
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