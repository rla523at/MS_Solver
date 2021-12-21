#include "../INC/Log.h"


Log& Log::operator<<(const Command& command)
{
	switch (command)
	{
	case Command::print:
	{
		this->print();
		break;
	}
	case Command::print_off:
	{
		this->is_print_on = false;
		break;
	}
	case Command::print_on:
	{
		this->is_print_on = true;
		break;
	}
	default:
	{
		EXCEPTION("not supported command");
		break;
	}
	}
	return *this;
}

void Log::clear(void)
{
	this->log_txt_.clear();
}

Log& Log::get_instance(void)
{
	static Log instance;
	return instance;
}

void Log::set_configuration(const Configuration& configuration)
{
	const auto& write_log_file = configuration.get_write_log_file();

	if (ms::compare_icase(write_log_file, "Yes"))
	{
		this->do_write = true;
	}
	else if (ms::compare_icase(write_log_file, "No"))
	{
		this->do_write = false;
	}
}

void Log::set_path(const std::string& grid_folder_path)
{
	this->log_folder_path_ = grid_folder_path;
	this->error_log_folder_path_ = this->make_error_file_path(grid_folder_path);
}

void Log::write(void) const
{
	if (!this->do_write)
	{
		return;
	}

	REQUIRE(!this->log_folder_path_.empty(), "path should be set before write");

	const auto file_path = this->log_folder_path_ + "_log.txt";
	this->log_txt_.write(file_path);
}

void Log::write_error_text(const std::string& grid_file_name, const std::vector<double>& error_values) const
{
	if (!this->do_write)
	{
		return;
	}

	Sentence error_sentence = grid_file_name + "      \t";
	error_sentence.insert_with_space(error_values);
	error_sentence << "\n";

	Text txt = { error_sentence };
	txt.add_write(this->error_log_folder_path_ + "error.txt");
}

void Log::print(void)
{
	auto log_str = this->content_.str();

	if (this->is_print_on)
		std::cout << log_str;

	this->log_txt_ << std::move(log_str);

	std::ostringstream tmp;
	this->content_.swap(tmp);
}

std::string Log::make_error_file_path(const std::string& path)
{
	auto result = path;

	const auto pos = ms::rfind_nth(result, "/", 2);
	result.erase(pos + 1);

	return result;
}