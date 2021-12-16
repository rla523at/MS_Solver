#pragma once
#include "Discrete_Solution.h"
#include "Exact_Solution.h"

class Error
{
public:
	virtual std::vector<double> calculate_error_norms(const Grid& grid, const Discrete_Solution_DG& discrete_soltuion, const double end_time) const abstract;

protected:
	std::unique_ptr<Exact_Solution> exact_solution_;
};

class No_Error_Calculation : public Error
{
public:
	std::vector<double> calculate_error_norms(const Grid& grid, const Discrete_Solution_DG& discrete_soltuion, const double end_time) const override;
};

class Global_Error : public Error
{
public:
	Global_Error(std::unique_ptr<Exact_Solution>&& exact_solution);

public:
	std::vector<double> calculate_error_norms(const Grid& grid, const Discrete_Solution_DG& discrete_soltuion, const double end_time) const override;

private:
	double calculate_global_error_L1_norm(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const;
	double calculate_global_error_L2_norm(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const;
	double calculate_global_error_Linf_norm(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const;
};

class Local_Average_Error : public Error
{
public:
	Local_Average_Error(std::unique_ptr<Exact_Solution>&& exact_solution);

public:
	std::vector<double> calculate_error_norms(const Grid& grid, const Discrete_Solution_DG& discrete_soltuion, const double end_time) const override;

private:
	double calculate_local_error_L1_norm(const uint cell_index, const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const;
	double calculate_local_error_L2_norm(const uint cell_index, const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const;
	double calculate_local_error_Linf_norm(const uint cell_index, const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const;
};

class Error_Factory
{
public:
	static std::unique_ptr<Error> make_unqiue(const Configuration& configuration);

private:
	Error_Factory(void) = delete;
};