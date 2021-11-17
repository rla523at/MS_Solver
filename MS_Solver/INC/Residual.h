#pragma once
#include "Exception.h"
#include "Matrix.h"

#include <vector>

using uint = unsigned int;

class Residual
{
public:
    Residual(const size_t num_values, const std::vector<size_t>& coefficient_start_indexes);         

public://Command
    void update_rhs(const uint cell_index, const Matrix& delta_rhs);
    std::vector<double>&& move_values(void);

private:
    double* get_icell_residual(const uint cell_index);

private:
    const std::vector<size_t>& coefficieint_start_indexes_;
	std::vector<double> values_;
};