#include "../INC/Residual.h"

Residual::Residual(const size_t num_values, const std::vector<size_t>& coefficient_start_indexes)
    : coefficieint_start_indexes_(coefficient_start_indexes)
    , num_cell_(coefficient_start_indexes.size())
    , values_(num_values) {};

void Residual::update_rhs(const uint cell_index, const double* delta_rhs_ptr)
{
    const auto n = static_cast<MKL_INT>(this->num_icell_values(cell_index));
    const auto a = 1.0;
    const auto incx = 1;
    const auto incy = 1;

    auto icell_rhs_ptr = this->get_icell_residual(cell_index);
    cblas_daxpy(n, a, delta_rhs_ptr, incx, icell_rhs_ptr, incy);
}

std::vector<double>&& Residual::move_values(void)
{
    return std::move(this->values_);
}

double* Residual::get_icell_residual(const uint cell_index)
{
    return this->values_.data() + coefficieint_start_indexes_[cell_index];
}

size_t Residual::num_icell_values(const uint cell_index)
{
    if (cell_index == this->num_cell_ - 1)
    {
        return values_.size() - coefficieint_start_indexes_[cell_index];
    }
    else
    {
        return coefficieint_start_indexes_[cell_index + 1] - coefficieint_start_indexes_[cell_index];
    }
}