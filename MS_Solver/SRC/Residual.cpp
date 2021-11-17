#include "../INC/Residual.h"

Residual::Residual(const size_t num_values, const std::vector<size_t>& coefficient_start_indexes)
    : coefficieint_start_indexes_(coefficient_start_indexes)
{
    this->values_.resize(num_values);
}

void Residual::update_rhs(const uint cell_index, const Matrix& delta_rhs)
{
    const auto n = static_cast<MKL_INT>(delta_rhs.num_values());
    const auto a = 1.0;
    const auto incx = 1;
    const auto incy = 1;

    auto icell_rhs_ptr = this->get_icell_residual(cell_index);
    cblas_daxpy(n, a, delta_rhs.data(), incx, icell_rhs_ptr, incy);
}

std::vector<double>&& Residual::move_values(void)
{
    return std::move(this->values_);
}


double* Residual::get_icell_residual(const uint cell_index)
{
    return this->values_.data() + coefficieint_start_indexes_[cell_index];
}