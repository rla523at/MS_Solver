#include "../INC/MBLAS.h"

namespace ms::BLAS
{
	void abs_x(const int n, double* x_ptr)
	{
		for (int i = 0; i < n; ++i)
		{
			x_ptr[i] = std::abs(x_ptr[i]);
		}
	}

	double abs_sum_x(const int n, const double* x_ptr)
	{
		double result = 0.0;

		if (n <= ms::BLAS::dasum_criteria)
		{
			for (size_t i = 0; i < n; ++i)
			{
				result += std::abs(x_ptr[i]);
			}
		}
		else
		{
			const auto incx = 1;

			result = cblas_dasum(n, x_ptr, incx);
		}

		return result;
	}

	void copy(const int n, const double* x_ptr, double* result_ptr)
	{
		if (n <= ms::BLAS::dcopy_criteria)
		{
			for (int i = 0; i < n; ++i)
			{
				result_ptr[i] = x_ptr[i];
			}
		}
		else
		{
			const auto incx = 1;
			const auto incy = 1;
			cblas_dcopy(n, x_ptr, incx, result_ptr, incy);
		}
	}

	void cx(const double c, const int n, double* x_ptr)
	{
		if (n <= ms::BLAS::dscal_criteria)
		{
			for (int i = 0; i < n; ++i)
			{
				x_ptr[i] *= c;
			}
		}
		else
		{
			const auto incx = 1;
			cblas_dscal(n, c, x_ptr, incx);
		}

	}

	double x_dot_y(const int n, const double* x_ptr, const double* y_ptr)
	{
		double result = 0.0;

		if (n <= ms::BLAS::dot_criteria)
		{
			for (int i = 0; i < n; ++i)
			{
				result += x_ptr[i] * y_ptr[i];
			}
		}
		else
		{
			const auto incx = 1;
			const auto incy = 1;

			result = cblas_ddot(n, x_ptr, incx, y_ptr, incy);
		}

		return result;
	}

	void x_plus_y(const int n, const double* x_ptr, const double* y_ptr, double* result_ptr)
	{
		ms::BLAS::copy(n, x_ptr, result_ptr);
		ms::BLAS::x_plus_assign_y(n, result_ptr, y_ptr);
	}

	void x_plus_assign_y(const int n, double* x_ptr, const double* y_ptr)
	{
		if (n <= ms::BLAS::axpy_criteria)
		{
			for (int i = 0; i < n; ++i)
			{
				x_ptr[i] += y_ptr[i];
			}
		}
		else
		{
			const auto a = 1.0;
			const auto incx = 1;
			const auto incy = 1;
			cblas_daxpy(n, a, y_ptr, incx, x_ptr, incy);
		}
	}

	void x_plus_assign_cy(const int n, double* x_ptr, const double c, const double* y_ptr)
	{
		if (n <= ms::BLAS::axpy_criteria)
		{
			for (int i = 0; i < n; ++i)
			{
				x_ptr[i] += c * y_ptr[i];
			}
		}
		else
		{			
			const auto incx = 1;
			const auto incy = 1;
			cblas_daxpy(n, c, y_ptr, incx, x_ptr, incy);
		}
	}

	void x_minus_y(const int n, const double* x_ptr, const double* y_ptr, double* result_ptr)
	{
		ms::BLAS::copy(n, x_ptr, result_ptr);
		ms::BLAS::x_minus_assign_y(n, result_ptr, y_ptr);
	}

	void x_minus_assign_y(const int n, double* x_ptr, const double* y_ptr)
	{
		if (n <= ms::BLAS::axpy_criteria)
		{
			for (int i = 0; i < n; ++i)
			{
				x_ptr[i] -= y_ptr[i];
			}
		}
		else
		{
			const auto a = -1.0;
			const auto incx = 1;
			const auto incy = 1;
			cblas_daxpy(n, a, y_ptr, incx, x_ptr, incy);
		}
	}

}