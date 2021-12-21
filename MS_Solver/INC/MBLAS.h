#pragma once
#include <cmath>
#include <mkl.h>

using ushort = unsigned short;

namespace ms::BLAS
{
	inline constexpr ushort dcopy_criteria = 50;
	inline constexpr ushort dscal_criteria = 10;
	inline constexpr ushort dasum_criteria = 10;
	inline constexpr ushort axpy_criteria = 20;
	inline constexpr ushort dot_criteria = 15;

	double abs_x(const int n, const double* x_ptr);
	void copy(const int n, const double* x_ptr, double* result_ptr);
	void cx(const double c, const int n, double* x_ptr);
	double x_dot_y(const int n, const double* x_ptr, const double* y_ptr);
	void x_plus_y(const int n, const double* x_ptr, const double* y_ptr, double* result_ptr);
	void x_plus_assign_y(const int n, double* x_ptr, const double* y_ptr);
	void x_plus_assign_cy(const int n, double* x_ptr, const double c, const double* y_ptr);
	void x_minus_y(const int n, const double* x_ptr, const double* y_ptr, double* result_ptr);
	void x_minus_assign_y(const int n, double* x_ptr, const double* y_ptr);
}