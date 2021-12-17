#pragma once
#include "Exception.h"
#include "MBLAS.h"

#include <iomanip>
#include <sstream>
#include <type_traits>
#include <vector>

namespace ms
{
	inline constexpr ushort blas_mv_criteria = 50;
	inline constexpr ushort blas_mm_criteria = 25;
}

class Matrix;
class Constant_Matrix_Wrapper
{
public:
	Constant_Matrix_Wrapper(void) = default;
	Constant_Matrix_Wrapper(const size_t num_row, const size_t num_column, const double* ptr);

public: //Query
	Matrix operator*(const double constant) const;	
	template <typename V, std::enable_if_t<!std::is_base_of_v<Constant_Matrix_Wrapper,V> && std::is_class_v<V>, bool> = true>  V operator*(const V& vec) const
	{
		REQUIRE(this->num_columns_ == vec.size(), "size should be mathced");
		std::vector<double> result(this->num_rows_);

		if (this->num_values() <= ms::blas_mv_criteria)
		{
			for (size_t i = 0; i < this->num_rows_; ++i)
			{
				for (size_t j = 0; j < this->num_columns_; ++j)
				{
					result[i] += this->at(i, j) * vec.at(j);
				}
			}
		}
		else
		{
			const auto Layout = CBLAS_LAYOUT::CblasRowMajor;
			const auto trans = this->transpose_type_;
			const auto m = static_cast<MKL_INT>(this->num_rows_);
			const auto n = static_cast<MKL_INT>(this->num_columns_);
			const auto alpha = 1.0;
			const auto lda = static_cast<MKL_INT>(this->leading_dimension());
			const auto incx = 1;
			const auto beta = 0;
			const auto incy = 1;

			cblas_dgemv(Layout, trans, m, n, alpha, this->const_data_ptr_, lda, vec.data(), incx, beta, result.data(), incy);
		}

		return result;
	}
	Matrix operator*(const Constant_Matrix_Wrapper& other) const;
	Matrix operator+(const Constant_Matrix_Wrapper& other) const;

	double at(const size_t row, const size_t column) const;
	std::vector<double> column(const size_t column_index) const;
	void column(const size_t column_index, double* value_ptr) const;
	const double* data(void) const;
	Matrix get_transpose(void) const;
	Matrix get_inverse(void) const;
	bool is_finite(void) const;
	size_t leading_dimension(void) const;
	size_t num_column(void) const;
	size_t num_row(void) const;
	size_t num_values(void) const;
	std::vector<double> row(const size_t row_index) const;
	std::pair<size_t, size_t> size(void) const;
	std::string to_string(void) const;
	CBLAS_TRANSPOSE transpose_type(void) const;
	bool is_transposed(void) const;

protected:
	bool is_square_matrix(void) const;
	bool is_in_range(const size_t irow, const size_t jcolumn) const;	

protected:
	CBLAS_TRANSPOSE transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	size_t num_rows_ = 0;
	size_t num_columns_ = 0;
	const double* const_data_ptr_ = nullptr;
};

class Matrix_Wrapper : public Constant_Matrix_Wrapper
{
public:
	Matrix_Wrapper(void) = default;
	Matrix_Wrapper(const size_t num_row, const size_t num_column, double* ptr)
		: Constant_Matrix_Wrapper(num_row, num_column, ptr)
		, data_ptr_(ptr) {};

public://Command
	void operator*=(const double constant);
	
	double& at(const size_t row, const size_t column);
	template <typename V>	void change_column(const size_t column_index, const V& vec)
	{
		REQUIRE(column_index < this->num_columns_, "column idnex can not exceed number of column");

		for (size_t i = 0; i < this->num_rows_; ++i)
			this->at(i, column_index) = vec.at(i);
	}
	template <typename V>	void change_row(const size_t start_row_index, const V& vec)
	{
		REQUIRE(start_row_index <= this->num_rows_, "index can not exceed given range");
		REQUIRE(!this->is_transposed(), "it should be not transposed for this routine");

		const auto jump_index = start_row_index * this->num_columns_;
		std::copy(vec.begin(), vec.end(), this->data_ptr_ + jump_index);

	}
	void change_columns(const size_t start_column_index, const Constant_Matrix_Wrapper& A);
	void change_columns(const size_t start_column_index, const size_t end_column_index, const double value);
	void change_rows(const size_t start_row_index, const Constant_Matrix_Wrapper& A);
	void change_rows(const size_t start_row_index, const size_t end_row_index, const double value);
	double* data(void);
	void inverse(void);
	void scalar_multiplcation_at_columns(const size_t start_column_index, const size_t end_column_index, const double scalar);
	void transpose(void);

public://Query
	double at(const size_t row, const size_t column) const;
	const double* data(void) const;

private:
	std::vector<int> PLU_decomposition(void);

protected:
	double* data_ptr_ = nullptr;
};

class Matrix : public Matrix_Wrapper
{
public:
	Matrix(void) = default;
	Matrix(const size_t matrix_order);
	Matrix(const size_t matrix_order, const std::vector<double>& value);
	Matrix(const size_t num_row, const size_t num_column);
	Matrix(const size_t num_row, const size_t num_column, const double* value_ptr);
	Matrix(const size_t num_row, const size_t num_column, std::vector<double>&& values);
	Matrix(const Matrix& other);
	Matrix(Matrix&& other) noexcept;

public://Command 
	void operator=(const Matrix& other);
	void operator=(Matrix&& other) noexcept;

public://Query
	bool operator==(const Matrix& other) const;

private:
	std::vector<double> values_;
};

template <typename Function>
class Matrix_Function
{
public:
	Matrix_Function(const size_t num_row, const size_t num_column)
		: num_row_(num_row),
		num_column_(num_column)
	{
		this->functions_.resize(this->num_row_ * this->num_column_);
	}
	Matrix_Function(const size_t num_row, const size_t num_column, const std::vector<Function>& functions)
		: num_row_(num_row),
		num_column_(num_column),
		functions_(functions) {};

public://Query
	template <typename V>	Matrix operator()(const V& space_vector) const
	{
		std::vector<double> values(this->num_row_ * this->num_column_);

		for (ushort i = 0; i < this->num_row_; ++i)
			for (ushort j = 0; j < this->num_column_; ++j)
				values[i * this->num_column_ + j] = this->at(i, j)(space_vector);

		return { this->num_row_,this->num_column_, std::move(values) };
	}
	bool operator==(const Matrix_Function& other) const
	{
		return this->functions_ == other.functions_;
	}
	const Function& at(const ushort row_index, const ushort column_index) const
	{
		REQUIRE(row_index < this->num_row_&& column_index < this->num_column_, "index can not exceed given range");
		return this->functions_[row_index * this->num_column_ + column_index];
	}
	template <typename VF> void change_column(const ushort column_index, const VF& vector_function)
	{
		REQUIRE(this->num_row_ == vector_function.size(), "column vector should have num row range dimension");

		for (ushort i = 0; i < this->num_row_; ++i)
			this->function_at(i, column_index) = vector_function[i];
	}
	std::string to_string(void) const
	{
		std::string str;
		for (ushort i = 0; i < this->num_row_; ++i)
		{
			for (ushort j = 0; j < this->num_column_; ++j)
			{
				str += this->at(i, j).to_string() + "\t";
			}
			str += "\n";
		}
		return str;
	}

private:
	Function& function_at(const ushort row_index, const ushort column_index)
	{
		REQUIRE(row_index < this->num_row_&& column_index < this->num_column_, "index can not exceed given range");
		return this->functions_[row_index * this->num_column_ + column_index];
	}

private:
	size_t num_row_ = 0;
	size_t num_column_ = 0;
	std::vector<Function> functions_;
};

namespace ms
{
	void gemm(const Constant_Matrix_Wrapper& A, const Constant_Matrix_Wrapper& B, double* output_ptr);
	void mpm(const Constant_Matrix_Wrapper& M1, const Constant_Matrix_Wrapper& M2, double* result);
	template <typename V> void mv(const Constant_Matrix_Wrapper& M, const V& v, double* result_ptr)
	{
		const auto m = static_cast<int>(M.num_row());
		const auto n = static_cast<int>(M.num_column());

		REQUIRE(n == v.size(), "size should be mathced");

		if (m * n <= ms::blas_mv_criteria)
		{
			for (size_t i = 0; i < m; ++i)
			{
				for (size_t j = 0; j < n; ++j)
				{
					result_ptr[i] += M.at(i, j) * v.at(j);
				}
			}
		}
		else
		{
			const auto Layout = CBLAS_LAYOUT::CblasRowMajor;
			const auto trans = M.transpose_type();
			const auto alpha = 1.0;
			const auto lda = static_cast<MKL_INT>(M.leading_dimension());
			const auto incx = 1;
			const auto beta = 0;
			const auto incy = 1;

			cblas_dgemv(Layout, trans, m, n, alpha, M.data(), lda, v.data(), incx, beta, result_ptr, incy);
		}
	}
}

Matrix operator*(const double constant, const Constant_Matrix_Wrapper& M);
std::ostream& operator<<(std::ostream& os, const Constant_Matrix_Wrapper& m);
std::ostream& operator<<(std::ostream& os, const Matrix& m);