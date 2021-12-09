#pragma once
#include "Euclidean_Vector.h"

using ushort = unsigned short;

namespace ms
{
	inline constexpr ushort blas_mv_criteria = 50;
	inline constexpr ushort blas_mm_criteria = 25;
}

class Matrix;
class Matrix_Base
{
public: //Command
	void transpose(void);

public: //Query
	Matrix operator*(const double constant) const;	
	Euclidean_Vector operator*(const Euclidean_Vector& vec) const;
	Matrix operator*(const Matrix_Base& other) const;
	Matrix operator+(const Matrix_Base& other) const;

	double at(const size_t row, const size_t column) const;
	Euclidean_Vector column(const size_t column_index) const;
	void column(const size_t column_index, double* value_ptr) const;
	const double* data(void) const;
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

class Matrix : public Matrix_Base
{
	friend class Matrix_Base;

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
	void operator*=(const double constant);

	double* data(void);
	Matrix& inverse(void);
	template <typename V>	void change_column(const size_t column_index, const V& vec) 
	{
		REQUIRE(column_index < this->num_columns_, "column idnex can not exceed number of column");

		for (size_t i = 0; i < this->num_rows_; ++i)
			this->value_at(i, column_index) = vec.at(i);
	}
	template <typename V>	void change_row(const size_t start_row_index, const V& vec) 
	{
		REQUIRE(start_row_index <= this->num_rows_, "index can not exceed given range");
		REQUIRE(!this->is_transposed(), "it should be not transposed for this routine");

		const auto jump_index = start_row_index * this->num_columns_;
		std::copy(vec.begin(), vec.end(), this->values_.data() + jump_index);

	}
	void change_rows(const size_t start_row_index, const Matrix& A);
	void change_columns(const size_t start_column_index, const Matrix& A);	

public://Query
	bool operator==(const Matrix& other) const;
	
	const double* data(void) const;
	Matrix get_transpose(void) const;
	Matrix get_inverse(void) const;
	double& value_at(const size_t row, const size_t column); // optimize
	double* value_ptr(void);

private:
	std::vector<int> PLU_decomposition(void);

private:
	std::vector<double> values_;
};

class Matrix_Constant_Wrapper : public Matrix_Base
{
public:
	Matrix_Constant_Wrapper(const size_t num_row, const size_t num_column, const double* ptr);
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
	void gemm(const Matrix_Base& A, const Matrix_Base& B, double* output_ptr);
	void mpm(const Matrix_Base& M1, const Matrix_Base& M2, double* result);
	void mv(const Matrix_Base& M, const Euclidean_Vector& v, double* result);
}

Matrix operator*(const double constant, const Matrix_Base& M);
std::ostream& operator<<(std::ostream& os, const Matrix_Base& m);
std::ostream& operator<<(std::ostream& os, const Matrix& m);