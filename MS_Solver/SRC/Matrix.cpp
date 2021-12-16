#include "../INC/Matrix.h"

Constant_Matrix_Wrapper::Constant_Matrix_Wrapper(const size_t num_row, const size_t num_column, const double* ptr)
	: num_rows_(num_row)
	, num_columns_(num_column)
	, const_data_ptr_(ptr) 
{
	REQUIRE(num_row * num_column != 0, "number of row or number of column can not be 0");
};


Matrix Constant_Matrix_Wrapper::operator*(const double constant) const
{
	Matrix result(this->num_rows_, this->num_columns_, this->const_data_ptr_);
	result *= constant;
	return result;
}

Matrix Constant_Matrix_Wrapper::operator+(const Constant_Matrix_Wrapper& other) const
{
	REQUIRE(this->size() == other.size(), "two matrix should be same size");
	REQUIRE(!this->is_transposed() && !other.is_transposed(), "both matrixes should not be transposed");

	Matrix result(this->num_rows_, this->num_columns_);
	const auto n = static_cast<MKL_INT>(this->num_values());	
	ms::BLAS::x_plus_y(n, this->const_data_ptr_, other.const_data_ptr_, result.data());

	return result;
}

Matrix Constant_Matrix_Wrapper::operator*(const Constant_Matrix_Wrapper& other) const 
{	
	Matrix result(this->num_rows_, other.num_columns_);
	ms::gemm(*this, other, result.data());
	
	return result;
}

double Constant_Matrix_Wrapper::at(const size_t row, const size_t column) const 
{
	REQUIRE(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
	{
		return this->const_data_ptr_[column * this->num_rows_ + row];
	}
	else
	{
		return this->const_data_ptr_[row * this->num_columns_ + column];
	}
}

std::vector<double> Constant_Matrix_Wrapper::column(const size_t column_index) const
{
	REQUIRE(column_index < this->num_columns_, "index can not exceed given range");

	std::vector<double> column_values(this->num_rows_);

	for (size_t i = 0; i < this->num_rows_; ++i)
	{
		column_values[i] = this->at(i, column_index);
	}

	return column_values;
}

void Constant_Matrix_Wrapper::column(const size_t column_index, double* value_ptr) const
{
	REQUIRE(column_index < this->num_columns_, "index can not exceed given range");

	for (size_t i = 0; i < this->num_rows_; ++i)
	{
		value_ptr[i] = this->at(i, column_index);
	}
}

const double* Constant_Matrix_Wrapper::data(void) const
{
	return this->const_data_ptr_;
}

Matrix Constant_Matrix_Wrapper::get_transpose(void) const
{
	Matrix result(this->num_rows_, this->num_columns_, this->const_data_ptr_);
	result.transpose();
	return result;
}

Matrix Constant_Matrix_Wrapper::get_inverse(void) const
{
	Matrix result(this->num_rows_, this->num_columns_, this->const_data_ptr_);
	result.inverse();
	return result;
}

bool Constant_Matrix_Wrapper::is_finite(void) const 
{
	for (size_t i = 0; i < this->num_values(); ++i) 
	{
		if (!std::isfinite(this->const_data_ptr_[i]))
			return false;
	}
	return true;
}

size_t Constant_Matrix_Wrapper::leading_dimension(void) const
{
	// num column before OP()
	if (this->is_transposed())
		return this->num_rows_;
	else
		return this->num_columns_;
}

size_t Constant_Matrix_Wrapper::num_column(void) const 
{
	return this->num_columns_;
}

size_t Constant_Matrix_Wrapper::num_row(void) const
{
	return this->num_rows_;
}

size_t Constant_Matrix_Wrapper::num_values(void) const
{
	return this->num_rows_ * this->num_columns_;
}

std::vector<double> Constant_Matrix_Wrapper::row(const size_t row_index) const
{
	REQUIRE(row_index < this->num_rows_, "index can not exceed given range");

	std::vector<double> row_values(this->num_columns_);

	for (size_t i = 0; i < this->num_columns_; ++i)
		row_values[i] = this->at(row_index, i);

	return row_values;
}

std::pair<size_t, size_t> Constant_Matrix_Wrapper::size(void) const
{
	return { this->num_rows_, this->num_columns_ };
}

std::string Constant_Matrix_Wrapper::to_string(void) const 
{
	std::ostringstream oss;
	oss << std::setprecision(16) << std::showpoint << std::left;
	for (size_t i = 0; i < this->num_rows_; ++i) 
	{
		for (size_t j = 0; j < this->num_columns_; ++j)
			oss << std::setw(25) << this->at(i, j);
		oss << "\n";
	}
	return oss.str();
}

CBLAS_TRANSPOSE Constant_Matrix_Wrapper::transpose_type(void) const
{
	return this->transpose_type_;
}

bool Constant_Matrix_Wrapper::is_transposed(void) const 
{
	return this->transpose_type_ == CBLAS_TRANSPOSE::CblasTrans;
}

bool Constant_Matrix_Wrapper::is_square_matrix(void) const 
{
	return this->num_rows_ == this->num_columns_;
}

bool Constant_Matrix_Wrapper::is_in_range(const size_t irow, const size_t jcolumn) const 
{
	return irow < this->num_rows_&& jcolumn < this->num_columns_;
}

void Matrix_Wrapper::operator*=(const double constant)
{
	ms::BLAS::cx(constant, static_cast<int>(this->num_values()), this->data_ptr_);
}

double& Matrix_Wrapper::at(const size_t row, const size_t column)
{
	REQUIRE(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
	{
		return this->data_ptr_[column * this->num_rows_ + row];
	}
	else
	{
		return this->data_ptr_[row * this->num_columns_ + column];
	}
}

void Matrix_Wrapper::change_columns(const size_t start_column_index, const Constant_Matrix_Wrapper& A)
{
	const auto [num_A_rows, num_A_columns] = A.size();

	REQUIRE(this->num_rows_ == num_A_rows, "size should be matched");
	REQUIRE(start_column_index + num_A_columns <= this->num_columns_, "index can not exceed given range");

	for (size_t i = 0; i < num_A_rows; ++i)
	{
		for (size_t j = 0; j < num_A_columns; ++j)
		{
			this->at(i, start_column_index + j) = A.at(i, j);
		}
	}
}

void Matrix_Wrapper::change_columns(const size_t start_column_index, const size_t end_column_index, const double value)
{
	//[start, end)
	REQUIRE(start_column_index <= end_column_index, "start index can not exceed end index");
	REQUIRE(end_column_index <= this->num_columns_, "index can not exceed given range");

	for (size_t i = 0; i < this->num_rows_; ++i)
	{
		for (size_t j = start_column_index; j < end_column_index; ++j)
		{
			this->at(i, j) = value;
		}
	}
}

void Matrix_Wrapper::change_rows(const size_t start_row_index, const Constant_Matrix_Wrapper& A)
{
	const auto [num_A_rows, num_A_columns] = A.size();
	const auto A_data_ptr = A.data();

	REQUIRE(this->num_columns_ == num_A_columns, "size should be matched");
	REQUIRE(start_row_index + num_A_rows <= this->num_rows_, "index can not exceed given range");
	REQUIRE(!this->is_transposed() && !A.is_transposed(), "it should be not transposed for this routine");

	const auto jump_index = start_row_index * this->num_columns_;
	std::copy(A_data_ptr, A_data_ptr + num_A_rows * num_A_columns, this->data_ptr_ + jump_index);
}

double* Matrix_Wrapper::data(void)
{
	return this->data_ptr_;
}

void Matrix_Wrapper::inverse(void)
{
	REQUIRE(this->is_square_matrix(), "invertable matrix should be square matrix");

	const auto ipiv = this->PLU_decomposition();

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int n = static_cast<int>(this->num_rows_);
	const lapack_int lda = n;
	const lapack_int info = LAPACKE_dgetri(matrix_layout, n, this->data_ptr_, lda, ipiv.data());

	REQUIRE(info == 0, "info should be 0 when success matrix get_inverse");
	//info > 0 "U is singular matrix in L-U decomposition"
	//info < 0 "fail to inverse the matrix"
}

void Matrix_Wrapper::scalar_multiplcation_at_columns(const size_t start_column_index, const size_t end_column_index, const double scalar)
{
	//[start, end)
	REQUIRE(start_column_index <= end_column_index, "start index can not exceed end index");
	REQUIRE(end_column_index <= this->num_columns_, "index can not exceed given range");

	for (size_t i = 0; i < this->num_rows_; ++i)
	{
		for (size_t j = start_column_index; j < end_column_index; ++j)
		{
			this->at(i, j) *= scalar;
		}
	}
}

void Matrix_Wrapper::transpose(void)
{
	std::swap(this->num_rows_, this->num_columns_);

	if (this->is_transposed())
	{
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	}
	else
	{
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasTrans;
	}
}

double Matrix_Wrapper::at(const size_t row, const size_t column) const
{
	REQUIRE(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
	{
		return this->const_data_ptr_[column * this->num_rows_ + row];
	}
	else
	{
		return this->const_data_ptr_[row * this->num_columns_ + column];
	}
}

const double* Matrix_Wrapper::data(void) const
{
	return this->const_data_ptr_;
}

std::vector<int> Matrix_Wrapper::PLU_decomposition(void)
{
	REQUIRE(!this->is_transposed(), "matrix should not be transposed");

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int m = static_cast<int>(this->num_rows_);
	const lapack_int n = static_cast<int>(this->num_columns_);
	const lapack_int lda = n;
	std::vector<int> ipiv(std::min(m, n));
	lapack_int info = LAPACKE_dgetrf(matrix_layout, m, n, this->data_ptr_, lda, ipiv.data());

	REQUIRE(0 <= info, "info should be greater than 0 when sucess PLU decomposition");
	return ipiv;
}


Matrix::Matrix(const size_t matrix_order)
	: Matrix_Wrapper(matrix_order, matrix_order, this->values_.data())
	, values_(matrix_order* matrix_order)
{
	//prevent dangling pointer caused by vector reallocation
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();

	for (size_t i = 0; i < matrix_order; ++i)
		this->at(i, i) = 1.0;
}

Matrix::Matrix(const size_t matrix_order, const std::vector<double>& diagonal_values)
	: Matrix_Wrapper(matrix_order, matrix_order, this->values_.data())
	, values_(matrix_order * matrix_order)
{
	REQUIRE(matrix_order == diagonal_values.size(), "number of diagonal value of square matrix should be same with matrix order");

	//prevent dangling pointer caused by vector reallocation
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();

	for (size_t i = 0; i < matrix_order; ++i)
		this->at(i, i) = diagonal_values[i];
}

Matrix::Matrix(const size_t num_row, const size_t num_column)
	: Matrix_Wrapper(num_row, num_column, this->values_.data())
	, values_(num_row * num_column)
{
	//prevent dangling pointer caused by vector reallocation
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Matrix::Matrix(const size_t num_row, const size_t num_column, const double* value_ptr)
	: Matrix_Wrapper(num_row, num_column, this->values_.data())
	, values_(value_ptr, value_ptr + num_row * num_column)
{
	//prevent dangling pointer caused by vector reallocation
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Matrix::Matrix(const size_t num_row, const size_t num_column, std::vector<double>&& values)
	: Matrix_Wrapper(num_row, num_column, this->values_.data())
	, values_(std::move(values))
{
	REQUIRE(num_row * num_column == this->values_.size(), "num value should be same with matrix size");

	//prevent dangling pointer caused by vector reallocation
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Matrix::Matrix(const Matrix& other) 
	: Matrix_Wrapper(other.num_rows_, other.num_columns_, this->values_.data())
	, values_(other.values_)
{
	this->transpose_type_ = other.transpose_type_;

	//prevent dangling pointer caused by vector reallocation
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Matrix::Matrix(Matrix&& other) noexcept
	: Matrix_Wrapper(other.num_rows_, other.num_columns_, this->values_.data())
	, values_(std::move(other.values_))
{
	this->transpose_type_ = other.transpose_type_;

	//prevent dangling pointer caused by vector reallocation
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
	
	other.const_data_ptr_ = nullptr;
	other.data_ptr_ = nullptr;
}

void Matrix::operator=(const Matrix& other) 
{
	this->transpose_type_ = other.transpose_type_;
	this->num_rows_ = other.num_rows_;
	this->num_columns_ = other.num_columns_;
	this->values_ = other.values_;
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

void Matrix::operator=(Matrix&& other) noexcept
{
	this->transpose_type_ = other.transpose_type_;
	this->num_rows_ = other.num_rows_;
	this->num_columns_ = other.num_columns_;
	this->values_ = std::move(other.values_);
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();

	other.const_data_ptr_ = nullptr;
	other.data_ptr_ = nullptr;

}



bool Matrix::operator==(const Matrix& other) const 
{
	if (this->num_rows_ != other.num_rows_)
		return false;

	if (this->num_columns_ != other.num_columns_)
		return false;

	if (this->transpose_type_ == other.transpose_type_)
	{
		return std::equal(this->const_data_ptr_, this->const_data_ptr_ + this->num_values(), other.const_data_ptr_);
	}
	else
	{
		for (size_t i = 0; i < this->num_rows_; ++i) 
		{
			for (size_t j = 0; j < this->num_columns_; ++j)
			{
				if (this->at(i,j) != other.at(i,j))
					return false;
			}
		}
		return true;
	}
}

//Constant_Matrix_Wrapper::Constant_Matrix_Wrapper(const size_t num_row, const size_t num_column, const double* ptr)
//{
//	this->num_rows_ = num_row;
//	this->num_columns_ = num_column;
//	this->const_data_ptr_ = ptr;
//}

namespace ms
{
	void gemm(const Constant_Matrix_Wrapper& A, const Constant_Matrix_Wrapper& B, double* output_ptr)
	{
		const auto m = static_cast<MKL_INT>(A.num_row());
		const auto n = static_cast<MKL_INT>(B.num_column());
		const auto k = static_cast<MKL_INT>(A.num_column());

		REQUIRE(A.num_column() == B.num_row(), "size should be matched for matrix multiplication");

		if (std::max(A.num_values(), B.num_values()) <= ms::blas_mm_criteria)
		{
			for (int i = 0; i < m; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					for (int l = 0; l < k; ++l)
					{
						output_ptr[i * n + j] += A.at(i, l) * B.at(l, j);
					}
				}
			}
		}
		else
		{
			const auto layout = CBLAS_LAYOUT::CblasRowMajor;
			const auto transA = A.transpose_type();
			const auto transB = B.transpose_type();
			const auto alpha = 1.0;
			const auto lda = static_cast<MKL_INT>(A.leading_dimension());
			const auto ldb = static_cast<MKL_INT>(B.leading_dimension());
			const auto beta = 0.0;
			const auto ldc = n;

			cblas_dgemm(layout, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, output_ptr, ldc);
		}
	}

	void mpm(const Constant_Matrix_Wrapper& M1, const Constant_Matrix_Wrapper& M2, double* result_ptr)
	{
		REQUIRE(M1.size() == M2.size(), "two matrix should be same size");
		REQUIRE(!M1.is_transposed() && !M2.is_transposed(), "both matrixes should not be transposed");
		
		const auto n = static_cast<int>(M1.num_values());
		ms::BLAS::x_plus_y(n, M1.data(), M2.data(), result_ptr);
	}
}

Matrix operator*(const double constant, const Constant_Matrix_Wrapper& M)
{
	return M * constant;
}

std::ostream& operator<<(std::ostream& os, const Constant_Matrix_Wrapper& m) 
{
	return os << m.to_string();
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) 
{
	return os << m.to_string();
}