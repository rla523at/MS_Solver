#include "../INC/Matrix.h"

void Matrix_Base::be_transpose(void) 
{
	std::swap(this->num_rows_, this->num_columns_);

	if (this->is_transposed())
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	else
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasTrans;
}

Matrix Matrix_Base::operator*(const double constant) const
{
	const auto num_values = this->num_values();
	std::vector<double> values(num_values);

	const auto n = static_cast<int>(num_values);
	const auto incx = 1;
	cblas_dscal(n, constant, values.data(), incx);

	return { this->num_rows_, this->num_columns_, std::move(values) };
}

Matrix Matrix_Base::operator+(const Matrix_Base& other) const
{
	REQUIRE(this->size() == other.size(), "two matrix should be same size");

	const auto n = static_cast<MKL_INT>(this->num_values());
	const auto a = 1.0;
	const auto incx = 1;
	const auto incy = 1;

	std::vector<double> values(n);
	cblas_daxpy(n, a, other.const_data_ptr_, incx, values.data(), incy);

	return { this->num_rows_,this->num_columns_,std::move(values) };
}


Matrix Matrix_Base::operator*(const Matrix_Base& other) const 
{
	const auto [num_row, num_column] = other.size();

	std::vector<double> values(this->num_rows_ * num_column);
	ms::gemm(*this, other, values.data());
	return { this->num_rows_, num_column, std::move(values) };
}

double Matrix_Base::at(const size_t row, const size_t column) const 
{
	REQUIRE(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
		return this->const_data_ptr_[column * this->num_rows_ + row];
	else
		return this->const_data_ptr_[row * this->num_columns_ + column];
}

std::vector<double> Matrix_Base::column(const size_t column_index) const
{
	REQUIRE(column_index < this->num_columns_, "index can not exceed given range");

	std::vector<double> column_values(this->num_rows_);

	for (size_t i = 0; i < this->num_rows_; ++i)
		column_values[i] = this->at(i, column_index);

	return column_values;
}

const double* Matrix_Base::data(void) const
{
	return this->const_data_ptr_;
}

bool Matrix_Base::is_finite(void) const 
{
	for (size_t i = 0; i < this->num_values(); ++i) 
	{
		if (!std::isfinite(this->const_data_ptr_[i]))
			return false;
	}
	return true;
}

size_t Matrix_Base::leading_dimension(void) const
{
	// num column before OP()
	if (this->is_transposed())
		return this->num_rows_;
	else
		return this->num_columns_;
}

size_t Matrix_Base::num_column(void) const 
{
	return this->num_columns_;
}

size_t Matrix_Base::num_row(void) const
{
	return this->num_rows_;
}

size_t Matrix_Base::num_values(void) const
{
	return this->num_rows_ * this->num_columns_;
}

std::vector<double> Matrix_Base::row(const size_t row_index) const
{
	REQUIRE(row_index < this->num_rows_, "index can not exceed given range");

	std::vector<double> row_values(this->num_columns_);

	for (size_t i = 0; i < this->num_columns_; ++i)
		row_values[i] = this->at(row_index, i);

	return row_values;
}

std::pair<size_t, size_t> Matrix_Base::size(void) const
{
	return { this->num_rows_, this->num_columns_ };
}

std::string Matrix_Base::to_string(void) const 
{
	std::ostringstream oss;
	oss << std::setprecision(16) << std::showpoint << std::left;
	for (size_t i = 0; i < this->num_rows_; ++i) {
		for (size_t j = 0; j < this->num_columns_; ++j)
			oss << std::setw(25) << this->at(i, j);
		oss << "\n";
	}
	return oss.str();
}

CBLAS_TRANSPOSE Matrix_Base::transpose_type(void) const
{
	return this->transpose_type_;
}

bool Matrix_Base::is_transposed(void) const 
{
	return this->transpose_type_ == CBLAS_TRANSPOSE::CblasTrans;
}

bool Matrix_Base::is_square_matrix(void) const 
{
	return this->num_rows_ == this->num_columns_;
}

bool Matrix_Base::is_in_range(const size_t irow, const size_t jcolumn) const 
{
	return irow < this->num_rows_&& jcolumn < this->num_columns_;
}



Matrix::Matrix(const size_t matrix_order) {
	REQUIRE(matrix_order != 0, "matrix order can not be 0");

	this->num_rows_ = matrix_order;
	this->num_columns_ = matrix_order;
	this->values_.resize(this->num_values());
	this->const_data_ptr_ = this->values_.data();

	for (size_t i = 0; i < matrix_order; ++i)
		this->value_at(i, i) = 1.0;
}

Matrix::Matrix(const size_t matrix_order, const std::vector<double>& value) {
	REQUIRE(matrix_order != 0, "matrix order can not be 0");
	REQUIRE(matrix_order == value.size(), "num value of square matrix should be same with matrix order");

	this->num_rows_ = matrix_order;
	this->num_columns_ = matrix_order;
	this->values_.resize(this->num_values());
	this->const_data_ptr_ = this->values_.data();

	for (size_t i = 0; i < matrix_order; ++i)
		this->value_at(i, i) = value[i];
}

Matrix::Matrix(const size_t num_row, const size_t num_column) {
	REQUIRE(num_row * num_column != 0, "number of row or number of column can not be 0");

	this->num_rows_ = num_row;
	this->num_columns_ = num_column;
	this->values_.resize(this->num_values());
	this->const_data_ptr_ = this->values_.data();
}

Matrix::Matrix(const size_t num_row, const size_t num_column, std::vector<double>&& value) {
	REQUIRE(num_row * num_column != 0, "number of row or number of column can not be 0");
	REQUIRE(num_row * num_column == value.size(), "num value should be same with matrix size");

	this->num_rows_ = num_row;
	this->num_columns_ = num_column;
	this->values_ = std::move(value);
	this->const_data_ptr_ = this->values_.data();
}

Matrix::Matrix(const Matrix& other) 
{
	this->transpose_type_ = other.transpose_type_;
	this->num_rows_ = other.num_rows_;
	this->num_columns_ = other.num_columns_;
	this->values_ = other.values_;
	this->const_data_ptr_ = this->values_.data();
}

Matrix::Matrix(Matrix&& other) noexcept
{
	this->transpose_type_ = other.transpose_type_;
	this->num_rows_ = other.num_rows_;
	this->num_columns_ = other.num_columns_;
	this->values_ = std::move(other.values_);
	this->const_data_ptr_ = this->values_.data();
	other.const_data_ptr_ = nullptr;
}

void Matrix::operator=(const Matrix& other) 
{
	this->transpose_type_ = other.transpose_type_;
	this->num_rows_ = other.num_rows_;
	this->num_columns_ = other.num_columns_;
	this->values_ = other.values_;
	this->const_data_ptr_ = this->values_.data();
}

void Matrix::operator=(Matrix&& other) noexcept
{
	this->transpose_type_ = other.transpose_type_;
	this->num_rows_ = other.num_rows_;
	this->num_columns_ = other.num_columns_;
	this->values_ = std::move(other.values_);
	this->const_data_ptr_ = this->values_.data();
	other.const_data_ptr_ = nullptr;
}

void Matrix::operator*=(const double constant)
{
	const auto n = this->num_values();
	const auto incx = 1;
	cblas_dscal(n, constant, this->values_.data(), incx);	
}

Matrix& Matrix::be_inverse(void) {
	REQUIRE(this->is_square_matrix(), "invertable matrix should be square matrix");

	const auto ipiv = this->PLU_decomposition();

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int n = static_cast<int>(this->num_rows_);
	const lapack_int lda = n;
	const lapack_int info = LAPACKE_dgetri(matrix_layout, n, this->values_.data(), lda, ipiv.data());

	REQUIRE(info == 0, "info should be 0 when success matrix inverse");
	//info > 0 "U is singular matrix in L-U decomposition"
	//info < 0 "fail to inverse the matrix"

	return *this;
}

void Matrix::change_columns(const size_t start_column_index, const Matrix& A)
{
	REQUIRE(this->num_rows_ == A.num_rows_, "size should be matched");
	REQUIRE(start_column_index + A.num_columns_ <= this->num_columns_, "index can not exceed given range");

	for (size_t i = 0; i < A.num_rows_; ++i)
	{
		for (size_t j = 0; j < A.num_columns_; ++j)
		{
			this->value_at(i, start_column_index + j) = A.at(i, j);
		}
	}
}

void Matrix::change_rows(const size_t start_row_index, const Matrix& A) 
{
	REQUIRE(this->num_columns_ == A.num_columns_, "size should be matched");
	REQUIRE(start_row_index + A.num_rows_ <= this->num_rows_, "index can not exceed given range");
	REQUIRE(!this->is_transposed() && !A.is_transposed(), "it should be not transposed for this routine");

	const auto jump_index = start_row_index * this->num_columns_;
	std::copy(A.values_.begin(), A.values_.end(), this->values_.begin() + jump_index);
}

//Matrix Matrix::operator*(const Matrix& other) const 
//{
//	REQUIRE(this->num_column_ == other.num_row_, "size should be matched to operate *");
//
//	Matrix result(this->num_row_, other.num_column_);
//	ms::gemm(*this, other, result.values_.data());
//	return result;
//}

bool  Matrix::operator==(const Matrix& other) const 
{
	if (this->num_rows_ != other.num_rows_)
		return false;

	if (this->num_columns_ != other.num_columns_)
		return false;

	for (size_t i = 0; i < this->num_values(); ++i) {
		if (this->const_data_ptr_[i] != other.const_data_ptr_[i])
			return false;
	}

	return true;
}

Matrix Matrix::transpose(void) const 
{
	auto result = *this;
	result.be_transpose();
	return result;
}

Matrix Matrix::inverse(void) const 
{
	auto result = *this;
	return result.be_inverse();
}

double& Matrix::value_at(const size_t row, const size_t column) 
{
	if (this->is_transposed())
		return this->values_[column * this->num_rows_ + row];
	else
		return this->values_[row * this->num_columns_ + column];
}

std::vector<int> Matrix::PLU_decomposition(void) 
{
	REQUIRE(!this->is_transposed(), "matrix should not be transposed");

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int m = static_cast<int>(this->num_rows_);
	const lapack_int n = static_cast<int>(this->num_columns_);
	const lapack_int lda = n;
	std::vector<int> ipiv(std::min(m, n));
	lapack_int info = LAPACKE_dgetrf(matrix_layout, m, n, this->values_.data(), lda, ipiv.data());

	REQUIRE(0 <= info, "info should be greater than 0 when sucess PLU decomposition");
	return ipiv;
}

Matrix_Wrapper::Matrix_Wrapper(const size_t num_row, const size_t num_column, const double* ptr)
{
	this->num_rows_ = num_row;
	this->num_columns_ = num_column;
	this->const_data_ptr_ = ptr;
}

namespace ms
{
	void gemm(const Matrix_Base& A, const Matrix_Base& B, double* output_ptr)
	{
		REQUIRE(A.num_column() == B.num_row(), "size should be matched for matrix multiplication");

		const auto layout = CBLAS_LAYOUT::CblasRowMajor;
		const auto transA = A.transpose_type();
		const auto transB = B.transpose_type();
		const auto m = static_cast<MKL_INT>(A.num_row());
		const auto n = static_cast<MKL_INT>(B.num_column());
		const auto k = static_cast<MKL_INT>(A.num_column());
		const auto alpha = 1.0;
		const auto lda = static_cast<MKL_INT>(A.leading_dimension());
		const auto ldb = static_cast<MKL_INT>(B.leading_dimension());
		const auto beta = 0.0;
		const auto ldc = n;

		cblas_dgemm(layout, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, output_ptr, ldc);
	}
}

std::ostream& operator<<(std::ostream& os, const Matrix_Base& m) 
{
	return os << m.to_string();
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) 
{
	return os << m.to_string();
}