#include "../INC/Matrix.h"

Matrix::Matrix(const size_t matrix_order) {
	REQUIRE(matrix_order != 0, "matrix order can not be 0");

	this->num_row_ = matrix_order;
	this->num_column_ = matrix_order;
	this->values_.resize(this->num_values());
	this->const_data_ptr_ = this->values_.data();
	
	for (size_t i = 0; i < matrix_order; ++i)
		this->value_at(i, i) = 1.0;
}
Matrix::Matrix(const size_t matrix_order, const std::vector<double>& value) {
	REQUIRE(matrix_order != 0, "matrix order can not be 0");
	REQUIRE(matrix_order == value.size(), "num value of square matrix should be same with matrix order");

	this->num_row_ = matrix_order;
	this->num_column_ = matrix_order;
	this->values_.resize(this->num_values());
	this->const_data_ptr_ = this->values_.data();

	for (size_t i = 0; i < matrix_order; ++i)
		this->value_at(i, i) = value[i];
}
Matrix::Matrix(const size_t num_row, const size_t num_column) {
	REQUIRE(num_row * num_column != 0, "number of row or number of column can not be 0");

	this->num_row_ = num_row;
	this->num_column_ = num_column;
	this->values_.resize(this->num_values());
	this->const_data_ptr_ = this->values_.data();
}
Matrix::Matrix(const size_t num_row, const size_t num_column, std::vector<double>&& value) {
	REQUIRE(num_row * num_column != 0, "number of row or number of column can not be 0");
	REQUIRE(num_row * num_column == value.size(), "num value should be same with matrix size");

	this->num_row_ = num_row;
	this->num_column_ = num_column;
	this->values_ = std::move(value);
	this->const_data_ptr_ = this->values_.data();
}

Matrix& Matrix::be_inverse(void) {
	REQUIRE(this->is_square_matrix(), "invertable matrix should be square matrix");

	const auto ipiv = this->PLU_decomposition();

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int n = static_cast<int>(this->num_row_);
	const lapack_int lda = n;
	const lapack_int info = LAPACKE_dgetri(matrix_layout, n, this->values_.data(), lda, ipiv.data());

	REQUIRE(info == 0, "info should be 0 when success matrix inverse");
	//info > 0 "U is singular matrix in L-U decomposition"
	//info < 0 "fail to inverse the matrix"

	return *this;
}
void Matrix::change_rows(const size_t start_row_index, const Matrix& A) {
	REQUIRE(this->num_column_ == A.num_column_, "dimension should be matched");
	REQUIRE(start_row_index + A.num_row_ <= this->num_row_, "index can not exceed given range");
	REQUIRE(!this->is_transposed() && !A.is_transposed(), "it should be not transposed for this routine");

	const auto jump_index = start_row_index * this->num_column_;
	std::copy(A.values_.begin(), A.values_.end(), this->values_.begin() + jump_index);
}

Matrix Matrix::operator*(const Matrix& other) const {
	REQUIRE(this->num_column_ == other.num_row_, "dimension should be matched to operate *");

	Matrix result(this->num_row_, other.num_column_);
	ms::gemm(*this, other, result.values_.data());
	return result;
}
bool  Matrix::operator==(const Matrix& other) const {
	if (this->num_row_ != other.num_row_)
		return false;

	if (this->num_column_ != other.num_column_)
		return false;
		
	for (size_t i = 0; i < this->num_values(); ++i) {
		if (this->const_data_ptr_[i] != other.const_data_ptr_[i])
			return false;
	}

	return true;
}

Matrix Matrix::transpose(void) const {
	auto result = *this;
	result.be_transpose();
	return result;
}
Matrix Matrix::inverse(void) const {
	auto result = *this;
	return result.be_inverse();
}

double& Matrix::value_at(const size_t row, const size_t column) {
	if (this->is_transposed())
		return this->values_[column * this->num_row_ + row];
	else
		return this->values_[row * this->num_column_ + column];
}
std::vector<int> Matrix::PLU_decomposition(void) {
	REQUIRE(!this->is_transposed(), "matrix should not be transposed");

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int m = static_cast<int>(this->num_row_);
	const lapack_int n = static_cast<int>(this->num_column_);
	const lapack_int lda = n;
	std::vector<int> ipiv(std::min(m, n));
	lapack_int info = LAPACKE_dgetrf(matrix_layout, m, n, this->values_.data(), lda, ipiv.data());

	REQUIRE(0 <= info, "info should be greater than 0 when sucess PLU decomposition");
	return ipiv;
}

Matrix Matrix_Wrapper::operator*(const Matrix& m) const {
	const auto [num_row, num_column] = m.size();

	std::vector<double> values(this->num_row_ * num_column);
	ms::gemm(*this, m, values.data());
	return { this->num_row_, num_column, std::move(values) };
}

std::ostream& operator<<(std::ostream& os, const Matrix_Base& m) {
	return os << m.to_string();
}
