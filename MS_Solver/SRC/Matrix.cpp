#include "../INC/Matrix.h"

Dynamic_Matrix::Dynamic_Matrix(const size_t matrix_order) {
	*this = Dynamic_Matrix(matrix_order, matrix_order);
}

Dynamic_Matrix::Dynamic_Matrix(const size_t matrix_order, const std::vector<double>& value) {
	dynamic_require(matrix_order == value.size(), "num value of square matrix should be same with matrix order");
	this->num_row_ = matrix_order;
	this->num_column_ = matrix_order;
	this->values_.resize(matrix_order * matrix_order);

	for (size_t i = 0; i < matrix_order; ++i)
		this->value_at(i, i) = value[i];
}

Dynamic_Matrix::Dynamic_Matrix(const size_t num_row, const size_t num_column) {
	this->num_row_ = num_row;
	this->num_column_ = num_column;
	this->values_.resize(num_row * num_column);
}

Dynamic_Matrix::Dynamic_Matrix(const size_t num_row, const size_t num_column, std::vector<double>&& value) {
	dynamic_require(num_row * num_column == value.size(), "num value should be same with matrix size");

	this->num_row_ = num_row;
	this->num_column_ = num_column;
	this->values_ = std::move(value);
}

Dynamic_Matrix Dynamic_Matrix::operator*(const Dynamic_Matrix& other) const {
	Dynamic_Matrix result(this->num_row_, other.num_column_);
	ms::gemm(*this, other, result.values_.data());
	return result;
}

bool  Dynamic_Matrix::operator==(const Dynamic_Matrix& other) const {
	return this->num_row_ == other.num_row_ && this->values_ == other.values_;
}

double Dynamic_Matrix::at(const size_t row, const size_t column) const {
	dynamic_require(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
		return this->values_[column * this->num_row_ + row];
	else
		return this->values_[row * this->num_column_ + column];
}

Dynamic_Matrix Dynamic_Matrix::transpose(void) const {
	auto result = *this;
	return result.be_transpose();
}

std::string Dynamic_Matrix::to_string(void) const {
	std::string result;
	for (size_t i = 0; i < this->num_row_; ++i) {
		for (size_t j = 0; j < this->num_column_; ++j)
			result += ms::double_to_string(this->at(i, j)) + "   \t";
		result += "\n";
	}
	return result;
}

Dynamic_Matrix Dynamic_Matrix::inverse(void) const {
	auto result = *this;
	return result.be_inverse();
}

std::pair<size_t, size_t> Dynamic_Matrix::size(void) const {
	return std::make_pair(this->num_row_, this->num_column_);
}

Dynamic_Matrix& Dynamic_Matrix::be_transpose(void) {
	std::swap(this->num_row_, this->num_column_);

	if (this->is_transposed())
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	else
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasTrans;

	return *this;
}

Dynamic_Matrix& Dynamic_Matrix::be_inverse(void) {
	dynamic_require(this->is_square_matrix(), "invertable matrix should be square matrix");

	const auto ipiv = this->PLU_decomposition();

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int n = static_cast<int>(this->num_row_);
	const lapack_int lda = n;
	const lapack_int info = LAPACKE_dgetri(matrix_layout, n, this->values_.data(), lda, ipiv.data());

	dynamic_require(info == 0, "info should be 0 when success matrix inverse");
	//info > 0 "U is singular matrix in L-U decomposition"
	//info < 0 "fail to inverse the matrix"

	return *this;
}

void Dynamic_Matrix::change_column(const size_t column_index, const Dynamic_Euclidean_Vector& vec) {
	dynamic_require(column_index < this->num_column_, "column idnex can not exceed number of column");
	dynamic_require(this->num_row_ == vec.dimension(), "vector dimension should be matched with number of column");

	for (size_t i = 0; i < this->num_row_; ++i)
		this->value_at(i, column_index) = vec.at(i);
}

void Dynamic_Matrix::change_row(const size_t start_row_index, const Dynamic_Euclidean_Vector& vec) {
	dynamic_require(this->num_column_ == vec.dimension(), "dimension should be matched");
	dynamic_require(start_row_index <= this->num_row_, "range can't not exceed given range");
	dynamic_require(!this->is_transposed(), "it should be not transposed for this routine");

	const auto jump_index = start_row_index * this->num_column_;
	std::copy(vec.begin(), vec.end(), this->values_.begin() + jump_index);
}

void Dynamic_Matrix::change_rows(const size_t start_row_index, const Dynamic_Matrix& A) {
	dynamic_require(this->num_column_ == A.num_column_, "dimension should be matched");
	dynamic_require(start_row_index + A.num_row_ <= this->num_row_, "range can't not exceed given range");
	dynamic_require(!this->is_transposed() && A.is_transposed(), "it should be not transposed for this routine");

	const auto jump_index = start_row_index * this->num_column_;
	std::copy(A.values_.begin(), A.values_.end(), this->values_.begin() + jump_index);
}

bool Dynamic_Matrix::is_square_matrix(void) const {
	return this->num_row_ == this->num_column_;
}

bool Dynamic_Matrix::is_transposed(void) const {
	return this->transpose_type_ != CBLAS_TRANSPOSE::CblasNoTrans;
}

bool Dynamic_Matrix::is_in_range(const size_t row, const size_t column) const {
	return row < this->num_row_&& column < this->num_column_;
}

size_t Dynamic_Matrix::leading_dimension(void) const {
	// num column before OP()
	if (this->is_transposed())
		return this->num_row_;
	else
		return this->num_column_;
}

double& Dynamic_Matrix::value_at(const size_t row, const size_t column) {
	//dynamic_require(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
		return this->values_[column * this->num_row_ + row];
	else
		return this->values_[row * this->num_column_ + column];
}

std::vector<int> Dynamic_Matrix::PLU_decomposition(void) {
	dynamic_require(!this->is_transposed(), "matrix should not be transposed");

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int m = static_cast<int>(this->num_row_);
	const lapack_int n = static_cast<int>(this->num_column_);
	const lapack_int lda = n;
	std::vector<int> ipiv(std::min(m, n));
	lapack_int info = LAPACKE_dgetrf(matrix_layout, m, n, this->values_.data(), lda, ipiv.data());

	dynamic_require(0 <= info, "info should be greater than 0 when sucess PLU decomposition");
	return ipiv;
}

std::ostream& operator<<(std::ostream& os, const Dynamic_Matrix& m) {
	return os << m.to_string();
}

namespace ms {
	void gemm(const Dynamic_Matrix& A, const Dynamic_Matrix& B, double* output_ptr) {
		dynamic_require(A.num_column_ == B.num_row_, "dimension should be matched for matrix multiplication");

		const auto layout = CBLAS_LAYOUT::CblasRowMajor;
		const auto transA = A.transpose_type_;
		const auto transB = B.transpose_type_;
		const auto m = static_cast<MKL_INT>(A.num_row_);
		const auto n = static_cast<MKL_INT>(B.num_column_);
		const auto k = static_cast<MKL_INT>(A.num_column_);
		const auto alpha = 1.0;
		const auto lda = static_cast<MKL_INT>(A.leading_dimension());
		const auto ldb = static_cast<MKL_INT>(B.leading_dimension());
		const auto beta = 0.0;
		const auto ldc = n;

		cblas_dgemm(layout, transA, transB, m, n, k, alpha, A.values_.data(), lda, B.values_.data(), ldb, beta, output_ptr, ldc);
	}

	void gemvpv(const Dynamic_Matrix& A, const Dynamic_Euclidean_Vector& v1, Dynamic_Euclidean_Vector& v2) {
		const auto layout = CBLAS_LAYOUT::CblasRowMajor;
		const auto transA = A.transpose_type_;
		const auto m = static_cast<MKL_INT>(A.num_row_);
		const auto n = static_cast<MKL_INT>(A.num_column_);
		const auto alpha = 1.0;
		const auto lda = static_cast<MKL_INT>(A.leading_dimension());		
		const auto incx = 1;
		const auto beta = 1.0;
		const auto incy = 1;

		cblas_dgemv(layout, transA, m, n, alpha, A.values_.data(), lda, v1.data(), incx, beta, v2.data(), incy);
	}

}