#include "../INC/Matrix.h"

Dynamic_Matrix_::Matrix(const size_t matrix_order) {
	*this = Matrix(matrix_order, matrix_order);
}

Dynamic_Matrix_::Matrix(const size_t matrix_order, const std::vector<double>& value){
	dynamic_require(matrix_order == value.size(), "num value of square matrix should be same with matrix order");
	this->num_row_ = matrix_order;
	this->num_column_ = matrix_order;
	this->value_.resize(matrix_order * matrix_order);

	for (size_t i = 0; i < matrix_order; ++i)
		this->at(i, i) = value[i];
}

Dynamic_Matrix_::Matrix(const size_t num_row, const size_t num_column){
	this->num_row_ = num_row;
	this->num_column_ = num_column;
	this->value_.resize(num_row * num_column);
}

Dynamic_Matrix_ Dynamic_Matrix_::operator*(const Dynamic_Matrix_& other) const {
	Dynamic_Matrix_ result(this->num_row_, other.num_column_);
	ms::gemm(*this, other, result.value_.data());
	return result;
}

bool  Dynamic_Matrix_::operator==(const Dynamic_Matrix_& other) const {
	return this->value_ == other.value_;
}

double Dynamic_Matrix_::at(const size_t row, const size_t column) const {
	dynamic_require(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
		return this->value_[column * this->num_row_ + row];
	else
		return this->value_[row * this->num_column_ + column];
}

Dynamic_Matrix_ Dynamic_Matrix_::transpose(void) const {
	auto result = *this;
	return result.be_transpose();
}

std::string Dynamic_Matrix_::to_string(void) const {
	std::string result;
	for (size_t i = 0; i < this->num_row_; ++i) {
		for (size_t j = 0; j < this->num_column_; ++j)
			result += ms::double_to_string(this->at(i, j)) + "   \t";
		result += "\n";
	}
	return result;
}

Dynamic_Matrix_ Dynamic_Matrix_::inverse(void) const {
	auto result = *this;
	return result.be_inverse();
}

std::pair<size_t, size_t> Dynamic_Matrix_::size(void) const {
	return std::make_pair(this->num_row_, this->num_column_);
}

Dynamic_Matrix_& Dynamic_Matrix_::be_transpose(void) {
	std::swap(this->num_row_, this->num_column_);

	if (this->is_transposed())
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	else
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasTrans;

	return *this;
}

Dynamic_Matrix_& Dynamic_Matrix_::be_inverse(void) {
	dynamic_require(this->is_square_matrix(), "invertable matrix should be square matrix");

	const auto ipiv = this->PLU_decomposition();

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int n = static_cast<int>(this->num_row_);
	const lapack_int lda = n;
	const lapack_int info = LAPACKE_dgetri(matrix_layout, n, this->value_.data(), lda, ipiv.data());

	dynamic_require(info == 0, "info should be 0 when success matrix inverse");
	//info > 0 "U is singular matrix in L-U decomposition"
	//info < 0 "fail to inverse the matrix"

	return *this;
}

void Dynamic_Matrix_::change_column(const size_t column_index, const Dynamic_Euclidean_Vector_& vec) {
	dynamic_require(column_index < this->num_column_, "column idnex can not exceed number of column");
	dynamic_require(this->num_row_ == vec.dimension(), "vector dimension should be matched with number of column");

	for (size_t i = 0; i < this->num_row_; ++i)
		this->at(i, column_index) = vec.at(i);
}

bool Dynamic_Matrix_::is_square_matrix(void) const {
	return this->num_row_ == this->num_column_;
}

bool Dynamic_Matrix_::is_transposed(void) const {
	return this->transpose_type_ != CBLAS_TRANSPOSE::CblasNoTrans;
}

bool Dynamic_Matrix_::is_in_range(const size_t row, const size_t column) const {
	return row < this->num_row_&& column < this->num_column_;
}

size_t Dynamic_Matrix_::leading_dimension(void) const {
	// num column before OP()
	if (this->is_transposed())
		return this->num_row_;
	else
		return this->num_column_;
}

double& Dynamic_Matrix_::at(const size_t row, const size_t column) {
	dynamic_require(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
		return this->value_[column * this->num_row_ + row];
	else
		return this->value_[row * this->num_column_ + column];
}

std::vector<int> Dynamic_Matrix_::PLU_decomposition(void) {
	dynamic_require(!this->is_transposed(), "matrix should not be transposed");

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int m = static_cast<int>(this->num_row_);
	const lapack_int n = static_cast<int>(this->num_column_);
	const lapack_int lda = n;
	std::vector<int> ipiv(std::min(m, n));
	lapack_int info = LAPACKE_dgetrf(matrix_layout, m, n, this->value_.data(), lda, ipiv.data());

	dynamic_require(0 <= info, "info should be greater than 0 when sucess PLU decomposition");
	return ipiv;
}


namespace ms {
	void gemm(const Dynamic_Matrix_& A, const Dynamic_Matrix_& B, double* output_ptr) {
		if (A.num_column_ != B.num_row_)
			throw std::length_error("length is not matched");

		const CBLAS_LAYOUT layout = CBLAS_LAYOUT::CblasRowMajor;
		const CBLAS_TRANSPOSE transA = A.transpose_type_;
		const CBLAS_TRANSPOSE transB = B.transpose_type_;
		const MKL_INT m = static_cast<MKL_INT>(A.num_row_);
		const MKL_INT n = static_cast<MKL_INT>(B.num_column_);
		const MKL_INT k = static_cast<MKL_INT>(A.num_column_);
		const double alpha = 1;
		const MKL_INT lda = static_cast<MKL_INT>(A.leading_dimension());
		const MKL_INT ldb = static_cast<MKL_INT>(B.leading_dimension());
		const double beta = 0;
		const MKL_INT ldc = n;

		cblas_dgemm(layout, transA, transB, m, n, k, alpha, A.value_.data(), lda, B.value_.data(), ldb, beta, output_ptr, ldc);
	}
}