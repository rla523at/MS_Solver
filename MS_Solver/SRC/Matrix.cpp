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


//double Matrix::at(const size_t row, const size_t column) const {
//	REQUIRE(this->is_in_range(row, column), "matrix indexes should not exceed given range");
//	if (this->is_transposed())
//		return this->values_[column * this->num_row_ + row];
//	else
//		return this->values_[row * this->num_column_ + column];
//}

Matrix Matrix::transpose(void) const {
	auto result = *this;
	result.be_transpose();
	return result;
}
//
//std::string Matrix::to_string(void) const {
//	std::ostringstream oss;
//	oss << std::setprecision(16) << std::showpoint << std::left;
//	for (size_t i = 0; i < this->num_row_; ++i) {
//		for (size_t j = 0; j < this->num_column_; ++j)
//			oss << std::setw(25) << this->at(i, j);
//		oss << "\n";
//	}
//	return oss.str();
//}
//
//Matrix Matrix::inverse(void) const {
//	auto result = *this;
//	return result.be_inverse();
//}
//
//std::pair<size_t, size_t> Matrix::size(void) const {
//	return std::make_pair(this->num_row_, this->num_column_);
//}
//
//bool Matrix::is_finite(void) const {
//	for (const auto value : this->values_) {
//		if (!std::isfinite(value))
//			return false;
//	}
//	return true;
//}
//


//
//void Matrix::change_column(const size_t column_index, const Dynamic_Euclidean_Vector& vec) {
//	REQUIRE(column_index < this->num_column_, "column idnex can not exceed number of column");
//	REQUIRE(this->num_row_ == vec.dimension(), "vector dimension should be matched with number of column");
//
//	for (size_t i = 0; i < this->num_row_; ++i)
//		this->value_at(i, column_index) = vec.at(i);
//}
//
//void Matrix::change_row(const size_t start_row_index, const Dynamic_Euclidean_Vector& vec) {
//	REQUIRE(this->num_column_ == vec.dimension(), "dimension should be matched");
//	REQUIRE(start_row_index <= this->num_row_, "index can not exceed given range");
//	REQUIRE(!this->is_transposed(), "it should be not transposed for this routine");
//
//	const auto jump_index = start_row_index * this->num_column_;
//	std::copy(vec.begin(), vec.end(), this->values_.begin() + jump_index);
//}
//
//void Matrix::change_rows(const size_t start_row_index, const Matrix& A) {
//	REQUIRE(this->num_column_ == A.num_column_, "dimension should be matched");
//	REQUIRE(start_row_index + A.num_row_ <= this->num_row_, "index can not exceed given range");
//	REQUIRE(!this->is_transposed() && !A.is_transposed(), "it should be not transposed for this routine");
//
//	const auto jump_index = start_row_index * this->num_column_;
//	std::copy(A.values_.begin(), A.values_.end(), this->values_.begin() + jump_index);
//}
//
//std::vector<double> Matrix::row(const size_t row_index) const {
//	REQUIRE(row_index < this->num_row_, "index can not exceed given range");
//
//	std::vector<double> row_values(this->num_column_);
//
//	if (this->is_transposed()) {
//		for (size_t i = 0; i < this->num_column_; ++i)
//			row_values[i] = this->at(row_index, i);
//	}
//	else {
//		const auto start_index = row_index * this->num_column_;
//		const auto start_iter = this->values_.begin() + start_index;		
//		std::copy(start_iter, start_iter + this->num_column_, row_values.begin());
//	}
//
//	return row_values;
//}
//
//bool Matrix::is_square_matrix(void) const {
//	return this->num_row_ == this->num_column_;
//}
//
//bool Matrix::is_transposed(void) const {
//	return this->transpose_type_ != CBLAS_TRANSPOSE::CblasNoTrans;
//}
//
//bool Matrix::is_in_range(const size_t row, const size_t column) const {
//	return row < this->num_row_&& column < this->num_column_;
//}
//
//size_t Matrix::leading_dimension(void) const {
//	// num column before OP()
//	if (this->is_transposed())
//		return this->num_row_;
//	else
//		return this->num_column_;
//}
//

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

//std::ostream& operator<<(std::ostream& os, const Matrix& m) {
//	return os << m.to_string();
//}
//
//namespace ms {
//	void BLAS::gemm(const Matrix& A, const Matrix& B, double* output_ptr) {
		//REQUIRE(A.num_column_ == B.num_row_, "dimension should be matched for matrix multiplication");

		//const auto layout = CBLAS_LAYOUT::CblasRowMajor;
		//const auto transA = A.transpose_type_;
		//const auto transB = B.transpose_type_;
		//const auto m = static_cast<MKL_INT>(A.num_row_);
		//const auto n = static_cast<MKL_INT>(B.num_column_);
		//const auto k = static_cast<MKL_INT>(A.num_column_);
		//const auto alpha = 1.0;
		//const auto lda = static_cast<MKL_INT>(A.leading_dimension());
		//const auto ldb = static_cast<MKL_INT>(B.leading_dimension());
		//const auto beta = 0.0;
		//const auto ldc = n;

		//cblas_dgemm(layout, transA, transB, m, n, k, alpha, A.values_.data(), lda, B.values_.data(), ldb, beta, output_ptr, ldc);
//	}
//
//	void BLAS::gemvpv(const Matrix& A, const Dynamic_Euclidean_Vector& v1, Dynamic_Euclidean_Vector& v2) {
//		REQUIRE(A.num_column_ == v1.dimension(), "dimension should be matched for matrix vector multiplication");
//
//		const auto layout = CBLAS_LAYOUT::CblasRowMajor;
//		const auto transA = A.transpose_type_;
//		const auto m = static_cast<MKL_INT>(A.num_row_);
//		const auto n = static_cast<MKL_INT>(A.num_column_);
//		const auto alpha = 1.0;
//		const auto lda = static_cast<MKL_INT>(A.leading_dimension());
//		const auto incx = 1;
//		const auto beta = 1.0;
//		const auto incy = 1;
//
//		cblas_dgemv(layout, transA, m, n, alpha, A.values_.data(), lda, v1.data(), incx, beta, v2.values_.data(), incy);
//	};
//}