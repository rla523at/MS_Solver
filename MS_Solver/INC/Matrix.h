#pragma once
#include "Euclidean_Vector.h"

using ushort = unsigned short;

template<size_t num_row, size_t num_column>
class Matrix;

using Dynamic_Matrix_ = Matrix<0, 0>;


namespace ms {
	inline constexpr ushort blas_axpy_criteria = 20;
	inline constexpr ushort blas_mv_criteria = 50;

	void gemm(const Dynamic_Matrix_& A, const Dynamic_Matrix_& B, double* output_ptr);

	template<size_t num_row, size_t num_column>
	void gemm(const Dynamic_Matrix_& A, const Dynamic_Matrix_& B, Matrix<num_row, num_column>& output_matrix);
	
}


// Matrix template class
template<size_t num_row, size_t num_column>
class Matrix
{
private:
	friend Dynamic_Matrix_;

	template<size_t num_row, size_t num_column>
	friend void ms::gemm(const Dynamic_Matrix_& A, const Dynamic_Matrix_& B, Matrix<num_row, num_column>& output_matrix);

	static constexpr size_t num_value_ = num_row * num_column;

private:
	std::array<double, num_value_> values_ = { 0 };

public:
	Matrix(void) = default;
	Matrix(const std::array<double, num_value_>& values) : values_{ values } {};
	Matrix(const Matrix<0, 0>& dynamic_matrix);
	template <typename... Args>
	Matrix(Args... args);

	Matrix& operator+=(const Matrix& A);
	Matrix operator+(const Matrix& A) const;
	Matrix operator*(const double scalar) const;
	Matrix<0, 0> operator*(const Matrix<0, 0>& dynamic_matrix) const;
	Euclidean_Vector<num_row> operator*(const Euclidean_Vector<num_column>& x) const;
	bool operator==(const Matrix & A) const;

	double at(const size_t row_index, const size_t column_index) const;
	Euclidean_Vector<num_row> column(const size_t column_index) const;
	std::string to_string(void) const;

	//void change_column(const size_t column_index, const Euclidean_Vector<num_row>& x);

//private:
	//double& at(const size_t row_index, const size_t column_index);
};


// Dynamic matrix := Matrix template class specialization
template<>
class Matrix<0, 0>
{
	template<size_t num_row, size_t num_column>
	friend class Matrix;

	friend void ms::gemm(const Dynamic_Matrix_& A, const Dynamic_Matrix_& B, double* ptr);


private:
	CBLAS_TRANSPOSE transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	size_t num_row_ = 0;
	size_t num_column_ = 0;
	std::vector<double> values_;

public:
	Matrix(const size_t matrix_order);
	Matrix(const size_t matrix_order, const std::vector<double>& value);
	Matrix(const size_t num_row, const size_t num_column);
	Matrix(const size_t num_row, const size_t num_column, std::vector<double>&& value);		

	Dynamic_Matrix_ operator*(const Dynamic_Matrix_& other) const;
	bool operator==(const Dynamic_Matrix_& other) const;
		
	double at(const size_t row, const size_t column) const;
	Dynamic_Matrix_ transpose(void) const;
	std::string to_string(void) const;
	Dynamic_Matrix_ inverse(void) const;
	std::pair<size_t, size_t> size(void) const;

	Dynamic_Matrix_& be_transpose(void);
	Dynamic_Matrix_& be_inverse(void);
	void change_column(const size_t row_index, const Dynamic_Euclidean_Vector_& vec);

	template <size_t num_row>
	Euclidean_Vector<num_row> column(const size_t column_index) const;

	template <size_t dimension>
	void change_column(const size_t column_index, const Euclidean_Vector<dimension>& vec);

	template <size_t num_row, size_t num_column>
	void change_rows(const size_t start_row_index, const Matrix<num_row, num_column>& A);

	template <size_t num_row, size_t num_column>
	void change_columns(const size_t start_column_index, const Matrix<num_row, num_column>& A);


private:
	bool is_square_matrix(void) const;
	bool is_transposed(void) const;
	bool is_in_range(const size_t irow, const size_t jcolumn) const;
	size_t leading_dimension(void) const;

	double& at(const size_t row, const size_t column);
	std::vector<int> PLU_decomposition(void);
};

template<size_t num_row, size_t num_column>
Matrix<num_row, num_column> operator*(const double scalar, const Matrix<num_row, num_column>& A);

template<size_t num_row, size_t num_column>
std::ostream& operator<<(std::ostream& os, const Matrix<num_row, num_column>& m);


//template definition part
namespace ms {
	template<size_t num_row, size_t num_column>
	void gemm(const Dynamic_Matrix_& A, const Dynamic_Matrix_& B, Matrix<num_row, num_column>& output_matrix) {
		ms::gemm(A, B, output_matrix.values_.data());
	}
}

template<size_t num_row, size_t num_column>
Matrix<num_row,num_column>::Matrix(const Matrix<0, 0>& dynamic_matrix) {
	dynamic_require(num_row == dynamic_matrix.num_row_ && num_column == dynamic_matrix.num_column_, "both matrix should be same size");
	std::copy(dynamic_matrix.values_.begin(), dynamic_matrix.values_.end(), this->values_.begin());
}

template<size_t num_row, size_t num_column>
template <typename... Args>
Matrix<num_row, num_column>::Matrix(Args... args) : values_{ static_cast<double>(args)... } {
	static_require(sizeof...(Args) <= num_row * num_column, "Number of arguments should be less then (num row * num column)");
	static_require(ms::are_arithmetics<Args...>, "every arguments should be arithmetics");
};

template<size_t num_row, size_t num_column>
Matrix<num_row, num_column>& Matrix<num_row, num_column>::operator+=(const Matrix& A) {
	if constexpr (this->num_value_ < ms::blas_axpy_criteria) {
		for (size_t i = 0; i < this->num_value_; ++i)
			this->values_[i] += A.values_[i];
	}
	else 
		cblas_daxpy(this->num_value_, 1.0, A.values_.data(), 1, this->values_.data(), 1);

	return *this;
}

template<size_t num_row, size_t num_column>
Matrix<num_row, num_column> Matrix<num_row, num_column>::operator+(const Matrix& A) const {	
	Matrix result = *this;
	return result += A;
	//for (size_t i = 0; i < num_row * num_column; ++i)
	//	result.values_[i] += A.values_[i];
	//return result;
}

template<size_t num_row, size_t num_column>
Matrix<num_row, num_column> Matrix<num_row, num_column>::operator*(const double scalar) const {
	Matrix result = *this;

	if constexpr (this->num_value_ < 10) {
		for (size_t i = 0; i < num_row * num_column; ++i)
			result.values_[i] *= scalar;
	}
	else 
		cblas_dscal(this->num_value_, scalar, result.values_.data(), 1);

	return result;
}

template<size_t num_row, size_t num_column>
Dynamic_Matrix_ Matrix<num_row, num_column>::operator*(const Dynamic_Matrix_& B) const {
	dynamic_require(num_column == B.num_row_, "dimension should be matched for matrix multiplication");

	const auto layout = CBLAS_LAYOUT::CblasRowMajor;
	const auto transA = CBLAS_TRANSPOSE::CblasNoTrans;
	const auto transB = B.transpose_type_;
	const auto m = static_cast<MKL_INT>(num_row);
	const auto n = static_cast<MKL_INT>(B.num_column_);
	const auto k = static_cast<MKL_INT>(num_column);
	const auto alpha = 1.0;
	const auto lda = static_cast<MKL_INT>(num_column);
	const auto ldb = static_cast<MKL_INT>(B.leading_dimension());
	const auto beta = 0.0;
	const auto ldc = n;

	Dynamic_Matrix_ result(m, n);
	cblas_dgemm(layout, transA, transB, m, n, k, alpha, this->values_.data(), lda, B.values_.data(), ldb, beta, result.values_.data(), ldc);

	return result;
}

template<size_t num_row, size_t num_column>
bool Matrix<num_row, num_column>::operator==(const Matrix& A) const {
	return this->values_ == A.values_;
}

template<size_t num_row, size_t num_column>
Euclidean_Vector<num_row> Matrix<num_row, num_column>::operator*(const Euclidean_Vector<num_column>& x) const {
	std::array<double, num_row> result = { 0 };
	if constexpr (num_row * num_column < ms::blas_mv_criteria) {
		for (size_t i = 0; i < num_row; ++i)
			for (size_t j = 0; j < num_column; ++j)
				result[i] += this->values_[i * num_column + j] * x.at(j);
	}
	return result;
}

template<size_t num_row, size_t num_column>
double Matrix<num_row, num_column>::at(const size_t row_index, const size_t column_index) const {
	dynamic_require(row_index < num_row && column_index < num_column, "index can not exceed given range");
	return this->values_[row_index * num_column + column_index];
}

template<size_t num_row, size_t num_column>
Euclidean_Vector<num_row> Matrix<num_row, num_column>::column(const size_t column_index) const {
	dynamic_require(column_index < num_column, "index can not exceed given range");

	std::array<double, num_row> column_value;	
	for (size_t i = 0; i < num_row; ++i)
		column_value[i] = this->at(i, column_index);

	return column_value;
}

template<size_t num_row, size_t num_column>
std::string Matrix<num_row, num_column>::to_string(void) const {
	std::string result;
	for (size_t i = 0; i < num_row; ++i) {
		for (size_t j = 0; j < num_column; ++j)
			result += ms::double_to_string(this->at(i, j)) + "   \t";
		result += "\n";
	}
	return result;
}

//template<size_t num_row, size_t num_column>
//void Matrix<num_row, num_column>::change_column(const size_t column_index, const Euclidean_Vector<num_row>& x) {
//	dynamic_require(column_index < num_column, "index can not exceed given range");
//	
//	for (size_t i = 0; i < num_row; ++i)
//		this->at(i, column_index) = x[i];
//}
//
//template<size_t num_row, size_t num_column>
//double& Matrix<num_row, num_column>::at(const size_t row_index, const size_t column_index) {
//	dynamic_require(row_index < num_row && column_index < num_column, "index can not exceed given range");
//	return this->values_[row_index * num_column + column_index];
//}

template<size_t num_row, size_t num_column>
std::ostream& operator<<(std::ostream& os, const Matrix<num_row, num_column>& m) {
	return os << m.to_string();
}

template<size_t num_row, size_t num_column>
Matrix<num_row, num_column> operator*(const double scalar, const Matrix<num_row, num_column>& A) {
	return A * scalar;
}

template <size_t dimension>
void Dynamic_Matrix_::change_column(const size_t column_index, const Euclidean_Vector<dimension>& vec) {
	dynamic_require(column_index < this->num_column_,	"column idnex can not exceed number of column");
	dynamic_require(this->num_row_ == dimension,		"vector dimension should be matched with number of row");

	for (size_t i = 0; i < this->num_row_; ++i)
		this->at(i, column_index) = vec.at(i);
}

template <size_t num_row, size_t num_column>
void Dynamic_Matrix_::change_rows(const size_t start_row_index, const Matrix<num_row, num_column>& A) {
	dynamic_require(this->num_column_ == num_column,				"dimension should be matched");
	dynamic_require(start_row_index + num_row <= this->num_row_,	"range can't not exceed given range");
	dynamic_require(!this->is_transposed(),							"it should be not transposed for this routine");

	const auto jump_index = start_row_index * this->num_column_;
	std::copy(A.values_.begin(), A.values_.end(), this->values_.begin() + jump_index);
}

template <size_t num_row, size_t num_column>
void Dynamic_Matrix_::change_columns(const size_t start_column_index, const Matrix<num_row, num_column>& A) {
	dynamic_require(this->num_row_ == num_row,								"dimension should be matched");
	dynamic_require(start_column_index + num_column <= this->num_column_,	"range can't not exceed given range");
	dynamic_require(!this->is_transposed(),									"it should be not transposed for this routine");

	for (size_t i = 0; i < this->num_row_; ++i) 
		for (size_t j = 0; j < num_column; ++j) 
			this->at(i, start_column_index + j) = A.at(i, j);	
}