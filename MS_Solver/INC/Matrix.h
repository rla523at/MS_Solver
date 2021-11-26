#pragma once
#include "Exception.h"

#include <mkl.h>
#include <sstream>
#include <iomanip>
#include <vector>

using ushort = unsigned short;

class Matrix;
class Matrix_Base
{
	template <typename T>
	using is_not_matrix = std::enable_if_t<!std::is_base_of_v<Matrix_Base, T>, bool>;

	template <typename T>
	using is_class = std::enable_if_t<std::is_class_v<T>, bool>;

public: //Command
	void transpose(void);

public: //Query
	Matrix operator*(const double constant) const;	
	template <typename V, is_class<V> = true, is_not_matrix<V> = true>	V operator*(const V& vec) const
	{
		REQUIRE(this->num_columns_ == vec.size(), "size should be mathced");

		const auto Layout = CBLAS_LAYOUT::CblasRowMajor;
		const auto trans = this->transpose_type_;
		const auto m = static_cast<MKL_INT>(this->num_rows_);
		const auto n = static_cast<MKL_INT>(this->num_columns_);
		const auto alpha = 1.0;
		const auto lda = static_cast<MKL_INT>(this->leading_dimension());
		const auto incx = 1;
		const auto beta = 0;
		const auto incy = 1;

		std::vector<double> values(this->num_rows_);
		cblas_dgemv(Layout, trans, m, n, alpha, this->const_data_ptr_, lda, vec.data(), incx, beta, values.data(), incy);

		return values;
	}
	Matrix operator*(const Matrix_Base& other) const;
	Matrix operator+(const Matrix_Base& other) const;

	double at(const size_t row, const size_t column) const;
	std::vector<double> column(const size_t column_index) const;
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

protected:
	bool is_transposed(void) const;
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
public:
	Matrix(void) = default;
	Matrix(const size_t matrix_order);
	Matrix(const size_t matrix_order, const std::vector<double>& value);
	Matrix(const size_t num_row, const size_t num_column);
	Matrix(const size_t num_row, const size_t num_column, std::vector<double>&& value);
	Matrix(const Matrix& other);
	Matrix(Matrix&& other) noexcept;

public://Command 
	void operator=(const Matrix& other);
	void operator=(Matrix&& other) noexcept;
	void operator*=(const double constant);

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
		std::copy(vec.begin(), vec.end(), this->values_.begin() + jump_index);

	}
	void change_rows(const size_t start_row_index, const Matrix& A);
	void change_columns(const size_t start_column_index, const Matrix& A);	

public://Query
	bool operator==(const Matrix& other) const;

	Matrix get_transpose(void) const;
	Matrix get_inverse(void) const;
	
private:
	double& value_at(const size_t row, const size_t column);
	std::vector<int> PLU_decomposition(void);

private:
	std::vector<double> values_;
};

class Matrix_Wrapper : public Matrix_Base
{
public:
	Matrix_Wrapper(const size_t num_row, const size_t num_column, const double* ptr);
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
}


std::ostream& operator<<(std::ostream& os, const Matrix_Base& m);
std::ostream& operator<<(std::ostream& os, const Matrix& m);


















//class Matrix;
//
//
//// Matrix template class
//template<size_t num_row, size_t num_column>
//class Static_Matrix
//{
//private:
//	template<size_t num_other_row, size_t num_other_column>
//	friend class Static_Matrix;
//	friend class Matrix;
//	friend class ms::BLAS;
//
//	static constexpr size_t num_value_ = num_row * num_column;
//
//private:
//	std::array<double, num_value_> values_ = { 0 };
//
//public:
//	Static_Matrix(void) = default;				
//	Static_Matrix(const std::array<double, num_value_>& values) : values_{ values } {};
//	Static_Matrix(const Matrix& dynamic_matrix);
//
//	template <typename... Args>
//	Static_Matrix(Args... args);
//
//	// diagonal matrix constructor
//	// SFINAE to resolve 1x1 matrix constructor ambiguity
//	template<size_t temp = num_row * num_column, std::enable_if_t<temp != 1, bool> = true>
//	Static_Matrix(const std::array<double, num_row>& values);
//
//public:
//	Static_Matrix& operator+=(const Static_Matrix& A);
//	Static_Matrix& operator-=(const Static_Matrix& A);
//	Static_Matrix& operator*=(const double scalar);
//	Static_Matrix& operator*=(const Static_Matrix<num_column, num_column>& A);
//
//	Static_Matrix operator+(const Static_Matrix& A) const;
//	Static_Matrix operator*(const double scalar) const;
//	Euclidean_Vector<num_row> operator*(const Euclidean_Vector<num_column>& x) const;
//	Matrix operator*(const Matrix& dynamic_matrix) const;
//
//	template <size_t other_num_column>
//	Static_Matrix<num_row, other_num_column> operator*(const Static_Matrix<num_column, other_num_column>& A) const;
//
//	bool operator==(const Static_Matrix & A) const;
//
//public:
//	void change_columns(const size_t column_index, const Matrix& A);
//
//	template <size_t num_other_column>
//	void change_columns(const size_t column_index, const Static_Matrix<num_row, num_other_column>& A);
//
//public:
//	double at(const size_t row_index, const size_t column_index) const;
//	Euclidean_Vector<num_row> column(const size_t column_index) const;
//	std::string to_string(void) const;
//
//private:
//	double& value_at(const size_t row_index, const size_t column_index);
//};
//
//template <size_t matrix_order>
//Static_Matrix(const std::array<double, matrix_order>& values)->Static_Matrix<matrix_order, matrix_order>;
//
//template<size_t num_row, size_t num_column>
//Static_Matrix<num_row, num_column> operator*(const double scalar, const Static_Matrix<num_row, num_column>& A);
//
//template<size_t num_row, size_t num_column>
//std::ostream& operator<<(std::ostream& os, const Static_Matrix<num_row, num_column>& m);









//std::ostream& operator<<(std::ostream& os, const Matrix& m);





////template definition part
//template<size_t num_row, size_t num_column>
//template<size_t temp, std::enable_if_t<temp != 1, bool>>
//Static_Matrix<num_row, num_column>::Static_Matrix(const std::array<double, num_row>& values) {
//	static_require(num_row == num_column, "It should be square matrix");
//	for (size_t i = 0; i < num_row; ++i)
//		this->value_at(i, i) = values[i];
//}
//
//template<size_t num_row, size_t num_column>
//Static_Matrix<num_row,num_column>::Static_Matrix(const Matrix& dynamic_matrix) {
//	dynamic_require(dynamic_matrix.size() == std::make_pair(num_row, num_column), "both matrix should be same size");
//	std::copy(dynamic_matrix.values_.begin(), dynamic_matrix.values_.end(), this->values_.begin());
//}
//
//template<size_t num_row, size_t num_column>
//template <typename... Args>
//Static_Matrix<num_row, num_column>::Static_Matrix(Args... args) : values_{ static_cast<double>(args)... } {
//	static_require(sizeof...(Args) == this->num_value_, "Number of arguments should be same with (num row * num column)");
//	static_require(ms::are_arithmetics<Args...>, "every arguments should be arithmetics");
//};
//
//template<size_t num_row, size_t num_column>
//Static_Matrix<num_row, num_column>& Static_Matrix<num_row, num_column>::operator+=(const Static_Matrix& A) {
//	if constexpr (this->num_value_ < ms::blas_axpy_criteria) {
//		for (size_t i = 0; i < this->num_value_; ++i)
//			this->values_[i] += A.values_[i];
//	}
//	else 
//		cblas_daxpy(this->num_value_, 1.0, A.values_.data(), 1, this->values_.data(), 1);
//
//	return *this;
//}
//
//template<size_t num_row, size_t num_column>
//Static_Matrix<num_row, num_column>& Static_Matrix<num_row, num_column>::operator-=(const Static_Matrix& A) {
//	if constexpr (this->num_value_ < ms::blas_axpy_criteria) {
//		for (size_t i = 0; i < this->num_value_; ++i)
//			this->values_[i] -= A.values_[i];
//	}
//	else
//		cblas_daxpy(this->num_value_, -1.0, A.values_.data(), 1, this->values_.data(), 1);
//
//	return *this;
//}
//
//template<size_t num_row, size_t num_column>
//Static_Matrix<num_row, num_column>& Static_Matrix<num_row, num_column>::operator*=(const double scalar) {
//	if constexpr (this->num_value_ < ms::blas_dscal_criteria) {
//		for (size_t i = 0; i < num_row * num_column; ++i)
//			this->values_[i] *= scalar;
//	}
//	else
//		cblas_dscal(this->num_value_, scalar, this->values_.data(), 1);
//
//	return *this;
//}
//
//template<size_t num_row, size_t num_column>
//Static_Matrix<num_row, num_column>& Static_Matrix<num_row, num_column>::operator*=(const Static_Matrix<num_column, num_column>& A) {
//	*this = (*this) * A;
//	return *this;
//}
//
//template<size_t num_row, size_t num_column>
//Static_Matrix<num_row, num_column> Static_Matrix<num_row, num_column>::operator+(const Static_Matrix& A) const {	
//	Static_Matrix result = *this;
//	return result += A;
//}
//
//template<size_t num_row, size_t num_column>
//Static_Matrix<num_row, num_column> Static_Matrix<num_row, num_column>::operator*(const double scalar) const {
//	Static_Matrix result = *this;
//
//	if constexpr (this->num_value_ < ms::blas_dscal_criteria) {
//		for (size_t i = 0; i < num_row * num_column; ++i)
//			result.values_[i] *= scalar;
//	}
//	else 
//		cblas_dscal(this->num_value_, scalar, result.values_.data(), 1);
//
//	return result;
//}
//
//template<size_t num_row, size_t num_column>
//Matrix Static_Matrix<num_row, num_column>::operator*(const Matrix& B) const {
//	dynamic_require(num_column == B.num_row_, "dimension should be matched for matrix multiplication");
//
//	const auto layout = CBLAS_LAYOUT::CblasRowMajor;
//	const auto transA = CBLAS_TRANSPOSE::CblasNoTrans;
//	const auto transB = B.transpose_type_;
//	const auto m = static_cast<MKL_INT>(num_row);
//	const auto n = static_cast<MKL_INT>(B.num_column_);
//	const auto k = static_cast<MKL_INT>(num_column);
//	const auto alpha = 1.0;
//	const auto lda = static_cast<MKL_INT>(num_column);
//	const auto ldb = static_cast<MKL_INT>(B.leading_dimension());
//	const auto beta = 0.0;
//	const auto ldc = n;
//
//	Matrix result(m, n);
//	cblas_dgemm(layout, transA, transB, m, n, k, alpha, this->values_.data(), lda, B.values_.data(), ldb, beta, result.values_.data(), ldc);
//
//	return result;
//}
//
//template <size_t num_row, size_t num_column>
//template <size_t other_num_column>
//Static_Matrix<num_row, other_num_column> Static_Matrix<num_row, num_column>::operator*(const Static_Matrix<num_column, other_num_column>& A) const {
//	const auto layout = CBLAS_LAYOUT::CblasRowMajor;
//	const auto transA = CBLAS_TRANSPOSE::CblasNoTrans;
//	const auto transB = CBLAS_TRANSPOSE::CblasNoTrans;
//	const auto m = static_cast<MKL_INT>(num_row);
//	const auto n = static_cast<MKL_INT>(other_num_column);
//	const auto k = static_cast<MKL_INT>(num_column);
//	const double alpha = 1;
//	const MKL_INT lda = static_cast<MKL_INT>(num_column);
//	const MKL_INT ldb = static_cast<MKL_INT>(other_num_column);
//	const double beta = 0;
//	const MKL_INT ldc = n;
//
//	Static_Matrix<num_row, other_num_column> result;
//	cblas_dgemm(layout, transA, transB, m, n, k, alpha, this->values_.data(), lda, A.values_.data(), ldb, beta, result.values_.data(), ldc);
//
//	return result;
//}
//
//template<size_t num_row, size_t num_column>
//Euclidean_Vector<num_row> Static_Matrix<num_row, num_column>::operator*(const Euclidean_Vector<num_column>& x) const {
//	std::array<double, num_row> result = { 0 };
//	if constexpr (this->num_value_ < ms::blas_mv_criteria) {
//		for (size_t i = 0; i < num_row; ++i)
//			for (size_t j = 0; j < num_column; ++j)
//				result[i] += this->values_[i * num_column + j] * x.at(j);
//	}
//	return result;
//}
//
//template<size_t num_row, size_t num_column>
//bool Static_Matrix<num_row, num_column>::operator==(const Static_Matrix& A) const {
//	return this->values_ == A.values_;
//}
//
//
//template<size_t num_row, size_t num_column>
//double Static_Matrix<num_row, num_column>::at(const size_t row_index, const size_t column_index) const {
//	dynamic_require(row_index < num_row && column_index < num_column, "index can not exceed given range");
//	return this->values_[row_index * num_column + column_index];
//}
//
//template<size_t num_row, size_t num_column>
//Euclidean_Vector<num_row> Static_Matrix<num_row, num_column>::column(const size_t column_index) const {
//	dynamic_require(column_index < num_column, "index can not exceed given range");
//
//	std::array<double, num_row> column_value;	
//	for (size_t i = 0; i < num_row; ++i)
//		column_value[i] = this->at(i, column_index);
//
//	return column_value;
//}
//
//template<size_t num_row, size_t num_column>
//std::string Static_Matrix<num_row, num_column>::to_string(void) const {
//	std::ostringstream oss;
//	oss << std::setprecision(16) << std::showpoint << std::left;
//	for (size_t i = 0; i < num_row; ++i) {
//		for (size_t j = 0; j < num_column; ++j)
//			oss << std::setw(25) << this->at(i, j);
//		oss << "\n";
//	}
//	return oss.str();
//}
//
//template<size_t num_row, size_t num_column>
//void Static_Matrix<num_row, num_column>::change_columns(const size_t start_column_index, const Matrix& A) {
//	const auto [num_A_row, num_A_column] = A.size();
//	dynamic_require(num_A_row == num_row, "number of row should be same");
//	dynamic_require(start_column_index + num_A_column <= num_column, "index can not exceed given range");
//
//	for (size_t i = 0; i < num_row; ++i)
//		for (size_t j = 0; j < num_A_column; ++j)
//			this->value_at(i, start_column_index + j) = A.at(i, j);
//}
//
//template<size_t num_row, size_t num_column>
//template <size_t num_other_column>
//void Static_Matrix<num_row, num_column>::change_columns(const size_t start_column_index, const Static_Matrix<num_row, num_other_column>& A) {
//	dynamic_require(start_column_index + num_other_column <= num_column, "index can not exceed given range");
//
//	for (size_t i = 0; i < num_row; ++i)
//		for (size_t j = 0; j < num_other_column; ++j)
//			this->value_at(i, start_column_index + j) = A.at(i, j);
//}
//
////template<size_t num_row, size_t num_column>
////void Matrix<num_row, num_column>::change_column(const size_t column_index, const Euclidean_Vector<num_row>& x) {
////	dynamic_require(column_index < num_column, "index can not exceed given range");
////	
////	for (size_t i = 0; i < num_row; ++i)
////		this->at(i, column_index) = x[i];
////}
////
//
//template<size_t num_row, size_t num_column>
//double& Static_Matrix<num_row, num_column>::value_at(const size_t row_index, const size_t column_index) {
//	// we should range check before call private at!
//	// dynamic_require(row_index < num_row && column_index < num_column, "index can not exceed given range"); 
//	return this->values_[row_index * num_column + column_index];
//}
//
//template<size_t num_row, size_t num_column>
//std::ostream& operator<<(std::ostream& os, const Static_Matrix<num_row, num_column>& m) {
//	return os << m.to_string();
//}
//
//template<size_t num_row, size_t num_column>
//Static_Matrix<num_row, num_column> operator*(const double scalar, const Static_Matrix<num_row, num_column>& A) {
//	return A * scalar;
//}
//
//template <size_t num_row>
//Euclidean_Vector<num_row> Matrix::column(const size_t column_index) const {
//	dynamic_require(this->num_row_ == num_row, "this dynamic matrix should have given row dimension");
//
//	std::array<double, num_row> column_values = { 0 };
//
//	for (size_t i = 0; i < num_row; ++i)
//		column_values[i] = this->at(i, column_index);
//
//	return column_values;
//}
//
//template <size_t dimension>
//void Matrix::change_row(const size_t row_index, const Euclidean_Vector<dimension>& vec) {
//	dynamic_require(row_index < this->num_row_,		"column idnex can not exceed number of column");
//	dynamic_require(this->num_column_ == dimension, "vector dimension should be matched with number of column");
//
//	const auto jump_index = row_index * this->num_column_;
//	std::copy(vec.cbegin(), vec.cend(), this->values_.begin() + jump_index);
//}
//
//template <size_t dimension>
//void Matrix::change_column(const size_t column_index, const Euclidean_Vector<dimension>& vec) {
//	dynamic_require(column_index < this->num_column_,	"column idnex can not exceed number of column");
//	dynamic_require(this->num_row_ == dimension,		"vector dimension should be matched with number of row");
//
//	for (size_t i = 0; i < this->num_row_; ++i)
//		this->value_at(i, column_index) = vec.at(i);
//}
//
//template <size_t num_row, size_t num_column>
//void Matrix::change_rows(const size_t start_row_index, const Static_Matrix<num_row, num_column>& A) {
//	dynamic_require(this->num_column_ == num_column,				"dimension should be matched");
//	dynamic_require(start_row_index + num_row <= this->num_row_,	"range can't not exceed given range");
//	dynamic_require(!this->is_transposed(),							"it should be not transposed for this routine");
//
//	const auto jump_index = start_row_index * this->num_column_;
//	std::copy(A.values_.begin(), A.values_.end(), this->values_.begin() + jump_index);
//}
//
//template <size_t num_row, size_t num_column>
//void Matrix::change_columns(const size_t start_column_index, const Static_Matrix<num_row, num_column>& A) {
//	dynamic_require(this->num_row_ == num_row,								"dimension should be matched");
//	dynamic_require(start_column_index + num_column <= this->num_column_,	"range can't not exceed given range");
//	dynamic_require(!this->is_transposed(),									"it should be not transposed for this routine");
//
//	for (size_t i = 0; i < this->num_row_; ++i) 
//		for (size_t j = 0; j < num_column; ++j) 
//			this->value_at(i, start_column_index + j) = A.at(i, j);
//}