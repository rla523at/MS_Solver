#pragma once
#include "Exception.h"

#include <initializer_list>
#include <vector>

using ushort = unsigned short;

template <typename Function>
class Vector_Function
{
public:
	Vector_Function(void) = default;
	Vector_Function(std::vector<Function>&& functions) : functions_(std::move(functions)) {};
	Vector_Function(const std::initializer_list<Function> list) : functions_(list) {};

public://Query
	template <typename V>	std::vector<double> operator()(const V& space_vector) const {
		const auto range_dimension = this->size();

		std::vector<double> result(range_dimension);
		for (size_t i = 0; i < range_dimension; ++i)
			result[i] = this->functions_[i](space_vector);

		return result;
	}
	const Function& operator[](const size_t index) const {
		REQUIRE(index < this->size(), "index can not exceed range size");
		return this->functions_[index];
	}
	bool operator==(const Vector_Function& other) const {
		return this->functions_ == other.functions_;
	}

	const Function& at(const size_t index) const {
		REQUIRE(index < this->size(), "index can not exceed range size");
		return this->functions_[index];
	}
	Vector_Function<Function> cross_product(const Vector_Function& other) const {
		constexpr auto result_range_dimension = 3;
		std::vector<Function> result(result_range_dimension);

		if (this->size() == 2)
		{
			result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);
		}
		else if (this->size() == 3)
		{
			result[0] = this->at(1) * other.at(2) - this->at(2) * other.at(1);
			result[1] = this->at(2) * other.at(0) - this->at(0) * other.at(2);
			result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);
		}
		else
			EXCEPTION("cross product only valid for R^3 dimension");

		return result;
	};
	auto L2_norm(void) const {
		Function result = 0.0;

		for (const auto& function : this->functions_)
			result += (function ^ 2);

		return result.root(0.5);
	}
	Vector_Function<Function> get_differentiate(const ushort variable_index) const {
		auto differentiate_functions = this->functions_;
		
		for (auto& differentiate_function : differentiate_functions)
			differentiate_function.differentiate(variable_index);

		return differentiate_functions;
	}
	size_t size(void) const {
		return this->functions_.size();
	}
	std::string to_string(void) const {
		std::string result;

		for (const auto& function : this->functions_)
			result += function.to_string() + "\t";

		return result;
	}

private:
	std::vector<Function> functions_;
};

template <typename Function>
std::ostream& operator<<(std::ostream& os, const Vector_Function<Function>& vf) {
	return os << vf.to_string();
}



//template<typename Function, ushort range_dimension_>
//class Vector_Function
//{
//private:
//	static constexpr size_t domain_dimension_ = Function::domain_dimension();
//
//private:
//	std::array<Function, range_dimension_> functions_ = { 0 };
//
//public:
//	Vector_Function(void) = default;
//	Vector_Function(const std::array<Function, range_dimension_>& functions) : functions_(functions) {};
//
//	template <typename... Args>
//	Vector_Function(const Args... args) {
//		static_require(sizeof...(Args) <= range_dimension_, "Number of arguments can not exceed dimension");
//
//		functions_ = { Function(args)... };
//	}
//
//	Euclidean_Vector<range_dimension_> operator()(const Euclidean_Vector<domain_dimension_>& space_vector) const {
//
//		std::array<double, range_dimension_> result = { 0 };
//		for (size_t i = 0; i < range_dimension_; ++i)
//			result[i] = functions_[i](space_vector);
//
//		return result;
//	}
//
//	std::vector<Euclidean_Vector<range_dimension_>> operator()(const std::vector<Euclidean_Vector<domain_dimension_>>& space_vectors) const {
//		const auto num_vector = space_vectors.size();
//
//		std::vector<Euclidean_Vector<range_dimension_>> result;
//		result.reserve(num_vector);
//
//		for (size_t i = 0; i < num_vector; ++i)
//			result.push_back((*this)(space_vectors[i]));
//
//		return result;
//	}
//
//	bool operator==(const Vector_Function& other) const {
//		return this->functions_ == other.functions_;
//	}
//
//	const Function& operator[](const size_t index) const {
//		dynamic_require(index < range_dimension_, "index can not exceed range dimension");
//		return this->functions_[index];
//	}
//
//	const Function& at(const size_t index) const {
//		dynamic_require(index < range_dimension_, "index can not exceed range dimension");
//		return this->functions_[index];
//	}
//
	//template <ushort temp = range_dimension_, std::enable_if_t<temp == 2, bool> = true >
	//Vector_Function<Function, range_dimension_ + 1> cross_product(const Vector_Function& other) const {
	//	std::array<Function, range_dimension_ + 1> result;
	//	result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);

	//	return result;
	//}

	//template <ushort temp = range_dimension_, std::enable_if_t<temp == 3, bool> = true >
	//Vector_Function<Function, range_dimension_> cross_product(const Vector_Function& other) const {
	//	static_require(range_dimension_ == 3, "cross product only valid for R^3 dimension");

	//	std::array<Function, range_dimension_> result;
	//	result[0] = this->at(1) * other.at(2) - this->at(2) * other.at(1);
	//	result[1] = this->at(2) * other.at(0) - this->at(0) * other.at(2);
	//	result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);

	//	return result;
	//}
//
//	Function* data(void) {
//		return this->functions_.data();
//	}
//
	//Vector_Function<Function, range_dimension_> differentiate(const ushort variable_index) const {
	//	dynamic_require(variable_index < domain_dimension_, "variable index can not exceed domain dimension");

	//	std::array<Function, range_dimension_> result;
	//	for (size_t i = 0; i < range_dimension_; ++i)
	//		result[i] = this->functions_[i].differentiate(variable_index);

	//	return result;
	//}
//
//	Function inner_product(const Vector_Function& other) const {
//		Function result;
//
//		for (ushort i = 0; i < range_dimension_; ++i)
//			result += this->functions_[i] * other.functions_[i];
//
//		return result;
//	}
//

//
//	std::string to_string(void) const {
//		std::string result;
//
//		for (const auto& function : this->functions_)
//			result += function.to_string() + "\n";
//
//		return result;
//	}
//
//	constexpr ushort range_dimension(void) const {
//		return range_dimension_;
//	}
//};
//
//

//
//
//template <typename Function, size_t range_dimension>
//std::ostream& operator<<(std::ostream& os, const Vector_Function<Function, range_dimension>& vf) {
//	return os << vf.to_string();
//}



// Vector_Function class template for Range dimension is not compile time constant
//
//
//template <typename Function, ushort range_num_row, ushort range_num_column>
//class Matrix_Function
//{
//private:
//	static constexpr ushort num_value_ = range_num_row * range_num_column;
//	static constexpr ushort domain_dimension_ = Function::domain_dimension();
//
//private:
//	std::array<Function, num_value_> functions_ = { 0 };
//
//public:
//	Matrix_Function(void) = default;
//	Matrix_Function(const std::array<Function, num_value_>& functions) : functions_(functions) {};
//
//	template <typename... Args>
//	Matrix_Function(const Args... args) {
//		static_require(sizeof...(Args) <= num_value_, "Number of arguments can not exceed dimension");
//
//		functions_ = { Function(args)... };
//	}
//
//	Static_Matrix<range_num_row, range_num_column> operator()(const Euclidean_Vector<domain_dimension_>& space_vector) const {
//		std::array<double, num_value_> values;
//
//		for (ushort i = 0; i < range_num_row; ++i)
//			for (ushort j = 0; j < range_num_column; ++j)
//				values[i * range_num_column + j] = this->at(i, j)(space_vector);
//
//		return values;
//	}
//
//	bool operator==(const Matrix_Function& other) const {
//		return this->functions_ == other.functions_;
//	}
//
//	const Function& at(const ushort row_index, const ushort column_index) const {
//		dynamic_require(row_index < range_num_row&& column_index < range_num_column, "index can not exceed given range");
//		return this->functions_[row_index * range_num_column + column_index];
//	}
//
//	void change_column(const ushort column_index, const Vector_Function<Function, range_num_row>& column_vector) {
//		for (ushort i = 0; i < range_num_row; ++i)
//			this->at(i, column_index) = column_vector[i];
//	}
//
//	std::string to_string(void) const {
//		std::string str;
//		for (ushort i = 0; i < range_num_row; ++i) {
//			for (ushort j = 0; j < range_num_column; ++j)
//				str += this->at(i, j).to_string() + "\t";
//			str += "\n";
//		}
//		return str;
//	}
//
//private:
//	Function& at(const ushort row_index, const ushort column_index) {
//		dynamic_require(row_index < range_num_row&& column_index < range_num_column, "index can not exceed given range");
//		return this->functions_[row_index * range_num_column + column_index];
//	}
//};
//
//template <typename Function, ushort range_num_row, ushort range_num_column>
//std::ostream& operator<<(std::ostream& os, const Matrix_Function<Function, range_num_row, range_num_column>& mf) {
//	return os << mf.to_string();
//}
//
//
//namespace ms {
//	template <typename Function>
//	void gemv(const Matrix& A, const Dynamic_Vector_Function<Function>& v, Function* ptr) {
//		//code for dynmaic matrix * dynmaic vector function => vector function
//		const auto [num_row, num_column] = A.size();
//		const auto range_dimension = v.range_dimension();
//		dynamic_require(num_column == range_dimension, "number of column should be same with range dimension");
//
//		for (size_t i = 0; i < num_row; ++i)
//			for (size_t j = 0; j < num_column; ++j)
//				ptr[i] += A.at(i, j) * v.at(j);
//	}
//
//	template <typename Function, ushort range_dimension>
//	Matrix_Function<Function, range_dimension, range_dimension> Jacobian(const Vector_Function<Function, range_dimension>& vector_function) {
//		constexpr auto num_value = range_dimension * range_dimension;
//		std::array<Function, num_value> functions;
//
//		for (ushort i = 0; i < range_dimension; ++i) {
//			const auto differential_vector_function = vector_function.differentiate(i);
//			for (ushort j = 0; j < range_dimension; ++j)
//				functions[j * range_dimension + i] = differential_vector_function[j];
//		}
//
//		return functions;
//	}
//
//	template <typename Function>
//	Function scalar_triple_product(const Vector_Function<Function, 3>& a, const Vector_Function<Function, 3>& b, const Vector_Function<Function, 3>& c) {
//		return a.inner_product(b.cross_product(c));
//	}
//}





























































//template<typename Function, ushort range_dimension_>
//class Vector_Function
//{
//private:
//	static constexpr size_t domain_dimension_ = Function::domain_dimension();
//
//private:
//	std::array<Function, range_dimension_> functions_ = { 0 };
//
//public:
//	Vector_Function(void) = default;
//	Vector_Function(const std::array<Function, range_dimension_>& functions) : functions_(functions) {};
//	
//	template <typename... Args>
//	Vector_Function(const Args... args) {
//		static_require(sizeof...(Args) <= range_dimension_, "Number of arguments can not exceed dimension");
//		
//		functions_ = { Function(args)... };		
//	}
//
//	Euclidean_Vector<range_dimension_> operator()(const Euclidean_Vector<domain_dimension_>& space_vector) const {
//
//		std::array<double, range_dimension_> result = { 0 };
//		for (size_t i = 0; i < range_dimension_; ++i)
//			result[i] = functions_[i](space_vector);
//
//		return result;
//	}
//
//	std::vector<Euclidean_Vector<range_dimension_>> operator()(const std::vector<Euclidean_Vector<domain_dimension_>>& space_vectors) const {
//		const auto num_vector = space_vectors.size();
//
//		std::vector<Euclidean_Vector<range_dimension_>> result;
//		result.reserve(num_vector);
//
//		for (size_t i = 0; i < num_vector; ++i)
//			result.push_back((*this)(space_vectors[i]));		
//
//		return result;
//	}
//
//	bool operator==(const Vector_Function& other) const {
//		return this->functions_ == other.functions_;
//	}
//
//	const Function& operator[](const size_t index) const {
//		dynamic_require(index < range_dimension_, "index can not exceed range dimension");
//		return this->functions_[index];
//	}
//
//	const Function& at(const size_t index) const {
//		dynamic_require(index < range_dimension_, "index can not exceed range dimension");
//		return this->functions_[index];
//	}
//
//	template <ushort temp = range_dimension_, std::enable_if_t<temp == 2, bool> = true >
//	Vector_Function<Function, range_dimension_ + 1> cross_product(const Vector_Function& other) const {
//		std::array<Function, range_dimension_ + 1> result;
//			result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);
//
//		return result;
//	}
//
//	template <ushort temp = range_dimension_, std::enable_if_t<temp == 3, bool> = true >
//	Vector_Function<Function, range_dimension_> cross_product(const Vector_Function& other) const {
//		static_require(range_dimension_ == 3, "cross product only valid for R^3 dimension");
//
//		std::array<Function, range_dimension_> result;
//		result[0] = this->at(1) * other.at(2) - this->at(2) * other.at(1);
//		result[1] = this->at(2) * other.at(0) - this->at(0) * other.at(2);
//		result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);
//
//		return result;
//	}
//
//	Function* data(void) {
//		return this->functions_.data();
//	}
//
//	Vector_Function<Function, range_dimension_> differentiate(const ushort variable_index) const {
//		dynamic_require(variable_index < domain_dimension_, "variable index can not exceed domain dimension");
//				
//		std::array<Function, range_dimension_> result;
//		for (size_t i = 0; i < range_dimension_; ++i)
//			result[i] = this->functions_[i].differentiate(variable_index);
//
//		return result;
//	}
//
//	Function inner_product(const Vector_Function& other) const {
//		Function result;
//
//		for (ushort i = 0; i < range_dimension_; ++i)
//			result += this->functions_[i] * other.functions_[i];
//
//		return result;
//	}
//
//	auto L2_norm(void) const {
//		Function result(0);
//
//		for (const auto& function : this->functions_)
//			result += (function ^ 2);
//
//		return result.root(0.5);
//	}
//
//	std::string to_string(void) const {
//		std::string result;
//
//		for (const auto& function : this->functions_)
//			result += function.to_string() + "\n";
//
//		return result;
//	}
//
//	constexpr ushort range_dimension(void) const {
//		return range_dimension_;
//	}
//};
//
//
//template <typename Function, size_t num_row, size_t num_column>
//Vector_Function<Function, num_row> operator*(const Static_Matrix<num_row, num_column>& matrix, const Vector_Function<Function, num_column>& vector_function) {
//	std::array<Function, num_row> functions;
//
//	for (size_t i = 0; i < num_row; ++i)
//		for (size_t j = 0; j < num_column; ++j)
//			functions[i] += matrix.at(i, j) * vector_function[j];
//
//	return functions;
//}
//
//
//template <typename Function, size_t range_dimension>
//std::ostream& operator<<(std::ostream& os, const Vector_Function<Function, range_dimension>& vf) {
//	return os << vf.to_string();
//}
//
//
//
//// Vector_Function class template for Range dimension is not compile time constant
//template <typename Function>
//class Dynamic_Vector_Function
//{
//private:
//	static constexpr size_t domain_dimension_ = Function::domain_dimension();
//
//private:
//	std::vector<Function> functions_;
//
//public:
//	Dynamic_Vector_Function(std::vector<Function>&& functions) : functions_(std::move(functions)) {};
//	Dynamic_Vector_Function(const std::initializer_list<Function> list) : functions_(list) {};
//
//	Dynamic_Euclidean_Vector operator()(const Euclidean_Vector<domain_dimension_>& space_vector) const {
//		const auto range_dimension_ = this->range_dimension();
//		
//		std::vector<double> result(range_dimension_);
//		for (size_t i = 0; i < range_dimension_; ++i)
//			result[i] = functions_.at(i)(space_vector);
//
//		return std::move(result);
//	}
//
//	const Function& operator[](const size_t index) const {
//		dynamic_require(index < this->range_dimension(), "index can not exceed range dimension");
//		return this->functions_[index];
//	}
//
//	bool operator==(const Dynamic_Vector_Function& other) const {
//		return this->functions_ == other.functions_;
//	}
//
//	const Function& at(const size_t index) const {
//		dynamic_require(index < this->range_dimension(), "index can not exceed range dimension");
//		return this->functions_[index];
//	}
//
//	size_t range_dimension(void) const {
//		return this->functions_.size();
//	}
//
//	std::string to_string(void) const {
//		std::string result;
//
//		for (const auto& function : this->functions_)
//			result += function.to_string() + "\n";
//
//		return result;
//	}
//};
//
//
//template <typename Function>
//std::ostream& operator<<(std::ostream& os, const Dynamic_Vector_Function<Function>& vf) {
//	return os << vf.to_string();
//}
//
//
//template <typename Function, ushort range_num_row, ushort range_num_column>
//class Matrix_Function
//{
//private:
//	static constexpr ushort num_value_			= range_num_row * range_num_column;
//	static constexpr ushort domain_dimension_	= Function::domain_dimension();
//
//private:
//	std::array<Function, num_value_> functions_ = { 0 };
//
//public:
//	Matrix_Function(void) = default;
//	Matrix_Function(const std::array<Function, num_value_>&functions) : functions_(functions) {};
//
//	template <typename... Args>
//	Matrix_Function(const Args... args) {
//		static_require(sizeof...(Args) <= num_value_, "Number of arguments can not exceed dimension");
//
//		functions_ = { Function(args)... };
//	}
//
//	Static_Matrix<range_num_row, range_num_column> operator()(const Euclidean_Vector<domain_dimension_>& space_vector) const {
//		std::array<double, num_value_> values;
//
//		for (ushort i = 0; i < range_num_row; ++i)
//			for (ushort j = 0; j < range_num_column; ++j)
//				values[i * range_num_column + j] = this->at(i, j)(space_vector);
//		
//		return values;
//	}
//
//	bool operator==(const Matrix_Function& other) const {
//		return this->functions_ == other.functions_;
//	}
//
//	const Function& at(const ushort row_index, const ushort column_index) const {
//		dynamic_require(row_index < range_num_row && column_index < range_num_column, "index can not exceed given range");
//		return this->functions_[row_index * range_num_column + column_index];
//	}
//
//	void change_column(const ushort column_index, const Vector_Function<Function, range_num_row>& column_vector) {
//		for (ushort i = 0; i < range_num_row; ++i) 
//			this->at(i, column_index) = column_vector[i];		
//	}
//
//	std::string to_string(void) const {
//		std::string str;
//		for (ushort i = 0; i < range_num_row; ++i) {
//			for (ushort j = 0; j < range_num_column; ++j)
//				str += this->at(i, j).to_string() + "\t";
//			str += "\n";
//		}
//		return str;
//	}
//
//private:
//	Function& at(const ushort row_index, const ushort column_index) {
//		dynamic_require(row_index < range_num_row&& column_index < range_num_column, "index can not exceed given range");
//		return this->functions_[row_index * range_num_column + column_index];
//	}
//};
//
//template <typename Function, ushort range_num_row, ushort range_num_column>
//std::ostream& operator<<(std::ostream& os, const Matrix_Function<Function, range_num_row, range_num_column>& mf) {
//	return os << mf.to_string();
//}
//
//
//namespace ms {
//	template <typename Function>
//	void gemv(const Matrix& A, const Dynamic_Vector_Function<Function>& v, Function* ptr) {
//		//code for dynmaic matrix * dynmaic vector function => vector function
//		const auto [num_row, num_column] = A.size();
//		const auto range_dimension = v.range_dimension();
//		dynamic_require(num_column == range_dimension, "number of column should be same with range dimension");
//
//		for (size_t i = 0; i < num_row; ++i)
//			for (size_t j = 0; j < num_column; ++j)
//				ptr[i] += A.at(i, j) * v.at(j);
//	}
//
//	template <typename Function, ushort range_dimension>
//	Matrix_Function<Function, range_dimension, range_dimension> Jacobian(const Vector_Function<Function, range_dimension>& vector_function) {
//		constexpr auto num_value = range_dimension * range_dimension;
//		std::array<Function, num_value> functions;
//
//		for (ushort i = 0; i < range_dimension; ++i) {
//			const auto differential_vector_function = vector_function.differentiate(i);
//			for (ushort j = 0; j < range_dimension; ++j)
//				functions[j * range_dimension + i] = differential_vector_function[j];
//		}
//		
//		return functions;
//	}
//
//	template <typename Function>
//	Function scalar_triple_product(const Vector_Function<Function, 3>& a, const Vector_Function<Function, 3>& b, const Vector_Function<Function, 3>& c) {
//		return a.inner_product(b.cross_product(c));
//	}
//}