#pragma once
#include "Euclidean_Vector.h"

#include <initializer_list>


using ushort = unsigned short;


template<typename Function, ushort range_dimension>
class Vector_Function
{
private:
	static constexpr size_t domain_dimension_ = Function::domain_dimension();

private:
	std::array<Function, range_dimension> functions_ = { 0 };

public:
	Vector_Function(void) = default;
	Vector_Function(const std::array<Function, range_dimension>& functions) : functions_(functions) {};
	
	template <typename... Args>
	Vector_Function(const Args... args) {
		static_require(sizeof...(Args) <= range_dimension, "Number of arguments can not exceed dimension");
		
		functions_ = { Function(args)... };		
	}

	Euclidean_Vector<range_dimension> operator()(const Euclidean_Vector<domain_dimension_>& space_vector) const {

		std::array<double, range_dimension> result = { 0 };
		for (size_t i = 0; i < range_dimension; ++i)
			result[i] = functions_[i](space_vector);

		return result;
	}

	std::vector<Euclidean_Vector<range_dimension>> operator()(const std::vector<Euclidean_Vector<domain_dimension_>>& space_vectors) const {
		const auto num_vector = space_vectors.size();

		std::vector<Euclidean_Vector<range_dimension>> result;
		result.reserve(num_vector);

		for (size_t i = 0; i < num_vector; ++i)
			result.push_back((*this)(space_vectors[i]));		

		return result;
	}

	bool operator==(const Vector_Function& other) const {
		return this->functions_ == other.functions_;
	}

	//Function& operator[](const size_t index) {
	//	dynamic_require(index < this->range_dimension(), "index can not exceed range dimension");
	//	return this->functions_[index];
	//}

	const Function& at(const size_t index) const {
		dynamic_require(index < range_dimension, "index can not exceed range dimension");
		return this->functions_[index];
	}

	template <ushort temp = range_dimension, std::enable_if_t<temp == 2, bool> = true >
	Vector_Function<Function, range_dimension + 1> cross_product(const Vector_Function& other) const {
		std::array<Function, range_dimension + 1> result;
			result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);

		return result;
	}

	template <ushort temp = range_dimension, std::enable_if_t<temp == 3, bool> = true >
	Vector_Function<Function, range_dimension> cross_product(const Vector_Function& other) const {
		static_require(range_dimension == 3, "cross product only valid for R^3 dimension");

		std::array<Function, range_dimension> result;
		result[0] = this->at(1) * other.at(2) - this->at(2) * other.at(1);
		result[1] = this->at(2) * other.at(0) - this->at(0) * other.at(2);
		result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);

		return result;
	}

	Function* data(void) {
		return this->functions_.data();
	}

	template <size_t variable_index>
	Vector_Function<Function, range_dimension> differentiate(void) const {
		static_require(variable_index < domain_dimension_, "variable index can not exceed domain dimension");
				
		std::array<Function, range_dimension> result;
		for (size_t i = 0; i < range_dimension; ++i)
			result[i] = this->functions_[i].differentiate<variable_index>();

		return result;
	}

	auto L2_norm(void) const {
		Function result(0);

		for (const auto& function : this->functions_)
			result += (function ^ 2);

		return result.root(0.5);
	}

	std::string to_string(void) const {
		std::string result;

		for (const auto& function : this->functions_)
			result += function.to_string() + "\n";

		return result;
	}
};


// Vector_Function class template for Range dimension is not compile time constant
template <typename Function>
using Dynamic_Vector_Function_ = Vector_Function<Function, 0>;


template <typename Function>
class Vector_Function<Function, 0>
{
private:
	static constexpr size_t domain_dimension_ = Function::domain_dimension();

private:
	std::vector<Function> functions_;

public:
	Vector_Function(std::vector<Function>&& functions) : functions_(std::move(functions)) {};
	Vector_Function(const std::initializer_list<Function> list) : functions_(list) {};

	Dynamic_Euclidean_Vector_ operator()(const Euclidean_Vector<domain_dimension_>& space_vector) const {
		const auto range_dimension = this->range_dimension();
		
		Dynamic_Euclidean_Vector_ result(range_dimension);
		for (size_t i = 0; i < range_dimension; ++i)
			result[i] = functions_.at(i)(space_vector);

		return result;
	}

	bool operator==(const Vector_Function& other) const {
		return this->functions_ == other.functions_;
	}

	const Function& at(const size_t index) const {
		dynamic_require(index < this->range_dimension(), "index can not exceed range dimension");
		return this->functions_[index];
	}

	size_t range_dimension(void) const {
		return this->functions_.size();
	}

	std::string to_string(void) const {
		std::string result;

		for (const auto& function : this->functions_)
			result += function.to_string() + "\n";

		return result;
	}

	//std::vector<Dynamic_Euclidean_Vector_> operator()(const std::vector<Euclidean_Vector<domain_dimension_>>& space_vectors) const {
	//	const auto num_vector = space_vectors.size();
	//	
	//	std::vector<Dynamic_Euclidean_Vector_> result;
	//	result.reserve(num_vector);
	//	for (size_t i = 0; i < num_vector; ++i)
	//		result.push_back( (*this)(space_vectors[i]));

	//	return result;
	//}

	//Function& operator[](const size_t index) {
	//	dynamic_require(index < this->range_dimension(), "index can not exceed range dimension");
	//	return this->functions_[index];
	//}

	//Vector_Function<Function> cross_product(const Vector_Function& other) const {
	//	dynamic_require(this->range_dimension() <= 3 && other.range_dimension() <= 3, "cross product only defined in range dimension 3 space");

	//	const auto range_dimension = this->range_dimension();

	//	Vector_Function<Function> result(3);
	//	if (range_dimension == 2)
	//		result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);
	//	else {
	//		result[0] = this->at(1) * other.at(2) - this->at(2) * other.at(1);
	//		result[1] = this->at(2) * other.at(0) - this->at(0) * other.at(2);
	//		result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);
	//	}


	//	return result;
	//}

	//template <size_t variable_index>
	//Vector_Function<Function> differentiate(void) const {
	//	static_require(variable_index < domain_dimension_, "variable index can not exceed domain dimension");

	//	const auto range_dimension = this->range_dimension();
	//	Vector_Function<Function> result(range_dimension);
	//	
	//	for (size_t i = 0; i < range_dimension; ++i) 
	//		result[i] = this->functions_[i].differentiate<variable_index>();
	//	
	//	return result;
	//}

	//auto L2_norm(void) const {		
	//	Function result(0);

	//	for (const auto& function : this->functions_)
	//		result += (function ^ 2);

	//	return result.root(0.5);
	//}
};


template <typename Function>
std::ostream& operator<<(std::ostream& os, const Dynamic_Vector_Function_<Function>& vf) {
	return os << vf.to_string();
}