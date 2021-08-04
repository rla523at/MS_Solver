#pragma once
#include "EuclideanVector.h"

#include <initializer_list>


//template<size_t domain_dimension, size_t range_dimension = 0>
//class Vector_Function
//{
//private:
//	std::array<Polynomial<domain_dimension>, range_dimension> functions_;
//};


// Vector_Function class template for Range dimension is not compile time constant
template<typename Function>
class Vector_Function
{
private:
	static constexpr size_t domain_dimension_ = Function::domain_dimension();

private:
	std::vector<Function> functions_;

public:
	Vector_Function(const size_t range_dimension) : functions_(range_dimension, 0) {};
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

	Function& operator[](const size_t index) {
		dynamic_require(index < this->range_dimension(), "index can not exceed range dimension");
		return this->functions_[index];
	}

	const Function& at(const size_t index) const {
		dynamic_require(index < this->range_dimension(), "index can not exceed range dimension");
		return this->functions_[index];
	}

	

	template <size_t variable_index>
	Vector_Function<Function> differentiate(void) const {
		static_require(variable_index < domain_dimension_, "variable index can not exceed domain dimension");

		const auto range_dimension = this->range_dimension();
		Vector_Function<Function> result(range_dimension);
		
		for (size_t i = 0; i < range_dimension; ++i) 
			result[i] = this->functions_[i].differentiate<variable_index>();
		
		return result;
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
};


template <typename Function>
std::ostream& operator<<(std::ostream& os, const Vector_Function<Function>& vf) {
	return os << vf.to_string();
}