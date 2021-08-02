#pragma once
#include "Polynomial.h"
#include <initializer_list>

template<size_t domain_dimension, size_t range_dimension = 0>
class Vector_Function
{
private:
	std::array<Polynomial<domain_dimension>, range_dimension> functions_;
};


template<size_t domain_dimension>
class Vector_Function<domain_dimension, 0>
{
private:
	std::vector<Polynomial<domain_dimension>> functions_;

public:
	Vector_Function(const size_t range_dimension) : functions_(range_dimension, 0) {};
	Vector_Function(const std::initializer_list<Polynomial<domain_dimension>> list) : functions_(list) {};

	Dynamic_Euclidean_Vector_ operator()(const Euclidean_Vector<domain_dimension>& space_vector) const {
		const auto range_dimension = this->range_dimension();
		
		Dynamic_Euclidean_Vector_ result(range_dimension);
		for (size_t i = 0; i < range_dimension; ++i)
			result[i] = functions_.at(i)(space_vector);

		return result;
	}

	bool operator==(const Vector_Function& other) const {
		return this->functions_ == other.functions_;
	}

	Polynomial<domain_dimension>& operator[](const size_t index) {
		dynamic_require(index < this->range_dimension(), "index can not exceed range dimension");
		return this->functions_[index];
	}

	const Polynomial<domain_dimension>& at(const size_t index) const {
		dynamic_require(index < this->range_dimension(), "index can not exceed range dimension");
		return this->functions_[index];
	}


	size_t range_dimension(void) const {
		return this->functions_.size();
	}
};