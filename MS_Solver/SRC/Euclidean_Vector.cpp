#include "../INC/EuclideanVector.h"

double& Dynamic_Euclidean_Vector_::operator[](const size_t position) {
	dynamic_require(position < this->dimension(), "position should be less then dimension");
	return this->vals_[position];
}

double Dynamic_Euclidean_Vector_::at(const size_t position) const {
	dynamic_require(position < this->dimension(), "position should be less then dimension");
	return this->vals_[position];
}

size_t Dynamic_Euclidean_Vector_::dimension(void) const {
	return vals_.size();
}

std::string Dynamic_Euclidean_Vector_::to_string(void) const {
	std::string result;
	for (const auto& element : this->vals_)
		result += ms::double_to_string(element) + " ";
	result.pop_back();
	return result;
}