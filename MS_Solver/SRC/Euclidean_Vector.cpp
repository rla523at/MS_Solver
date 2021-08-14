#include "../INC/Euclidean_Vector.h"

double Dynamic_Euclidean_Vector::operator[](const size_t position) const {
	dynamic_require(position < this->dimension(), "position should be less then dimension");
	return this->values_[position];
}

bool Dynamic_Euclidean_Vector::operator==(const Dynamic_Euclidean_Vector& other) const {
	return this->values_ == other.values_;
}

double Dynamic_Euclidean_Vector::at(const size_t position) const {
	dynamic_require(position < this->dimension(), "position should be less then dimension");
	return this->values_[position];
}

std::vector<double>::const_iterator Dynamic_Euclidean_Vector::begin(void) const {
	return this->values_.cbegin();
}

std::vector<double>::const_iterator Dynamic_Euclidean_Vector::end(void) const {
	return this->values_.cend();
}

const double* Dynamic_Euclidean_Vector::data(void) const {
	return this->values_.data();
};

size_t Dynamic_Euclidean_Vector::dimension(void) const {
	return values_.size();
}

std::string Dynamic_Euclidean_Vector::to_string(void) const {
	std::string result;
	for (const auto& element : this->values_)
		result += ms::double_to_string(element) + " ";
	result.pop_back();
	return result;
}

double* Dynamic_Euclidean_Vector::data(void) {
	return this->values_.data(); 
};

std::ostream& operator<<(std::ostream& os, const Dynamic_Euclidean_Vector& x) {
	return os << x.to_string();
}
