#include "../INC/Euclidean_Vector.h"

double Euclidean_Vector::operator[](const size_t position) const {
	REQUIRE(position < this->dimension(), "position should be less then dimension");
	return this->values_[position];
}
size_t Euclidean_Vector::dimension(void) const {
	return values_.size();
}

//
//Dynamic_Euclidean_Vector& Dynamic_Euclidean_Vector::operator-=(const Dynamic_Euclidean_Vector& other) {
//	dynamic_require(this->dimension() == other.dimension(), "dimension should be matched");
//
//	cblas_daxpy(static_cast<MKL_INT>(this->dimension()), -1.0, other.values_.data(), 1, this->values_.data(), 1);
//	return *this;
//}
//
//Dynamic_Euclidean_Vector& Dynamic_Euclidean_Vector::operator*=(const double constant) {
//	cblas_dscal(static_cast<MKL_INT>(this->dimension()), constant, this->values_.data(), 1);
//	return *this;
//}
//
//Dynamic_Euclidean_Vector Dynamic_Euclidean_Vector::operator-(const Dynamic_Euclidean_Vector& other) const {
//	auto result = *this;
//	return result -= other;
//}
//
//double Dynamic_Euclidean_Vector::operator[](const size_t position) const {
//	dynamic_require(position < this->dimension(), "position should be less then dimension");
//	return this->values_[position];
//}
//
//bool Dynamic_Euclidean_Vector::operator==(const Dynamic_Euclidean_Vector& other) const {
//	return this->values_ == other.values_;
//}
//
//double Dynamic_Euclidean_Vector::at(const size_t position) const {
//	dynamic_require(position < this->dimension(), "position should be less then dimension");
//	return this->values_[position];
//}
//
//std::vector<double>::const_iterator Dynamic_Euclidean_Vector::begin(void) const {
//	return this->values_.cbegin();
//}
//
//std::vector<double>::const_iterator Dynamic_Euclidean_Vector::end(void) const {
//	return this->values_.cend();
//}
//
//const double* Dynamic_Euclidean_Vector::data(void) const {
//	return this->values_.data();
//};
//
//size_t Dynamic_Euclidean_Vector::dimension(void) const {
//	return values_.size();
//}
//
//std::string Dynamic_Euclidean_Vector::to_string(void) const {
//	std::string result;
//	for (const auto& element : this->values_)
//		result += ms::double_to_string(element) + " ";
//	result.pop_back();
//	return result;
//}
//
//double Dynamic_Euclidean_Vector::inner_product(const Dynamic_Euclidean_Vector& other) const {
//	dynamic_require(this->dimension() == other.dimension(), "dimension should be matched");
//
//	return cblas_ddot(static_cast<MKL_INT>(this->dimension()), this->values_.data(), 1, other.values_.data(), 1);
//}
//
//
//void Dynamic_Euclidean_Vector::be_absolute(void){
//	for (auto& value : this->values_)
//		value = std::abs(value);	
//};
//
//std::ostream& operator<<(std::ostream& os, const Dynamic_Euclidean_Vector& x) {
//	return os << x.to_string();
//}
