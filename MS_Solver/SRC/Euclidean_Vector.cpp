#include "../INC/Euclidean_Vector.h"

double Euclidean_Vector::operator[](const size_t position) const {
	REQUIRE(position < this->dimension(), "position should be less then dimension");
	return this->values_[position];
}
bool Euclidean_Vector::operator==(const Euclidean_Vector& other) const {
	return this->values_ == other.values_;
}
size_t Euclidean_Vector::dimension(void) const {
	return values_.size();
}
std::string Euclidean_Vector::to_string(void) const {
	std::ostringstream oss;
	oss << std::setprecision(16) << std::showpoint << std::left;
	for (const auto value : this->values_)
		oss << std::setw(25) << value;
	return oss.str();
}




std::ostream& operator<<(std::ostream& os, const Euclidean_Vector& x) {
	return os << x.to_string();
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

