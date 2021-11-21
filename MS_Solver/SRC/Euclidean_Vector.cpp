#include "../INC/Euclidean_Vector.h"

Euclidean_Vector::Euclidean_Vector(const size_t size)
	: values_(size)
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(const std::initializer_list<double> list)
	: values_(list)
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(std::vector<double>&& values)
	: values_(std::move(values))
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
};

Euclidean_Vector::Euclidean_Vector(const std::vector<double>& values)
	: values_(values)
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(const Euclidean_Vector& other)
{
	this->num_values_ = other.num_values_;
	this->values_ = other.values_;
	this->const_data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(Euclidean_Vector&& other) noexcept
{
	this->num_values_ = other.num_values_;
	this->values_ = std::move(other.values_);
	this->const_data_ptr_ = this->values_.data();
}

void Euclidean_Vector::operator=(const Euclidean_Vector& other)
{
	this->num_values_ = other.num_values_;
	this->values_ = other.values_;
	this->const_data_ptr_ = this->values_.data();
}

void Euclidean_Vector::operator=(Euclidean_Vector&& other) noexcept
{
	this->num_values_ = other.num_values_;
	this->values_ = std::move(other.values_);
	this->const_data_ptr_ = this->values_.data();
}

Euclidean_Vector Euclidean_Vector_Base::operator*(const double constant) const
{
	std::vector<double> values = { this->const_data_ptr_, this->const_data_ptr_ + this->num_values_ };
	
	const auto n = this->num_values_;
	const auto incx = 1;
	cblas_dscal(n, constant, values.data(), incx);
	return values;
}

Euclidean_Vector Euclidean_Vector_Base::operator+(const Euclidean_Vector_Base& other) const
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");

	std::vector<double> values = { this->const_data_ptr_, this->const_data_ptr_ + this->num_values_ };

	const auto n = this->num_values_;
	const auto a = 1.0;
	const auto incx = 1;
	const auto incy = 1;

	cblas_daxpy(n, a, other.const_data_ptr_, incx, values.data(), incy);
	
	return values;
}

Euclidean_Vector Euclidean_Vector_Base::operator-(const Euclidean_Vector_Base& other) const
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");

	std::vector<double> values = { this->const_data_ptr_, this->const_data_ptr_ + this->num_values_ };

	const auto n = this->num_values_;
	const auto a = -1.0;
	const auto incx = 1;
	const auto incy = 1;

	cblas_daxpy(n, a, other.const_data_ptr_, incx, values.data(), incy);

	return values;
}

double Euclidean_Vector_Base::operator[](const size_t position) const
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->const_data_ptr_[position];
}

bool Euclidean_Vector_Base::operator==(const Euclidean_Vector_Base& other) const
{
	if (this->num_values_ != other.num_values_)
		return false;
	
	for (int i = 0; i < this->num_values_; ++i)
	{
		if (this->const_data_ptr_[i] != other.const_data_ptr_[i])
			return false;
	}

	return true;
}

double Euclidean_Vector_Base::at(const size_t position) const
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->const_data_ptr_[position];
}

const double* Euclidean_Vector_Base::begin(void) const
{
	return this->const_data_ptr_;
}

const double* Euclidean_Vector_Base::data(void) const
{
	return this->const_data_ptr_;
}

double Euclidean_Vector_Base::L1_norm(void) const
{
	const auto n = this->num_values_;
	const auto incx = 1;

	return cblas_dasum(n, this->const_data_ptr_, incx);
}

double Euclidean_Vector_Base::L2_norm(void) const
{
	return std::sqrt(this->inner_product(*this));
}

double Euclidean_Vector_Base::inner_product(const Euclidean_Vector_Base& other) const
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");

	const auto n = this->num_values_;
	const auto incx = 1;
	const auto incy = 1;

	return cblas_ddot(n, this->const_data_ptr_, incx, other.const_data_ptr_, incy);
}

bool Euclidean_Vector_Base::is_axis_translation(const Euclidean_Vector_Base& other) const
{
	const auto line_vector = *this - other;
	const auto L1_norm = line_vector.L1_norm();
	const auto L2_norm = line_vector.L2_norm();

	constexpr auto epsilon = 1.0E-10;
	if (std::abs(L1_norm - L2_norm) <= epsilon)
		return true;
	else
		return false;
}

size_t Euclidean_Vector_Base::size(void) const
{
	return this->num_values_;
}

std::string Euclidean_Vector_Base::to_string(void) const
{
	std::ostringstream oss;
	oss << std::setprecision(16) << std::showpoint << std::left;
	for (int i = 0; i < this->num_values_; ++i)
		oss << std::setw(25) << this->const_data_ptr_[i];
	return oss.str();
}

Euclidean_Vector& Euclidean_Vector::operator*=(const double constant)
{
	cblas_dscal(static_cast<MKL_INT>(this->size()), constant, this->values_.data(), 1);
	return *this;
}

Euclidean_Vector& Euclidean_Vector::operator+=(const Euclidean_Vector& other)
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");

	const auto n = this->num_values_;
	const auto a = 1.0;
	const auto incx = 1;
	const auto incy = 1;

	cblas_daxpy(n, a, other.const_data_ptr_, incx, this->values_.data(), incy);

	return *this;
}

//Euclidean_Vector& Euclidean_Vector::operator+=(const Euclidean_Vector& other)
//{
//	REQUIRE(this->size() == other.size(), "other vector should be same size");
//
//	const auto n = static_cast<MKL_INT>(this->size());
//	const auto a = 1.0;
//	const auto incx = 1;
//	const auto incy = 1;
//
//	cblas_daxpy(n, a, other.values_.data(), incx, this->values_.data(), incy);
//	
//	return *this;
//}
//
//Euclidean_Vector& Euclidean_Vector::operator-=(const Euclidean_Vector& other)
//{
//	REQUIRE(this->size() == other.size(), "other vector should be same size");
//
//	const auto n = static_cast<MKL_INT>(this->size());
//	const auto a = -1.0;
//	const auto incx = 1;
//	const auto incy = 1;
//
//	cblas_daxpy(n, a, other.values_.data(), incx, this->values_.data(), incy);
//
//	return *this;
//}

Euclidean_Vector& Euclidean_Vector::normalize(void)
{
	return *this *= 1.0 / this->L2_norm();
}

//Euclidean_Vector Euclidean_Vector::operator+(const Euclidean_Vector& other) const
//{
//	auto result = *this;
//	return result += other;
//}
//
//Euclidean_Vector Euclidean_Vector::operator-(const Euclidean_Vector& other) const
//{
//	auto result = *this;
//	return result -= other;
//}
//
//Euclidean_Vector Euclidean_Vector::operator*(const double constant) const
//{
//	auto result = *this;
//	return result *= constant;	
//}

//double Euclidean_Vector::operator[](const size_t position) const 
//{
//	REQUIRE(position < this->size(), "position should be less then size");
//	return this->values_[position];
//}
//
//bool Euclidean_Vector::operator==(const Euclidean_Vector& other) const 
//{
//	return this->values_ == other.values_;
//}
//
//double Euclidean_Vector::at(const size_t position) const 
//{
//	REQUIRE(position < this->size(), "position should be less then size");
//	return this->values_[position];
//}
//
//const double* Euclidean_Vector::begin(void) const 
//{
//	return this->values_.data();
//}
//
//double Euclidean_Vector::L2_norm(void) const
//{
//	return std::sqrt(this->inner_product(*this));
//}
//
//double Euclidean_Vector::inner_product(const Euclidean_Vector& other) const
//{
//	REQUIRE(this->size() == other.size(), "other vector should be same size");
//
//	const auto n = static_cast<MKL_INT>(this->values_.size());
//	const auto incx = 1;
//	const auto incy = 1;
//
//	return cblas_ddot(n, this->values_.data(), incx, other.values_.data(), incy);
//}
//
//size_t Euclidean_Vector::size(void) const 
//{
//	return values_.size();
//}
//
//std::string Euclidean_Vector::to_string(void) const 
//{
//	std::ostringstream oss;
//	oss << std::setprecision(16) << std::showpoint << std::left;
//	for (const auto value : this->values_)
//		oss << std::setw(25) << value;
//	return oss.str();
//}


Euclidean_Vector_Wrapper::Euclidean_Vector_Wrapper(const size_t num_value, const double* ptr)
{
	this->num_values_ = static_cast<int>(num_value);
	this->const_data_ptr_ = ptr;
}

Euclidean_Vector_Wrapper::Euclidean_Vector_Wrapper(const std::vector<double>& values)
{
	this->num_values_ = static_cast<int>(values.size());
	this->const_data_ptr_ = values.data();
}



Euclidean_Vector operator*(const double constant, const Euclidean_Vector_Base& x)
{
	return x * constant;
}

std::ostream& operator<<(std::ostream& os, const Euclidean_Vector& x) 
{
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

