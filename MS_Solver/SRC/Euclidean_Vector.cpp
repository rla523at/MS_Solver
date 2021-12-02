#include "../INC/Euclidean_Vector.h"


Euclidean_Vector Euclidean_Vector_Constant_Base::operator*(const double constant) const
{
	std::vector<double> values = { this->const_data_ptr_, this->const_data_ptr_ + this->num_values_ };

	const auto n = this->num_values_;
	const auto incx = 1;
	cblas_dscal(n, constant, values.data(), incx);
	return values;
}

Euclidean_Vector Euclidean_Vector_Constant_Base::operator+(const Euclidean_Vector_Constant_Base& other) const
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

Euclidean_Vector Euclidean_Vector_Constant_Base::operator-(const Euclidean_Vector_Constant_Base& other) const
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

double Euclidean_Vector_Constant_Base::operator[](const size_t position) const
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->const_data_ptr_[position];
}

bool Euclidean_Vector_Constant_Base::operator==(const Euclidean_Vector_Constant_Base& other) const
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

double Euclidean_Vector_Constant_Base::at(const size_t position) const
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->const_data_ptr_[position];
}

const double* Euclidean_Vector_Constant_Base::begin(void) const
{
	return this->const_data_ptr_;
}

std::vector<double> Euclidean_Vector_Constant_Base::copy_values(void) const
{
	return { this->const_data_ptr_, this->const_data_ptr_ + this->num_values_ };
}

const double* Euclidean_Vector_Constant_Base::data(void) const
{
	return this->const_data_ptr_;
}

const double* Euclidean_Vector_Constant_Base::end(void) const
{
	return this->const_data_ptr_ + this->num_values_;
}


double Euclidean_Vector_Constant_Base::L1_norm(void) const
{
	const auto n = this->num_values_;
	const auto incx = 1;

	return cblas_dasum(n, this->const_data_ptr_, incx);
}

double Euclidean_Vector_Constant_Base::L2_norm(void) const
{
	return std::sqrt(this->inner_product(*this));
}

double Euclidean_Vector_Constant_Base::Linf_norm(void) const
{
	const auto n = this->num_values_;
	const auto incx = 1;
	const auto pos = cblas_idamax(n, this->const_data_ptr_, incx);

	return this->const_data_ptr_[pos];
}

double Euclidean_Vector_Constant_Base::inner_product(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");

	const auto n = this->num_values_;
	const auto incx = 1;
	const auto incy = 1;

	return cblas_ddot(n, this->const_data_ptr_, incx, other.const_data_ptr_, incy);
}

bool Euclidean_Vector_Constant_Base::is_axis_translation(const Euclidean_Vector_Constant_Base& other) const
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

size_t Euclidean_Vector_Constant_Base::size(void) const
{
	return this->num_values_;
}

std::string Euclidean_Vector_Constant_Base::to_string(void) const
{
	std::ostringstream oss;
	oss << std::setprecision(16) << std::showpoint << std::left;
	for (int i = 0; i < this->num_values_; ++i)
	{
		oss << std::setw(25) << this->const_data_ptr_[i];
	}
	return oss.str();
}

void Euclidean_Vector_Base::operator*=(const double constant)
{
	const auto n = this->num_values_;
	const auto incx = 1;

	cblas_dscal(n, constant, this->data_ptr_, incx);
}

void Euclidean_Vector_Base::operator+=(const Euclidean_Vector_Base& other)
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");

	const auto n = this->num_values_;
	const auto a = 1.0;
	const auto incx = 1;
	const auto incy = 1;

	cblas_daxpy(n, a, other.const_data_ptr_, incx, this->data_ptr_, incy);
}

void Euclidean_Vector_Base::normalize(void)
{
	*this *= 1.0 / this->L2_norm();
}

Euclidean_Vector::Euclidean_Vector(const size_t size)
	: values_(size)
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(const std::initializer_list<double> list)
	: values_(list)
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(std::vector<double>&& values)
	: values_(std::move(values))
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
};

Euclidean_Vector::Euclidean_Vector(const std::vector<double>& values)
	: values_(values)
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(const Euclidean_Vector& other)
{
	this->num_values_ = other.num_values_;
	this->values_ = other.values_;
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(Euclidean_Vector&& other) noexcept
{
	this->num_values_ = other.num_values_;
	this->values_ = std::move(other.values_);
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

void Euclidean_Vector::operator=(const Euclidean_Vector& other)
{
	this->num_values_ = other.num_values_;
	this->values_ = other.values_;
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

void Euclidean_Vector::operator=(Euclidean_Vector&& other) noexcept
{
	this->num_values_ = other.num_values_;
	this->values_ = std::move(other.values_);
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

std::vector<double>&& Euclidean_Vector::move_values(void)
{
	this->num_values_ = 0;
	this->const_data_ptr_ = nullptr;
	this->data_ptr_ = nullptr;
	return std::move(this->values_);
}

Euclidean_Vector_Constant_Wrapper::Euclidean_Vector_Constant_Wrapper(const std::vector<double>& values)
	:values_constant_wrapper_(values)
{
	this->num_values_ = static_cast<int>(values.size());
	this->const_data_ptr_ = values.data();
}

Euclidean_Vector_Wrapper::Euclidean_Vector_Wrapper(std::vector<double>& values)
	:values_wrapper_(values)
{
	this->num_values_ = static_cast<int>(this->values_wrapper_.size());
	this->const_data_ptr_ = this->values_wrapper_.data();
	this->data_ptr_ = this->values_wrapper_.data();
}

void Euclidean_Vector_Wrapper::operator=(Euclidean_Vector&& other) noexcept
{
	this->values_wrapper_ = std::move(other.move_values());
	this->num_values_ = static_cast<int>(this->values_wrapper_.size());
	this->const_data_ptr_ = this->values_wrapper_.data();
	this->data_ptr_ = this->values_wrapper_.data();
}

Euclidean_Vector operator*(const double constant, const Euclidean_Vector_Constant_Base& x)
{
	return x * constant;
}

std::ostream& operator<<(std::ostream& os, const Euclidean_Vector_Constant_Base& x)
{
	return os << x.to_string();
}
