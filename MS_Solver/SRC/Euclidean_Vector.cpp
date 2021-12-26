#include "../INC/Euclidean_Vector.h"

Constant_Euclidean_Vector_Wrapper::Constant_Euclidean_Vector_Wrapper(const size_t num_values, const double* const_data_ptr)
{
	REQUIRE((0 < num_values) && (num_values <= std::numeric_limits<int>::max()), "number of values should be in range");
	REQUIRE(const_data_ptr != nullptr, "data ptr should not be nullptr");
	this->num_values_ = static_cast<int>(num_values);
	this->const_data_ptr_ = const_data_ptr;
}

Constant_Euclidean_Vector_Wrapper::Constant_Euclidean_Vector_Wrapper(const std::vector<double>& vec)
{
	const auto num_values = vec.size();
	
	REQUIRE((0 < num_values) && (num_values <= std::numeric_limits<int>::max()), "number of values should be in range");
	this->num_values_ = static_cast<int>(num_values);
	this->const_data_ptr_ = vec.data();
}


Euclidean_Vector Constant_Euclidean_Vector_Wrapper::operator*(const double constant) const
{
	Euclidean_Vector result(this->const_data_ptr_, this->const_data_ptr_ + this->num_values_);
	
	ms::BLAS::cx(constant, this->num_values_, result.data());

	return result;
}

Euclidean_Vector Constant_Euclidean_Vector_Wrapper::operator+(const Constant_Euclidean_Vector_Wrapper& other) const
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");
	
	Euclidean_Vector result(this->const_data_ptr_, this->const_data_ptr_ + this->num_values_);
	ms::vpv(*this, other, result.data());

	return result;
}

Euclidean_Vector Constant_Euclidean_Vector_Wrapper::operator-(const Constant_Euclidean_Vector_Wrapper& other) const
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");

	Euclidean_Vector result(this->const_data_ptr_, this->const_data_ptr_ + this->num_values_);
	ms::vmv(*this, other, result.data());

	return result;
}

double Constant_Euclidean_Vector_Wrapper::operator[](const size_t position) const
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->const_data_ptr_[position];
}

bool Constant_Euclidean_Vector_Wrapper::operator==(const Constant_Euclidean_Vector_Wrapper& other) const
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

double Constant_Euclidean_Vector_Wrapper::at(const size_t position) const
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->const_data_ptr_[position];
}

const double* Constant_Euclidean_Vector_Wrapper::begin(void) const
{
	return this->const_data_ptr_;
}

std::vector<double> Constant_Euclidean_Vector_Wrapper::copy_values(void) const
{
	return { this->const_data_ptr_, this->const_data_ptr_ + this->num_values_ };
}

const double* Constant_Euclidean_Vector_Wrapper::data(void) const
{
	return this->const_data_ptr_;
}

const double* Constant_Euclidean_Vector_Wrapper::end(void) const
{
	return this->const_data_ptr_ + this->num_values_;
}


double Constant_Euclidean_Vector_Wrapper::L1_norm(void) const
{
	return ms::BLAS::abs_x(this->num_values_, this->const_data_ptr_);
}

double Constant_Euclidean_Vector_Wrapper::L2_norm(void) const
{
	return std::sqrt(this->inner_product(*this));
}

double Constant_Euclidean_Vector_Wrapper::Linf_norm(void) const
{
	const auto n = this->num_values_;
	const auto incx = 1;
	const auto pos = cblas_idamax(n, this->const_data_ptr_, incx);

	return this->const_data_ptr_[pos];
}

double Constant_Euclidean_Vector_Wrapper::inner_product(const Constant_Euclidean_Vector_Wrapper& other) const
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");
	
	return ms::BLAS::x_dot_y(this->num_values_, this->const_data_ptr_, other.const_data_ptr_);
}

bool Constant_Euclidean_Vector_Wrapper::is_axis_translation(const Constant_Euclidean_Vector_Wrapper& other) const
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

size_t Constant_Euclidean_Vector_Wrapper::size(void) const
{
	return this->num_values_;
}

std::string Constant_Euclidean_Vector_Wrapper::to_string(void) const
{
	std::ostringstream oss;
	oss << std::setprecision(16) << std::showpoint << std::left;
	for (int i = 0; i < this->num_values_; ++i)
	{
		oss << std::setw(25) << this->const_data_ptr_[i];
	}
	return oss.str();
}

void Euclidean_Vector_Wrapper::operator*=(const double constant)
{
	ms::BLAS::cx(constant, this->num_values_, this->data_ptr_);
}

void Euclidean_Vector_Wrapper::operator+=(const Constant_Euclidean_Vector_Wrapper& other)
{
	REQUIRE(this->num_values_ == other.size(), "other vector should be same size");

	ms::BLAS::x_plus_assign_y(this->num_values_, this->data_ptr_, other.data());
}

void Euclidean_Vector_Wrapper::operator-=(const Constant_Euclidean_Vector_Wrapper& other)
{
	REQUIRE(this->num_values_ == other.size(), "other vector should be same size");

	ms::BLAS::x_minus_assign_y(this->num_values_, this->data_ptr_, other.data());
}

double& Euclidean_Vector_Wrapper::operator[](const size_t position)
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->data_ptr_[position];
}

double& Euclidean_Vector_Wrapper::at(const size_t position)
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->data_ptr_[position];
}

double* Euclidean_Vector_Wrapper::data(void)
{
	return this->data_ptr_;
}

void Euclidean_Vector_Wrapper::initalize(void)
{
	for (int i = 0; i < this->num_values_; ++i)
	{
		this->data_ptr_[i] = 0.0;
	}
}

void Euclidean_Vector_Wrapper::normalize(void)
{
	*this *= 1.0 / this->L2_norm();
}

double Euclidean_Vector_Wrapper::operator[](const size_t position) const
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->const_data_ptr_[position];
}

double Euclidean_Vector_Wrapper::at(const size_t position) const
{
	REQUIRE(position < this->num_values_, "position should be less then size");
	return this->const_data_ptr_[position];
}

const double* Euclidean_Vector_Wrapper::data(void) const
{
	return this->const_data_ptr_;
}

Euclidean_Vector::Euclidean_Vector(const size_t size)
	:values_(size)
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(const std::initializer_list<double> list)
	:values_(list)
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(std::vector<double>&& values)
	:values_(std::move(values))
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
};

Euclidean_Vector::Euclidean_Vector(const std::vector<double>& values)
	:values_(values)
{
	this->num_values_ = static_cast<int>(this->values_.size());
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(const Euclidean_Vector& other)
	:values_(other.values_)
{
	this->num_values_ = other.num_values_;
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

Euclidean_Vector::Euclidean_Vector(Euclidean_Vector&& other) noexcept
	:values_(std::move(other.values_))
{
	this->num_values_ = other.num_values_;
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();

	other.num_values_ = 0;
	other.const_data_ptr_ = nullptr;
	other.data_ptr_ = nullptr;
}

void Euclidean_Vector::operator=(const Euclidean_Vector& other)
{
	this->values_ = other.values_;
	this->num_values_ = other.num_values_;
	this->const_data_ptr_ = this->values_.data();
	this->data_ptr_ = this->values_.data();
}

void Euclidean_Vector::operator=(Euclidean_Vector&& other) noexcept
{
	this->values_ = std::move(other.values_);
	this->num_values_ = other.num_values_;
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

Euclidean_Vector operator*(const double constant, const Constant_Euclidean_Vector_Wrapper& x)
{
	return x * constant;
}

std::ostream& operator<<(std::ostream& os, const Constant_Euclidean_Vector_Wrapper& x)
{
	return os << x.to_string();
}

namespace ms
{
	void vpv(const Constant_Euclidean_Vector_Wrapper& v1, const Constant_Euclidean_Vector_Wrapper& v2, double* result_ptr)
	{
		const auto n = static_cast<int>(v1.size());
		REQUIRE(n == v2.size(), "size should be same");
		ms::BLAS::x_plus_y(n, v1.data(), v2.data(), result_ptr);
	}

	void vmv(const Constant_Euclidean_Vector_Wrapper& v1, const Constant_Euclidean_Vector_Wrapper& v2, double* result_ptr)
	{
		const auto n = static_cast<int>(v1.size());
		REQUIRE(n == v2.size(), "size should be same");
		ms::BLAS::x_minus_y(n, v1.data(), v2.data(), result_ptr);
	}
}