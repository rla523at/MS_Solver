#include "../INC/Euclidean_Vector.h"

Euclidean_Vector_Constant_Base::Euclidean_Vector_Constant_Base(const size_t num_values, const double* const_data_ptr)
{
	REQUIRE((0 < num_values) && (num_values <= std::numeric_limits<int>::max()), "number of values should be in range");
	REQUIRE(const_data_ptr != nullptr, "data ptr should not be nullptr");
	this->num_values_ = static_cast<int>(num_values);
	this->const_data_ptr_ = const_data_ptr;
}

Euclidean_Vector Euclidean_Vector_Constant_Base::operator*(const double constant) const
{
	Euclidean_Vector result(this->const_data_ptr_, this->const_data_ptr_ + this->num_values_);
	
	if (this->num_values_ <= ms::blas_dscal_criteria)
	{
		for (int i = 0; i < this->num_values_; ++i)
		{
			result.data_ptr_[i] *= constant;
		}
	}
	else
	{
		const auto n = this->num_values_;
		const auto incx = 1;
		cblas_dscal(n, constant, result.data_ptr_, incx);
	}

	return result;
}

Euclidean_Vector Euclidean_Vector_Constant_Base::operator+(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");
	
	Euclidean_Vector result(this->const_data_ptr_, this->const_data_ptr_ + this->num_values_);

	if (this->num_values_ <= ms::blas_axpy_criteria)
	{
		for (int i = 0; i < this->num_values_; ++i)
		{
			result.data_ptr_[i] += other.const_data_ptr_[i];
		}
	}
	else
	{
		const auto n = this->num_values_;
		const auto a = 1.0;
		const auto incx = 1;
		const auto incy = 1;

		cblas_daxpy(n, a, other.const_data_ptr_, incx, result.data_ptr_, incy);
	}

	return result;
}

Euclidean_Vector Euclidean_Vector_Constant_Base::operator-(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");

	Euclidean_Vector result(this->const_data_ptr_, this->const_data_ptr_ + this->num_values_);

	if (this->num_values_ <= ms::blas_axpy_criteria)
	{
		for (int i = 0; i < this->num_values_; ++i)
		{
			result.data_ptr_[i] -= other.const_data_ptr_[i];
		}
	}
	else
	{
		const auto n = this->num_values_;
		const auto a = -1.0;
		const auto incx = 1;
		const auto incy = 1;

		cblas_daxpy(n, a, other.const_data_ptr_, incx, result.data_ptr_, incy);
	}

	return result;
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
	double result = 0.0;

	if (this->num_values_ <= ms::blas_dasum_criteria) 
	{
		for (size_t i = 0; i < this->num_values_; ++i)
		{
			result += std::abs(this->const_data_ptr_[i]);
		}
	}
	else
	{
		const auto n = this->num_values_;
		const auto incx = 1;

		result = cblas_dasum(n, this->const_data_ptr_, incx);
	}

	return result;
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
	
	double result = 0.0;

	if (this->num_values_ <= ms::blas_dot_criteria)
	{
		for (int i = 0; i < this->num_values_; ++i)
		{
			result += this->const_data_ptr_[i] * other.const_data_ptr_[i];
		}
	}
	else
	{
		const auto n = this->num_values_;
		const auto incx = 1;
		const auto incy = 1;

		result = cblas_ddot(n, this->const_data_ptr_, incx, other.const_data_ptr_, incy);
	}

	return result;
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
	if (this->num_values_ <= ms::blas_dscal_criteria)
	{
		for (int i = 0; i < this->num_values_; ++i)
		{
			this->data_ptr_[i] *= constant;
		}
	}
	else
	{
		const auto n = this->num_values_;
		const auto incx = 1;

		cblas_dscal(n, constant, this->data_ptr_, incx);
	}
}

void Euclidean_Vector_Base::operator+=(const Euclidean_Vector_Constant_Base& other)
{
	REQUIRE(this->num_values_ == other.num_values_, "other vector should be same size");

	const double* other_data_ptr = other.data();
	if (this->num_values_ <= ms::blas_axpy_criteria)
	{
		for (int i = 0; i < this->num_values_; ++i)
		{
			this->data_ptr_[i] += other_data_ptr[i];
		}
	}
	else
	{
		const auto n = this->num_values_;
		const auto a = 1.0;
		const auto incx = 1;
		const auto incy = 1;

		cblas_daxpy(n, a, other_data_ptr, incx, this->data_ptr_, incy);
	}
}

void Euclidean_Vector_Base::normalize(void)
{
	*this *= 1.0 / this->L2_norm();
}

Euclidean_Vector::Euclidean_Vector(const size_t size)
{
	this->num_values_ = static_cast<int>(size);
	if (this->num_values_ <= this->small_criterion_)
	{
		this->const_data_ptr_ = this->small_buffer_.data();
		this->data_ptr_ = this->small_buffer_.data();
	}
	else
	{
		this->values_.resize(num_values_);
		this->const_data_ptr_ = this->values_.data();
		this->data_ptr_ = this->values_.data();
	}
}

Euclidean_Vector::Euclidean_Vector(const std::initializer_list<double> list)
{
	this->num_values_ = static_cast<int>(list.size());

	if (this->num_values_ <= this->small_criterion_)
	{
		std::copy(list.begin(), list.end(), this->small_buffer_.begin());
		this->const_data_ptr_ = this->small_buffer_.data();
		this->data_ptr_ = this->small_buffer_.data();
	}
	else
	{
		this->values_ = list;
		this->const_data_ptr_ = this->values_.data();
		this->data_ptr_ = this->values_.data();
	}
}

Euclidean_Vector::Euclidean_Vector(std::vector<double>&& values)
{
	this->num_values_ = static_cast<int>(values.size());

	if (this->num_values_ <= this->small_criterion_)
	{
		std::copy(values.begin(), values.end(), this->small_buffer_.begin());
		this->const_data_ptr_ = this->small_buffer_.data();
		this->data_ptr_ = this->small_buffer_.data();
	}
	else
	{
		this->values_ = std::move(values);
		this->const_data_ptr_ = this->values_.data();
		this->data_ptr_ = this->values_.data();
	}
};

Euclidean_Vector::Euclidean_Vector(const std::vector<double>& values)
{
	this->num_values_ = static_cast<int>(values.size());

	if (this->num_values_ <= this->small_criterion_)
	{
		std::copy(values.begin(), values.end(), this->small_buffer_.begin());
		this->const_data_ptr_ = this->small_buffer_.data();
		this->data_ptr_ = this->small_buffer_.data();
	}
	else
	{
		this->values_ = values;
		this->const_data_ptr_ = this->values_.data();
		this->data_ptr_ = this->values_.data();
	}
}

Euclidean_Vector::Euclidean_Vector(const Euclidean_Vector& other)
{
	this->num_values_ = other.num_values_;
	
	if (this->num_values_ <= this->small_criterion_)
	{
		this->small_buffer_ = other.small_buffer_;
		this->const_data_ptr_ = this->small_buffer_.data();
		this->data_ptr_ = this->small_buffer_.data();
	}
	else
	{
		this->values_ = other.values_;
		this->const_data_ptr_ = this->values_.data();
		this->data_ptr_ = this->values_.data();
	}
}

Euclidean_Vector::Euclidean_Vector(Euclidean_Vector&& other) noexcept
{
	this->num_values_ = other.num_values_;

	if (this->num_values_ <= this->small_criterion_)
	{
		this->small_buffer_ = other.small_buffer_;
		this->const_data_ptr_ = this->small_buffer_.data();
		this->data_ptr_ = this->small_buffer_.data();
	}
	else
	{
		this->values_ = std::move(other.values_);
		this->const_data_ptr_ = this->values_.data();
		this->data_ptr_ = this->values_.data();
	}

	other.num_values_ = 0;
	other.const_data_ptr_ = nullptr;
	other.data_ptr_ = nullptr;
}

void Euclidean_Vector::operator=(const Euclidean_Vector& other)
{
	this->num_values_ = other.num_values_;

	if (this->num_values_ <= this->small_criterion_)
	{
		std::copy(other.small_buffer_.begin(), other.small_buffer_.begin() + this->num_values_, this->small_buffer_.begin());

		this->const_data_ptr_ = this->small_buffer_.data();
		this->data_ptr_ = this->small_buffer_.data();
	}
	else
	{
		this->values_ = other.values_;
		this->const_data_ptr_ = this->values_.data();
		this->data_ptr_ = this->values_.data();
	}
}

void Euclidean_Vector::operator=(Euclidean_Vector&& other) noexcept
{
	this->num_values_ = other.num_values_;

	if (this->num_values_ <= this->small_criterion_)
	{
		std::copy(other.small_buffer_.begin(), other.small_buffer_.begin() + this->num_values_, this->small_buffer_.begin());

		this->const_data_ptr_ = this->small_buffer_.data();
		this->data_ptr_ = this->small_buffer_.data();
	}
	else
	{
		this->values_ = std::move(other.values_);
		this->const_data_ptr_ = this->values_.data();
		this->data_ptr_ = this->values_.data();
	}
}

std::vector<double>&& Euclidean_Vector::move_values(void)
{
	this->num_values_ = 0;
	this->const_data_ptr_ = nullptr;
	this->data_ptr_ = nullptr;

	//if (this->num_values_ <= )
	return std::move(this->values_);
}

bool Euclidean_Vector::is_small(void) const
{
	return this->values_.empty();
}



Euclidean_Vector Euclidean_Vector_Constant_Wrapper::operator-(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_ - other;
}
Euclidean_Vector Euclidean_Vector_Constant_Wrapper::operator+(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_ + other;
}
Euclidean_Vector Euclidean_Vector_Constant_Wrapper::operator*(const double constant) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_ * constant;
}
double Euclidean_Vector_Constant_Wrapper::operator[](const size_t position) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_[position];
}
bool Euclidean_Vector_Constant_Wrapper::operator==(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_ == other;
}

double Euclidean_Vector_Constant_Wrapper::at(const size_t position) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.at(position);
}

const double* Euclidean_Vector_Constant_Wrapper::begin(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.begin();
}

std::vector<double> Euclidean_Vector_Constant_Wrapper::copy_values(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.copy_values();
}

const double* Euclidean_Vector_Constant_Wrapper::data(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.data();
}

const double* Euclidean_Vector_Constant_Wrapper::end(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.data();
}

double Euclidean_Vector_Constant_Wrapper::L1_norm(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.L1_norm();
}

double Euclidean_Vector_Constant_Wrapper::L2_norm(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.L2_norm();
}

double Euclidean_Vector_Constant_Wrapper::Linf_norm(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.Linf_norm();
}

double Euclidean_Vector_Constant_Wrapper::inner_product(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.inner_product(other);
}

bool Euclidean_Vector_Constant_Wrapper::is_axis_translation(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.is_axis_translation(other);
}

size_t Euclidean_Vector_Constant_Wrapper::size(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.size();
}

std::string Euclidean_Vector_Constant_Wrapper::to_string(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.to_string();
}

bool Euclidean_Vector_Constant_Wrapper::is_sync(void) const
{
	return this->base_.data() == this->values_wrapper_.data() && this->base_.size() == this->values_wrapper_.size();
}

void Euclidean_Vector_Wrapper::operator*=(const double constant)
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	this->base_ *= constant;
}

void Euclidean_Vector_Wrapper::operator+=(const Euclidean_Vector_Constant_Base& other)
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	this->base_ += other;
}

void Euclidean_Vector_Wrapper::operator=(Euclidean_Vector&& other) noexcept
{
	this->num_values_ = static_cast<int>(other.size());

	if (other.is_small())
	{
		std::copy(other.begin(), other.end(), this->values_wrapper_.begin());
	}
	else
	{
		this->values_wrapper_ = std::move(other.move_values());
	}

	this->base_ = Euclidean_Vector_Base(this->values_wrapper_.size(), this->values_wrapper_.data());

	this->const_data_ptr_ = this->values_wrapper_.data();
	this->data_ptr_ = this->values_wrapper_.data();
}

void Euclidean_Vector_Wrapper::normalize(void)
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	this->base_.normalize();
}

Euclidean_Vector Euclidean_Vector_Wrapper::operator-(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_ - other;
}
Euclidean_Vector Euclidean_Vector_Wrapper::operator+(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_ + other;
}
Euclidean_Vector Euclidean_Vector_Wrapper::operator*(const double constant) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_ * constant;
}
double Euclidean_Vector_Wrapper::operator[](const size_t position) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_[position];
}
bool Euclidean_Vector_Wrapper::operator==(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_ == other;
}

double Euclidean_Vector_Wrapper::at(const size_t position) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.at(position);
}

const double* Euclidean_Vector_Wrapper::begin(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.begin();
}

std::vector<double> Euclidean_Vector_Wrapper::copy_values(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.copy_values();
}

const double* Euclidean_Vector_Wrapper::data(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.data();
}

const double* Euclidean_Vector_Wrapper::end(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.data();
}

double Euclidean_Vector_Wrapper::L1_norm(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.L1_norm();
}

double Euclidean_Vector_Wrapper::L2_norm(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.L2_norm();
}

double Euclidean_Vector_Wrapper::Linf_norm(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.Linf_norm();
}

double Euclidean_Vector_Wrapper::inner_product(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.inner_product(other);
}

bool Euclidean_Vector_Wrapper::is_axis_translation(const Euclidean_Vector_Constant_Base& other) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.is_axis_translation(other);
}

size_t Euclidean_Vector_Wrapper::size(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.size();
}

std::string Euclidean_Vector_Wrapper::to_string(void) const
{
	REQUIRE(this->is_sync(), "wrapper should be sync");
	return this->base_.to_string();
}

bool Euclidean_Vector_Wrapper::is_sync(void) const
{
	return this->base_.data() == this->values_wrapper_.data() && this->base_.size() == this->values_wrapper_.size();
}









Euclidean_Vector operator*(const double constant, const Euclidean_Vector_Constant_Base& x)
{
	return x * constant;
}

std::ostream& operator<<(std::ostream& os, const Euclidean_Vector_Constant_Base& x)
{
	return os << x.to_string();
}
