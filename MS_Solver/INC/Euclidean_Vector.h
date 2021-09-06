#pragma once
#include <algorithm>
#include <array>
#include <iomanip>
#include <initializer_list>
#include <mkl.h>
#include <string>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>


#define static_require static_assert
#define dynamic_require(requirement, state) if (!(requirement)) throw std::runtime_error(state)


using ushort = unsigned short;


namespace ms {
	class BLAS;

	inline constexpr ushort blas_dasum_criteria = 10;
	inline constexpr ushort blas_axpy_criteria = 20;
	inline constexpr ushort blas_dot_criteria = 15;
	
	template <typename... Args>
	inline constexpr bool are_arithmetics = (... && std::is_arithmetic_v<Args>);
}


//Template Class Euclidean_Vector
template <size_t dimension_>
class Euclidean_Vector
{
private:
	std::array<double, dimension_> values_ = { 0 };

public:
	Euclidean_Vector(void) = default;
	Euclidean_Vector(const std::array<double, dimension_>& other) : values_(other) {};
	template <typename... Args>
	Euclidean_Vector(Args... args);

public:
	Euclidean_Vector& operator+=(const Euclidean_Vector& y);
	Euclidean_Vector& operator-=(const Euclidean_Vector& y);
	Euclidean_Vector& operator*=(const double scalar);

public:
	Euclidean_Vector operator+(const Euclidean_Vector& y) const;	
	Euclidean_Vector operator-(const Euclidean_Vector& y) const;
	Euclidean_Vector operator*(const double scalar) const;
	bool operator==(const Euclidean_Vector& y) const;
	double operator[](const size_t position) const;

public:
	Euclidean_Vector& be_normalize(void);

public:
	double at(const size_t position) const;	
	std::array<double, dimension_>::const_iterator cbegin(void) const;
	std::array<double, dimension_>::const_iterator cend(void) const;
	double inner_product(const Euclidean_Vector& y) const;
	double L1_norm(void) const;
	double L2_norm(void) const;
	bool is_axis_translation(const Euclidean_Vector& other) const;
	std::string to_string(void) const;

public:
	static constexpr size_t dimension(void);
};


//user-defined deduction guides
template <typename... Args>
Euclidean_Vector(Args... args)->Euclidean_Vector<sizeof...(Args)>;

template <size_t dimension_> 
std::ostream& operator<<(std::ostream& os, const Euclidean_Vector<dimension_>& x);

template <size_t dimension_>
Euclidean_Vector<dimension_> operator*(const double constant, const Euclidean_Vector<dimension_>& x);


class Dynamic_Euclidean_Vector
{
private:
	template <size_t dimension_>
	friend class Euclidean_Vector;
	friend class ms::BLAS;

private:
	std::vector<double> values_;

public:
	Dynamic_Euclidean_Vector(const size_t dimension) : values_(dimension) {};
	Dynamic_Euclidean_Vector(const std::initializer_list<double> list) : values_(list) {};
	Dynamic_Euclidean_Vector(std::vector<double>&& values) : values_(std::move(values)) {};
	Dynamic_Euclidean_Vector(const std::vector<double>& values) : values_(values) {};

public:
	Dynamic_Euclidean_Vector& operator-=(const Dynamic_Euclidean_Vector& other);

public:
	Dynamic_Euclidean_Vector operator-(const Dynamic_Euclidean_Vector& other) const;
	double operator[](const size_t position) const;
	bool operator==(const Dynamic_Euclidean_Vector& other) const;

public:
	template <typename Function>
	void apply(const Function& f);
	void be_absolute(void);

public:
	double at(const size_t position) const;
	std::vector<double>::const_iterator begin(void) const;
	std::vector<double>::const_iterator end(void) const;
	const double* data(void) const;
	size_t dimension(void) const;
	double inner_product(const Dynamic_Euclidean_Vector& other) const;
	std::string to_string(void) const;
};

std::ostream& operator<<(std::ostream& os, const Dynamic_Euclidean_Vector& x);


namespace ms {
	inline std::string double_to_string(const double val) {
		constexpr size_t precision = 16;
		std::stringstream stream;
		stream << std::setprecision(precision) << std::noshowpoint << val;
		return stream.str();
	}

	inline bool compare_double(const double d1, const double d2, const size_t ULP_precision = 4) {
		const auto lower_ULP = d1 - std::nextafter(d1, std::numeric_limits<double>::lowest());
		const auto upper_ULP = std::nextafter(d1, std::numeric_limits<double>::max()) - d1;

		return d1 - ULP_precision * lower_ULP <= d2 && d2 <= d1 + ULP_precision * upper_ULP;
	}
}


// Template Definition Part
template <size_t dimension_>
template <typename... Args>
Euclidean_Vector<dimension_>::Euclidean_Vector(Args... args) : values_{ static_cast<double>(args)... } {
	static_require(sizeof...(Args) <= dimension_, "Number of arguments can not exceed dimension");
	static_require(ms::are_arithmetics<Args...>, "every arguments should be arithmetics");
};


template <size_t dimension_>
Euclidean_Vector<dimension_>& Euclidean_Vector<dimension_>::operator+=(const Euclidean_Vector& y) {
	for (size_t i = 0; i < dimension_; ++i)
		this->values_[i] += y.values_[i];
	return *this;
}

template <size_t dimension_>
Euclidean_Vector<dimension_>& Euclidean_Vector<dimension_>::operator-=(const Euclidean_Vector& y) {
	for (size_t i = 0; i < dimension_; ++i)
		this->values_[i] -= y.values_[i];
	return *this;
}

template <size_t dimension_>
Euclidean_Vector<dimension_>& Euclidean_Vector<dimension_>::operator*=(const double scalar) {
	if constexpr (dimension_ < 10) {
		for (size_t i = 0; i < dimension_; ++i)
			this->values_[i] *= scalar;
	}
	else {
		cblas_dscal(dimension_, scalar, this->values_.data(), 1);
	}
	return *this;
}

template <size_t dimension_> 
Euclidean_Vector<dimension_> Euclidean_Vector<dimension_>::operator+(const Euclidean_Vector& y) const {
	auto result = *this;
	return result += y;
}

template <size_t dimension_> Euclidean_Vector<dimension_> Euclidean_Vector<dimension_>::operator-(const Euclidean_Vector& y) const {
	auto result = *this;
	return result -= y;
}

template <size_t dimension_> Euclidean_Vector<dimension_> Euclidean_Vector<dimension_>::operator*(const double scalar) const {
	auto result = *this;
	return result *= scalar;
}

template <size_t dimension_> 
bool Euclidean_Vector<dimension_>::operator==(const Euclidean_Vector& y) const {
	for (size_t i = 0; i < dimension_; ++i) {
		if (this->values_[i] != y.values_[i])
			return false;
	}
	return true;
}

//template <size_t dimension_>
//bool Euclidean_Vector<dimension_>::operator==(const Dynamic_Euclidean_Vector& y) const {
//	if (dimension_ != y.dimension())
//		return false;
//	
//	for (size_t i = 0; i < dimension_; ++i) {
//		if (this->values_[i] != y.values_[i])
//			return false;
//	}
//	return true;
//}

template <size_t dimension_> 
double Euclidean_Vector<dimension_>::operator[](const size_t position) const {
	dynamic_require(position <= dimension_, "Position should be less than dimension");
	return this->values_[position];
}

template <size_t dimension_>
double Euclidean_Vector<dimension_>::at(const size_t position) const {
	dynamic_require(position <= dimension_, "Position should be less than dimension");
	return this->values_[position];
}

template <size_t dimension_>
std::array<double, dimension_>::const_iterator Euclidean_Vector<dimension_>::cbegin(void) const {
	return this->values_.cbegin();
}

template <size_t dimension_>
std::array<double, dimension_>::const_iterator Euclidean_Vector<dimension_>::cend(void) const {
	return this->values_.cend();
}

template <size_t dimension_>
Euclidean_Vector<dimension_>& Euclidean_Vector<dimension_>::be_normalize(void) {
	const auto scale_factor = 1.0 / this->L2_norm();
	return (*this) *= scale_factor;
}

template <size_t dimension_>
constexpr size_t Euclidean_Vector<dimension_>::dimension(void) {
	return dimension_;
}

template <size_t dimension_>
double Euclidean_Vector<dimension_>::inner_product(const Euclidean_Vector& y) const {
	if constexpr (dimension_ < ms::blas_dot_criteria) {
		double result = 0;
		for (size_t i = 0; i < dimension_; ++i)
			result += this->values_[i] * y.values_[i];
		return result;
	}
	else
		return cblas_ddot(static_cast<MKL_INT>(dimension_), this->values_.data(), 1, y.values_.data(), 1);
}

template <size_t dimension_>
double Euclidean_Vector<dimension_>::L1_norm(void) const {
	if constexpr (dimension_ < ms::blas_dasum_criteria) {
		double L1_norm = 0.0;
		for (size_t i = 0; i < dimension_; ++i)
			L1_norm += std::abs(this->values_[i]);
		return L1_norm;
	}
	else 
		return cblas_dasum(static_cast<MKL_INT>(dimension_), this->values_.data(), 1);
}

template <size_t dimension_>
double Euclidean_Vector<dimension_>::L2_norm(void) const {
	return std::sqrt(this->inner_product(*this));
}

template <size_t dimension_>
bool Euclidean_Vector<dimension_>::is_axis_translation(const Euclidean_Vector& other) const {
	const auto line_vector = *this - other;
	const auto L1_norm = line_vector.L1_norm();
	const auto L2_norm = line_vector.L2_norm();

	constexpr auto epsilon = 1.0E-10;
	if (std::abs(L1_norm - L2_norm) <= epsilon)
		return true;
	else
		return false;
}

template <size_t dimension_> 
std::string Euclidean_Vector<dimension_>::to_string(void) const {
	std::string result;
	for (const auto& element : this->values_)
		result += ms::double_to_string(element) + " ";
	result.pop_back();	
	return result;
}

template <size_t dimension_> std::ostream& operator<<(std::ostream& os, const Euclidean_Vector<dimension_>& x) {
	return os << x.to_string();
};

template <size_t dimension_>
Euclidean_Vector<dimension_> operator*(const double constant, const Euclidean_Vector<dimension_>& x) {
	return x * constant;
}

template <typename Function>
void Dynamic_Euclidean_Vector::apply(const Function& f) {
	std::transform(this->values_.begin(), this->values_.end(), this->values_.begin(), f);
}