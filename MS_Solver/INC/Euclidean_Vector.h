#pragma once
#include "Exception.h"

#include <iomanip>
#include <mkl.h>
#include <sstream>
#include <vector>

class Euclidean_Vector;
class Euclidean_Vector_Base
{
public:
	Euclidean_Vector operator-(const Euclidean_Vector_Base& other) const;
	Euclidean_Vector operator+(const Euclidean_Vector_Base& other) const;
	Euclidean_Vector operator*(const double constant) const;
	double operator[](const size_t position) const;
	bool operator==(const Euclidean_Vector_Base& other) const;

	double at(const size_t position) const;
	const double* begin(void) const;
	const double* data(void) const;
	double L1_norm(void) const;
	double L2_norm(void) const;
	double inner_product(const Euclidean_Vector_Base& other) const;
	bool is_axis_translation(const Euclidean_Vector_Base& other) const;
	size_t size(void) const;
	std::string to_string(void) const;

protected:
	int num_values_ = 0;
	const double* const_data_ptr_ = nullptr;
};

class Euclidean_Vector : public Euclidean_Vector_Base
{
public:
	Euclidean_Vector(void) = default;
	Euclidean_Vector(const size_t size);
	Euclidean_Vector(const std::initializer_list<double> list);
	Euclidean_Vector(std::vector<double>&& values);
	Euclidean_Vector(const std::vector<double>& values);
	template <typename Iter>	Euclidean_Vector(Iter first, Iter last) : values_(first, last) {
		this->num_values_ = static_cast<int>(this->values_.size());
		this->const_data_ptr_ = this->values_.data();
	};
	Euclidean_Vector(const Euclidean_Vector& other);
	Euclidean_Vector(Euclidean_Vector&& other) noexcept;

public://Command	
	void operator=(const Euclidean_Vector& other);
	void operator=(Euclidean_Vector&& other) noexcept;
	Euclidean_Vector& operator*=(const double constant);

	Euclidean_Vector& normalize(void);

private:
	std::vector<double> values_;
};

class Euclidean_Vector_Wrapper : public Euclidean_Vector_Base
{
public:
	Euclidean_Vector_Wrapper(const size_t num_value, const double* ptr);
	Euclidean_Vector_Wrapper(const std::vector<double>& values);
};


//class Euclidean_Vector
//{
//public:
//	Euclidean_Vector(const size_t size) : values_(size) {};
//	Euclidean_Vector(const std::initializer_list<double> list) : values_(list) {};
//	Euclidean_Vector(std::vector<double>&& values) : values_(std::move(values)) {};
//	Euclidean_Vector(const std::vector<double>& values) : values_(values) {};
//	template <typename Iter>	Euclidean_Vector(Iter first, Iter last) : values_(first, last) {};
//
//public://Command
//	Euclidean_Vector& operator*=(const double constant);
//	Euclidean_Vector& operator+=(const Euclidean_Vector& other);
//	Euclidean_Vector& operator-=(const Euclidean_Vector& other);
//
//	Euclidean_Vector& normalize(void);
//
//public://Query
//	//Euclidean_Vector operator-(const Euclidean_Vector& other) const;
//	//Euclidean_Vector operator+(const Euclidean_Vector& other) const;
//	//Euclidean_Vector operator*(const double constant) const;
//	double operator[](const size_t position) const;
//	bool operator==(const Euclidean_Vector& other) const;
//
//	double at(const size_t position) const;
//	const double* begin(void) const;	
//	double L2_norm(void) const;
//	double inner_product(const Euclidean_Vector& other) const;
//	size_t size(void) const;
//	std::string to_string(void) const;
//
//private:
//	std::vector<double> values_;
//
//
////public:
////	Euclidean_Vector& operator-=(const Euclidean_Vector& other);
////	Euclidean_Vector& operator*=(const double constant);
////
////
////public:
////	Euclidean_Vector operator-(const Euclidean_Vector& other) const;
////
////public:
////	template <typename Function>
////	void apply(const Function& f);
////	void be_absolute(void);
////
////public:
////	std::vector<double>::const_iterator end(void) const;
////	const double* data(void) const;
////	double inner_product(const Euclidean_Vector& other) const;
//};

Euclidean_Vector operator*(const double constant, const Euclidean_Vector_Base& x);
std::ostream& operator<<(std::ostream& os, const Euclidean_Vector& x);







































//#include <algorithm>
//#include <array>
//#include <iomanip>
//#include <initializer_list>
//#include <mkl.h>
//#include <string>
//#include <sstream>
//#include <stdexcept>
//#include <type_traits>
//#include <vector>




//
//
//#define static_require static_assert
//#define dynamic_require(requirement, state) if (!(requirement)) throw std::runtime_error(state)
//
//
//using ushort = unsigned short;
//
//
//namespace ms {
//	class BLAS;
//
//	inline constexpr ushort blas_dasum_criteria = 10;
//	inline constexpr ushort blas_axpy_criteria = 20;
//	inline constexpr ushort blas_dot_criteria = 15;
//	
//	template <typename... Args>
//	inline constexpr bool are_arithmetics = (... && std::is_arithmetic_v<Args>);
//}
//
//
////Template Class Euclidean_Vector

//class Euclidean_Vector
//{
//private:
//	std::array<double, dimension_> values_ = { 0 };
//
//public:
//	Euclidean_Vector(void) = default;
//	Euclidean_Vector(const std::array<double, dimension_>& other) : values_(other) {};
//	template <typename... Args>
//	Euclidean_Vector(Args... args);
//
//public:
//	template <size_t temp = dimension_, std::enable_if_t<temp == 1, bool> = true>
//	operator double(void) const { return values_[0]; };
//	operator std::array<double, dimension_>(void) const { return values_; };
//
//public:
//	Euclidean_Vector& operator+=(const Euclidean_Vector& y);
//	Euclidean_Vector& operator-=(const Euclidean_Vector& y);
//	Euclidean_Vector& operator*=(const double scalar);
//
//public:
//	template <size_t temp = dimension_, std::enable_if_t<temp != 1, bool> = true>
//	Euclidean_Vector operator+(const Euclidean_Vector& y) const;	
//
//	Euclidean_Vector operator-(const Euclidean_Vector& y) const;
//
//	template <size_t temp = dimension_, std::enable_if_t<temp != 1, bool> = true>
//	Euclidean_Vector operator*(const double scalar) const;
//
//	bool operator==(const Euclidean_Vector& y) const;
//	double operator[](const size_t position) const;
//
//public:
//	Euclidean_Vector& be_normalize(void);
//
//public:
//	double at(const size_t position) const;	
//	std::array<double, dimension_>::const_iterator cbegin(void) const;
//	std::array<double, dimension_>::const_iterator cend(void) const;
//	double inner_product(const Euclidean_Vector& y) const;
//	double L1_norm(void) const;
//	double L2_norm(void) const;
//	bool is_axis_translation(const Euclidean_Vector& other) const;
//	std::string to_string(void) const;
//
//public:
//	static constexpr size_t dimension(void);
//};
//
//
////user-defined deduction guides
//template <typename... Args>
//Euclidean_Vector(Args... args)->Euclidean_Vector<sizeof...(Args)>;
//
 
//std::ostream& operator<<(std::ostream& os, const Euclidean_Vector& x);
//
//template <size_t dimension_, std::enable_if_t<dimension_ != 1, bool> = true>
//Euclidean_Vector operator*(const double constant, const Euclidean_Vector& x);
//
//
//class Dynamic_Euclidean_Vector
//{
//private:
//	template <size_t dimension_>
//	friend class Euclidean_Vector;
//	friend class ms::BLAS;
//
//private:
//	std::vector<double> values_;
//
//public:
//	Dynamic_Euclidean_Vector(const size_t dimension) : values_(dimension) {};
//	Dynamic_Euclidean_Vector(const std::initializer_list<double> list) : values_(list) {};
//	Dynamic_Euclidean_Vector(std::vector<double>&& values) : values_(std::move(values)) {};
//	Dynamic_Euclidean_Vector(const std::vector<double>& values) : values_(values) {};
//
//public:
//	Dynamic_Euclidean_Vector& operator-=(const Dynamic_Euclidean_Vector& other);
//	Dynamic_Euclidean_Vector& operator*=(const double constant);
//
//
//public:
//	Dynamic_Euclidean_Vector operator-(const Dynamic_Euclidean_Vector& other) const;
//	double operator[](const size_t position) const;
//	bool operator==(const Dynamic_Euclidean_Vector& other) const;
//
//public:
//	template <typename Function>
//	void apply(const Function& f);
//	void be_absolute(void);
//
//public:
//	double at(const size_t position) const;
//	std::vector<double>::const_iterator begin(void) const;
//	std::vector<double>::const_iterator end(void) const;
//	const double* data(void) const;
//	size_t dimension(void) const;
//	double inner_product(const Dynamic_Euclidean_Vector& other) const;
//	std::string to_string(void) const;
//};
//
//std::ostream& operator<<(std::ostream& os, const Dynamic_Euclidean_Vector& x);
//
//
//namespace ms {
//	inline std::string double_to_string(const double val) {
//		constexpr size_t precision = 16;
//		std::stringstream stream;
//		stream << std::setprecision(precision) << std::noshowpoint << val;
//		return stream.str();
//	}
//
//	inline bool compare_double(const double d1, const double d2, const size_t ULP_precision = 4) {
//		const auto lower_ULP = d1 - std::nextafter(d1, std::numeric_limits<double>::lowest());
//		const auto upper_ULP = std::nextafter(d1, std::numeric_limits<double>::max()) - d1;
//
//		return d1 - ULP_precision * lower_ULP <= d2 && d2 <= d1 + ULP_precision * upper_ULP;
//	}
//	template <size_t dimension>
//	Euclidean_Vector<dimension> arithmetic_mean(const std::vector<Euclidean_Vector<dimension>>& euclidean_vectors) {
//		Euclidean_Vector<dimension> arithmetic_mean;
//
//		for (const auto& vector : euclidean_vectors)
//			arithmetic_mean += vector;
//
//		const auto num_vector = euclidean_vectors.size();
//		arithmetic_mean *= 1.0 / num_vector;
//
//		return arithmetic_mean;
//	}
//	template <size_t dimension>
//	Euclidean_Vector<dimension> gather_min_value(const std::vector<Euclidean_Vector<dimension>>& euclideean_vectors) {
//		dynamic_require(!euclideean_vectors.empty(), "min value can not collected from empty");
//		
//		const auto num_vector = euclideean_vectors.size();
//
//		std::array<double, dimension> result = euclideean_vectors.front();
//		for (size_t i = 1; i < num_vector; ++i) {
//			for (size_t j = 0; j < dimension; ++j) 
//				result[j] = (std::min)(result[j], euclideean_vectors[i][j]);
//		}
//
//		return result;
//	}
//	template <size_t dimension>
//	Euclidean_Vector<dimension> gather_max_value(const std::vector<Euclidean_Vector<dimension>>& euclideean_vectors) {
//		dynamic_require(!euclideean_vectors.empty(), "min value can not collected from empty");
//
//		const auto num_vector = euclideean_vectors.size();
//
//		std::array<double, dimension> result = euclideean_vectors.front();
//		for (size_t i = 1; i < num_vector; ++i) {
//			for (size_t j = 0; j < dimension; ++j)
//				result[j] = (std::max)(result[j], euclideean_vectors[i][j]);
//		}
//
//		return result;
//	}
//	template <size_t dimension>
//	std::vector<std::vector<double>> gather_by_element_order(const std::vector<Euclidean_Vector<dimension>>& euclideean_vectors) {
//		dynamic_require(!euclideean_vectors.empty(), "can not gather from empty");
//
//		std::vector<std::vector<double>> gathered_vector(dimension);
//
//		const auto num_vector = euclideean_vectors.size();
//
//		for (size_t j = 0; j < dimension; ++j)
//			for (size_t i = 0; i < num_vector; ++i) 
//				gathered_vector[j][i] = euclideean_vectors[i][j];		
//
//		return result;
//	}
//
//
//}
//
//
//// Template Definition Part

//template <typename... Args>
//Euclidean_Vector::Euclidean_Vector(Args... args) : values_{ static_cast<double>(args)... } {
//	static_require(sizeof...(Args) <= dimension_, "Number of arguments can not exceed dimension");
//	static_require(ms::are_arithmetics<Args...>, "every arguments should be arithmetics");
//};
//
//

//Euclidean_Vector& Euclidean_Vector::operator+=(const Euclidean_Vector& y) {
//	for (size_t i = 0; i < dimension_; ++i)
//		this->values_[i] += y.values_[i];
//	return *this;
//}
//

//Euclidean_Vector& Euclidean_Vector::operator-=(const Euclidean_Vector& y) {
//	for (size_t i = 0; i < dimension_; ++i)
//		this->values_[i] -= y.values_[i];
//	return *this;
//}
//

//Euclidean_Vector& Euclidean_Vector::operator*=(const double scalar) {
//	if constexpr (dimension_ < 10) {
//		for (size_t i = 0; i < dimension_; ++i)
//			this->values_[i] *= scalar;
//	}
//	else {
//		cblas_dscal(dimension_, scalar, this->values_.data(), 1);
//	}
//	return *this;
//}
//
 
//template <size_t temp, std::enable_if_t<temp != 1, bool>>
//Euclidean_Vector Euclidean_Vector::operator+(const Euclidean_Vector& y) const {
//	auto result = *this;
//	return result += y;
//}
//
 //Euclidean_Vector Euclidean_Vector::operator-(const Euclidean_Vector& y) const {
//	auto result = *this;
//	return result -= y;
//}
//
 
//template <size_t temp, std::enable_if_t<temp != 1, bool>>
//Euclidean_Vector Euclidean_Vector::operator*(const double scalar) const {
//	auto result = *this;
//	return result *= scalar;
//}
//
 
//bool Euclidean_Vector::operator==(const Euclidean_Vector& y) const {
//	for (size_t i = 0; i < dimension_; ++i) {
//		if (this->values_[i] != y.values_[i])
//			return false;
//	}
//	return true;
//}
//
//
////bool Euclidean_Vector::operator==(const Dynamic_Euclidean_Vector& y) const {
////	if (dimension_ != y.dimension())
////		return false;
////	
////	for (size_t i = 0; i < dimension_; ++i) {
////		if (this->values_[i] != y.values_[i])+
////			return false;
////	}
////	return true;
////}
//
 
//double Euclidean_Vector::operator[](const size_t position) const {
//	dynamic_require(position <= dimension_, "Position should be less than dimension");
//	return this->values_[position];
//}
//

//double Euclidean_Vector::at(const size_t position) const {
//	dynamic_require(position <= dimension_, "Position should be less than dimension");
//	return this->values_[position];
//}
//

//std::array<double, dimension_>::const_iterator Euclidean_Vector::cbegin(void) const {
//	return this->values_.cbegin();
//}
//

//std::array<double, dimension_>::const_iterator Euclidean_Vector::cend(void) const {
//	return this->values_.cend();
//}
//

//Euclidean_Vector& Euclidean_Vector::be_normalize(void) {
//	const auto scale_factor = 1.0 / this->L2_norm();
//	return (*this) *= scale_factor;
//}
//

//constexpr size_t Euclidean_Vector::dimension(void) {
//	return dimension_;
//}
//

//double Euclidean_Vector::inner_product(const Euclidean_Vector& y) const {
//	if constexpr (dimension_ < ms::blas_dot_criteria) {
//		double result = 0;
//		for (size_t i = 0; i < dimension_; ++i)
//			result += this->values_[i] * y.values_[i];
//		return result;
//	}
//	else
//		return cblas_ddot(static_cast<MKL_INT>(dimension_), this->values_.data(), 1, y.values_.data(), 1);
//}
//

//double Euclidean_Vector::L1_norm(void) const {
//	if constexpr (dimension_ < ms::blas_dasum_criteria) {
//		double L1_norm = 0.0;
//		for (size_t i = 0; i < dimension_; ++i)
//			L1_norm += std::abs(this->values_[i]);
//		return L1_norm;
//	}
//	else 
//		return cblas_dasum(static_cast<MKL_INT>(dimension_), this->values_.data(), 1);
//}
//

//double Euclidean_Vector::L2_norm(void) const {
//	return std::sqrt(this->inner_product(*this));
//}
//


//
 
//std::string Euclidean_Vector::to_string(void) const {
//	std::string result;
//	for (const auto& element : this->values_)
//		result += ms::double_to_string(element) + " ";
//	result.pop_back();	
//	return result;
//}
//
//  std::ostream& operator<<(std::ostream& os, const Euclidean_Vector& x) {
//	return os << x.to_string();
//};
//
//template <size_t dimension_, std::enable_if_t<dimension_ != 1, bool>>
//Euclidean_Vector operator*(const double constant, const Euclidean_Vector& x) {
//	return x * constant;
//}
//
//template <typename Function>
//void Dynamic_Euclidean_Vector::apply(const Function& f) {
//	std::transform(this->values_.begin(), this->values_.end(), this->values_.begin(), f);
//}