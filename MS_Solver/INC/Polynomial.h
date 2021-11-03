#pragma once
//#include "../INC/Exception.h"
#include "../INC/Vector_Function.h"

#include <array>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>

using ushort = unsigned short;


//Linear combination of monomial degree 1
class PolyTerm;
class Simple_Poly_Term
{
public:
	Simple_Poly_Term(const double constant) : constant_(constant) {};
	Simple_Poly_Term(const std::string& variable);
	Simple_Poly_Term(const std::vector<double>& coefficients, const double constant = 0);
	Simple_Poly_Term(const Simple_Poly_Term& other);

public://Command
	Simple_Poly_Term& operator+=(const Simple_Poly_Term& other);
	Simple_Poly_Term& operator-=(const Simple_Poly_Term& other);
	Simple_Poly_Term& operator*=(const double constant);
	void operator=(const Simple_Poly_Term& other);

public://Query
	Simple_Poly_Term operator+(const Simple_Poly_Term& other) const;
	Simple_Poly_Term operator*(const double constant) const;	
	PolyTerm operator*(const Simple_Poly_Term& other) const;
	template <typename V> double operator()(const V& domain_vector) const {
		REQUIRE(this->domain_dimension_ <= domain_vector.size(), "dimension of domain vector should be greater than domain dimension for clarity");

		auto coef_iter = this->coefficient_ptr_;
		auto dv_iter = domain_vector.begin();

		auto result = this->constant_;
		for (ushort i = 0; i < this->domain_dimension_; ++i, ++coef_iter, ++dv_iter)
			result += (*coef_iter) * (*dv_iter);
		
		return result;
	}	
	bool operator==(const Simple_Poly_Term& other) const;
	bool operator!=(const Simple_Poly_Term& other) const;
	bool operator<(const Simple_Poly_Term& other) const;
	//bool operator>(const Simple_Poly_Term& other) const;
	
	double to_constant(void) const;
	double differentiate(const ushort variable_index) const;
	ushort degree(void) const;
	ushort domain_dimension(void) const;
	bool is_constant(void) const;
	std::string to_string(void) const;

private:
	void change_domain_dimension(const ushort new_domain_dimension);
	bool is_small(void) const;

private:
	static constexpr ushort small_criterion_ = 3;

	ushort domain_dimension_ = 0;
	double constant_ = 0.0;
	double* coefficient_ptr_ = this->small_buffer_.data();	
	std::array<double, small_criterion_> small_buffer_ = { 0 };
	std::vector<double> coefficients_;
};


class PoweredPolyTerm
{
public:
	PoweredPolyTerm(void) = default;
	PoweredPolyTerm(const double constant) : base_(constant) {};
	PoweredPolyTerm(const Simple_Poly_Term& simple_poly_term) : base_(simple_poly_term) {};
	PoweredPolyTerm(const Simple_Poly_Term& simple_poly_term, const ushort exponent) : base_(simple_poly_term), exponent_(exponent) {};

public://Command
	void multiply_assign_with_same_base(const PoweredPolyTerm& other);

public://Query
	PolyTerm operator*(const double constant) const;
	template <typename V> double operator()(const V& domain_vector) const {
		return std::pow(this->base_(domain_vector), this->exponent_);
	}
	bool operator==(const PoweredPolyTerm& other) const;
	bool operator<(const PoweredPolyTerm& other) const;
	//bool operator>(const PoweredPolyTerm& other) const;
	
	PolyTerm differentiate(const ushort variable_index) const;
	ushort degree(void) const;
	ushort domain_dimension(void) const {
		return this->base_.domain_dimension();
	}
	double to_constant(void) const;
	Simple_Poly_Term be_simple(void) const;
	bool has_same_base(const PoweredPolyTerm& other) const;
	bool is_constant(void) const;
	bool is_simple(void) const;
	std::string to_string(void) const;

private:
	Simple_Poly_Term base_ = 0.0;
	ushort exponent_ = 1;
};


// constant * powered poly term1 * powered poly term2 * ...
class Polynomial;
class PolyTerm	// SSO
{
public:
	//PolyTerm(const double coeeficient);
	PolyTerm(const Simple_Poly_Term& simple_poly_term);
	PolyTerm(const PoweredPolyTerm& powered_poly_term);
	PolyTerm(const double constant, const PoweredPolyTerm& powered_poly_term = 1.0);
	PolyTerm(const PolyTerm& other);

public://command
	PolyTerm& operator*=(const double constant);
	PolyTerm& operator*=(const PolyTerm& other);
	PolyTerm& operator=(const PolyTerm& other);

	void add_assign_with_same_form(const PolyTerm& other);
	void minus_assign_with_same_form(const PolyTerm& other);

public://Query
	PolyTerm operator*(const double constant) const;
	PolyTerm operator*(const PolyTerm& other) const;
	template <typename V> double operator()(const V& domain_vector) const {
		auto result = this->constant_;
		for (ushort i = 0; i < this->num_term_; ++i)
			result *= this->term_ptr_[i](domain_vector);
		return result;
	}
	bool operator==(const PolyTerm& other) const;
	bool operator!=(const PolyTerm& other) const;

	Simple_Poly_Term be_simple(void) const;
	ushort degree(void) const;
	ushort domain_dimension(void) const;
	Polynomial differentiate(const ushort variable_index) const;
	double to_constant(void) const;
	bool has_same_form(const PolyTerm& other) const;
	bool is_simple(void) const;
	bool is_zero(void) const;
	std::string to_string(void) const;

private:
	void add_term(const PoweredPolyTerm& powered_poly_term);
	void multiply_assign_powered_poly_term(const PoweredPolyTerm& power_poly_term);
	bool is_constant(void) const;
	bool is_small(void) const;

private:
	static constexpr ushort small_criterion_ = 4;

	ushort num_term_ = 0;
	double constant_ = 0.0;
	PoweredPolyTerm* term_ptr_ = small_buffer_.data();
	std::array<PoweredPolyTerm, small_criterion_> small_buffer_ = { 0 };
	std::vector<PoweredPolyTerm> terms;
};


class Polynomial
{
public:
	Polynomial(void) = default;
	Polynomial(const double coeeficient) : simple_poly_term_(coeeficient) {};
	Polynomial(const std::string& variable) : simple_poly_term_(variable) {};
	Polynomial(const Simple_Poly_Term& simple_poly_term);
	Polynomial(const PolyTerm& poly_term);

public://Command
	Polynomial& operator+=(const Polynomial& other);
	Polynomial& operator-=(const Polynomial& other);
	Polynomial& operator*=(const double constant);

	Polynomial& be_absolute(void);
	Polynomial& be_derivative(const ushort variable_index);

public://Query
	Polynomial operator+(const Polynomial& other) const;
	Polynomial operator-(const Polynomial& other) const;
	Polynomial operator*(const Polynomial& other) const;
	Polynomial operator*(const double constant) const;
	Polynomial operator^(const ushort power_index) const;
	template<typename V> double operator()(const V& domain_vector) const {
		auto result = this->simple_poly_term_(domain_vector);
		for (const auto& poly_term : this->poly_terms_)
			result += poly_term(domain_vector);

		if (this->is_absolute_)
			return std::abs(result);
		else
			return result;
	}
	bool operator==(const Polynomial& other) const;

	ushort degree(void) const;
	ushort domain_dimension(void) const;
	Polynomial differentiate(const ushort variable_index) const;
	Vector_Function<Polynomial> gradient(void) const;
	size_t num_term(void) const;
	double to_constant(void) const;
	std::string to_string(void) const;
	//Irrational_Function<domain_dimension_> root(const double root_index) const;

private:
	void add_assign_poly_term(const PolyTerm& term);
	void minus_assign_poly_term(const PolyTerm& term);

	bool is_operable(void) const;

private:
	std::vector<PolyTerm> poly_terms_;
	Simple_Poly_Term simple_poly_term_ = 0.0;
	bool is_absolute_ = false;
};


Simple_Poly_Term operator*(const double constant, const Simple_Poly_Term& simple_poly_term);
PolyTerm operator*(const double constant, const PoweredPolyTerm& powered_poly_term);
PolyTerm operator*(const double constant, const PolyTerm& poly_term);
Polynomial operator*(const double constant, const Polynomial& polynomial);

std::ostream& operator<<(std::ostream& ostream, const Simple_Poly_Term& simple_poly_term);
std::ostream& operator<<(std::ostream& ostream, const PolyTerm& poly_term);
std::ostream& operator<<(std::ostream& ostream, const Polynomial& polynomial);







//#pragma once
////#include "Vector_Function.h"
//
//#include <algorithm>
//#include <iostream>
//
//
//
//template <ushort domain_dimension_>
//class Irrational_Function;
//
//
//template <ushort domain_dimension_>
//class Polynomial
//{
////private: // for test
//public: 
//	class SimplePolyTerm;
//	class PoweredPolyTerm;
//	class PolyTerm;
//
//private:
//	std::vector<PolyTerm> added_poly_term_set_;
//	SimplePolyTerm simple_poly_term_ = 0.0;
//	bool is_absolute_ = false;
//
//public:
//	static constexpr ushort domain_dimension(void);
//
//public:
//	Polynomial(void) = default;
//	Polynomial(const double coeeficient) : simple_poly_term_(coeeficient) {};
//	Polynomial(const std::string& variable) : simple_poly_term_(variable) {};
//
//public:
//	Polynomial& operator+=(const Polynomial& other);
//	Polynomial& operator-=(const Polynomial& other);
//	Polynomial& operator*=(const double constant);
//
//public:
//	Polynomial operator+(const Polynomial& other) const;
//	Polynomial operator-(const Polynomial& other) const;
//	Polynomial operator*(const Polynomial& other) const;
//	Polynomial operator*(const double constant) const;
//	Polynomial operator^(const ushort power_index) const;
//	double operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const;
//	bool operator==(const Polynomial& other) const;
//
//public:
//	Polynomial& be_absolute(void);
//	Polynomial& be_derivative(const ushort variable_index);	
//
//public:
//	size_t num_term(void) const { return this->added_poly_term_set_.size(); };
//	Polynomial differentiate(const ushort variable_index) const;
//	ushort degree(void) const;
//	double to_constant(void) const;
//	std::string to_string(void) const;
//	Irrational_Function<domain_dimension_> root(const double root_index) const;
//	Vector_Function<Polynomial<domain_dimension_>, domain_dimension_> gradient(void) const;
//
//private: 
//	void add_assign_poly_term(const PolyTerm& term);
//
//private:
//	bool is_operable(void) const;
//
////Inner classes defintion
//private: 
//	class SimplePolyTerm
//	{
//	private:
//		std::array<double, domain_dimension_> coefficients_ = { 0 };
//		double constant_ = 0.0;
//
//	public:
//		SimplePolyTerm(const double constant) : constant_(constant) {};
//		SimplePolyTerm(const std::string& variable);
//
//	public:
//		SimplePolyTerm& operator+=(const SimplePolyTerm& other);
//		SimplePolyTerm& operator-=(const SimplePolyTerm& other);
//		SimplePolyTerm& operator*=(const double constant);
//		SimplePolyTerm operator+(const SimplePolyTerm& other) const;
//		SimplePolyTerm operator*(const double constant) const;
//		double operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const;
//		bool operator==(const SimplePolyTerm& other) const;
//		bool operator!=(const SimplePolyTerm& other) const;
//		bool operator<(const SimplePolyTerm& other) const;
//		bool operator>(const SimplePolyTerm& other) const;
//
//	public:
//		double to_constant(void) const;
//		double differentiate(const ushort variable_index) const;
//		ushort degree(void) const;
//		bool is_constant(void) const;
//		std::string to_string(void) const;
//	};
//
//
	//class PoweredPolyTerm
	//{
	//private:
	//	SimplePolyTerm base_ = 0.0;
	//	ushort exponent_ = 1;

	//public:
	//	PoweredPolyTerm(void) = default;
	//	PoweredPolyTerm(const double constant) : base_(constant) {};
	//	PoweredPolyTerm(const SimplePolyTerm& simple_poly_term) : base_(simple_poly_term) {};
	//	PoweredPolyTerm(const SimplePolyTerm& simple_poly_term, const ushort exponent) : base_(simple_poly_term), exponent_(exponent) {};

	//public:
	//	void multiply_assign_with_same_base(const PoweredPolyTerm& other);
	//	double operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const;
	//	bool operator==(const PoweredPolyTerm& other) const;
	//	bool operator<(const PoweredPolyTerm& other) const;
	//	bool operator>(const PoweredPolyTerm& other) const;

	//public:
	//	double to_constant(void) const;
	//	SimplePolyTerm be_simple(void) const;
	//	PolyTerm differentiate(const ushort variable_index) const;
	//	bool has_same_base(const PoweredPolyTerm& other) const;
	//	bool is_constant(void) const;
	//	bool is_simple(void) const;
	//	ushort degree(void) const;
	//	std::string to_string(void) const;		
	//};
//
//
//	class PolyTerm	// SSO
//	{
//	private:
//		double coefficient_ = 1.0;
//		std::vector<PoweredPolyTerm> multiplied_powered_poly_term_set_;
//
//		ushort num_term_ = 0;
//		std::array<PoweredPolyTerm, 4> small_buffer_ = { 0 };
//		PoweredPolyTerm* data_ptr_ = small_buffer_.data();
//
//	public:
//		PolyTerm(const double coeeficient) : coefficient_(coeeficient) {};
//		PolyTerm(const SimplePolyTerm& simple_poly_term);
//		PolyTerm(const PoweredPolyTerm& powered_poly_term);
//		PolyTerm(const PolyTerm& other);
//
//	public:
//		void add_assign_with_same_form(const PolyTerm& other);
//		PolyTerm& operator*=(const PolyTerm& other);
//		PolyTerm& operator=(const PolyTerm& other);
//
//	public:
//		PolyTerm operator*(const PolyTerm& other) const;
//		double operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const;
//		bool operator==(const PolyTerm& other) const; 
//		bool operator!=(const PolyTerm& other) const;
//
//	public:
//		double to_constant(void) const;
//		SimplePolyTerm be_simple(void) const;
//		Polynomial differentiate(const ushort variable_index) const;
//		ushort degree(void) const;
//		bool has_same_form(const PolyTerm& other) const;
//		bool is_simple(void) const;
//		bool is_zero(void) const;
//		std::string to_string(void) const;
//
//	private:
//		void multiply_assign_powered_poly_term(const PoweredPolyTerm& power_poly_term);
//		bool is_small(void) const;
//	};
//};
//
//
//template <ushort domain_dimension_> 
//std::ostream& operator<<(std::ostream& ostream, const Polynomial<domain_dimension_>& polynomial);
//
//template <ushort domain_dimension_> 
//Polynomial<domain_dimension_> operator+(const double constant, const Polynomial<domain_dimension_>& polynomial);
//
//template <ushort domain_dimension_>
//Polynomial<domain_dimension_> operator-(const double constant, const Polynomial<domain_dimension_>& polynomial);
//
//template <ushort domain_dimension_> 
//Polynomial<domain_dimension_> operator*(const double constant, const Polynomial<domain_dimension_>& polynomial);
//
//
//template <ushort domain_dimension_>
//class Irrational_Function
//{
//private:
//	Polynomial<domain_dimension_> base_ = 0.0;
//	double exponent_ = 1.0;
//
//public:
//	Irrational_Function(void) = default;
//	Irrational_Function(const Polynomial<domain_dimension_>& polynomial, const double root_index = 1.0) : base_(polynomial), exponent_(root_index) {};
//
//	double operator()(const Euclidean_Vector<domain_dimension_>& value_vector) const;
//	bool operator==(const Irrational_Function& other) const;		
//	
//	std::string to_string(void) const;
//	//ushort order(void) const;
//};
//
//template <ushort domain_dimension>
//std::ostream& operator<<(std::ostream& ostream, const Irrational_Function<domain_dimension>& irrational_function);
//
//namespace ms {
//	template <ushort domain_dimension_> 
//	std::vector<Euclidean_Vector<domain_dimension_>> polynomial_compare_node_set(const ushort polynomial_order);
//
//	inline bool is_positive_odd_number(const double val);
//	inline bool is_natural_number(const double val);
//	constexpr ushort combination(const ushort n, const ushort k);
//	constexpr ushort combination_with_repetition(const ushort n, const ushort k);
//}
//
//
////template definition part
//template <ushort domain_dimension_> Polynomial<domain_dimension_>& Polynomial<domain_dimension_>::operator+=(const Polynomial& other) {
//	dynamic_require(this->is_operable() && other.is_operable(), "polynomials should be operable");
//
//	this->simple_poly_term_ += other.simple_poly_term_;
//	for (const auto& poly_term : other.added_poly_term_set_)
//		this->add_assign_poly_term(poly_term);
//
//	return *this;
//}
//
//template <ushort domain_dimension_> Polynomial<domain_dimension_>& Polynomial<domain_dimension_>::operator-=(const Polynomial& other) {
//	return *this += (-1 * other);
//}
//
//template <ushort domain_dimension_> Polynomial<domain_dimension_>& Polynomial<domain_dimension_>::operator*=(const double constant) {
//	if (constant == 0.0)
//		return *this = 0.0;
//
//	dynamic_require(this->is_operable(), "polynomial should be operable");
//
//	for (auto& poly_term : this->added_poly_term_set_)
//		poly_term *= constant;
//
//	this->simple_poly_term_ *= constant;
//
//	return *this;
//}
//
//template <ushort domain_dimension_> Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator+(const Polynomial& other) const {
//	Polynomial result(*this);
//	return result += other;
//}
//
//template <ushort domain_dimension_> Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator-(const Polynomial& other) const {
//	Polynomial result(*this);
//	return result += -1 * other;
//}
//
//template <ushort domain_dimension_> Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator*(const Polynomial& other) const {
//	dynamic_require(this->is_operable() && other.is_operable(), "polynomials should be operable");
//
//	Polynomial result = 0.0;
//
//	const auto num_this_term = this->added_poly_term_set_.size();
//	const auto num_other_term = other.added_poly_term_set_.size();
//	for (ushort i = 0; i < num_this_term; ++i)
//		for (ushort j = 0; j < num_other_term; ++j)
//			result.add_assign_poly_term(this->added_poly_term_set_[i] * other.added_poly_term_set_[j]);
//
//	if (this->simple_poly_term_ != 0.0)
//		for (ushort j = 0; j < num_other_term; ++j)
//			result.add_assign_poly_term(other.added_poly_term_set_[j] * this->simple_poly_term_);
//
//	if (other.simple_poly_term_ != 0.0)
//		for (ushort i = 0; i < num_this_term; ++i) 
//			result.add_assign_poly_term(this->added_poly_term_set_[i] * other.simple_poly_term_);
//
//	if (this->simple_poly_term_ != 0.0 && other.simple_poly_term_ != 0.0) {
//		if (this->simple_poly_term_.is_constant())
//			result.simple_poly_term_ = other.simple_poly_term_ * this->simple_poly_term_.to_constant();
//		else if (other.simple_poly_term_.is_constant())
//			result.simple_poly_term_ = this->simple_poly_term_ * other.simple_poly_term_.to_constant();
//		else
//			result.add_assign_poly_term(PolyTerm(this->simple_poly_term_) * other.simple_poly_term_);
//	}
//
//	return result;
//}
//
//template <ushort domain_dimension_> Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator*(const double constant) const {
//	Polynomial result(*this);
//	return result *= constant;
//}
//
//template <ushort domain_dimension_> 
//Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator^(const ushort power_index) const {
//	if (power_index == 0)
//		return 1;
//
//	dynamic_require(this->is_operable(), "polynomials should be operable");
//
//	auto result = *this;
//	for (ushort i = 1; i < power_index; ++i)
//		result = std::move(result * *this);
//	return result;
//}
//
//template <ushort domain_dimension_> 
//double Polynomial<domain_dimension_>::operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const {
//	auto result = this->simple_poly_term_(domain_vector);
//	for (const auto& poly_term : this->added_poly_term_set_)
//		result += poly_term(domain_vector);
//
//	if (this->is_absolute_)
//		return std::abs(result);
//	else
//		return result;
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::operator==(const Polynomial& other) const {
//	const auto max_order = std::max(this->degree(), other.degree());
//	const auto compare_node_set = ms::polynomial_compare_node_set<domain_dimension_>(max_order);
//
//	for (const auto& compare_node : compare_node_set) {
//		if (!ms::compare_double((*this)(compare_node), other(compare_node)))
//			return false;
//	}
//
//	return true;
//}
//
//template <ushort domain_dimension_>
//Polynomial<domain_dimension_>& Polynomial<domain_dimension_>::be_absolute(void) {
//	this->is_absolute_ = true;
//	return *this;
//}
//
//template <ushort domain_dimension_>
//Polynomial<domain_dimension_>& Polynomial<domain_dimension_>::be_derivative(const ushort variable_index) {
//	auto result = this->differentiate(variable_index);
//	return *this = std::move(result);
//};
//
//template <ushort domain_dimension_>
//Polynomial<domain_dimension_> Polynomial<domain_dimension_>::differentiate(const ushort variable_index) const {
//	if (domain_dimension_ <= variable_index)
//		return 0.0;
//
//	dynamic_require(this->is_operable(), "polynomials should be operable");
//
//	Polynomial result = 0.0;
//	result.simple_poly_term_ = this->simple_poly_term_.differentiate(variable_index);
//	for (const auto& poly_term : this->added_poly_term_set_)
//		result += poly_term.differentiate(variable_index);
//
//	return result;
//}
//
//template <ushort domain_dimension_>
//constexpr ushort Polynomial<domain_dimension_>::domain_dimension(void) {
//	return domain_dimension_;
//}
//
//template <ushort domain_dimension_> 
//ushort Polynomial<domain_dimension_>::degree(void) const {
//	ushort result = this->simple_poly_term_.degree();
//	for (const auto& term : this->added_poly_term_set_)
//		result = std::max(result, term.degree());
//	return result;
//}
//
//template <ushort domain_dimension_>
//double Polynomial<domain_dimension_>::to_constant(void) const {
//	dynamic_require(this->added_poly_term_set_.empty(), "constant polynomial should have empty poly term set");
//
//	if (this->is_absolute_)
//		return std::abs(this->simple_poly_term_.to_constant());
//	else
//		return this->simple_poly_term_.to_constant();
//}
//
//
//template <ushort domain_dimension_> 
//std::string Polynomial<domain_dimension_>::to_string(void) const {
//	if (this->added_poly_term_set_.empty())
//		return this->simple_poly_term_.to_string();
//
//	std::string str;
//	for (const auto& poly_term : this->added_poly_term_set_)
//		str += poly_term.to_string();
//
//	if (!(this->simple_poly_term_ == 0))
//		str += "+" + this->simple_poly_term_.to_string();
//
//	if (str.front() == '+')
//		str.erase(str.begin());
//
//	return str;
//}
//
//template <ushort domain_dimension_>
//Irrational_Function<domain_dimension_> Polynomial<domain_dimension_>::root(const double root_index) const {
//	return Irrational_Function<domain_dimension_>(*this, root_index);
//}
//
//template <ushort domain_dimension_>
//Vector_Function<Polynomial<domain_dimension_>, domain_dimension_> Polynomial<domain_dimension_>::gradient(void) const {
//	std::array<Polynomial<domain_dimension_>, domain_dimension_> gradient;
//
//	for (ushort i = 0; i < domain_dimension_; ++i)
//		gradient[i] = this->differentiate(i);	
//
//	return gradient;
//}
//
//template <ushort domain_dimension_> 
//void Polynomial<domain_dimension_>::add_assign_poly_term(const PolyTerm& term) {
//	dynamic_require(this->is_operable(), "polynomials should be operable");
//
//	for (auto iter = this->added_poly_term_set_.begin(); iter != this->added_poly_term_set_.end(); ++iter) {
//		if (iter->has_same_form(term)) {
//			iter->add_assign_with_same_form(term);
//			if (iter->is_zero())
//				this->added_poly_term_set_.erase(iter);
//			return;
//		}
//	}
//	this->added_poly_term_set_.push_back(term);
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::is_operable(void) const {
//	return !this->is_absolute_;
//}
//
//template <ushort domain_dimension_> 
//std::ostream& operator<<(std::ostream& ostream, const Polynomial<domain_dimension_>& polynomial) {
//	return ostream << polynomial.to_string();
//};
//
//template <ushort domain_dimension_> 
//Polynomial<domain_dimension_> operator+(const double constant, const Polynomial<domain_dimension_>& polynomial) {
//	return polynomial + constant;
//};
//
//template <ushort domain_dimension_>
//Polynomial<domain_dimension_> operator-(const double constant, const Polynomial<domain_dimension_>& polynomial) {
//	return -1 * polynomial + constant;
//}
//
//template <ushort domain_dimension_> 
//Polynomial<domain_dimension_> operator*(const double constant, const Polynomial<domain_dimension_>& polynomial) {
//	return polynomial * constant;
//};
//
//
//template <ushort domain_dimension_> 
//Polynomial<domain_dimension_>::SimplePolyTerm::SimplePolyTerm(const std::string& variable) {
//	dynamic_require(variable.front() == 'x', "variable should be start with 'x'");
//
//	constexpr ushort index_pos = 1;
//	const auto variable_index = std::stoull(variable.substr(index_pos));
//	dynamic_require(variable_index < domain_dimension_, "variable index should be less than space dimension");
//
//	this->coefficients_[variable_index] = 1.0;
//}
//
//template <ushort domain_dimension_> 
//typename Polynomial<domain_dimension_>::SimplePolyTerm& Polynomial<domain_dimension_>::SimplePolyTerm::operator+=(const SimplePolyTerm& other) {
//	this->constant_ += other.constant_;
//	for (ushort i = 0; i < domain_dimension_; ++i)
//		this->coefficients_[i] += other.coefficients_[i];
//	return *this;
//}
//
//template <ushort domain_dimension_> 
//typename Polynomial<domain_dimension_>::SimplePolyTerm& Polynomial<domain_dimension_>::SimplePolyTerm::operator-=(const SimplePolyTerm& other) {
//	this->constant_ -= other.constant_;
//	for (ushort i = 0; i < domain_dimension_; ++i)
//		this->coefficients_[i] -= other.coefficients_[i];
//	return *this;
//}
//
//template <ushort domain_dimension_> 
//typename Polynomial<domain_dimension_>::SimplePolyTerm& Polynomial<domain_dimension_>::SimplePolyTerm::operator*=(const double constant) {
//	this->constant_ *= constant;
//	for (ushort i = 0; i < domain_dimension_; ++i)
//		this->coefficients_[i] *= constant;
//	return *this;
//}
//
//template <ushort domain_dimension_>
//typename Polynomial<domain_dimension_>::SimplePolyTerm Polynomial<domain_dimension_>::SimplePolyTerm::operator+(const SimplePolyTerm& other) const {
//	auto result = *this;
//	return result += other;
//}
//
//template <ushort domain_dimension_>
//typename Polynomial<domain_dimension_>::SimplePolyTerm Polynomial<domain_dimension_>::SimplePolyTerm::operator*(const double constant) const {
//	auto result = *this;
//	return result *= constant;
//}
//
//template <ushort domain_dimension_> 
//double Polynomial<domain_dimension_>::SimplePolyTerm::operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const {
//	auto result = this->constant_;
//	for (ushort i = 0; i < domain_dimension_; ++i)
//		result += this->coefficients_[i] * domain_vector.at(i);
//	return result;
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::SimplePolyTerm::operator==(const SimplePolyTerm& other) const {
//	return this->constant_ == other.constant_ && this->coefficients_ == other.coefficients_;
//}
//
//template <ushort domain_dimension_> 
//bool Polynomial<domain_dimension_>::SimplePolyTerm::operator!=(const SimplePolyTerm& other) const {
//	return !(*this == other);
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::SimplePolyTerm::operator<(const SimplePolyTerm& other) const {
//	if (this->coefficients_ == other.coefficients_)
//		return this->constant_ < other.constant_;
//	else
//		return this->coefficients_ < other.coefficients_;
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::SimplePolyTerm::operator>(const SimplePolyTerm& other) const {
//	if (this->coefficients_ == other.coefficients_)
//		return this->constant_ > other.constant_;
//	else
//		return this->coefficients_ > other.coefficients_;
//}
//
//template <ushort domain_dimension_> 
//double Polynomial<domain_dimension_>::SimplePolyTerm::to_constant(void) const {
//	dynamic_require(this->is_constant(), "to be a constant, it should be constant");
//	return this->constant_;
//}
//
//template <ushort domain_dimension_> 
//double Polynomial<domain_dimension_>::SimplePolyTerm::differentiate(const ushort variable_index) const {
//	return this->coefficients_[variable_index];
//}
//
//template <ushort domain_dimension_> ushort Polynomial<domain_dimension_>::SimplePolyTerm::degree(void) const {
//	if (this->is_constant())
//		return 0;
//	else
//		return 1;
//}
//
//template <ushort domain_dimension_> bool Polynomial<domain_dimension_>::SimplePolyTerm::is_constant(void) const {
//	return this->coefficients_ == std::array<double, domain_dimension_>();
//}
//
//template <ushort domain_dimension_> std::string Polynomial<domain_dimension_>::SimplePolyTerm::to_string(void) const {
//	if (this->is_constant())
//		return +"[" + ms::double_to_string(this->constant_) + "]";
//
//	std::string str = "[";
//	for (ushort i = 0; i < domain_dimension_; ++i){
//		if (this->coefficients_[i] == 0.0)
//			continue;
//		else if (this->coefficients_[i] == 1.0)
//			str += "+x" + std::to_string(i);
//		else if (this->coefficients_[i] > 0)
//			str += "+" + ms::double_to_string(this->coefficients_[i]) + "(x" + std::to_string(i) + ")";
//		else if (this->coefficients_[i] == -1.0)
//			str += "-x" + std::to_string(i);
//		else
//			str += ms::double_to_string(this->coefficients_[i]) + "(x" + std::to_string(i) + ")";
//	}
//
//	constexpr ushort position = 1;
//	constexpr ushort size = 1;
//	if (str.at(position) == '+')
//		str.erase(position, size);
//
//	if (this->constant_ == 0)
//		return str += "]";
//	else if (this->constant_ > 0)
//		return str += '+' + ms::double_to_string(this->constant_) + "]";
//	else
//		return str += ms::double_to_string(this->constant_) + "]";
//}
//
//template <ushort domain_dimension_>
//void Polynomial<domain_dimension_>::PoweredPolyTerm::multiply_assign_with_same_base(const PoweredPolyTerm& other) {
//	this->exponent_ += other.exponent_;
//}
//
//template <ushort domain_dimension_> 
//double Polynomial<domain_dimension_>::PoweredPolyTerm::operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const {
//	return std::pow(this->base_(domain_vector), this->exponent_);
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::PoweredPolyTerm::operator==(const PoweredPolyTerm& other) const {
//	return this->base_ == other.base_ && this->exponent_ == other.exponent_;
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::PoweredPolyTerm::operator<(const PoweredPolyTerm& other) const {
//	if (this->exponent_ == other.exponent_)
//		return this->base_ < other.base_;
//	else
//		return this->exponent_ < other.exponent_;
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::PoweredPolyTerm::operator>(const PoweredPolyTerm& other) const {
//	if (this->exponent_ == other.exponent_)
//		return this->base_ > other.base_;
//	else
//		return this->exponent_ > other.exponent_;
//}
//
//template <ushort domain_dimension_>
//double Polynomial<domain_dimension_>::PoweredPolyTerm::to_constant(void) const {
//	return std::pow(this->base_.to_constant(), this->exponent_);
//}
//
//template <ushort domain_dimension_>
//
//typename Polynomial<domain_dimension_>::SimplePolyTerm Polynomial<domain_dimension_>::PoweredPolyTerm::be_simple(void) const {
//	return this->base_;
//}
//
//template <ushort domain_dimension_>
//typename Polynomial<domain_dimension_>::PolyTerm Polynomial<domain_dimension_>::PoweredPolyTerm::differentiate(const ushort variable_index) const {
//	const auto base_derivative = this->base_.differentiate(variable_index);
//	if (base_derivative == 0.0)
//		return 0.0;
//
//	if (this->exponent_ == 1)
//		return base_derivative;
//	else {
//		auto tmp = *this;
//		tmp.exponent_--;
//		PolyTerm result = tmp;
//		return result * (base_derivative * this->exponent_);
//	}
//}
//
//template <ushort domain_dimension_> 
//bool Polynomial<domain_dimension_>::PoweredPolyTerm::has_same_base(const PoweredPolyTerm& other) const {
//	return this->base_ == other.base_;
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::PoweredPolyTerm::is_constant(void) const {
//	return this->base_.is_constant();
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::PoweredPolyTerm::is_simple(void) const {
//	return this->exponent_ == 1;
//}
//
//template <ushort domain_dimension_> ushort Polynomial<domain_dimension_>::PoweredPolyTerm::degree(void) const {
//	if (this->is_constant())
//		return 0;
//	else
//		return this->exponent_;
//}
//
//template <ushort domain_dimension_> 
//std::string Polynomial<domain_dimension_>::PoweredPolyTerm::to_string(void) const {
//	auto str = this->base_.to_string();
//	if (this->exponent_ != 1)
//		return  str + "^" + std::to_string(this->exponent_);
//	else
//		return str;
//}
//
//
//template <ushort domain_dimension_> 
//Polynomial<domain_dimension_>::PolyTerm::PolyTerm(const SimplePolyTerm& simple_poly_term) {
//	if (simple_poly_term.is_constant())
//		this->coefficient_ = simple_poly_term.to_constant();
//	else
//		this->data_ptr_[this->num_term_++] = simple_poly_term;
//}
//
//template <ushort domain_dimension_> Polynomial<domain_dimension_>::PolyTerm::PolyTerm(const PoweredPolyTerm& powered_poly_term) {
//	if (powered_poly_term.is_constant())
//		this->coefficient_ = powered_poly_term.to_constant();
//	else
//		this->data_ptr_[this->num_term_++] = powered_poly_term;
//}
//
//template <ushort domain_dimension_> Polynomial<domain_dimension_>::PolyTerm::PolyTerm(const PolyTerm& other) {
//	this->coefficient_ = other.coefficient_;
//	this->num_term_ = other.num_term_;
//	if (other.is_small()) {
//		this->small_buffer_ = other.small_buffer_;
//		this->data_ptr_ = this->small_buffer_.data();
//	}
//	else {
//		this->multiplied_powered_poly_term_set_ = other.multiplied_powered_poly_term_set_;
//		this->data_ptr_ = this->multiplied_powered_poly_term_set_.data();
//	}
//}
//
//template <ushort domain_dimension_> void Polynomial<domain_dimension_>::PolyTerm::add_assign_with_same_form(const PolyTerm& other) {
//	this->coefficient_ += other.coefficient_;
//	if (this->coefficient_ == 0.0)
//		*this = 0;
//}
//
//template <ushort domain_dimension_> 
//typename Polynomial<domain_dimension_>::PolyTerm& Polynomial<domain_dimension_>::PolyTerm::operator*=(const PolyTerm& other) {
//	this->coefficient_ *= other.coefficient_;
//	if (this->coefficient_ == 0.0)
//		return *this = 0.0;
//
//	for (ushort i = 0; i < other.num_term_; ++i)
//		this->multiply_assign_powered_poly_term(other.data_ptr_[i]);
//
//	if(this->is_small())
//		std::sort(this->small_buffer_.begin(), this->small_buffer_.begin() + this->num_term_);
//	else
//		std::sort(this->multiplied_powered_poly_term_set_.begin(), this->multiplied_powered_poly_term_set_.end());
//
//	return *this;
//}
//
//template <ushort domain_dimension_> 
//typename Polynomial<domain_dimension_>::PolyTerm Polynomial<domain_dimension_>::PolyTerm::operator*(const PolyTerm& other) const {
//	auto result = *this;
//	return result *= other;
//}
//
//template <ushort domain_dimension_> 
//double Polynomial<domain_dimension_>::PolyTerm::operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const {
//	auto result = this->coefficient_;
//	for (ushort i = 0; i < this->num_term_; ++i)
//		result *= this->data_ptr_[i](domain_vector);
//	return result;
//}
//
//template <ushort domain_dimension_> 
//bool Polynomial<domain_dimension_>::PolyTerm::operator==(const PolyTerm& other) const {
//	return this->coefficient_ == other.coefficient_ &&  this->small_buffer_ == other.small_buffer_;
//}
//
//template <ushort domain_dimension_> 
//bool Polynomial<domain_dimension_>::PolyTerm::operator!=(const PolyTerm& other) const {
//	return !(*this == other);
//}
//
//template <ushort domain_dimension_> 
//typename Polynomial<domain_dimension_>::PolyTerm& Polynomial<domain_dimension_>::PolyTerm::operator=(const PolyTerm& other) {
//	this->coefficient_ = other.coefficient_;
//	this->num_term_ = other.num_term_;
//	if (other.is_small()) {
//		this->small_buffer_ = other.small_buffer_;
//		this->data_ptr_ = this->small_buffer_.data();
//	}
//	else {
//		this->multiplied_powered_poly_term_set_ = other.multiplied_powered_poly_term_set_;
//		this->data_ptr_ = this->multiplied_powered_poly_term_set_.data();
//	}
//
//	return *this;
//}
//
//template <ushort domain_dimension_> 
//double Polynomial<domain_dimension_>::PolyTerm::to_constant(void) const {
//	return this->coefficient_;
//}
//
//template <ushort domain_dimension_>
//typename Polynomial<domain_dimension_>::SimplePolyTerm Polynomial<domain_dimension_>::PolyTerm::be_simple(void) const {
//	return this->data_ptr_[0].be_simple() * this->coefficient_;
//}
//
//template <ushort domain_dimension_>
//Polynomial<domain_dimension_> Polynomial<domain_dimension_>::PolyTerm::differentiate(const ushort variable_index) const {
//	Polynomial<domain_dimension_> result = 0.0;
//
//	for (ushort i = 0; i < this->num_term_; ++i) {
//		const auto diff_term = this->data_ptr_[i].differentiate(variable_index);
//
//		if (diff_term.is_zero())
//			continue;
//
//		PolyTerm derivative = this->coefficient_;
//		if (this->is_small()) {
//			ushort index = 0;
//			for (ushort j = 0; j < this->num_term_; ++j) {
//				if (j == i)
//					continue;
//				else 
//					derivative.small_buffer_[index++] = this->data_ptr_[j];					
//			}
//			derivative.num_term_ = index;
//		}
//		else {
//			derivative = *this;
//			derivative.multiplied_powered_poly_term_set_.erase(derivative.multiplied_powered_poly_term_set_.begin() + i);
//		}
//
//		derivative *= diff_term;
//
//		if (derivative.is_simple())
//			result.simple_poly_term_ += derivative.be_simple();
//		else
//			result.add_assign_poly_term(derivative);
//	}
//
//	return result;
//}
//
//template <ushort domain_dimension_> 
//ushort Polynomial<domain_dimension_>::PolyTerm::degree(void) const {
//	ushort result = 0;
//
//	for (ushort i = 0; i < this->num_term_; ++i)
//		result += this->data_ptr_[i].degree();
//
//	return result;
//}
//
//template <ushort domain_dimension_> 
//bool Polynomial<domain_dimension_>::PolyTerm::has_same_form(const PolyTerm& other) const {
//	if (this->num_term_ != other.num_term_)
//		return false;
//
//	for (ushort i = 0; i < this->num_term_; ++i) {
//		if (this->data_ptr_[i] != other.data_ptr_[i])
//			return false;
//	}
//
//	return true;
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::PolyTerm::is_simple(void) const {
//	if (this->num_term_ != 1)
//		return false;
//
//	return this->data_ptr_[0].is_simple();
//}
//
//template <ushort domain_dimension_>
//bool Polynomial<domain_dimension_>::PolyTerm::is_zero(void) const {
//	return this->coefficient_ == 0.0;
//}
//
//template <ushort domain_dimension_>
//std::string Polynomial<domain_dimension_>::PolyTerm::to_string(void) const {
//	std::string str;
//	if (std::abs(this->coefficient_) != 1.0) {
//		if (this->coefficient_ > 0)
//			str += "+";
//
//		str += ms::double_to_string(this->coefficient_);
//	}
//	else if (this->coefficient_ == 1.0)
//		str += "+";
//	else
//		str += "-";
//
//	for (ushort i = 0; i < this->num_term_; ++i)
//		str += this->data_ptr_[i].to_string();
//
//	return str;
//}
//
//template <ushort domain_dimension_> 
//void Polynomial<domain_dimension_>::PolyTerm::multiply_assign_powered_poly_term(const PoweredPolyTerm& power_poly_term) {
//	for (ushort i = 0; i < this->num_term_; ++i) {
//		if (this->data_ptr_[i].has_same_base(power_poly_term)) {
//			this->data_ptr_[i].multiply_assign_with_same_base(power_poly_term);
//			return;
//		}
//	}
//
//	this->num_term_++;
//
//	if (this->small_buffer_.size() < this->num_term_) {
//		std::cout << "\n poly term exceed small buffer \n";
//
//		this->multiplied_powered_poly_term_set_.resize(this->num_term_, 0.0);
//
//		if (this->num_term_ == this->small_buffer_.size() + 1) {
//			std::copy(this->small_buffer_.begin(), this->small_buffer_.end(), this->multiplied_powered_poly_term_set_.begin());
//			this->data_ptr_ = this->multiplied_powered_poly_term_set_.data();
//		}
//	}
//
//	const auto position = this->num_term_ - 1;
//	this->data_ptr_[position] = power_poly_term;
//}
//
//template <ushort domain_dimension_> 
//bool Polynomial<domain_dimension_>::PolyTerm::is_small(void) const {
//	return this->multiplied_powered_poly_term_set_.empty();
//}
//
//template <ushort domain_dimension_>
//double Irrational_Function<domain_dimension_>::operator()(const Euclidean_Vector<domain_dimension_>& value_vector) const {
//	return std::pow(this->base_(value_vector), this->exponent_);
//}
//
//template <ushort domain_dimension_>
//bool Irrational_Function<domain_dimension_>::operator==(const Irrational_Function& other) const {
//	return this->base_ == other.base_ && this->exponent_ == other.exponent_;
//}
//
//template <ushort domain_dimension_>
//std::string Irrational_Function<domain_dimension_>::to_string(void) const {
//	if (this->exponent_ != 1.0)
//		return "[" + this->base_.to_string() + "]^" + ms::double_to_string(this->exponent_);
//	else
//		return this->base_.to_string();
//}	
//
//template <ushort domain_dimension>
//std::ostream& operator<<(std::ostream& ostream, const Irrational_Function<domain_dimension>& irrational_function) {
//	return ostream << irrational_function.to_string();
//}
//
////
////
////class IrrationalFunction::PoweredPolynomial
////{
////public:
////	PoweredPolynomial(const Polynomial& base, const double exponent = 1.0) : base_(base), exponent_(exponent) {};
////
////	double operator()(const MathVector& value_vector) const;
////
////	ushort domain_dimension_(void) const;
////	ushort order(void) const;
////	std::string to_string(void) const;
////
////private:
////	Polynomial base_ = 0.0;
////	double exponent_ = 1.0;
////};
////
////
////class IrrationalFunction::Term
////{
////public:
////	Term(const PoweredPolynomial& powered_polynomial);
////
////	double operator()(const MathVector& value_vector) const;
////
////	ushort domain_dimension_(void) const;
////	ushort order(void) const;
////	std::string to_string(void) const;
////
////private:
////	double coefficient_ = 1.0;
////	std::vector<PoweredPolynomial> multiplied_powered_polynomial_set_;
////};
////
////
////template <ushort DomainDim, ushort RangeDim> class PolynomialVectorFunction
////{
////public:
////	template <typename... Args> VectorFunction(Args... args) : elements_{ args... } {};
////
////	Euclidean_Vector<RangeDim> operator()(const Euclidean_Vector < )
////		const Function& operator[](const ushort position) const;
////
////private:
////	std::array<Function, RangeDim> elements_;
////};
////template <typename... Args> struct ms::Enforce_polynomial_type;
////template <ushort DomainDim, ushort RangeDim, typename... Args> PolynomialVectorFunction(Args...)->PolynomialVectorFunction<typename ms::Enforce_polynomial_type<Args>::type, 1 + sizeof...(Rest)>; // user defined deduction guides
////
////
////
//
//

namespace ms {


	constexpr ushort combination(const ushort n, const ushort k) {
		//calculate nCk
		//the combination of n things taken k at a time without repetition.
		if (n == k || k == 0)
			return 1;
		else
			return combination(n - 1, k - 1) + combination(n - 1, k);
	}

	constexpr ushort combination_with_repetition(const ushort n, const ushort k) {
		//calculate nHk
		//the combination of n things taken k at a time with repetition.
		return combination(n + k - 1, k);
	}

	inline bool is_positive_odd_number(const double val) {
		if (val < 0)
			return false;

		if (val - std::floor(val) == 0)
			return static_cast<size_t>(val) % 2 == 0;
		else
			return false;
	}

	inline bool is_natural_number(const double val) {
		if (val < 0)
			return false;

		if (val - std::floor(val) == 0)
			return true;
		else
			return false;
	}

	inline bool compare_double(const double d1, const double d2, const size_t ULP_precision = 4) {
		const auto lower_ULP = d1 - std::nextafter(d1, std::numeric_limits<double>::lowest());
		const auto upper_ULP = std::nextafter(d1, std::numeric_limits<double>::max()) - d1;

		return d1 - ULP_precision * lower_ULP <= d2 && d2 <= d1 + ULP_precision * upper_ULP;
	}

	inline std::vector<std::vector<double>> polynomial_compare_node_set(const ushort polynomial_order, const ushort domain_dimension) {
		const auto num_node = ms::combination_with_repetition(polynomial_order + 1, domain_dimension);

		std::vector<std::vector<double>> compare_node_set;
		compare_node_set.reserve(num_node);

		std::vector<double> compare_node(domain_dimension);
		if (domain_dimension == 0) {
			compare_node_set.push_back(compare_node);
			return compare_node_set;
		}

		while (true) {
			auto iter = std::find(compare_node.begin(), compare_node.end(), polynomial_order);
			if (iter != compare_node.end()) {
				compare_node_set.push_back(compare_node);

				if (iter == compare_node.begin())
					break;

				std::fill(compare_node.begin(), compare_node.end(), 0);
				(*(--iter))++;

				if (compare_node.front() == polynomial_order) {
					compare_node_set.push_back(compare_node);
					break;
				}

			}

			double component_sum = 0;
			for (const auto& val : compare_node)
				component_sum += val;

			if (component_sum == polynomial_order) {
				compare_node_set.push_back(compare_node);
				const auto is_zero = [](const double i) {return i == 0; };
				auto iter = std::find_if_not(compare_node.rbegin(), compare_node.rend(), is_zero);
				*iter = 0;
				(*(++iter))++;
				continue;
			}

			compare_node_set.push_back(compare_node);
			compare_node.back()++;
		}

		return compare_node_set;
	}
}