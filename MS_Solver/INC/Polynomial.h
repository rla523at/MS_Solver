#pragma once
#include "Vector_Function.h"

#include <algorithm>
#include <iostream>

template <size_t domain_dimension_>
class Irrational_Function;

template <size_t domain_dimension_>
class Polynomial
{
//private: // for test
public: 
	class SimplePolyTerm;
	class PoweredPolyTerm;
	class PolyTerm;

private:
	std::vector<PolyTerm> added_poly_term_set_;
	SimplePolyTerm simple_poly_term_ = 0.0;

public:
	Polynomial(const double coeeficient) : simple_poly_term_(coeeficient) {};
	Polynomial(const std::string& variable) : simple_poly_term_(variable) {};

	Polynomial& operator+=(const Polynomial& other);
	Polynomial& operator-=(const Polynomial& other);
	Polynomial& operator*=(const double constant);
	Polynomial operator+(const Polynomial& other) const;
	Polynomial operator-(const Polynomial& other) const;
	Polynomial operator*(const Polynomial& other) const;
	Polynomial operator*(const double constant) const;
	Polynomial operator^(const size_t power_index) const;
	double operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const;
	bool operator==(const Polynomial& other) const;

	static constexpr size_t domain_dimension(void);
	size_t order(void) const;
	std::string to_string(void) const;
	Irrational_Function<domain_dimension_> root(const double root_index) const;

	template <size_t variable_index>
	Polynomial& be_derivative(void);

	template <size_t variable_index>
	Polynomial differentiate(void) const;

	//Vector_Function<Polynomial<domain_dimension_>> gradient(void) const;
	//Vector_Function<Polynomial<domain_dimension_>> gradient(const size_t domain_dimension_) const;


private: 
	void add_assign_poly_term(const PolyTerm& term);


//Inner classes defintion
private: 
	class SimplePolyTerm
	{
	private:
		std::array<double, domain_dimension_> coefficients_ = { 0 };
		double constant_ = 0.0;

	public:
		SimplePolyTerm(const double constant) : constant_(constant) {};
		SimplePolyTerm(const std::string& variable);

		SimplePolyTerm& operator+=(const SimplePolyTerm& other);
		SimplePolyTerm& operator-=(const SimplePolyTerm& other);
		SimplePolyTerm& operator*=(const double constant);
		SimplePolyTerm operator+(const SimplePolyTerm& other) const;
		SimplePolyTerm operator*(const double constant) const;
		double operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const;
		bool operator==(const SimplePolyTerm& other) const;
		bool operator!=(const SimplePolyTerm& other) const;
		bool operator<(const SimplePolyTerm& other) const;
		bool operator>(const SimplePolyTerm& other) const;

		double be_constant(void) const;
		size_t order(void) const;
		bool is_constant(void) const;
		std::string to_string(void) const;

		template <size_t VariableIndex>
		double differentiate(void) const;
	};


	class PoweredPolyTerm
	{
	private:
		SimplePolyTerm base_ = 0.0;
		size_t exponent_ = 1;

	public:
		PoweredPolyTerm(void) = default;
		PoweredPolyTerm(const double constant) : base_(constant) {};
		PoweredPolyTerm(const SimplePolyTerm& simple_poly_term) : base_(simple_poly_term) {};
		PoweredPolyTerm(const SimplePolyTerm& simple_poly_term, const size_t exponent) : base_(simple_poly_term), exponent_(exponent) {};

		void multiply_assign_with_same_base(const PoweredPolyTerm& other);
		double operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const;
		bool operator==(const PoweredPolyTerm& other) const;
		bool operator<(const PoweredPolyTerm& other) const;
		bool operator>(const PoweredPolyTerm& other) const;

		double be_constant(void) const;
		SimplePolyTerm be_simple(void) const;
		bool has_same_base(const PoweredPolyTerm& other) const;
		bool is_constant(void) const;
		bool is_simple(void) const;
		size_t order(void) const;
		std::string to_string(void) const;

		template <size_t VariableIndex>
		PolyTerm differentiate(void) const;
	};


	class PolyTerm	// SSO
	{
	private:
		double coefficient_ = 1.0;
		std::vector<PoweredPolyTerm> multiplied_powered_poly_term_set_;

		size_t num_term_ = 0;
		std::array<PoweredPolyTerm, 4> small_buffer_ = { 0 };
		PoweredPolyTerm* data_ptr_ = small_buffer_.data();

	public:
		PolyTerm(const double coeeficient) : coefficient_(coeeficient) {};
		PolyTerm(const SimplePolyTerm& simple_poly_term);
		PolyTerm(const PoweredPolyTerm& powered_poly_term);
		PolyTerm(const PolyTerm& other);

		void add_assign_with_same_form(const PolyTerm& other);
		PolyTerm& operator*=(const PolyTerm& other);
		PolyTerm operator*(const PolyTerm& other) const;
		double operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const;
		bool operator==(const PolyTerm& other) const; 
		bool operator!=(const PolyTerm& other) const;
		PolyTerm& operator=(const PolyTerm& other);

		double be_constant(void) const;
		SimplePolyTerm be_simple(void) const;
		size_t order(void) const;
		bool has_same_form(const PolyTerm& other) const;
		bool is_simple(void) const;
		bool is_zero(void) const;
		std::string to_string(void) const;

		template <size_t variable_index>
		Polynomial differentiate(void) const;

	private:
		void multiply_assign_powered_poly_term(const PoweredPolyTerm& power_poly_term);
		bool is_small(void) const;
	};
};


template <size_t domain_dimension_> 
std::ostream& operator<<(std::ostream& ostream, const Polynomial<domain_dimension_>& polynomial);

template <size_t domain_dimension_> 
Polynomial<domain_dimension_> operator+(const double constant, const Polynomial<domain_dimension_>& polynomial);

template <size_t domain_dimension_>
Polynomial<domain_dimension_> operator-(const double constant, const Polynomial<domain_dimension_>& polynomial);

template <size_t domain_dimension_> 
Polynomial<domain_dimension_> operator*(const double constant, const Polynomial<domain_dimension_>& polynomial);


template <size_t domain_dimension_>
class Irrational_Function
{
private:
	Polynomial<domain_dimension_> base_ = 0.0;
	double exponent_ = 1.0;

public:
	Irrational_Function(void) = default;
	Irrational_Function(const Polynomial<domain_dimension_>& polynomial, const double root_index = 1.0) : base_(polynomial), exponent_(root_index) {};

	double operator()(const Euclidean_Vector<domain_dimension_>& value_vector) const;
	bool operator==(const Irrational_Function& other) const;		
	
	std::string to_string(void) const;
	//size_t order(void) const;
};

template <size_t domain_dimension>
std::ostream& operator<<(std::ostream& ostream, const Irrational_Function<domain_dimension>& irrational_function);

namespace ms {
	template <size_t domain_dimension_> std::vector<Euclidean_Vector<domain_dimension_>> polynomial_compare_node_set(const size_t polynomial_order);
	size_t combination(const size_t n, const size_t k);
	size_t combination_with_repetition(const size_t n, const size_t k);
	bool is_positive_odd_number(const double val);
	bool is_natural_number(const double val);
}


//template definition part
template <size_t domain_dimension_> Polynomial<domain_dimension_>& Polynomial<domain_dimension_>::operator+=(const Polynomial& other) {
	this->simple_poly_term_ += other.simple_poly_term_;
	for (const auto& poly_term : other.added_poly_term_set_)
		this->add_assign_poly_term(poly_term);

	return *this;
}

template <size_t domain_dimension_> Polynomial<domain_dimension_>& Polynomial<domain_dimension_>::operator-=(const Polynomial& other) {
	return *this += (-1 * other);
}

template <size_t domain_dimension_> Polynomial<domain_dimension_>& Polynomial<domain_dimension_>::operator*=(const double constant) {
	if (constant == 0.0)
		return *this = 0.0;

	for (auto& poly_term : this->added_poly_term_set_)
		poly_term *= constant;

	this->simple_poly_term_ *= constant;

	return *this;
}

template <size_t domain_dimension_> Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator+(const Polynomial& other) const {
	Polynomial result(*this);
	return result += other;
}

template <size_t domain_dimension_> Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator-(const Polynomial& other) const {
	Polynomial result(*this);
	return result += -1 * other;
}

template <size_t domain_dimension_> Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator*(const Polynomial& other) const {
	Polynomial result = 0.0;

	const auto num_this_term = this->added_poly_term_set_.size();
	const auto num_other_term = other.added_poly_term_set_.size();
	for (size_t i = 0; i < num_this_term; ++i)
		for (size_t j = 0; j < num_other_term; ++j)
			result.add_assign_poly_term(this->added_poly_term_set_[i] * other.added_poly_term_set_[j]);

	if (this->simple_poly_term_ != 0.0)
		for (size_t j = 0; j < num_other_term; ++j)
			result.add_assign_poly_term(other.added_poly_term_set_[j] * this->simple_poly_term_);

	if (other.simple_poly_term_ != 0.0)
		for (size_t i = 0; i < num_this_term; ++i) 
			result.add_assign_poly_term(this->added_poly_term_set_[i] * other.simple_poly_term_);

	if (this->simple_poly_term_ != 0.0 && other.simple_poly_term_ != 0.0) {
		if (this->simple_poly_term_.is_constant())
			result.simple_poly_term_ = other.simple_poly_term_ * this->simple_poly_term_.be_constant();
		else if (other.simple_poly_term_.is_constant())
			result.simple_poly_term_ = this->simple_poly_term_ * other.simple_poly_term_.be_constant();
		else
			result.add_assign_poly_term(PolyTerm(this->simple_poly_term_) * other.simple_poly_term_);
	}

	return result;
}

template <size_t domain_dimension_> Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator*(const double constant) const {
	Polynomial result(*this);
	return result *= constant;
}

template <size_t domain_dimension_> 
Polynomial<domain_dimension_> Polynomial<domain_dimension_>::operator^(const size_t power_index) const {
	if (power_index == 0)
		return 1;

	auto result = *this;
	for (size_t i = 1; i < power_index; ++i)
		result = std::move(result * *this);
	return result;
}

template <size_t domain_dimension_> 
double Polynomial<domain_dimension_>::operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const {
	auto result = this->simple_poly_term_(domain_vector);
	for (const auto& poly_term : this->added_poly_term_set_)
		result += poly_term(domain_vector);
	return result;
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::operator==(const Polynomial& other) const {
	const auto max_order = std::max(this->order(), other.order());
	const auto compare_node_set = ms::polynomial_compare_node_set<domain_dimension_>(max_order);

	for (const auto& compare_node : compare_node_set) {
		if (!ms::compare_double((*this)(compare_node), other(compare_node)))
			return false;
	}

	return true;
}

template <size_t domain_dimension_>
template <size_t variable_index>
Polynomial<domain_dimension_>& Polynomial<domain_dimension_>::be_derivative(void) {
	auto result = this->differentiate<variable_index>();
	return *this = std::move(result);
};

template <size_t domain_dimension_>
template <size_t variable_index>
Polynomial<domain_dimension_> Polynomial<domain_dimension_>::differentiate(void) const {
	if constexpr (domain_dimension_ <= variable_index)
		return 0.0;

	Polynomial result = 0.0;
	result.simple_poly_term_ = this->simple_poly_term_.differentiate<variable_index>();
	for (const auto& poly_term : this->added_poly_term_set_)
		result += poly_term.differentiate<variable_index>();

	return result;
}

template <size_t domain_dimension_>
constexpr size_t Polynomial<domain_dimension_>::domain_dimension(void) {
	return domain_dimension_;
}

//template <size_t domain_dimension_>
//Vector_Function<Polynomial<domain_dimension_>> Polynomial<domain_dimension_>::gradient(void) const {
//	return this->gradient(domain_dimension_);
//}

//template <size_t domain_dimension_>
//Vector_Function<Polynomial<domain_dimension_>> Polynomial<domain_dimension_>::gradient(const size_t domain_dimension_) const {
//	Vector_Function<Polynomial<domain_dimension_>> result(domain_dimension_);
//	for (size_t i = 0; i < domain_dimension_; ++i)
//		result.push_back(this->differentiate(i));
//
//	return result;
//}
//


template <size_t domain_dimension_> 
size_t Polynomial<domain_dimension_>::order(void) const {
	size_t result = this->simple_poly_term_.order();
	for (const auto& term : this->added_poly_term_set_)
		result = std::max(result, term.order());
	return result;
}

template <size_t domain_dimension_> std::string Polynomial<domain_dimension_>::to_string(void) const {
	if (this->added_poly_term_set_.empty())
		return this->simple_poly_term_.to_string();

	std::string str;
	for (const auto& poly_term : this->added_poly_term_set_)
		str += poly_term.to_string();

	if (!(this->simple_poly_term_ == 0))
		str += "+" + this->simple_poly_term_.to_string();

	if (str.front() == '+')
		str.erase(str.begin());

	return str;
}

template <size_t domain_dimension_>
Irrational_Function<domain_dimension_> Polynomial<domain_dimension_>::root(const double root_index) const {
	return Irrational_Function<domain_dimension_>(*this, root_index);
}

template <size_t domain_dimension_> 
void Polynomial<domain_dimension_>::add_assign_poly_term(const PolyTerm& term) {
	for (auto iter = this->added_poly_term_set_.begin(); iter != this->added_poly_term_set_.end(); ++iter) {
		if (iter->has_same_form(term)) {
			iter->add_assign_with_same_form(term);
			if (iter->is_zero())
				this->added_poly_term_set_.erase(iter);
			return;
		}
	}
	this->added_poly_term_set_.push_back(term);
}

template <size_t domain_dimension_> 
std::ostream& operator<<(std::ostream& ostream, const Polynomial<domain_dimension_>& polynomial) {
	return ostream << polynomial.to_string();
};

template <size_t domain_dimension_> 
Polynomial<domain_dimension_> operator+(const double constant, const Polynomial<domain_dimension_>& polynomial) {
	return polynomial + constant;
};

template <size_t domain_dimension_>
Polynomial<domain_dimension_> operator-(const double constant, const Polynomial<domain_dimension_>& polynomial) {
	return -1 * polynomial + constant;
}

template <size_t domain_dimension_> 
Polynomial<domain_dimension_> operator*(const double constant, const Polynomial<domain_dimension_>& polynomial) {
	return polynomial * constant;
};


template <size_t domain_dimension_> 
Polynomial<domain_dimension_>::SimplePolyTerm::SimplePolyTerm(const std::string& variable) {
	dynamic_require(variable.front() == 'x', "variable should be start with 'x'");

	constexpr size_t index_pos = 1;
	const auto variable_index = std::stoull(variable.substr(index_pos));
	dynamic_require(variable_index < domain_dimension_, "variable index should be less than space dimension");

	this->coefficients_[variable_index] = 1.0;
}

template <size_t domain_dimension_> 
typename Polynomial<domain_dimension_>::SimplePolyTerm& Polynomial<domain_dimension_>::SimplePolyTerm::operator+=(const SimplePolyTerm& other) {
	this->constant_ += other.constant_;
	for (size_t i = 0; i < domain_dimension_; ++i)
		this->coefficients_[i] += other.coefficients_[i];
	return *this;
}

template <size_t domain_dimension_> 
typename Polynomial<domain_dimension_>::SimplePolyTerm& Polynomial<domain_dimension_>::SimplePolyTerm::operator-=(const SimplePolyTerm& other) {
	this->constant_ -= other.constant_;
	for (size_t i = 0; i < domain_dimension_; ++i)
		this->coefficients_[i] -= other.coefficients_[i];
	return *this;
}

template <size_t domain_dimension_> 
typename Polynomial<domain_dimension_>::SimplePolyTerm& Polynomial<domain_dimension_>::SimplePolyTerm::operator*=(const double constant) {
	this->constant_ *= constant;
	for (size_t i = 0; i < domain_dimension_; ++i)
		this->coefficients_[i] *= constant;
	return *this;
}

template <size_t domain_dimension_>
typename Polynomial<domain_dimension_>::SimplePolyTerm Polynomial<domain_dimension_>::SimplePolyTerm::operator+(const SimplePolyTerm& other) const {
	auto result = *this;
	return result += other;
}

template <size_t domain_dimension_>
typename Polynomial<domain_dimension_>::SimplePolyTerm Polynomial<domain_dimension_>::SimplePolyTerm::operator*(const double constant) const {
	auto result = *this;
	return result *= constant;
}

template <size_t domain_dimension_> 
double Polynomial<domain_dimension_>::SimplePolyTerm::operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const {
	auto result = this->constant_;
	for (size_t i = 0; i < domain_dimension_; ++i)
		result += this->coefficients_[i] * domain_vector.at(i);
	return result;
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::SimplePolyTerm::operator==(const SimplePolyTerm& other) const {
	return this->constant_ == other.constant_ && this->coefficients_ == other.coefficients_;
}

template <size_t domain_dimension_> 
bool Polynomial<domain_dimension_>::SimplePolyTerm::operator!=(const SimplePolyTerm& other) const {
	return !(*this == other);
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::SimplePolyTerm::operator<(const SimplePolyTerm& other) const {
	if (this->coefficients_ == other.coefficients_)
		return this->constant_ < other.constant_;
	else
		return this->coefficients_ < other.coefficients_;
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::SimplePolyTerm::operator>(const SimplePolyTerm& other) const {
	if (this->coefficients_ == other.coefficients_)
		return this->constant_ > other.constant_;
	else
		return this->coefficients_ > other.coefficients_;
}

template <size_t domain_dimension_> 
double Polynomial<domain_dimension_>::SimplePolyTerm::be_constant(void) const {
	return this->constant_;
}

template <size_t domain_dimension_> 
template <size_t variable_index> 
double Polynomial<domain_dimension_>::SimplePolyTerm::differentiate() const {
	return this->coefficients_[variable_index];
}

//size_t Polynomial::SimplePolyTerm::domain_dimension_(void) const {
//	return this->domain_dimension__;
//}

template <size_t domain_dimension_> size_t Polynomial<domain_dimension_>::SimplePolyTerm::order(void) const {
	if (this->is_constant())
		return 0;
	else
		return 1;
}

template <size_t domain_dimension_> bool Polynomial<domain_dimension_>::SimplePolyTerm::is_constant(void) const {
	return this->coefficients_ == std::array<double, domain_dimension_>();

/*	for (const auto& coefficient : this->coefficients_) {
		if (coefficient != 0)
			return false;
	}
	return true*/;
}

template <size_t domain_dimension_> std::string Polynomial<domain_dimension_>::SimplePolyTerm::to_string(void) const {
	if (this->is_constant())
		return +"[" + ms::double_to_string(this->constant_) + "]";

	std::string str = "[";
	for (size_t i = 0; i < domain_dimension_; ++i){
		if (this->coefficients_[i] == 0.0)
			continue;
		else if (this->coefficients_[i] == 1.0)
			str += "+x" + std::to_string(i);
		else if (this->coefficients_[i] > 0)
			str += "+" + ms::double_to_string(this->coefficients_[i]) + "(x" + std::to_string(i) + ")";
		else if (this->coefficients_[i] == -1.0)
			str += "-x" + std::to_string(i);
		else
			str += ms::double_to_string(this->coefficients_[i]) + "(x" + std::to_string(i) + ")";
	}

	constexpr size_t position = 1;
	constexpr size_t size = 1;
	if (str.at(position) == '+')
		str.erase(position, size);

	if (this->constant_ == 0)
		return str += "]";
	else if (this->constant_ > 0)
		return str += '+' + ms::double_to_string(this->constant_) + "]";
	else
		return str += ms::double_to_string(this->constant_) + "]";
}

template <size_t domain_dimension_>
void Polynomial<domain_dimension_>::PoweredPolyTerm::multiply_assign_with_same_base(const PoweredPolyTerm& other) {
	this->exponent_ += other.exponent_;
}

template <size_t domain_dimension_> 
double Polynomial<domain_dimension_>::PoweredPolyTerm::operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const {
	return std::pow(this->base_(domain_vector), this->exponent_);
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::PoweredPolyTerm::operator==(const PoweredPolyTerm& other) const {
	return this->base_ == other.base_ && this->exponent_ == other.exponent_;
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::PoweredPolyTerm::operator<(const PoweredPolyTerm& other) const {
	if (this->exponent_ == other.exponent_)
		return this->base_ < other.base_;
	else
		return this->exponent_ < other.exponent_;
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::PoweredPolyTerm::operator>(const PoweredPolyTerm& other) const {
	if (this->exponent_ == other.exponent_)
		return this->base_ > other.base_;
	else
		return this->exponent_ > other.exponent_;
}

template <size_t domain_dimension_>
double Polynomial<domain_dimension_>::PoweredPolyTerm::be_constant(void) const {
	return std::pow(this->base_.be_constant(), this->exponent_);
}

template <size_t domain_dimension_>

typename Polynomial<domain_dimension_>::SimplePolyTerm Polynomial<domain_dimension_>::PoweredPolyTerm::be_simple(void) const {
	return this->base_;
}

template <size_t domain_dimension_> 
bool Polynomial<domain_dimension_>::PoweredPolyTerm::has_same_base(const PoweredPolyTerm& other) const {
	return this->base_ == other.base_;
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::PoweredPolyTerm::is_constant(void) const {
	return this->base_.is_constant();
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::PoweredPolyTerm::is_simple(void) const {
	return this->exponent_ == 1;
}

template <size_t domain_dimension_> size_t Polynomial<domain_dimension_>::PoweredPolyTerm::order(void) const {
	if (this->is_constant())
		return 0;
	else
		return this->exponent_;
}

template <size_t domain_dimension_> 
std::string Polynomial<domain_dimension_>::PoweredPolyTerm::to_string(void) const {
	auto str = this->base_.to_string();
	if (this->exponent_ != 1)
		return  str + "^" + std::to_string(this->exponent_);
	else
		return str;
}

template <size_t domain_dimension_>
template <size_t variable_index>
typename Polynomial<domain_dimension_>::PolyTerm Polynomial<domain_dimension_>::PoweredPolyTerm::differentiate(void) const {
	const auto base_derivative = this->base_.differentiate<variable_index>();
	if (base_derivative == 0.0)
		return 0.0;

	if (this->exponent_ == 1)
		return base_derivative;
	else {
		auto tmp = *this;
		tmp.exponent_--;
		PolyTerm result = tmp;
		return result * (base_derivative * this->exponent_);
	}
}

template <size_t domain_dimension_> 
Polynomial<domain_dimension_>::PolyTerm::PolyTerm(const SimplePolyTerm& simple_poly_term) {
	if (simple_poly_term.is_constant())
		this->coefficient_ = simple_poly_term.be_constant();
	else
		this->data_ptr_[this->num_term_++] = simple_poly_term;
}

template <size_t domain_dimension_> Polynomial<domain_dimension_>::PolyTerm::PolyTerm(const PoweredPolyTerm& powered_poly_term) {
	if (powered_poly_term.is_constant())
		this->coefficient_ = powered_poly_term.be_constant();
	else
		this->data_ptr_[this->num_term_++] = powered_poly_term;
}

template <size_t domain_dimension_> Polynomial<domain_dimension_>::PolyTerm::PolyTerm(const PolyTerm& other) {
	this->coefficient_ = other.coefficient_;
	this->num_term_ = other.num_term_;
	if (other.is_small()) {
		this->small_buffer_ = other.small_buffer_;
		this->data_ptr_ = this->small_buffer_.data();
	}
	else {
		this->multiplied_powered_poly_term_set_ = other.multiplied_powered_poly_term_set_;
		this->data_ptr_ = this->multiplied_powered_poly_term_set_.data();
	}
}

template <size_t domain_dimension_> void Polynomial<domain_dimension_>::PolyTerm::add_assign_with_same_form(const PolyTerm& other) {
	this->coefficient_ += other.coefficient_;
	if (this->coefficient_ == 0.0)
		*this = 0;
}

template <size_t domain_dimension_> 
typename Polynomial<domain_dimension_>::PolyTerm& Polynomial<domain_dimension_>::PolyTerm::operator*=(const PolyTerm& other) {
	this->coefficient_ *= other.coefficient_;
	if (this->coefficient_ == 0.0)
		return *this = 0.0;

	for (size_t i = 0; i < other.num_term_; ++i)
		this->multiply_assign_powered_poly_term(other.data_ptr_[i]);

	if(this->is_small())
		std::sort(this->small_buffer_.begin(), this->small_buffer_.begin() + this->num_term_);
	else
		std::sort(this->multiplied_powered_poly_term_set_.begin(), this->multiplied_powered_poly_term_set_.end());

	return *this;
}

template <size_t domain_dimension_> 
typename Polynomial<domain_dimension_>::PolyTerm Polynomial<domain_dimension_>::PolyTerm::operator*(const PolyTerm& other) const {
	auto result = *this;
	return result *= other;
}

template <size_t domain_dimension_> 
double Polynomial<domain_dimension_>::PolyTerm::operator()(const Euclidean_Vector<domain_dimension_>& domain_vector) const {
	auto result = this->coefficient_;
	for (size_t i = 0; i < this->num_term_; ++i)
		result *= this->data_ptr_[i](domain_vector);
	return result;
}

template <size_t domain_dimension_> 
bool Polynomial<domain_dimension_>::PolyTerm::operator==(const PolyTerm& other) const {
	return this->coefficient_ == other.coefficient_ &&  this->small_buffer_ == other.small_buffer_;
}

template <size_t domain_dimension_> 
bool Polynomial<domain_dimension_>::PolyTerm::operator!=(const PolyTerm& other) const {
	return !(*this == other);
}

template <size_t domain_dimension_> 
typename Polynomial<domain_dimension_>::PolyTerm& Polynomial<domain_dimension_>::PolyTerm::operator=(const PolyTerm& other) {
	this->coefficient_ = other.coefficient_;
	this->num_term_ = other.num_term_;
	if (other.is_small()) {
		this->small_buffer_ = other.small_buffer_;
		this->data_ptr_ = this->small_buffer_.data();
	}
	else {
		this->multiplied_powered_poly_term_set_ = other.multiplied_powered_poly_term_set_;
		this->data_ptr_ = this->multiplied_powered_poly_term_set_.data();
	}

	return *this;
}

template <size_t domain_dimension_> 
double Polynomial<domain_dimension_>::PolyTerm::be_constant(void) const {
	return this->coefficient_;
}

template <size_t domain_dimension_>
typename Polynomial<domain_dimension_>::SimplePolyTerm Polynomial<domain_dimension_>::PolyTerm::be_simple(void) const {
	return this->data_ptr_[0].be_simple() * this->coefficient_;
}

template <size_t domain_dimension_>
template <size_t variable_index>
Polynomial<domain_dimension_> Polynomial<domain_dimension_>::PolyTerm::differentiate(void) const {
	Polynomial<domain_dimension_> result = 0.0;

	for (size_t i = 0; i < this->num_term_; ++i) {
		const auto diff_term = this->data_ptr_[i].differentiate<variable_index>();

		if (diff_term.is_zero())
			continue;

		PolyTerm derivative = this->coefficient_;
		if (this->is_small()) {
			size_t index = 0;
			for (size_t j = 0; j < this->num_term_; ++j) {
				if (j == i)
					continue;
				else 
					derivative.small_buffer_[index++] = this->data_ptr_[j];					
			}
			derivative.num_term_ = index;
		}
		else {
			derivative = *this;
			derivative.multiplied_powered_poly_term_set_.erase(derivative.multiplied_powered_poly_term_set_.begin() + i);
		}

		derivative *= diff_term;

		if (derivative.is_simple())
			result.simple_poly_term_ += derivative.be_simple();
		else
			result.add_assign_poly_term(derivative);
	}

	return result;
}

template <size_t domain_dimension_> 
size_t Polynomial<domain_dimension_>::PolyTerm::order(void) const {
	size_t result = 0;

	for (size_t i = 0; i < this->num_term_; ++i)
		result += this->data_ptr_[i].order();

	return result;
}

template <size_t domain_dimension_> 
bool Polynomial<domain_dimension_>::PolyTerm::has_same_form(const PolyTerm& other) const {
	if (this->num_term_ != other.num_term_)
		return false;

	for (size_t i = 0; i < this->num_term_; ++i) {
		if (this->data_ptr_[i] != other.data_ptr_[i])
			return false;
	}

	return true;
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::PolyTerm::is_simple(void) const {
	if (this->num_term_ != 1)
		return false;

	return this->data_ptr_[0].is_simple();
}

template <size_t domain_dimension_>
bool Polynomial<domain_dimension_>::PolyTerm::is_zero(void) const {
	return this->coefficient_ == 0.0;
}

template <size_t domain_dimension_>
std::string Polynomial<domain_dimension_>::PolyTerm::to_string(void) const {
	std::string str;
	if (std::abs(this->coefficient_) != 1.0) {
		if (this->coefficient_ > 0)
			str += "+";

		str += ms::double_to_string(this->coefficient_);
	}
	else if (this->coefficient_ == 1.0)
		str += "+";
	else
		str += "-";

	for (size_t i = 0; i < this->num_term_; ++i)
		str += this->data_ptr_[i].to_string();

	return str;
}

template <size_t domain_dimension_> 
void Polynomial<domain_dimension_>::PolyTerm::multiply_assign_powered_poly_term(const PoweredPolyTerm& power_poly_term) {
	for (size_t i = 0; i < this->num_term_; ++i) {
		if (this->data_ptr_[i].has_same_base(power_poly_term)) {
			this->data_ptr_[i].multiply_assign_with_same_base(power_poly_term);
			return;
		}
	}

	this->num_term_++;

	if (this->small_buffer_.size() < this->num_term_) {
		std::cout << "\n poly term exceed small buffer \n";

		this->multiplied_powered_poly_term_set_.resize(this->num_term_, 0.0);

		if (this->num_term_ == this->small_buffer_.size() + 1) {
			std::copy(this->small_buffer_.begin(), this->small_buffer_.end(), this->multiplied_powered_poly_term_set_.begin());
			this->data_ptr_ = this->multiplied_powered_poly_term_set_.data();
		}
	}

	const auto position = this->num_term_ - 1;
	this->data_ptr_[position] = power_poly_term;
}

template <size_t domain_dimension_> 
bool Polynomial<domain_dimension_>::PolyTerm::is_small(void) const {
	return this->multiplied_powered_poly_term_set_.empty();
}

template <size_t domain_dimension_>
double Irrational_Function<domain_dimension_>::operator()(const Euclidean_Vector<domain_dimension_>& value_vector) const {
	return std::pow(this->base_(value_vector), this->exponent_);
}

template <size_t domain_dimension_>
bool Irrational_Function<domain_dimension_>::operator==(const Irrational_Function& other) const {
	return this->base_ == other.base_ && this->exponent_ == other.exponent_;
}

template <size_t domain_dimension_>
std::string Irrational_Function<domain_dimension_>::to_string(void) const {
	if (this->exponent_ != 1.0)
		return "[" + this->base_.to_string() + "]^" + ms::double_to_string(this->exponent_);
	else
		return this->base_.to_string();
}	

template <size_t domain_dimension>
std::ostream& operator<<(std::ostream& ostream, const Irrational_Function<domain_dimension>& irrational_function) {
	return ostream << irrational_function.to_string();
}

//
//
//class IrrationalFunction::PoweredPolynomial
//{
//public:
//	PoweredPolynomial(const Polynomial& base, const double exponent = 1.0) : base_(base), exponent_(exponent) {};
//
//	double operator()(const MathVector& value_vector) const;
//
//	size_t domain_dimension_(void) const;
//	size_t order(void) const;
//	std::string to_string(void) const;
//
//private:
//	Polynomial base_ = 0.0;
//	double exponent_ = 1.0;
//};
//
//
//class IrrationalFunction::Term
//{
//public:
//	Term(const PoweredPolynomial& powered_polynomial);
//
//	double operator()(const MathVector& value_vector) const;
//
//	size_t domain_dimension_(void) const;
//	size_t order(void) const;
//	std::string to_string(void) const;
//
//private:
//	double coefficient_ = 1.0;
//	std::vector<PoweredPolynomial> multiplied_powered_polynomial_set_;
//};
//
//
//template <size_t DomainDim, size_t RangeDim> class PolynomialVectorFunction
//{
//public:
//	template <typename... Args> VectorFunction(Args... args) : elements_{ args... } {};
//
//	Euclidean_Vector<RangeDim> operator()(const Euclidean_Vector < )
//		const Function& operator[](const size_t position) const;
//
//private:
//	std::array<Function, RangeDim> elements_;
//};
//template <typename... Args> struct ms::Enforce_polynomial_type;
//template <size_t DomainDim, size_t RangeDim, typename... Args> PolynomialVectorFunction(Args...)->PolynomialVectorFunction<typename ms::Enforce_polynomial_type<Args>::type, 1 + sizeof...(Rest)>; // user defined deduction guides
//
//
//


namespace ms {
	template <size_t domain_dimension_> 
	std::vector<Euclidean_Vector<domain_dimension_>> polynomial_compare_node_set(const size_t polynomial_order) {
		const auto num_node = ms::combination_with_repetition(polynomial_order + 1, domain_dimension_);

		std::vector<Euclidean_Vector<domain_dimension_>> compare_node_set;
		compare_node_set.reserve(num_node);

		std::array<double, domain_dimension_> compare_node = { 0 };
		if constexpr (domain_dimension_ == 0) {
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