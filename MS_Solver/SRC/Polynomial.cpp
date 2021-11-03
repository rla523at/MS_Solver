#include "../INC/Polynomial.h"

Simple_Poly_Term::Simple_Poly_Term(const std::string& variable) {
	REQUIRE(variable.front() == 'x', "variable should be start with 'x'");

	constexpr ushort index_pos = 1;
	const auto variable_index = static_cast<ushort>(std::stoul(variable.substr(index_pos)));

	this->change_domain_dimension(variable_index + 1);
	this->coefficient_ptr_[variable_index] = 1.0;
}

Simple_Poly_Term::Simple_Poly_Term(const std::vector<double>& coefficients, const double constant) {
	this->change_domain_dimension(static_cast<ushort>(coefficients.size()));

	this->constant_ = constant;

	for (ushort i = 0; i < this->domain_dimension_; ++i)
		this->coefficient_ptr_[i] += coefficients[i];	
}

Simple_Poly_Term::Simple_Poly_Term(const Simple_Poly_Term& other) {
	*this = other;
}

void Simple_Poly_Term::operator=(const Simple_Poly_Term& other) {
	this->domain_dimension_ = other.domain_dimension_;
	this->constant_ = other.constant_;

	if (other.is_small()) {
		this->small_buffer_ = other.small_buffer_;
		this->coefficient_ptr_ = this->small_buffer_.data();
	}
	else {
		this->coefficients_ = other.coefficients_;
		this->coefficient_ptr_ = this->coefficients_.data();
	}
}

Simple_Poly_Term& Simple_Poly_Term::operator+=(const Simple_Poly_Term& other) {
	if (this->domain_dimension_ < other.domain_dimension_)
		this->change_domain_dimension(other.domain_dimension_);
	
	this->constant_ += other.constant_;

	for (ushort i = 0; i < other.domain_dimension_; ++i)
		this->coefficient_ptr_[i] += other.coefficient_ptr_[i];

	if (this->is_constant())
		*this = this->constant_;

	return *this;
}

Simple_Poly_Term& Simple_Poly_Term::operator-=(const Simple_Poly_Term& other) {
	if (this->domain_dimension_ < other.domain_dimension_)
		this->change_domain_dimension(other.domain_dimension_);

	this->constant_ -= other.constant_;

	for (ushort i = 0; i < other.domain_dimension_; ++i)
		this->coefficient_ptr_[i] -= other.coefficient_ptr_[i];

	return *this;
}

Simple_Poly_Term& Simple_Poly_Term::operator*=(const double constant) {
	this->constant_ *= constant;

	for (ushort i = 0; i < this->domain_dimension_; ++i)
		this->coefficient_ptr_[i] *= constant;
	
	return *this;
}

Simple_Poly_Term Simple_Poly_Term::operator+(const Simple_Poly_Term& other) const {
	auto result = *this;
	return result += other;
}

Simple_Poly_Term Simple_Poly_Term::operator*(const double constant) const {
	auto result = *this;
	return result *= constant;
}

PolyTerm Simple_Poly_Term::operator*(const Simple_Poly_Term& other) const {
	PolyTerm result = *this;
	return result *= other;
}


bool Simple_Poly_Term::operator==(const Simple_Poly_Term& other) const {
	if (this->domain_dimension_ != other.domain_dimension_)
		return false;
	
	if (this->constant_ != other.constant_)
		return false;

	for (ushort i = 0; i < this->domain_dimension_; ++i) {
		if (this->coefficient_ptr_[i] != other.coefficient_ptr_[i])
			return false;
	}

	return true;
}

bool Simple_Poly_Term::operator!=(const Simple_Poly_Term& other) const {
	return !(*this == other);
}

bool Simple_Poly_Term::operator<(const Simple_Poly_Term& other) const {
	if (this->domain_dimension_ != other.domain_dimension_)
		return this->domain_dimension_ < other.domain_dimension_;

	for (ushort i = 0; i < this->domain_dimension_; ++i) {
		if (this->coefficient_ptr_[i] != other.coefficient_ptr_[i])
			return this->coefficient_ptr_[i] < other.coefficient_ptr_[i];
	}

	return this->constant_ < other.constant_;
}

//bool SimplePolyTerm::operator>(const SimplePolyTerm& other) const {
//	if (this->coefficients_ == other.coefficients_)
//		return this->constant_ > other.constant_;
//	else
//		return this->coefficients_ > other.coefficients_;
//}

double Simple_Poly_Term::to_constant(void) const {
	REQUIRE(this->is_constant(), "it should be constant");
	return this->constant_;
}

double Simple_Poly_Term::differentiate(const ushort variable_index) const {
	if (this->domain_dimension_ <= variable_index)
		return 0;
	else
		return this->coefficient_ptr_[variable_index];
}

ushort Simple_Poly_Term::degree(void) const {
	if (this->is_constant())
		return 0;
	else
		return 1;
}

ushort Simple_Poly_Term::domain_dimension(void) const {
	return this->domain_dimension_;
}

bool Simple_Poly_Term::is_constant(void) const {	
	for (ushort i = 0; i < this->domain_dimension_; ++i) {
		if (this->coefficient_ptr_[i] != 0)
			return false;
	}

	return true;
}

std::string Simple_Poly_Term::to_string(void) const {
	std::ostringstream os;
	os << std::setprecision(16) << std::showpos;

	os << "[";
	for (ushort i = 0; i < this->domain_dimension_; ++i){
		if (this->coefficient_ptr_[i] == 0.0)
			continue;
		else if (this->coefficient_ptr_[i] == 1.0)
			os << "+x" << i;
		else if (this->coefficient_ptr_[i] == -1.0)
			os << "-x" << i;
		else
			os << this->coefficient_ptr_[i] << "(x" << i << ")";
	}

	if (this->constant_ == 0)
		os << "]";
	else
		os << this->constant_ << "]";

	auto str = os.str();

	constexpr ushort position = 1;
	constexpr ushort size = 1;
	if (str.at(position) == '+')
		str.erase(position, size);

	return str;
}

void Simple_Poly_Term::change_domain_dimension(const ushort new_domain_dimension) {
	if (this->is_small()) {
		if (new_domain_dimension > this->small_criterion_) {
			this->coefficients_.resize(new_domain_dimension);
			this->coefficient_ptr_ = this->coefficients_.data(); // resize can cuase reallcoation -> ptr should be updated

			std::copy(this->small_buffer_.begin(), small_buffer_.end(), this->coefficients_.begin());
			this->small_buffer_.fill(0);
		}
	}
	else {
		this->coefficients_.resize(new_domain_dimension);
		this->coefficient_ptr_ = this->coefficients_.data(); // resize can cuase reallcoation -> ptr should be updated
	}

	this->domain_dimension_ = new_domain_dimension;
}

bool Simple_Poly_Term::is_small(void) const {
	return this->coefficients_.empty();
}




void PoweredPolyTerm::multiply_assign_with_same_base(const PoweredPolyTerm& other) {
	this->exponent_ += other.exponent_;
}

PolyTerm PoweredPolyTerm::operator*(const double constant) const {
	return { constant, *this };
}


bool PoweredPolyTerm::operator==(const PoweredPolyTerm& other) const {
	return this->base_ == other.base_ && this->exponent_ == other.exponent_;
}

bool PoweredPolyTerm::operator<(const PoweredPolyTerm& other) const {
	if (this->exponent_ == other.exponent_)
		return this->base_ < other.base_;
	else
		return this->exponent_ < other.exponent_;
}

//bool PoweredPolyTerm::operator>(const PoweredPolyTerm& other) const {
//	if (this->exponent_ == other.exponent_)
//		return this->base_ > other.base_;
//	else
//		return this->exponent_ > other.exponent_;
//}

double PoweredPolyTerm::to_constant(void) const {
	return std::pow(this->base_.to_constant(), this->exponent_);
}

Simple_Poly_Term PoweredPolyTerm::be_simple(void) const {
	return this->base_;
}

PolyTerm PoweredPolyTerm::differentiate(const ushort variable_index) const {
	const auto base_derivative = this->base_.differentiate(variable_index);

	if (base_derivative == 0.0)
		return 0.0;

	if (this->exponent_ == 1)
		return base_derivative;
	else 
		return this->exponent_ * PoweredPolyTerm(this->base_, this->exponent_ - 1) * base_derivative;
}

bool PoweredPolyTerm::has_same_base(const PoweredPolyTerm& other) const {
	return this->base_ == other.base_;
}

bool PoweredPolyTerm::is_constant(void) const {
	return this->base_.is_constant();
}

bool PoweredPolyTerm::is_simple(void) const {
	return this->exponent_ == 1;
}

 ushort PoweredPolyTerm::degree(void) const {
	if (this->is_constant())
		return 0;
	else
		return this->exponent_;
}

std::string PoweredPolyTerm::to_string(void) const {
	auto str = this->base_.to_string();
	if (this->exponent_ != 1)
		return str + "^" + std::to_string(this->exponent_);
	else
		return str;
}

//PolyTerm::PolyTerm(const double constant) {
//	this->constant_ = constant;
//	this->add_term(1.0);
//}

PolyTerm::PolyTerm(const Simple_Poly_Term& simple_poly_term) {
	if (simple_poly_term.is_constant()) {
		this->constant_ = simple_poly_term.to_constant();
		this->add_term(1.0);
	}
	else {
		this->constant_ = 1.0;
		this->add_term(simple_poly_term);
	}
}

PolyTerm::PolyTerm(const PoweredPolyTerm& powered_poly_term) {
	 if (powered_poly_term.is_constant()) {
		 this->constant_ = powered_poly_term.to_constant();
		 this->add_term(1.0);
	 }
	 else {
		 this->constant_ = 1.0;
		 this->add_term(powered_poly_term);
	 }
}

PolyTerm::PolyTerm(const double constant, const PoweredPolyTerm& powered_poly_term) {
	if (powered_poly_term.is_constant()) {
		this->constant_ = constant * powered_poly_term.to_constant();
		this->add_term(1.0);
	}
	else {
		this->constant_ = constant;
		this->add_term(powered_poly_term);
	}
}


 PolyTerm::PolyTerm(const PolyTerm& other) {
	this->constant_ = other.constant_;
	this->num_term_ = other.num_term_;
	if (other.is_small()) {
		this->small_buffer_ = other.small_buffer_;
		this->term_ptr_ = this->small_buffer_.data();
	}
	else {
		this->terms = other.terms;
		this->term_ptr_ = this->terms.data();
	}
}

 void PolyTerm::add_assign_with_same_form(const PolyTerm& other) {
	this->constant_ += other.constant_;
	if (this->constant_ == 0.0)
		*this = 0.0;
}

void PolyTerm::minus_assign_with_same_form(const PolyTerm& other) {
	this->constant_ -= other.constant_;
	if (this->constant_ == 0.0)
		*this = 0.0;
}

PolyTerm& PolyTerm::operator*=(const double constant) {
	this->constant_ *= constant;

	if (this->constant_ == 0.0)
		*this = 0.0;

	return *this;
}

PolyTerm& PolyTerm::operator*=(const PolyTerm& other) {
	*this *= other.constant_;

	for (ushort i = 0; i < other.num_term_; ++i)
		this->multiply_assign_powered_poly_term(other.term_ptr_[i]);

	if (this->is_small())
		std::sort(this->small_buffer_.begin(), this->small_buffer_.begin() + this->num_term_);
	else
		std::sort(this->terms.begin(), this->terms.end());

	return *this;
}

PolyTerm PolyTerm::operator*(const double constant) const {
	auto result = *this;
	return result *= constant;
}

PolyTerm PolyTerm::operator*(const PolyTerm& other) const {
	auto result = *this;
	return result *= other;
}

bool PolyTerm::operator==(const PolyTerm& other) const {
	if (this->constant_ != other.constant_)
		return false;

	return this->has_same_form(other);
}

bool PolyTerm::operator!=(const PolyTerm& other) const {
	return !(*this == other);
}

PolyTerm& PolyTerm::operator=(const PolyTerm& other) {
	this->constant_ = other.constant_;
	this->num_term_ = other.num_term_;

	if (other.is_small()) {
		this->small_buffer_ = other.small_buffer_;
		this->term_ptr_ = this->small_buffer_.data();
	}
	else {
		this->terms = other.terms;
		this->term_ptr_ = this->terms.data();
	}

	return *this;
}

Simple_Poly_Term PolyTerm::be_simple(void) const {
	return this->term_ptr_[0].be_simple() * this->constant_;
}

Polynomial PolyTerm::differentiate(const ushort variable_index) const {
	Polynomial result = 0.0;

	for (ushort i = 0; i < this->num_term_; ++i) {
		const auto diff_term = this->term_ptr_[i].differentiate(variable_index);

		if (diff_term.is_zero())
			continue;

		PolyTerm derivative = this->constant_;
		for (ushort j = 0; j < this->num_term_; ++j) {
			if (j == i)
				continue;
			else
				derivative.add_term(this->term_ptr_[j]);
		}
		


		//if (this->is_small()) {
		//	ushort index = 0;
		//	for (ushort j = 0; j < this->num_term_; ++j) {
		//		if (j == i)
		//			continue;
		//		else
		//			derivative.small_buffer_[index++] = this->term_ptr_[j];
		//	}
		//	derivative.num_term_ = index;
		//}
		//else {
		//	derivative = *this;
		//	derivative.multiplied_powered_poly_term_set_.erase(derivative.multiplied_powered_poly_term_set_.begin() + i);
		//}

		derivative *= diff_term;

		if (derivative.is_simple())
			result += derivative.be_simple();
		else
			result += derivative;
	}

	return result;
}


double PolyTerm::to_constant(void) const {
	return this->constant_;
}



ushort PolyTerm::degree(void) const {
	ushort result = 0;

	for (ushort i = 0; i < this->num_term_; ++i)
		result += this->term_ptr_[i].degree();

	return result;
}

ushort PolyTerm::domain_dimension(void) const {
	ushort domain_dimension = 0;

	for (ushort i = 0; i < this->num_term_; ++i)
		domain_dimension = std::max(domain_dimension, this->term_ptr_[i].domain_dimension());

	return domain_dimension;
}

bool PolyTerm::has_same_form(const PolyTerm& other) const {
	if (this->num_term_ != other.num_term_)
		return false;

	if (this->is_small())
		return this->small_buffer_ == other.small_buffer_;
	else
		return this->terms == other.terms;
}

bool PolyTerm::is_simple(void) const {
	if (this->num_term_ != 1)
		return false;

	return this->term_ptr_[0].is_simple();
}

bool PolyTerm::is_zero(void) const {
	return this->constant_ == 0.0;
}


std::string PolyTerm::to_string(void) const {
	std::ostringstream oss;
	oss << std::setprecision(16) << std::showpos;

	if (std::abs(this->constant_) != 1.0) 
		oss << this->constant_;
	else if (this->constant_ == 1.0)
		oss << "+";
	else
		oss << "-";

	for (ushort i = 0; i < this->num_term_; ++i)
		oss << this->term_ptr_[i].to_string();

	return oss.str();
}

void PolyTerm::add_term(const PoweredPolyTerm& powered_poly_term) {
	const auto new_term_pos = this->num_term_;
	this->num_term_++;

	if (this->is_small()) {
		if (this->small_criterion_ < this->num_term_) {
			this->terms.resize(this->num_term_); // resize can cuase reallcoation -> ptr should be updated
			this->term_ptr_ = this->terms.data();

			std::copy(this->small_buffer_.begin(), this->small_buffer_.end(), this->terms.begin());
			this->small_buffer_.fill(0.0);
		}
	}
	else {
		this->terms.resize(this->num_term_); // resize can cuase reallcoation -> ptr should be updated
		this->term_ptr_ = this->terms.data();
	}

	this->term_ptr_[new_term_pos] = powered_poly_term;
}

void PolyTerm::multiply_assign_powered_poly_term(const PoweredPolyTerm& power_poly_term) {
	if (power_poly_term == 0.0)
		*this = 0.0;

	if (power_poly_term == 1.0)
		return;

	for (ushort i = 0; i < this->num_term_; ++i) {
		if (this->term_ptr_[i].has_same_base(power_poly_term)) {
			this->term_ptr_[i].multiply_assign_with_same_base(power_poly_term);
			return;
		}
	}

	this->add_term(power_poly_term);
}

bool PolyTerm::is_small(void) const {
	return this->terms.empty();
}



















Polynomial::Polynomial(const Simple_Poly_Term& simple_poly_term) {
	this->simple_poly_term_ += simple_poly_term;
}

Polynomial::Polynomial(const PolyTerm& poly_term) {
	this->add_assign_poly_term(poly_term);
}

Polynomial& Polynomial::operator+=(const Polynomial& other) {
	REQUIRE(this->is_operable() && other.is_operable(), "polynomials should be operable");

	this->simple_poly_term_ += other.simple_poly_term_;
	for (const auto& poly_term : other.poly_terms_)
		this->add_assign_poly_term(poly_term);

	return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& other) {
	REQUIRE(this->is_operable() && other.is_operable(), "polynomials should be operable");

	this->simple_poly_term_ -= other.simple_poly_term_;
	for (const auto& poly_term : other.poly_terms_)
		this->minus_assign_poly_term(poly_term);

	return *this;
}

Polynomial& Polynomial::operator*=(const double constant) {
	if (constant == 0.0)
		return *this = 0.0;

	REQUIRE(this->is_operable(), "polynomial should be operable");

	this->simple_poly_term_ *= constant;
	for (auto& poly_term : this->poly_terms_)
		poly_term *= constant;

	return *this;
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
	Polynomial result(*this);
	return result += other;
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
	Polynomial result(*this);
	return result -= other;
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
	REQUIRE(this->is_operable() && other.is_operable(), "polynomials should be operable");

	Polynomial result = 0.0;

	const auto num_this_term = this->poly_terms_.size();
	const auto num_other_term = other.poly_terms_.size();
	for (ushort i = 0; i < num_this_term; ++i)
		for (ushort j = 0; j < num_other_term; ++j)
			result.add_assign_poly_term(this->poly_terms_[i] * other.poly_terms_[j]);

	if (this->simple_poly_term_ != 0.0)
		for (ushort j = 0; j < num_other_term; ++j)
			result.add_assign_poly_term(other.poly_terms_[j] * this->simple_poly_term_);

	if (other.simple_poly_term_ != 0.0)
		for (ushort i = 0; i < num_this_term; ++i)
			result.add_assign_poly_term(this->poly_terms_[i] * other.simple_poly_term_);

	if (this->simple_poly_term_ != 0.0 && other.simple_poly_term_ != 0.0) {
		if (this->simple_poly_term_.is_constant())
			result.simple_poly_term_ = other.simple_poly_term_ * this->simple_poly_term_.to_constant();
		else if (other.simple_poly_term_.is_constant())
			result.simple_poly_term_ = this->simple_poly_term_ * other.simple_poly_term_.to_constant();
		else
			result.add_assign_poly_term(PolyTerm(this->simple_poly_term_) * other.simple_poly_term_);
	}

	return result;
}

Polynomial Polynomial::operator*(const double constant) const {
	Polynomial result = *this;
	return result *= constant;
}

Polynomial Polynomial::operator^(const ushort power_index) const {
	if (power_index == 0)
		return 1.0;

	REQUIRE(this->is_operable(), "polynomials should be operable");

	auto result = *this;
	for (ushort i = 1; i < power_index; ++i)
		result = std::move(result * *this);
		//result *= *this;

	return result;
}

bool Polynomial::operator==(const Polynomial& other) const {
	const auto max_degree = std::max(this->degree(), other.degree());
	const auto max_domain_dimension = std::max(this->domain_dimension(), other.domain_dimension());
	const auto compare_node_set = ms::polynomial_compare_node_set(max_degree, max_domain_dimension);

	for (const auto& compare_node : compare_node_set) {
		if (!ms::compare_double((*this)(compare_node), other(compare_node)))
			return false;
	}

	return true;
}

Polynomial& Polynomial::be_absolute(void) {
	this->is_absolute_ = true;
	return *this;
}

//Polynomial& Polynomial::be_derivative(const ushort variable_index) {
//	auto result = this->differentiate(variable_index);
//	return *this = std::move(result);
//};
//
Polynomial Polynomial::differentiate(const ushort variable_index) const {
	if (this->domain_dimension() <= variable_index)
		return 0.0;

	REQUIRE(this->is_operable(), "polynomials should be operable");

	Polynomial result = 0.0;
	result.simple_poly_term_ = this->simple_poly_term_.differentiate(variable_index);
	for (const auto& poly_term : this->poly_terms_)
		result += poly_term.differentiate(variable_index);

	return result;
}

ushort Polynomial::degree(void) const {
	ushort result = this->simple_poly_term_.degree();
	for (const auto& term : this->poly_terms_)
		result = std::max(result, term.degree());
	return result;
}

ushort Polynomial::domain_dimension(void) const {
	ushort domain_dimension = this->simple_poly_term_.domain_dimension();

	for (const auto& term : poly_terms_)
		domain_dimension = std::max(domain_dimension, term.domain_dimension());

	return domain_dimension;
}

size_t Polynomial::num_term(void) const {
	return this->poly_terms_.size();
};

//double Polynomial::to_constant(void) const {
//	REQUIRE(this->poly_terms_.empty(), "constant polynomial should have empty poly term set");
//
//	if (this->is_absolute_)
//		return std::abs(this->simple_poly_term_.to_constant());
//	else
//		return this->simple_poly_term_.to_constant();
//}
//
//
std::string Polynomial::to_string(void) const {
	if (this->poly_terms_.empty())
		return this->simple_poly_term_.to_string();

	std::string str;
	for (const auto& poly_term : this->poly_terms_)
		str += poly_term.to_string();

	if (!(this->simple_poly_term_ == 0))
		str += "+" + this->simple_poly_term_.to_string();

	if (str.front() == '+')
		str.erase(str.begin());

	return str;
}
//
//Irrational_Function<domain_dimension_> Polynomial::root(const double root_index) const {
//	return Irrational_Function<domain_dimension_>(*this, root_index);
//}
//
//Vector_Function<Polynomial, domain_dimension_> Polynomial::gradient(void) const {
//	std::array<Polynomial, domain_dimension_> gradient;
//
//	for (ushort i = 0; i < domain_dimension_; ++i)
//		gradient[i] = this->differentiate(i);
//
//	return gradient;
//}
//
void Polynomial::add_assign_poly_term(const PolyTerm& term) {
	REQUIRE(this->is_operable(), "polynomials should be operable");

	for (auto iter = this->poly_terms_.begin(); iter != this->poly_terms_.end(); ++iter) {
		if (iter->has_same_form(term)) {
			iter->add_assign_with_same_form(term);
			if (iter->is_zero())
				this->poly_terms_.erase(iter);
			return;
		}
	}
	this->poly_terms_.push_back(term);
}

void Polynomial::minus_assign_poly_term(const PolyTerm& term) {
	REQUIRE(this->is_operable(), "polynomials should be operable");

	for (auto iter = this->poly_terms_.begin(); iter != this->poly_terms_.end(); ++iter) {
		if (iter->has_same_form(term)) {
			iter->minus_assign_with_same_form(term);
			if (iter->is_zero())
				this->poly_terms_.erase(iter);
			return;
		}
	}
	this->poly_terms_.push_back(-1 * term);
}

bool Polynomial::is_operable(void) const {
	return !this->is_absolute_;
}






std::ostream& operator<<(std::ostream& ostream, const Simple_Poly_Term& simple_poly_term) {
	return ostream << simple_poly_term.to_string();
}
std::ostream& operator<<(std::ostream& ostream, const PolyTerm& poly_term) {
	return ostream << poly_term.to_string();
}
std::ostream& operator<<(std::ostream& ostream, const Polynomial& polynomial) {
	return ostream << polynomial.to_string();
}


PolyTerm operator*(const double constant, const PoweredPolyTerm& powered_poly_term) {
	return powered_poly_term * constant;
}
PolyTerm operator*(const double constant, const PolyTerm& poly_term) {
	return poly_term * constant;
}
Polynomial operator*(const double constant, const Polynomial& polynomial) {
	return polynomial * constant;
}
