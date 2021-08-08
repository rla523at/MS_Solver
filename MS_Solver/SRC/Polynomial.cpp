#include "../INC/Polynomial.h"

//Polynomial& Polynomial::operator+=(const Polynomial& other) {
//	this->simple_poly_term_ += other.simple_poly_term_;
//	for (const auto& poly_term : other.added_poly_term_set_)
//		this->add_assign_poly_term(poly_term);
//
//	return *this;
//}
//
//Polynomial& Polynomial::operator-=(const Polynomial& other) {
//	return *this += (-1 * other);
//}
//
//Polynomial& Polynomial::operator*=(const double constant) {
//	if (constant == 0.0)
//		return *this = 0.0;
//
//	for (auto& poly_term : this->added_poly_term_set_)
//		poly_term *= constant;
//
//	this->simple_poly_term_ *= constant;
//
//	return *this;
//}
//
//Polynomial Polynomial::operator+(const Polynomial& other) const {
//	Polynomial result(*this);
//	return result += other;
//}
//
//Polynomial Polynomial::operator-(const Polynomial& other) const {
//	Polynomial result(*this);
//	return result += -1 * other;
//}
//
//Polynomial Polynomial::operator*(const Polynomial& other) const {
//	Polynomial result = 0.0;
//
//	const auto num_this_term = this->added_poly_term_set_.size();
//	const auto num_other_term = other.added_poly_term_set_.size();
//	for (size_t i = 0; i < num_this_term; ++i)
//		for (size_t j = 0; j < num_other_term; ++j)
//			result.add_assign_poly_term(this->added_poly_term_set_[i] * other.added_poly_term_set_[j]);
//
//	if (this->simple_poly_term_ != 0.0)
//		for (size_t j = 0; j < num_other_term; ++j)
//			result.add_assign_poly_term(other.added_poly_term_set_[j] * this->simple_poly_term_);
//
//	if (other.simple_poly_term_ != 0.0)
//		for (size_t i = 0; i < num_this_term; ++i)
//			result.add_assign_poly_term(this->added_poly_term_set_[i] * other.simple_poly_term_);
//
//	//if (this->simple_poly_term_ != 0.0 && other.simple_poly_term_ != 0.0) {
//	//	const auto poly_term = PolyTerm(this->simple_poly_term_) * other.simple_poly_term_;
//	//	if (poly_term.is_constant())
//	//		result.simple_poly_term_ = poly_term.be_constant();
//	//	else
//	//		result.add_assign_poly_term(poly_term);
//	//}
//
//	//»ó¼ö * simple poly term = simple poly term
//	if (this->simple_poly_term_ != 0.0 && other.simple_poly_term_ != 0.0) {
//		if (this->simple_poly_term_.is_constant())
//			result.simple_poly_term_ = other.simple_poly_term_ * this->simple_poly_term_.be_constant();
//		else if (other.simple_poly_term_.is_constant())
//			result.simple_poly_term_ = this->simple_poly_term_ * other.simple_poly_term_.be_constant();
//		else
//			result.add_assign_poly_term(PolyTerm(this->simple_poly_term_) * other.simple_poly_term_);
//	}
//
//
//	return result;
//}
//
//Polynomial Polynomial::operator*(const double constant) const {
//	Polynomial result(*this);
//	return result *= constant;
//}
//
//Polynomial Polynomial::operator^(const size_t power_index) const {
//	if (power_index == 0)
//		return 1;
//
//	auto result = *this;
//	for (size_t i = 1; i < power_index; ++i)
//		result = std::move(result * *this);
//	return result;
//}
//
//double Polynomial::operator()(const MathVector& variable_vector) const {
//	auto result = this->simple_poly_term_(variable_vector);
//	for (const auto& poly_term : this->added_poly_term_set_) 
//		result += poly_term(variable_vector);
//	return result;
//}
//
//bool Polynomial::operator==(const Polynomial& other) const {
//	const auto max_order = std::max(this->order(), other.order());
//	const auto max_domain_dimension = std::max(this->domain_dimension(), other.domain_dimension());
//	const auto compare_node_set = ms::polynomial_compare_node_set(max_order, max_domain_dimension);
//
//	for (const auto& compare_node : compare_node_set) {
//		if (!ms::compare_double((*this)(compare_node), other(compare_node)))
//			return false;
//	}
//
//	return true;
//}
//
//size_t Polynomial::domain_dimension(void) const {
//	size_t result = this->simple_poly_term_.domain_dimension();
//	for (const auto& term : this->added_poly_term_set_)
//		result = std::max(result, term.domain_dimension());
//	return result;
//}
//
//Polynomial& Polynomial::be_derivative(const size_t variable_index) {
//	auto result = this->differentiate(variable_index);
//	return *this = std::move(result);
//}
//
//Polynomial Polynomial::differentiate(const size_t variable_index) const {
//	Polynomial result = 0.0;
//	result.simple_poly_term_ = this->simple_poly_term_.differentiate(variable_index);
//	for (const auto& poly_term : this->added_poly_term_set_)
//		result += poly_term.differentiate(variable_index);
//
//	return result;
//}
//
//VectorFunction<Polynomial> Polynomial::gradient(void) const {
//	const auto domain_dimension = this->domain_dimension();
//	return this->gradient(domain_dimension);
//}
//
//VectorFunction<Polynomial> Polynomial::gradient(const size_t domain_dimension) const {
//	VectorFunction<Polynomial> result;
//	result.reserve(domain_dimension);
//	for (size_t i = 0; i < domain_dimension; ++i)
//		result.push_back(this->differentiate(i));
//
//	return result;
//}
//
//IrrationalFunction Polynomial::root(const double root_index) const {
//	return IrrationalFunction(*this, root_index);
//}
//
//size_t Polynomial::order(void) const {
//	size_t result = this->simple_poly_term_.order();
//	for (const auto& term : this->added_poly_term_set_)
//		result = std::max(result, term.order());
//	return result;
//}
//
//std::string Polynomial::to_string(void) const {
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
//void Polynomial::add_assign_poly_term(const PolyTerm& term) {
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
//std::ostream& operator<<(std::ostream& ostream, const Polynomial& polynomial) {
//	return ostream << polynomial.to_string();
//}
//
//Polynomial operator+(const double constant, const Polynomial& polynomial) {
//	return polynomial + constant;
//}
//
//Polynomial operator*(const double constant, const Polynomial& polynomial) {
//	return polynomial * constant;
//}


//Polynomial::SimplePolyTerm::SimplePolyTerm(const std::string& variable) {
//	dynamic_require(variable.front() == 'x', "variable should be start with 'x'");
//
//	constexpr size_t position = 1;
//	const auto variable_index = std::stoull(variable.substr(position));
//	this->domain_dimension_ = variable_index + 1;
//		
//	if (this->small_buffer_.size() < this->domain_dimension_) {
//		this->coefficient_vector_.resize(this->domain_dimension_);
//		this->data_ptr_ = this->coefficient_vector_.data();
//	}
//
//	this->data_ptr_[variable_index] = 1.0;
//}

//Polynomial::SimplePolyTerm::SimplePolyTerm(const SimplePolyTerm& other) {
//	this->constant_ = other.constant_;
//	this->domain_dimension_ = other.domain_dimension_;
//	if (other.is_small()) 
//		this->small_buffer_ = other.small_buffer_;
//	else {
//		this->coefficient_vector_ = other.coefficient_vector_;
//		this->data_ptr_ = this->coefficient_vector_.data();
//	}
//}

//Polynomial::SimplePolyTerm& Polynomial::SimplePolyTerm::operator+=(const SimplePolyTerm& other) {
//	this->match_dimension(other.domain_dimension_);
//	this->constant_ += other.constant_;
//	for (size_t i = 0; i < other.domain_dimension_; ++i)
//		this->data_ptr_[i] += other.data_ptr_[i];
//	return *this;
//}
//
//Polynomial::SimplePolyTerm& Polynomial::SimplePolyTerm::operator-=(const SimplePolyTerm& other) {
//	this->match_dimension(other.domain_dimension_);
//	this->constant_ -= other.constant_;
//	for (size_t i = 0; i < other.domain_dimension_; ++i)
//		this->data_ptr_[i] -= other.data_ptr_[i];
//	return *this;
//}
//
//Polynomial::SimplePolyTerm& Polynomial::SimplePolyTerm::operator*=(const double constant) {
//	this->constant_ *= constant;
//	for (size_t i = 0; i < this->domain_dimension_; ++i)
//		this->data_ptr_[i] *= constant;
//	return *this;
//}

//Polynomial::SimplePolyTerm Polynomial::SimplePolyTerm::operator*(const double constant) const {
//	auto result = *this;
//	return result *= constant;
//}

//bool Polynomial::SimplePolyTerm::operator==(const SimplePolyTerm& other) const {
//	if (this->constant_ != other.constant_ || this->domain_dimension_ != other.domain_dimension_)
//		return false;
//
//	for (size_t i = 0; i < this->domain_dimension_; ++i) {
//		if (this->data_ptr_[i] != other.data_ptr_[i])
//			return false;
//	}	
//	return true;
//}
//
//bool Polynomial::SimplePolyTerm::operator!=(const SimplePolyTerm& other) const {
//	return !(*this == other);
//}
//
//Polynomial::SimplePolyTerm& Polynomial::SimplePolyTerm::operator=(const SimplePolyTerm& other) {
//	this->constant_ = other.constant_;
//	this->domain_dimension_ = other.domain_dimension_;
//	if (other.is_small()) {
//		this->small_buffer_ = other.small_buffer_;
//		this->data_ptr_ = this->small_buffer_.data();
//	}
//	else {
//		this->coefficient_vector_ = other.coefficient_vector_;
//		this->data_ptr_ = this->coefficient_vector_.data();
//	}
//	return *this;
//}
//
//double Polynomial::SimplePolyTerm::be_constant(void) const {
//	return this->constant_;
//}
//
//double Polynomial::SimplePolyTerm::differentiate(const size_t variable_index) const {
//	if (this->domain_dimension_ <= variable_index) 
//		return 0.0;
//	else 
//		return this->data_ptr_[variable_index];
//}
//
//size_t Polynomial::SimplePolyTerm::domain_dimension(void) const {
//	return this->domain_dimension_;
//}
//
//size_t Polynomial::SimplePolyTerm::order(void) const {
//	if (this->is_constant())
//		return 0;
//	else
//		return 1;
//}
//
//bool Polynomial::SimplePolyTerm::is_constant(void) const {
//	return this->domain_dimension_ == 0;
//}
//
//std::string Polynomial::SimplePolyTerm::to_string(void) const {
//	if (this->is_constant())
//		return +"[" + ms::double_to_string(this->constant_) + "]";
//
//	std::string str = "[";
//	for (size_t i = 0; i < this->domain_dimension_; ++i) {
//		if (this->data_ptr_[i] == 0.0)
//			continue;
//		else if (this->data_ptr_[i] == 1.0)
//			str += "+x" + std::to_string(i);
//		else if (this->data_ptr_[i] > 0)
//			str += "+" + ms::double_to_string(this->data_ptr_[i]) + "x" + std::to_string(i);
//		else if (this ->data_ptr_[i] == -1.0)
//			str += "-x" + std::to_string(i);
//		else
//			str += ms::double_to_string(this->data_ptr_[i]) + "x" + std::to_string(i);
//	}
//
//	constexpr size_t position = 1;
//	constexpr size_t size = 1;
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

//void Polynomial::SimplePolyTerm::match_dimension(const size_t other_domain_dimension) {
//	this->domain_dimension_ = std::max(this->domain_dimension_, other_domain_dimension);
//	if (this->small_buffer_.size() < this->domain_dimension_) {
//		this->coefficient_vector_.resize(this->domain_dimension_);
//		if (this->is_small()) {
//			std::copy(std::begin(this->small_buffer_), std::end(this->small_buffer_), this->coefficient_vector_.begin());
//			this->data_ptr_ = this->coefficient_vector_.data();
//		}
//	}
//}
//
//bool Polynomial::SimplePolyTerm::is_small(void) const {
//	return this->coefficient_vector_.empty();
//}


//void Polynomial::PoweredPolyTerm::multiply_assign_with_same_base(const PoweredPolyTerm& other) {
//	this->exponent_ += other.exponent_;
//}
//
//bool Polynomial::PoweredPolyTerm::operator==(const PoweredPolyTerm& other) const {
//	return this->base_ == other.base_ && this->exponent_ == other.exponent_;
//}
//
//double Polynomial::PoweredPolyTerm::be_constant(void) const {
//	return std::pow(this->base_.be_constant(), this->exponent_);
//}
//
//size_t Polynomial::PoweredPolyTerm::domain_dimension(void) const {
//	return this->base_.domain_dimension();
//}
//
//Polynomial::PolyTerm Polynomial::PoweredPolyTerm::differentiate(const size_t variable_index) const {
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
//bool Polynomial::PoweredPolyTerm::has_same_base(const PoweredPolyTerm& other) const {
//	return this->base_ == other.base_;
//}
//
//bool Polynomial::PoweredPolyTerm::is_constant(void) const {
//	return this->base_.is_constant();
//}
//
//size_t Polynomial::PoweredPolyTerm::order(void) const {
//	if (this->is_constant())
//		return 0;
//	else
//		return this->exponent_;
//}
//
//std::string Polynomial::PoweredPolyTerm::to_string(void) const {
//	auto str = this->base_.to_string();
//	if (this->exponent_ != 1)
//		return  str + "^" + std::to_string(this->exponent_);
//	else
//		return str;
//}
//
//
//Polynomial::PolyTerm::PolyTerm(const SimplePolyTerm& simple_poly_term) {	
//	if (simple_poly_term.is_constant())
//		this->coefficient_ = simple_poly_term.be_constant();
//	else
//		this->data_ptr_[this->num_term_++] = simple_poly_term;
//}
//
//Polynomial::PolyTerm::PolyTerm(const PoweredPolyTerm& powered_poly_term) {
//	if (powered_poly_term.is_constant())
//		this->coefficient_ = powered_poly_term.be_constant();
//	else
//		this->data_ptr_[this->num_term_++] = powered_poly_term;
//}
//
//Polynomial::PolyTerm::PolyTerm(const PolyTerm& other) {
//	this->coefficient_ = other.coefficient_;
//	this->num_term_ = other.num_term_;
//	if (this->is_small()) {
//		this->small_buffer_ = other.small_buffer_;
//		this->data_ptr_ = this->small_buffer_.data();
//	}
//	else {
//		this->multiplied_powered_poly_term_set_ = other.multiplied_powered_poly_term_set_;
//		this->data_ptr_ = this->multiplied_powered_poly_term_set_.data();
//	}
//}
//
//void Polynomial::PolyTerm::add_assign_with_same_form(const PolyTerm& other) {
//	this->coefficient_ += other.coefficient_;
//	if (this->coefficient_ == 0.0)
//		*this = 0;
//}
//
//Polynomial::PolyTerm& Polynomial::PolyTerm::operator*=(const PolyTerm& other) {
//	this->coefficient_ *= other.coefficient_;
//	if (this->coefficient_ == 0.0)
//		return *this = 0.0;
//
//	for (size_t i = 0; i < other.num_term_; ++i)
//		this->multiply_assign_powered_poly_term(other.data_ptr_[i]);
//	return *this;
//}
//
//Polynomial::PolyTerm Polynomial::PolyTerm::operator*(const PolyTerm& other) const {
//	auto result = *this;
//	return result *= other;
//}
//
//bool Polynomial::PolyTerm::operator==(const PolyTerm& other) const {
//	const auto max_order = std::max(this->order(), other.order());
//	const auto max_domain_dimension = std::max(this->domain_dimension(), other.domain_dimension());
//	const auto compare_node_set = ms::polynomial_compare_node_set(max_order, max_domain_dimension);
//
//	for (const auto& compare_node : compare_node_set) {
//		if (!ms::compare_double((*this)(compare_node), other(compare_node)))
//			return false;
//	}
//
//	return true;	
//}
//
//bool Polynomial::PolyTerm::operator!=(const PolyTerm& other) const {
//	return !(*this == other);
//}
//
//Polynomial::PolyTerm& Polynomial::PolyTerm::operator=(const PolyTerm& other) {
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
//double Polynomial::PolyTerm::be_constant(void) const {
//	return this->coefficient_;
//}
//
//size_t Polynomial::PolyTerm::domain_dimension(void) const {
//	size_t result = 0;
//	for (size_t i = 0; i < this->num_term_; ++i)
//		result = std::max(result, this->data_ptr_[i].domain_dimension());
//	return result;
//}
//
//Polynomial Polynomial::PolyTerm::differentiate(const size_t variable_index) const {
//	Polynomial result = 0.0;
//	
//	for (size_t i = 0; i < this->num_term_; ++i) {
//		const auto diff_term = this->data_ptr_[i].differentiate(variable_index);
//		
//		if (!diff_term.is_zero()) {
//			PolyTerm derivative = this->coefficient_;
//			if (this->is_small()) {
//				size_t index = 0;
//				for (size_t j = 0; j < this->num_term_; ++j) {
//					if (j == i)
//						continue;
//					else {
//						derivative.small_buffer_[index++] = this->data_ptr_[j];
//						derivative.num_term_++;
//					}
//				}
//			}
//			else {
//				derivative = *this;
//				derivative.multiplied_powered_poly_term_set_.erase(derivative.multiplied_powered_poly_term_set_.begin() + i);
//			}
//
//			derivative *= diff_term;
//
//			if (derivative.is_constant())
//				result.simple_poly_term_ += derivative.be_constant();
//			else
//				result.add_assign_poly_term(derivative);
//		}
//	}
//
//	return result;
//}
//
//size_t Polynomial::PolyTerm::order(void) const {
//	size_t result = 0;
//	for (size_t i = 0; i < this->num_term_; ++i)
//		result += this->data_ptr_[i].order();
//	return result;
//}
//
//bool Polynomial::PolyTerm::has_same_form(const PolyTerm& other) const {
//	if (this->num_term_ != other.num_term_)
//		return false;
//
//	for (size_t i = 0; i < this->num_term_; ++i) {
//		bool is_here = false;
//		for (size_t j = 0; j < other.num_term_; ++j) {
//			if (this->data_ptr_[i] == other.data_ptr_[j])
//				is_here = true;
//		}
//		if (!is_here)
//			return false;
//	}
//	return true;
//}
//
//bool Polynomial::PolyTerm::is_constant(void) const {
//	return this->num_term_ == 0;
//}
//
//bool Polynomial::PolyTerm::is_zero(void) const {
//	return this->coefficient_ == 0.0;
//}
//
//std::string Polynomial::PolyTerm::to_string(void) const {
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
//	for (size_t i = 0; i < this->num_term_; ++i)
//		str += this->data_ptr_[i].to_string();
//
//	return str;
//}
//
//void Polynomial::PolyTerm::multiply_assign_powered_poly_term(const PoweredPolyTerm& power_poly_term) {
//	for (size_t i = 0; i < this->num_term_; ++i) {
//		if (this->data_ptr_[i].has_same_base(power_poly_term)) {
//			this->data_ptr_[i].multiply_assign_with_same_base(power_poly_term);
//			return;
//		}
//	}
//
//	if (this->small_buffer_.size() < ++this->num_term_) {
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
//bool Polynomial::PolyTerm::is_small(void) const {
//	return this->multiplied_powered_poly_term_set_.empty();
//}
//
//
//IrrationalFunction::IrrationalFunction(const Polynomial& polynomial, const double root_index) {
//	if (root_index == 1.0)
//		this->polynomial_term_ = polynomial;
//	else
//		this->added_irrational_term_set_.push_back(PoweredPolynomial(polynomial, root_index));	
//}
//
//double IrrationalFunction::operator()(const MathVector& value_vector) const {
//	double result = this->polynomial_term_(value_vector);
//	for (const auto& term : this->added_irrational_term_set_)
//		result += term(value_vector);
//	return result;
//}
//
//bool IrrationalFunction::operator==(const IrrationalFunction& other) const {
//	//temporary
//	const auto max_order = std::max(this->order(), other.order());
//	const auto max_domain_dimension = std::max(this->domain_dimension(), other.domain_dimension());
//	const auto compare_node_set = ms::polynomial_compare_node_set(max_order, max_domain_dimension);
//
//	for (const auto& compare_node : compare_node_set) {
//		if (!ms::compare_double((*this)(compare_node), other(compare_node)))
//			return false;
//	}
//
//	return true;
//}
//
//size_t IrrationalFunction::domain_dimension(void) const {
//	size_t result = this->polynomial_term_.domain_dimension();
//	for (const auto& term : this->added_irrational_term_set_)
//		result = std::max(result, term.domain_dimension());
//	return result;
//}
//
//size_t IrrationalFunction::order(void) const {
//	size_t result = this->polynomial_term_.order();
//	for (const auto& term : this->added_irrational_term_set_)
//		result = std::max(result, term.order());
//	return result;
//}
//
//std::string IrrationalFunction::to_string(void) const {
//	if (this->added_irrational_term_set_.empty())
//		return this->polynomial_term_.to_string();
//
//	std::string str;
//	for (const auto& irrational_term : this->added_irrational_term_set_)
//		str += irrational_term.to_string();
//
//	if (!(this->polynomial_term_ == 0))
//		str += "+" + this->polynomial_term_.to_string();
//
//	if (str.front() == '+')
//		str.erase(str.begin());
//
//	return str;
//}
//
//std::ostream& operator<<(std::ostream& ostream, const IrrationalFunction& irrational_function) {
//	return ostream << irrational_function.to_string();
//}
//
//double IrrationalFunction::PoweredPolynomial::operator()(const MathVector& value_vector) const {
//	return std::pow(this->base_(value_vector), this->exponent_);
//}
//
//size_t IrrationalFunction::PoweredPolynomial::domain_dimension(void) const {
//	return this->base_.domain_dimension();
//}
//
//size_t IrrationalFunction::PoweredPolynomial::order(void) const {
//	return this->base_.order() * static_cast<size_t>(std::ceil(this->exponent_));
//}
//
//std::string IrrationalFunction::PoweredPolynomial::to_string(void) const {
//	auto str = this->base_.to_string();
//	if (this->exponent_ != 1)
//		return  str + "^" + ms::double_to_string(this->exponent_);
//	else
//		return str;
//}
//
//IrrationalFunction::Term::Term(const PoweredPolynomial& powered_polynomial) {
//	this->multiplied_powered_polynomial_set_.push_back(powered_polynomial);
//}
//
//double IrrationalFunction::Term::operator()(const MathVector& value_vector) const {
//	double result = this->coefficient_;
//	for (const auto& powered_polynomial : this->multiplied_powered_polynomial_set_)
//		result *= powered_polynomial(value_vector);
//	return result;
//}
//
//size_t IrrationalFunction::Term::domain_dimension(void) const {
//	size_t result = 0;
//	for (const auto& term : this->multiplied_powered_polynomial_set_)
//		result = std::max(result, term.domain_dimension());
//	return result;
//}
//
//size_t IrrationalFunction::Term::order(void) const {
//	size_t result = 0;
//	for (const auto& term : this->multiplied_powered_polynomial_set_)
//		result += term.order();
//	return result;
//}
//
//std::string IrrationalFunction::Term::to_string(void) const {
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
//	for (const auto& term : this->multiplied_powered_polynomial_set_)
//		str += term.to_string();
//
//	return str;
//}
//
//

namespace ms {


	bool is_positive_odd_number(const double val) {
		if (val < 0)
			return false;

		if (val - std::floor(val) == 0)
			return static_cast<size_t>(val) % 2 == 0;
		else
			return false;
	}

	bool is_natural_number(const double val) {
		if (val < 0)
			return false;

		if (val - std::floor(val) == 0)
			return true;
		else
			return false;
	}

}

