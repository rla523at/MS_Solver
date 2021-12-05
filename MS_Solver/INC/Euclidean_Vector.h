#pragma once
#include "Exception.h"

#include <array>
#include <iomanip>
#include <mkl.h>
#include <sstream>
#include <vector>

using ushort = unsigned short;


class Euclidean_Vector;
class Euclidean_Vector_Constant_Base
{
public:
	Euclidean_Vector_Constant_Base(void) = default;
	Euclidean_Vector_Constant_Base(const size_t num_value, const double* const_data_ptr);

public://Query
	Euclidean_Vector operator-(const Euclidean_Vector_Constant_Base& other) const;
	Euclidean_Vector operator+(const Euclidean_Vector_Constant_Base& other) const;
	Euclidean_Vector operator*(const double constant) const;
	double operator[](const size_t position) const;
	bool operator==(const Euclidean_Vector_Constant_Base& other) const;

	double at(const size_t position) const;
	const double* begin(void) const;
	std::vector<double> copy_values(void) const;
	const double* data(void) const;
	const double* end(void) const;
	double L1_norm(void) const;
	double L2_norm(void) const;
	double Linf_norm(void) const;
	double inner_product(const Euclidean_Vector_Constant_Base& other) const;
	bool is_axis_translation(const Euclidean_Vector_Constant_Base& other) const;
	size_t size(void) const;
	std::string to_string(void) const;

protected:
	int num_values_ = 0;
	const double* const_data_ptr_ = nullptr;
};

class Euclidean_Vector_Base : public Euclidean_Vector_Constant_Base
{
public:
	Euclidean_Vector_Base(void) = default;
	Euclidean_Vector_Base(const size_t num_values, double* data_ptr)
		: Euclidean_Vector_Constant_Base(num_values, data_ptr)
		, data_ptr_(data_ptr) {};
		
public://Command
	void operator*=(const double constant);
	void operator+=(const Euclidean_Vector_Constant_Base& other);

	void normalize(void);

protected:
	double* data_ptr_ = nullptr;
};

class Euclidean_Vector : public Euclidean_Vector_Base
{
public:
	Euclidean_Vector(void) = default;
	explicit Euclidean_Vector(const size_t size);
	Euclidean_Vector(const std::initializer_list<double> list);
	Euclidean_Vector(std::vector<double>&& values);
	Euclidean_Vector(const std::vector<double>& values);
	template <typename Iter>	Euclidean_Vector(Iter first, Iter last) 
	{
		this->num_values_ = static_cast<int>(last - first);

		if (this->num_values_ <= this->small_criterion_)
		{
			std::copy(first, last, this->small_buffer_.begin());
			this->const_data_ptr_ = this->small_buffer_.data();
			this->data_ptr_ = this->small_buffer_.data();
		}
		else
		{
			std::copy(first, last, this->values_.begin());
			this->const_data_ptr_ = this->values_.data();
			this->data_ptr_ = this->values_.data();
		}
	};
	Euclidean_Vector(const Euclidean_Vector& other);
	Euclidean_Vector(Euclidean_Vector&& other) noexcept;

public://Command	
	void operator=(const Euclidean_Vector& other);
	void operator=(Euclidean_Vector&& other) noexcept;

	std::vector<double>&& move_values(void);

private:
	//bool is_small(void) const;

private:
	static constexpr ushort small_criterion_ = 10;
	std::array<double, small_criterion_> small_buffer_ = { 0 };
	std::vector<double> values_;
};

//Wrapper class 안전하게 사용 하는 방법
//1. 매번 새로 생성해서 사용하기
//2. 상속 말고 합성으로 만들고 매번 검사하기 또는 반영하기 (V)
//3. 상속 안할시 upcasting이 안되서 인자로 Wrapper를 못 받는 상황이 발생
//4. 상속과 합성을 동시에 사용

class Euclidean_Vector_Constant_Wrapper : public Euclidean_Vector_Constant_Base
{
public:
	Euclidean_Vector_Constant_Wrapper(const std::vector<double>& values)
		: values_wrapper_(values)
		, base_(this->values_wrapper_.size(), this->values_wrapper_.data()) 
	{
		this->num_values_ = static_cast<int>(this->values_wrapper_.size());
		this->const_data_ptr_ = this->values_wrapper_.data();
	};

public://Query
	Euclidean_Vector operator-(const Euclidean_Vector_Constant_Base& other) const;
	Euclidean_Vector operator+(const Euclidean_Vector_Constant_Base& other) const;
	Euclidean_Vector operator*(const double constant) const;
	double operator[](const size_t position) const;
	bool operator==(const Euclidean_Vector_Constant_Base& other) const;

	double at(const size_t position) const;
	const double* begin(void) const;
	std::vector<double> copy_values(void) const;
	const double* data(void) const;
	const double* end(void) const;
	double L1_norm(void) const;
	double L2_norm(void) const;
	double Linf_norm(void) const;
	double inner_product(const Euclidean_Vector_Constant_Base& other) const;
	bool is_axis_translation(const Euclidean_Vector_Constant_Base& other) const;
	size_t size(void) const;
	std::string to_string(void) const;

private:
	bool is_sync(void) const;

private:
	const std::vector<double>& values_wrapper_;
	Euclidean_Vector_Constant_Base base_;
};

class Euclidean_Vector_Wrapper : public Euclidean_Vector_Base
{
public:
	Euclidean_Vector_Wrapper(std::vector<double>& values)
		: values_wrapper_(values)
		, base_(this->values_wrapper_.size(), this->values_wrapper_.data())
	{
		this->num_values_ = static_cast<int>(this->values_wrapper_.size());
		this->const_data_ptr_ = this->values_wrapper_.data();
		this->data_ptr_ = this->values_wrapper_.data();
	};

public://Command
	void operator*=(const double constant);
	void operator+=(const Euclidean_Vector_Constant_Base& other);
	void operator=(Euclidean_Vector&& other) noexcept;

	void normalize(void);

public://Query
	Euclidean_Vector operator-(const Euclidean_Vector_Constant_Base& other) const;
	Euclidean_Vector operator+(const Euclidean_Vector_Constant_Base& other) const;
	Euclidean_Vector operator*(const double constant) const;
	double operator[](const size_t position) const;
	bool operator==(const Euclidean_Vector_Constant_Base& other) const;

	double at(const size_t position) const;
	const double* begin(void) const;
	std::vector<double> copy_values(void) const;
	const double* data(void) const;
	const double* end(void) const;
	double L1_norm(void) const;
	double L2_norm(void) const;
	double Linf_norm(void) const;
	double inner_product(const Euclidean_Vector_Constant_Base& other) const;
	bool is_axis_translation(const Euclidean_Vector_Constant_Base& other) const;
	size_t size(void) const;
	std::string to_string(void) const;

private:
	bool is_sync(void) const;

private:
	std::vector<double>& values_wrapper_;
	Euclidean_Vector_Base base_;
};

template <typename Function>
class Vector_Function
{
public:
	Vector_Function(void) = default;
	Vector_Function(std::vector<Function>&& functions) : functions_(std::move(functions)) {};
	Vector_Function(const std::initializer_list<Function> list) : functions_(list) {};

public://Query
	template <typename V>	V operator()(const V& point_v) const
	{
		const auto range_dimension = this->size();

		std::vector<double> result(range_dimension);
		for (size_t i = 0; i < range_dimension; ++i)
		{
			result[i] = this->functions_[i](point_v);
		}

		return result;
	}
	const Function& operator[](const size_t index) const
	{
		REQUIRE(index < this->size(), "index can not exceed range size");
		return this->functions_[index];
	}
	bool operator==(const Vector_Function& other) const
	{
		return this->functions_ == other.functions_;
	}

	const Function& at(const size_t index) const
	{
		REQUIRE(index < this->size(), "index can not exceed range size");
		return this->functions_[index];
	}
	Vector_Function<Function> cross_product(const Vector_Function& other) const
	{
		constexpr auto result_range_dimension = 3;
		std::vector<Function> result(result_range_dimension);

		if (this->size() == 2)
		{
			result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);
		}
		else if (this->size() == 3)
		{
			result[0] = this->at(1) * other.at(2) - this->at(2) * other.at(1);
			result[1] = this->at(2) * other.at(0) - this->at(0) * other.at(2);
			result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);
		}
		else
			EXCEPTION("cross product only valid for R^3 dimension");

		return result;
	};
	auto L2_norm(void) const
	{
		Function result = 0.0;

		for (const auto& function : this->functions_)
		{
			result += (function ^ 2);
		}

		return result.root(0.5);
	}
	Vector_Function<Function> get_differentiate(const ushort variable_index) const
	{
		auto differentiate_functions = this->functions_;

		for (auto& differentiate_function : differentiate_functions)
		{
			differentiate_function.differentiate(variable_index);
		}
		return differentiate_functions;
	}
	size_t size(void) const
	{
		return this->functions_.size();
	}
	std::string to_string(void) const
	{
		std::string result;

		for (const auto& function : this->functions_)
		{
			result += function.to_string() + "\t";
		}

		return result;
	}

private:
	std::vector<Function> functions_;
};

Euclidean_Vector operator*(const double constant, const Euclidean_Vector_Constant_Base& x);
std::ostream& operator<<(std::ostream& os, const Euclidean_Vector_Constant_Base& x);
template <typename Function>	std::ostream& operator<<(std::ostream& os, const Vector_Function<Function>& vf)
{
	return os << vf.to_string();
}