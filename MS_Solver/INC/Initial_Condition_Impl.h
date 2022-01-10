#pragma once
#include "Configuration.h"
#include "Initial_Condition.h"

class Constant1 : public Initial_Condition
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Sine_Wave : public Initial_Condition
{
public:
	Sine_Wave(const std::vector<double>& wave_lengths);

public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;

private:
	ushort space_dimension_ = 0;
	std::vector<double> wave_numbers_;// sin(kx)sin(ky)sin(kz)
};

class Sine_Wave_2D : public Sine_Wave
{
public:
	Sine_Wave_2D(const double x_wave_length, const double y_wave_length)
		: Sine_Wave({ x_wave_length,y_wave_length }) {};
	Sine_Wave_2D(const double* wave_lengths_ptr)
		: Sine_Wave({ wave_lengths_ptr,wave_lengths_ptr + 2 }) {};
};

class Sine_Wave_3D : public Sine_Wave
{
public:
	Sine_Wave_3D(const double x_wave_length, const double y_wave_length, const double z_wave_length)
		: Sine_Wave({ x_wave_length,y_wave_length,z_wave_length }) {};
	Sine_Wave_3D(const double* wave_lengths_ptr)
		: Sine_Wave({ wave_lengths_ptr,wave_lengths_ptr + 3 }) {};
};

class Square_Wave : public Initial_Condition
{
public:
	Square_Wave(const ushort space_dimension)
		:space_dimension_(space_dimension) {};

public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;

protected:
	ushort space_dimension_ = 0;
};

class Square_Wave_2D : public Square_Wave
{
public:
	Square_Wave_2D(void)
		:Square_Wave(2) {};
};

class Square_Wave_3D : public Square_Wave
{
public:
	Square_Wave_3D(void)
		:Square_Wave(3) {};
};

class Circle_Wave : public Initial_Condition
{
public:
	Circle_Wave(const ushort space_dimension);

public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;

private:
	ushort space_dimension_ = 0;
	Euclidean_Vector center_v_;
};

class Circle_Wave_2D : public Circle_Wave
{
public:
	Circle_Wave_2D(void)
		:Circle_Wave(2) {};
};

class Circle_Wave_3D : public Circle_Wave
{
public:
	Circle_Wave_3D(void)
		:Circle_Wave(3) {};
};

class Gaussian_Wave : public Initial_Condition
{
public:
	Gaussian_Wave(const ushort space_dimension);

public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;

private:
	ushort space_dimension_ = 0;
	Euclidean_Vector center_v_;
};

class Gaussian_Wave_2D : public Gaussian_Wave
{
public:
	Gaussian_Wave_2D(void)
		:Gaussian_Wave(2) {};
};

class Gaussian_Wave_3D : public Gaussian_Wave
{
public:
	Gaussian_Wave_3D(void)
		:Gaussian_Wave(3) {};
};

class Euler_Shocktube_2D
{
protected:
	Euclidean_Vector primitive_to_conservative(const Euclidean_Vector& primitive_variable) const;
};

class Euler_Shocktube_3D
{
protected:
	Euclidean_Vector primitive_to_conservative(const Euclidean_Vector& primitive_variable) const;
};

class SOD_2D : public Initial_Condition, Euler_Shocktube_2D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class SOD_3D : public Initial_Condition, Euler_Shocktube_3D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Modified_SOD_2D : public Initial_Condition, Euler_Shocktube_2D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.2; };
};

class Modified_SOD_3D : public Initial_Condition, Euler_Shocktube_3D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.2; };
};

class Supersonic_Expansion_2D : public Initial_Condition
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Supersonic_Expansion_3D : public Initial_Condition
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Left_Half_Blast_Wave_2D : public Initial_Condition, Euler_Shocktube_2D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.012; };
};

class Left_Half_Blast_Wave_3D : public Initial_Condition, Euler_Shocktube_3D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.012; };
};

class Double_Strong_Shock_2D : public Initial_Condition, Euler_Shocktube_2D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.035; };
};

class Double_Strong_Shock_3D : public Initial_Condition, Euler_Shocktube_3D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.035; };
};

class Slowly_Moving_Contact_2D : public Initial_Condition, Euler_Shocktube_2D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.012; };
};

class Slowly_Moving_Contact_3D : public Initial_Condition, Euler_Shocktube_3D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.012; };
};

class Harten_Lax_2D : public Initial_Condition, Euler_Shocktube_2D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Harten_Lax_3D : public Initial_Condition, Euler_Shocktube_3D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Shu_Osher_2D : public Initial_Condition, Euler_Shocktube_2D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.178; };
};

class Shu_Osher_3D : public Initial_Condition, Euler_Shocktube_3D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
	constexpr double target_end_time(void) const override { return 0.178; };
};

class Blast_Wave_2D : public Initial_Condition, Euler_Shocktube_2D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Blast_Wave_3D : public Initial_Condition, Euler_Shocktube_3D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Explosion_2D : public Initial_Condition, Euler_Shocktube_2D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Explosion_3D : public Initial_Condition, Euler_Shocktube_3D
{
public:
	Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const override;
};

class Initial_Condition_Factory//static class
{
public:
	static std::unique_ptr<Initial_Condition> make_unique(const Configuration& configuration);
};