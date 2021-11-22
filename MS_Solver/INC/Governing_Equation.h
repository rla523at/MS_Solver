#pragma once
#include "Configuration.h"
#include "Euclidean_Vector.h"
#include "Matrix.h"

using ushort = unsigned short;

class Governing_Equation
{
public://Query    
    const std::vector<std::string>& get_solution_names(void) const;
    ushort num_equations(void) const;
    ushort num_solutions(void) const;
    ushort space_dimension(void) const;
    
    virtual std::vector<std::vector<double>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const abstract;
    virtual double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const abstract;
    virtual Matrix calculate_physical_flux(const Euclidean_Vector& solution) const abstract;
    virtual void extend_to_solution(Euclidean_Vector& governing_equation_solution) const abstract;
		
protected:
	ushort num_equations_;
    ushort num_solutions_;
    ushort space_dimension_;
    std::vector<std::string> solution_names_;
};

class Scalar_Equation : public Governing_Equation
{
public:
    Scalar_Equation(void);

public:
    void extend_to_solution(Euclidean_Vector& governing_equation_solution) const override {};
};

class Linear_Advection : public Scalar_Equation
{
public://Query
    std::vector<std::vector<double>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;

protected:
    Euclidean_Vector advection_speeds_;
};

class Linear_Advection_2D : public Linear_Advection
{
public:
    Linear_Advection_2D(const double x_advection_speed, const double y_advection_speed);
};

class Linear_Advection_3D : public Linear_Advection
{
public:
    Linear_Advection_3D(const double x_advection_speed, const double y_advection_speed, const double z_advection_speed);
};

class Burgers : public Scalar_Equation
{
public:
    std::vector<std::vector<double>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
    Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
};

class Burgers_2D : public Burgers
{
public:
    Burgers_2D(void);
};

class Burgers_3D : public Burgers
{
public:
    Burgers_3D(void);
};

class Euler_2D : public Governing_Equation
{
public:
	Euler_2D(void);

public://Query
    std::vector<std::vector<double>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Euclidean_Vector>& P0_solutions) const override;
    double calculate_inner_face_maximum_lambda(const Euclidean_Vector& oc_solution, const Euclidean_Vector& nc_solution, const Euclidean_Vector& nomal_vector) const override;
	Matrix calculate_physical_flux(const Euclidean_Vector& solution) const override;
    void extend_to_solution(Euclidean_Vector& governing_equation_solution) const override;

private:
    static constexpr auto gamma_ = 1.4;
};

class Governing_Equation_Factory//static class
{
public:
    static std::shared_ptr<Governing_Equation> make_shared(const Configuration& config);
    static std::unique_ptr<Governing_Equation> make_unique(const Configuration& config);

private:
	Governing_Equation_Factory(void) = delete;
};








//template <ushort space_dimension_>
//class Euler : public GE
//{
//private:
//    Euler(void) = delete;
//
//private:
//    static constexpr ushort num_equation_ = 2 + space_dimension_;
//
//    using This_                 = Euler<space_dimension_>;
//    using Space_Vector_         = Euclidean_Vector<space_dimension_>;
//    using Solution_             = Euclidean_Vector<num_equation_>;
//    using Physical_Flux_        = Matrix<num_equation_, space_dimension_>;
//
//public:
//    //static auto physical_flux(const Solution_& conservative_variable);
//    static auto physical_flux(const Solution_& cvariable) {
//        const auto pvariable = conservative_to_primitive(cvariable);
//        return physical_flux(cvariable, pvariable);
//    }
//    static auto physical_flux(const Solution_& conservative_variable, const Solution_& primitivie_variable);
//    static auto physical_fluxes(const std::vector<Solution_>& conservative_variables, const std::vector<Solution_>& primitive_variables);
//    static std::vector<std::array<double, space_dimension_>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& conservative_variables);
//    static double inner_face_maximum_lambda(const Solution_& oc_primitive_variable, const Solution_& nc_primitive_variable, const Space_Vector_& nomal_vector);
//    
//public:
//    static std::vector<double> pressures(const std::vector<Solution_>& cvariables);
//    static Solution_ conservative_to_primitive(const Solution_& conservative_variable);
//    static constexpr ushort space_dimension(void) { return space_dimension_; };
//    static constexpr ushort num_equation(void) { return num_equation_; };
//    static std::string name(void) { return "Euler" + std::to_string(space_dimension_) + "D"; };
//};
//
//
//namespace ms {
//    template <typename T>
//    inline constexpr bool is_governing_equation = std::is_base_of_v<GE, T>;
//    template <typename T>
//    inline constexpr bool is_SCL = std::is_base_of_v<SCL<T::space_dimension()>, T>;
//    template <typename T>
//    inline constexpr bool is_linear_advection = std::is_same_v<Linear_Advection<T::space_dimension()>, T>;
//    template <typename T>
//    inline constexpr bool is_Euler = std::is_same_v<Euler<T::space_dimension()>, T>;
//}
//
//
////Template Definition Part

//auto Linear_Advection<space_dimension>::physical_fluxes(const std::vector<Solution_>& solutions) {
//    const auto num_solution = solutions.size();    
//    std::vector<This_::Physical_Flux_> physical_fluxes(num_solution);
//
//    for (size_t i = 0; i < num_solution; ++i)
//        physical_fluxes[i] = This_::physical_flux(solutions[i]);
//
//    return physical_fluxes;
//}
//





//auto Burgers<space_dimension>::physical_fluxes(const std::vector<Solution_>& solutions) {
//    const auto num_solution = solutions.size();
//    std::vector<This_::Physical_Flux_> physical_fluxes(num_solution);
//
//    for (size_t i = 0; i < num_solution; ++i)
//        physical_fluxes[i] = This_::physical_flux(solutions[i]);
//
//    return physical_fluxes;
//}
//


////template <ushort space_dimension_>
////auto Euler<space_dimension_>::physical_flux(const Solution_& cvariable) {
////    const auto pvariable = conservative_to_primitive(cvariable);
////    return physical_flux(cvariable, pvariable);
////}
//
//template <ushort space_dimension_>
//auto Euler<space_dimension_>::physical_flux(const Solution_& conservative_variable, const Solution_& primitivie_variable) {
// //    else {
//        const auto rho = conservative_variable[0];
//        const auto rhou = conservative_variable[1];
//        const auto rhov = conservative_variable[2];
//        const auto rhow = conservative_variable[3];
//        const auto rhoE = conservative_variable[4];
//        const auto u = primitivie_variable[0];
//        const auto v = primitivie_variable[1];
//        const auto w = primitivie_variable[2];
//        const auto p = primitivie_variable[3];
//        const auto a = primitivie_variable[4];
//
//        const auto rhouv = rhou * v;
//        const auto rhouw = rhou * w;
//        const auto rhovw = rhov * w;
//
//        dynamic_require(rho >= 0 && p >= 0, "density and pressure shold be positive");
//
//        return Matrix<num_equation_, space_dimension_>({
//            rhou,				rhov,               rhow,
//            rhou * u + p,		rhouv,              rhouw,
//            rhouv,				rhov * v + p,       rhovw,
//            rhouw,              rhovw,              rhow * w + p,
//            (rhoE + p) * u,		(rhoE + p) * v,      (rhoE + p) * w
//        });
//    }
//}
//
//template <ushort space_dimension_>
//auto Euler<space_dimension_>::physical_fluxes(const std::vector<Solution_>& conservative_variables, const std::vector<Solution_>& primitive_variables) {
//    const size_t num_solution = conservative_variables.size();
//
//    std::vector<Physical_Flux_> physical_fluxes;
//    physical_fluxes.reserve(num_solution);
//
//    for (size_t i = 0; i < num_solution; ++i)
//        physical_fluxes.push_back(physical_flux(conservative_variables[i], primitive_variables[i]));
//
//    return physical_fluxes;
//}
//
//template <ushort space_dimension_>
//std::vector<std::array<double, space_dimension_>> Euler<space_dimension_>::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution_>& conservative_variables) {
    //auto num_solution = conservative_variables.size();

    //std::vector<Solution_> primitive_variables;
    //primitive_variables.reserve(num_solution);

    //for (size_t i = 0; i < num_solution; ++i)
    //    primitive_variables.push_back(Euler<space_dimension_>::conservative_to_primitive(conservative_variables[i]));

    //std::vector<std::array<double, space_dimension_>> coordinate_projected_maximum_lambdas(num_solution);

    //for (size_t i = 0; i < num_solution; ++i) {
    //    if constexpr (space_dimension_ == 2) {
    //        const auto u = primitive_variables[i].at(0);
    //        const auto v = primitive_variables[i].at(1);
    //        const auto a = primitive_variables[i].at(3);

    //        const auto x_projected_maximum_lambda = std::abs(u) + a;
    //        const auto y_projected_maximum_lambda = std::abs(v) + a;

    //        coordinate_projected_maximum_lambdas[i] = { x_projected_maximum_lambda, y_projected_maximum_lambda };
    //    }
    //    else {
    //        const auto u = primitive_variables[i].at(0);
    //        const auto v = primitive_variables[i].at(1);
    //        const auto w = primitive_variables[i].at(2);
    //        const auto a = primitive_variables[i].at(4);

    //        const auto x_projected_maximum_lambda = std::abs(u) + a;
    //        const auto y_projected_maximum_lambda = std::abs(v) + a;
    //        const auto z_projected_maximum_lambda = std::abs(w) + a;

    //        coordinate_projected_maximum_lambdas[i] = { x_projected_maximum_lambda, y_projected_maximum_lambda, z_projected_maximum_lambda };
    //    }
    //}

    //return coordinate_projected_maximum_lambdas;
//}
//
//
//template <ushort space_dimension_>
//double Euler<space_dimension_>::inner_face_maximum_lambda(const Solution_& oc_primitive_variable, const Solution_& nc_primitive_variable, const Space_Vector_& normal_vector) {
//    double oc_side_face_maximum_lambda = 0;
//    double nc_side_face_maximum_lambda = 0;

//    else {
//        const auto oc_u = oc_primitive_variable[0];
//        const auto oc_v = oc_primitive_variable[1];
//        const auto oc_w = oc_primitive_variable[2];
//        const auto oc_a = oc_primitive_variable[4];
//        oc_side_face_maximum_lambda = std::abs(oc_u * normal_vector[0] + oc_v * normal_vector[1] + oc_w * normal_vector[2]) + oc_a;
//
//        const auto nc_u = nc_primitive_variable[0];
//        const auto nc_v = nc_primitive_variable[1];
//        const auto nc_w = nc_primitive_variable[2];
//        const auto nc_a = nc_primitive_variable[4];
//        nc_side_face_maximum_lambda = std::abs(nc_u * normal_vector[0] + nc_v * normal_vector[1] + nc_w * normal_vector[2]) + nc_a;
//    }
//
//    return (std::max)(oc_side_face_maximum_lambda, nc_side_face_maximum_lambda);
//}
//
