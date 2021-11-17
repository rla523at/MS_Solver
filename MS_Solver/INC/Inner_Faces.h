#pragma once
#include "Discrete_Solution.h"
#include "Numerical_Flux_Function.h"
#include "Residual.h"

class Inner_Faces_DG
{
    Inner_Faces_DG(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function);

public:
	void calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_soltuion) const;

protected:
    std::shared_ptr<Numerical_Flux_Function> numerical_flux_function_;
    uint num_inner_faces_;
	ushort num_equations_;
	std::vector<std::pair<uint,uint>> oc_nc_index_pairs_;
	std::vector<std::pair<Matrix, Matrix>> oc_nc_side_basis_QPs_m_pairs_;
	std::vector<std::pair<Matrix, Matrix>> oc_nc_side_QWs_basis_m_pairs_;
	std::vector<std::vector<Euclidean_Vector>> set_of_normals_;
};


//#pragma once
//#include "Inner_Faces_FVM.h"
//#include "Inner_Faces_HOM.h"
//#include "Spatial_Discrete_Method.h"
//
//
//template <typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
//class Inner_Faces;
//
//
//template <typename Numerical_Flux_Function>
//class Inner_Faces<FVM, Constant_Reconstruction, Numerical_Flux_Function> : public Inner_Faces_FVM_Constant<Numerical_Flux_Function>
//{
//private:
//    static constexpr size_t space_dimension_ = Numerical_Flux_Function::space_dimension();
//
//public:
//    Inner_Faces(const Grid<space_dimension_>& grid, const Constant_Reconstruction& reconstruction_method)
//        : Inner_Faces_FVM_Constant<Numerical_Flux_Function>(grid) {};
//};
//
//
//template<typename Reconstruction_Method, typename Numerical_Flux_Function>
//class Inner_Faces<FVM, Reconstruction_Method, Numerical_Flux_Function> : public Inner_Faces_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>
//{
//private:
//    static constexpr size_t space_dimension_ = Numerical_Flux_Function::space_dimension();
//
//public:
//    Inner_Faces(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
//        : Inner_Faces_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>(grid, reconstruction_method) {};
//};
//
//
//template<typename Reconstruction_Method, typename Numerical_Flux_Function>
//class Inner_Faces<HOM, Reconstruction_Method, Numerical_Flux_Function> : public Inner_Faces_HOM<Reconstruction_Method, Numerical_Flux_Function>
//{
//private:
//    static constexpr size_t space_dimension_ = Numerical_Flux_Function::space_dimension();
//
//public:
//    Inner_Faces(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
//        : Inner_Faces_HOM<Reconstruction_Method, Numerical_Flux_Function>(grid, reconstruction_method) {};
//};