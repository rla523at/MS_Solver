#pragma once
#include "Discrete_Solution.h"
#include "Numerical_Flux_Function.h"
#include "Residual.h"

class Inner_Faces_DG
{
public:
    Inner_Faces_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function);

public:
	void calculate_RHS(Residual& residual, Discrete_Solution_DG& discrete_soltuion) const;

protected:
    std::shared_ptr<Numerical_Flux_Function> numerical_flux_function_;
    uint num_inner_faces_;
	std::vector<std::pair<uint,uint>> oc_nc_index_pairs_;
	std::vector<std::pair<Matrix, Matrix>> oc_nc_side_QWs_basis_m_pairs_;
	std::vector<std::vector<Euclidean_Vector>> set_of_normals_;

	//for construction optimization
	static constexpr ushort max_num_equation = 5;
	static constexpr ushort max_num_basis = 200;
	static constexpr ushort max_num_QPs = 200;
	
	mutable Euclidean_Vector numerical_flux_;
	mutable std::vector<Matrix> set_of_numerical_flux_QPs_m_;
	mutable std::array<double, max_num_equation* max_num_basis> residual_values_ = { 0 };
	mutable std::array<Euclidean_Vector, max_num_QPs> ocs_solution_v_at_QPs_ = {};
	mutable std::array<Euclidean_Vector, max_num_QPs> ncs_solution_v_at_QPs_ = {};

};