#pragma once
#include "Discrete_Solution.h"
#include "Numerical_Flux_Function.h"
#include "Residual.h"

class Inner_Faces_DG
{
public:
    Inner_Faces_DG(const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function, const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
	void calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_soltuion) const;

protected:
    std::shared_ptr<Numerical_Flux_Function> numerical_flux_function_;
    uint num_inner_faces_;
	ushort num_equations_;
	std::vector<std::pair<uint,uint>> oc_nc_index_pairs_;
	std::vector<std::pair<Matrix, Matrix>> oc_nc_side_QWs_basis_m_pairs_;
	std::vector<std::vector<Euclidean_Vector>> set_of_normals_;

	//for construction optimization
	static constexpr ushort max_num_equation = 5;
	static constexpr ushort max_num_basis = 200;
	
	std::vector<ushort> set_of_num_QPs_;
	mutable std::array<double, max_num_equation* max_num_basis> residual_values_;
};