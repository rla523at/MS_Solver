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