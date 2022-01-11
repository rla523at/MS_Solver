#include "../INC/Inner_Faces.h"

Inner_Faces_DG::Inner_Faces_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function)
	: numerical_flux_function_(numerical_flux_function)
{
	Profiler::set_time_point();

	this->num_inner_faces_ = static_cast<uint>(grid.num_inner_faces());

	this->oc_nc_index_pairs_.resize(this->num_inner_faces_);
	this->oc_nc_side_QWs_basis_m_pairs_.resize(this->num_inner_faces_);
	this->set_of_normals_.resize(this->num_inner_faces_);

	//for precalculation
	std::vector<uint> oc_indexes(this->num_inner_faces_);
	std::vector<uint> nc_indexes(this->num_inner_faces_);
	std::vector<std::vector<Euclidean_Vector>> set_of_ocs_QPs(this->num_inner_faces_);
	std::vector<std::vector<Euclidean_Vector>> set_of_ncs_QPs(this->num_inner_faces_);

	//for construct optimization
	const auto num_equations = discrete_solution.num_equations();
	const auto num_solutions = discrete_solution.num_solutions();

	this->numerical_flux_ = Euclidean_Vector(num_equations);
	this->set_of_numerical_flux_QPs_m_.resize(this->num_inner_faces_);
	this->ocs_solution_v_at_QPs_.fill(Euclidean_Vector(num_solutions));
	this->ncs_solution_v_at_QPs_.fill(Euclidean_Vector(num_solutions));

	// consider inner face 
	for (uint infc_index = 0; infc_index < this->num_inner_faces_; ++infc_index)
	{
		const auto [oc_index, nc_index] = grid.inner_face_oc_nc_index_pair(infc_index);
		this->oc_nc_index_pairs_[infc_index] = { oc_index,nc_index };

		const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
		const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
		const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);
		const auto integrand_degree = 2 * max_solution_degree + 1;

		const auto& [ocs_quadrature_rule, ncs_quadrature_rule] = grid.inner_face_quadrature_rule(infc_index, integrand_degree);

		const auto& ocs_QPs = ocs_quadrature_rule.points;
		const auto& ocs_QWs = ocs_quadrature_rule.weights;
		const auto& ncs_QPs = ncs_quadrature_rule.points;
		const auto& ncs_QWs = ncs_quadrature_rule.weights;
		const auto num_QPs = ocs_QPs.size();

		const auto oc_num_basis = discrete_solution.num_basis(oc_index);
		const auto nc_num_basis = discrete_solution.num_basis(nc_index);

		Matrix oc_side_QWs_basis(num_QPs, oc_num_basis);
		Matrix nc_side_QWs_basis(num_QPs, nc_num_basis);

		for (uint q = 0; q < num_QPs; ++q)
		{
			oc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(oc_index, ocs_QPs[q]) * ocs_QWs[q]);
			nc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(nc_index, ncs_QPs[q]) * ncs_QWs[q]);
		}
		oc_side_QWs_basis *= -1.0;//owner side

		this->oc_nc_side_QWs_basis_m_pairs_[infc_index] = { std::move(oc_side_QWs_basis), std::move(nc_side_QWs_basis) };
		this->set_of_normals_[infc_index] = grid.inner_face_normals(infc_index, oc_index, ocs_QPs);

		//for construct optimization
		this->set_of_numerical_flux_QPs_m_[infc_index] = Matrix(num_equations, num_QPs);

		//for precalculation
		oc_indexes[infc_index] = oc_index;
		nc_indexes[infc_index] = nc_index;
		set_of_ocs_QPs[infc_index] = ocs_quadrature_rule.points;
		set_of_ncs_QPs[infc_index] = ncs_quadrature_rule.points;
	}

	//for precaclulation
	discrete_solution.precalculate_infs_ocs_RHS_QPs_basis_values(oc_indexes, set_of_ocs_QPs);
	discrete_solution.precalculate_infs_ncs_RHS_QPs_basis_values(nc_indexes, set_of_ncs_QPs);

	LOG << std::left << std::setw(50) << "@ Inner faces DG precalculation" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << LOG.print_;
}

void Inner_Faces_DG::calculate_RHS(Residual& residual, Discrete_Solution_DG& discrete_solution) const
{
	for (uint infc_index = 0; infc_index < this->num_inner_faces_; ++infc_index)
	{
		const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[infc_index];
		discrete_solution.calculate_solution_at_infc_ocs_RHS_QPs(this->ocs_solution_v_at_QPs_.data(), infc_index, oc_index);
		discrete_solution.calculate_solution_at_infc_ncs_RHS_QPs(this->ncs_solution_v_at_QPs_.data(), infc_index, nc_index);

		const auto& normals = this->set_of_normals_[infc_index];
		const auto num_QPs = normals.size();

		auto& numerical_flux_QPs_m = this->set_of_numerical_flux_QPs_m_[infc_index];

		for (ushort q = 0; q < num_QPs; ++q)
		{
			this->numerical_flux_function_->calculate(this->numerical_flux_.data(), this->ocs_solution_v_at_QPs_[q], this->ncs_solution_v_at_QPs_[q], normals[q]);
			numerical_flux_QPs_m.change_column(q, this->numerical_flux_);
		}

		const auto& [oc_side_QWs_basis_m, nc_side_QWs_basis_m] = this->oc_nc_side_QWs_basis_m_pairs_[infc_index];

		std::fill(this->residual_values_.begin(), this->residual_values_.begin() + discrete_solution.num_values(oc_index), 0.0);
		ms::gemm(numerical_flux_QPs_m, oc_side_QWs_basis_m, this->residual_values_.data());
		residual.update_rhs(oc_index, this->residual_values_.data());

		std::fill(this->residual_values_.begin(), this->residual_values_.begin() + discrete_solution.num_values(nc_index), 0.0);
		ms::gemm(numerical_flux_QPs_m, nc_side_QWs_basis_m, this->residual_values_.data());
		residual.update_rhs(nc_index, this->residual_values_.data());
	}
}