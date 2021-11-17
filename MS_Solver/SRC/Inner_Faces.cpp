#include "../INC/Inner_Faces.h"

Inner_Faces_DG::Inner_Faces_DG(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function)
    : numerical_flux_function_(numerical_flux_function)
{
    SET_TIME_POINT;

    this->num_equations_ = discrete_solution.num_equations();

    const auto num_inner_faces = grid.num_inner_faces();
    const auto num_pbdry_pairs = grid.num_periodic_boundary_pairs();    
    this->num_inner_faces_ = num_inner_faces + num_pbdry_pairs; //periodic boundary can be seen as inner face

    this->oc_nc_index_pairs_.reserve(this->num_inner_faces_);
    this->oc_nc_side_basis_QPs_m_pairs_.reserve(this->num_inner_faces_);
    this->oc_nc_side_QWs_basis_m_pairs_.reserve(this->num_inner_faces_);
    this->set_of_normals_.reserve(this->num_inner_faces_);

    // consider inner face 
    for (int infc_index = 0; infc_index < num_inner_faces; ++infc_index)
    {
        //set oc nc index pair
        const auto [oc_index, nc_index] = grid.inner_face_oc_nc_index_pair(infc_index);
        this->oc_nc_index_pairs_.push_back({ oc_index,nc_index });

        //set oc nc side basis QPs m & oc nc side QWs basis m
        const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
        const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
        const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

        const auto quadrature_rule = grid.inner_face_quadrature_rule(infc_index, max_solution_degree);
        const auto& QPs = quadrature_rule.points;
        const auto& QWs = quadrature_rule.weights;
        const auto num_QPs = QPs.size();

        auto oc_side_basis_QPs_m = discrete_solution.calculate_basis_points_m(oc_index, QPs);
        auto nc_side_basis_QPs_m = discrete_solution.calculate_basis_points_m(nc_index, QPs);
        const auto oc_num_basis = oc_side_basis_QPs_m.num_row();
        const auto nc_num_basis = nc_side_basis_QPs_m.num_row();

        Matrix oc_side_QWs_basis(num_QPs, oc_num_basis);
        Matrix nc_side_QWs_basis(num_QPs, nc_num_basis);

        for (int q = 0; q < num_QPs; ++q)
        {
            oc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(oc_index, QPs[q]) * QWs[q]);
            nc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(nc_index, QPs[q]) * QWs[q]);
        }
        oc_side_QWs_basis *= -1.0;//owner side

        this->oc_nc_side_basis_QPs_m_pairs_.push_back({ std::move(oc_side_basis_QPs_m), std::move(nc_side_basis_QPs_m) });
        this->oc_nc_side_QWs_basis_m_pairs_.push_back({ std::move(oc_side_QWs_basis), std::move(nc_side_QWs_basis) });

        //set normals
        this->set_of_normals_.push_back(grid.inner_face_normals(infc_index, oc_index, QPs));
    }

    // consider periodic boundary
    for (int pbdry_pair_index = 0; pbdry_pair_index < num_pbdry_pairs; ++pbdry_pair_index)
    {
        //set oc nc index pair
        const auto [oc_index, nc_index] = grid.periodic_boundary_oc_nc_index_pair(pbdry_pair_index);
        this->oc_nc_index_pairs_.push_back({ oc_index,nc_index });

        //set oc nc side basis QPs m & oc nc side QWs basis m
        const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
        const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
        const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

        const auto [oc_quadrature_rule, nc_quadrature_rule] = grid.periodic_boundary_quadrature_rule_pair(pbdry_pair_index, max_solution_degree);
        const auto& oc_QPs = oc_quadrature_rule.points;
        const auto& oc_QWs = oc_quadrature_rule.weights;
        const auto& nc_QPs = nc_quadrature_rule.points;
        const auto& nc_QWs = nc_quadrature_rule.weights;
        const auto num_QPs = oc_QPs.size();

        auto oc_side_basis_QPs_m = discrete_solution.calculate_basis_points_m(oc_index, oc_QPs);
        auto nc_side_basis_QPs_m = discrete_solution.calculate_basis_points_m(nc_index, nc_QPs);
        const auto oc_num_basis = oc_side_basis_QPs_m.num_row();
        const auto nc_num_basis = nc_side_basis_QPs_m.num_row();

        Matrix oc_side_QWs_basis(num_QPs, oc_num_basis);
        Matrix nc_side_QWs_basis(num_QPs, nc_num_basis);

        for (int q = 0; q < num_QPs; ++q)
        {
            oc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(oc_index, oc_QPs[q]) * oc_QWs[q]);
            nc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(nc_index, nc_QPs[q]) * nc_QWs[q]);
        }
        oc_side_QWs_basis *= -1.0;//owner side

        this->oc_nc_side_basis_QPs_m_pairs_.push_back({ std::move(oc_side_basis_QPs_m), std::move(nc_side_basis_QPs_m) });
        this->oc_nc_side_QWs_basis_m_pairs_.push_back({ std::move(oc_side_QWs_basis), std::move(nc_side_QWs_basis) });

        //set normals
        this->set_of_normals_.push_back(grid.inner_face_normals(pbdry_pair_index, oc_index, oc_QPs));
    }

    LOG << std::left << std::setw(50) << "@ Inner faces DG precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n" << LOG.print_;
}

void Inner_Faces_DG::calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_soltuion) const
{
    for (int i = 0; i < this->num_inner_faces_; ++i)
    {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& [oc_side_basis_QPs_m, nc_side_basis_QPs_m] = this->oc_nc_side_basis_QPs_m_pairs_[i];

        const auto oc_side_solution_at_QPs = discrete_soltuion.calculate_solution_at_points(oc_index, oc_side_basis_QPs_m);
        const auto nc_side_solution_at_QPs = discrete_soltuion.calculate_solution_at_points(nc_index, nc_side_basis_QPs_m);
        const auto num_QPs = oc_side_solution_at_QPs.size();

        const auto& normals = this->set_of_normals_[i];

        Matrix numerical_flux_quadrature(this->num_equations_, num_QPs);

        for (ushort q = 0; q < num_QPs; ++q)
        {
            numerical_flux_quadrature.change_column(q, this->numerical_flux_function_->calculate(oc_side_solution_at_QPs[q], nc_side_solution_at_QPs[q], normals[q]));
        }

        const auto& [oc_side_QWs_basis_m, nc_side_QWs_basis_m] = this->oc_nc_side_QWs_basis_m_pairs_[i];
        const auto oc_side_delta_rhs = numerical_flux_quadrature * oc_side_QWs_basis_m;
        const auto nc_side_delta_rhs = numerical_flux_quadrature * nc_side_QWs_basis_m;

        residual.update_rhs(oc_index, oc_side_delta_rhs);
        residual.update_rhs(nc_index, nc_side_delta_rhs);
    }
}