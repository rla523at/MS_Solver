#include "../INC/Inner_Faces.h"

Inner_Faces_DG::Inner_Faces_DG(const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function, const Grid& grid, Discrete_Solution_DG& discrete_solution)
    : numerical_flux_function_(numerical_flux_function)
{
    Profiler::set_time_point();


    const auto num_inner_faces = grid.num_inner_faces();
    const auto num_pbdry_pairs = grid.num_periodic_boundary_pairs();    
    this->num_inner_faces_ = static_cast<uint>(num_inner_faces + num_pbdry_pairs); //periodic boundary can be seen as inner face

    this->oc_nc_index_pairs_.reserve(this->num_inner_faces_);
    this->oc_nc_side_QWs_basis_m_pairs_.reserve(this->num_inner_faces_);
    this->set_of_normals_.reserve(this->num_inner_faces_);

    //for precalculation
    std::vector<uint> oc_indexes;
    std::vector<uint> nc_indexes;
    std::vector<std::vector<Euclidean_Vector>> set_of_ocs_QPs;
    std::vector<std::vector<Euclidean_Vector>> set_of_ncs_QPs;

    oc_indexes.reserve(this->num_inner_faces_);
    nc_indexes.reserve(this->num_inner_faces_);
    set_of_ocs_QPs.reserve(this->num_inner_faces_);
    set_of_ncs_QPs.reserve(this->num_inner_faces_);

    //for construct optimization
    const auto num_equations = discrete_solution.num_equations();
    const auto num_solutions = discrete_solution.num_solutions();

    this->set_of_numerical_flux_QPs_m_.reserve(this->num_inner_faces_);    
    this->ocs_solution_v_at_QPs_.fill(Euclidean_Vector(num_solutions));
    this->ncs_solution_v_at_QPs_.fill(Euclidean_Vector(num_solutions));
    //this->ocs_solution_v_at_QPs_.resize(this->max_num_QPs, Euclidean_Vector(num_solutions));
    //this->ncs_solution_v_at_QPs_.resize(this->max_num_QPs, Euclidean_Vector(num_solutions));
    //

    // consider inner face 
    for (uint infc_index = 0; infc_index < num_inner_faces; ++infc_index)
    {        
        const auto [oc_index, nc_index] = grid.inner_face_oc_nc_index_pair(infc_index);
        this->oc_nc_index_pairs_.push_back({ oc_index,nc_index });
                
        const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
        const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
        const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);
        const auto integrand_degree = 2 * max_solution_degree + 1;

        auto quadrature_rule = grid.inner_face_quadrature_rule(infc_index, integrand_degree);
        const auto& QPs = quadrature_rule.points;
        const auto& QWs = quadrature_rule.weights;
        const auto num_QPs = QPs.size();

        //for construct optimization
        this->set_of_numerical_flux_QPs_m_.push_back({ num_equations,num_QPs });

        const auto oc_num_basis = discrete_solution.num_basis(oc_index);
        const auto nc_num_basis = discrete_solution.num_basis(nc_index);

        Matrix oc_side_QWs_basis(num_QPs, oc_num_basis);
        Matrix nc_side_QWs_basis(num_QPs, nc_num_basis);

        for (uint q = 0; q < num_QPs; ++q)
        {
            oc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(oc_index, QPs[q]) * QWs[q]);
            nc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(nc_index, QPs[q]) * QWs[q]);
        }
        oc_side_QWs_basis *= -1.0;//owner side

        this->oc_nc_side_QWs_basis_m_pairs_.push_back({ std::move(oc_side_QWs_basis), std::move(nc_side_QWs_basis) });
                
        this->set_of_normals_.push_back(grid.inner_face_normals(infc_index, oc_index, QPs));

        //for precalculation
        oc_indexes.push_back(oc_index);
        nc_indexes.push_back(nc_index);
        set_of_ocs_QPs.push_back(quadrature_rule.points);
        set_of_ncs_QPs.push_back(std::move(quadrature_rule.points));
    }

    // consider periodic boundary
    for (uint pbdry_pair_index = 0; pbdry_pair_index < num_pbdry_pairs; ++pbdry_pair_index)
    {        
        const auto [oc_index, nc_index] = grid.periodic_boundary_oc_nc_index_pair(pbdry_pair_index);
        this->oc_nc_index_pairs_.push_back({ oc_index,nc_index });

        const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
        const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
        const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);
        const auto integrand_degree = 2 * max_solution_degree + 1;

        auto [oc_quadrature_rule, nc_quadrature_rule] = grid.periodic_boundary_quadrature_rule_pair(pbdry_pair_index, integrand_degree);
        const auto& oc_QPs = oc_quadrature_rule.points;
        const auto& oc_QWs = oc_quadrature_rule.weights;
        const auto& nc_QPs = nc_quadrature_rule.points;
        const auto& nc_QWs = nc_quadrature_rule.weights;
        const auto num_QPs = oc_QPs.size();

        //for construct optimization
        this->set_of_numerical_flux_QPs_m_.push_back({ num_equations,num_QPs });

        const auto oc_num_basis = discrete_solution.num_basis(oc_index);
        const auto nc_num_basis = discrete_solution.num_basis(nc_index);

        Matrix oc_side_QWs_basis(num_QPs, oc_num_basis);
        Matrix nc_side_QWs_basis(num_QPs, nc_num_basis);

        for (uint q = 0; q < num_QPs; ++q)
        {
            oc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(oc_index, oc_QPs[q]) * oc_QWs[q]);
            nc_side_QWs_basis.change_row(q, discrete_solution.calculate_basis_point_v(nc_index, nc_QPs[q]) * nc_QWs[q]);
        }
        oc_side_QWs_basis *= -1.0;//owner side

        this->oc_nc_side_QWs_basis_m_pairs_.push_back({ std::move(oc_side_QWs_basis), std::move(nc_side_QWs_basis) });

        this->set_of_normals_.push_back(grid.periodic_boundary_normals(pbdry_pair_index, oc_index, oc_QPs));

        //for precalculation
        oc_indexes.push_back(oc_index);
        nc_indexes.push_back(nc_index);
        set_of_ocs_QPs.push_back(std::move(oc_quadrature_rule.points));
        set_of_ncs_QPs.push_back(std::move(nc_quadrature_rule.points));
    }

    //for precaclulation
    discrete_solution.precalculate_infs_ocs_QPs_basis_values(oc_indexes, set_of_ocs_QPs);
    discrete_solution.precalculate_infs_ncs_QPs_basis_values(nc_indexes, set_of_ncs_QPs);

    LOG << std::left << std::setw(50) << "@ Inner faces DG precalculation" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << LOG.print_;
}

void Inner_Faces_DG::calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_solution) const
{
    for (uint infc_index = 0; infc_index < this->num_inner_faces_; ++infc_index)
    {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[infc_index];
        discrete_solution.calculate_solution_at_infc_ocs_QPs(this->ocs_solution_v_at_QPs_.data(), infc_index, oc_index);
        discrete_solution.calculate_solution_at_infc_ncs_QPs(this->ncs_solution_v_at_QPs_.data(), infc_index, nc_index);
        //discrete_solution.calculate_solution_at_infc_ocs_QPs(this->ocs_solution_v_at_QPs_, infc_index, oc_index);
        //discrete_solution.calculate_solution_at_infc_ncs_QPs(this->ncs_solution_v_at_QPs_, infc_index, nc_index);

        const auto& normals = this->set_of_normals_[infc_index];
        const auto num_QPs = normals.size();

        auto& numerical_flux_QPs_m = this->set_of_numerical_flux_QPs_m_[infc_index];

        for (ushort q = 0; q < num_QPs; ++q)
        {            
            numerical_flux_QPs_m.change_column(q, this->numerical_flux_function_->calculate(this->ocs_solution_v_at_QPs_[q], this->ncs_solution_v_at_QPs_[q], normals[q]));
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