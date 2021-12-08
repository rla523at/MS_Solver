#include "../INC/Boundaries.h"

Boundaries_DG::Boundaries_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux)
{
    Profiler::set_time_point();

    this->num_boundaries_ = static_cast<uint>(grid.num_boundaries());

    this->oc_indexes_.resize(this->num_boundaries_);
    this->boundary_flux_functions_.resize(this->num_boundaries_);
    this->set_of_oc_side_QWs_basis_m_.resize(this->num_boundaries_);
    this->set_of_normals_.resize(this->num_boundaries_);

    //for construct optimization
    const auto num_equations = discrete_solution.num_equations();
    const auto num_solutions = discrete_solution.num_solutions();

    this->set_of_boundary_flux_QPs_m_.resize(this->num_boundaries_);
    this->solution_v_at_QPs_.resize(this->max_num_QPs, Euclidean_Vector(num_solutions));
    //

    //for precalculation
    std::vector<Quadrature_Rule> quadrature_rules(this->num_boundaries_);
    //

    for (uint bdry_index = 0; bdry_index < this->num_boundaries_; ++bdry_index)
    {
        const auto oc_index = grid.boundary_owner_cell_index(bdry_index);
        this->oc_indexes_[bdry_index] = oc_index;

        const auto boundary_type = grid.boundary_type(bdry_index);
        this->boundary_flux_functions_[bdry_index] = std::move(Boundary_Flux_Function_Factory::make_unique(boundary_type, numerical_flux));

        const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
        const auto integrand_degree = 2 * oc_solution_degree + 1;
        auto quadrature_rule = grid.boundary_quadrature_rule(bdry_index, integrand_degree);

        const auto& QPs = quadrature_rule.points;
        const auto& QWs = quadrature_rule.weights;
        const auto num_QPs = QPs.size();

        const auto num_basis = discrete_solution.num_basis(oc_index);

        Matrix QWs_basis_m(num_QPs, num_basis);
        for (ushort q = 0; q < num_QPs; ++q)
        {
            QWs_basis_m.change_row(q, discrete_solution.calculate_basis_point_v(oc_index, QPs[q]) * QWs[q]);
        }
        QWs_basis_m *= -1.0; //owner side

        this->set_of_oc_side_QWs_basis_m_[bdry_index] = std::move(QWs_basis_m);
        this->set_of_normals_[bdry_index] = std::move(grid.boundary_normals(bdry_index, oc_index, QPs));

        //for construct optimization
        REQUIRE(num_QPs <= this->max_num_QPs, "number of quadrature points should be less then expected maximum");
        this->set_of_boundary_flux_QPs_m_[bdry_index] = Matrix(num_equations, num_QPs);
        //

        //for precalculation
        quadrature_rules[bdry_index] = std::move(quadrature_rule);
        //
    }

    //for precalculation
    discrete_solution.precalculate_basis_bdry_QPs_basis_values(this->oc_indexes_, quadrature_rules);
    //

    LOG << std::left << std::setw(50) << "@ Boundaries DG precalculation" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << LOG.print_;
}

void Boundaries_DG::calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_solution) const
{
    for (uint bdry_index = 0; bdry_index < this->num_boundaries_; ++bdry_index)
    {
        const auto oc_index = this->oc_indexes_[bdry_index];
        discrete_solution.calculate_solution_at_bdry_QPs(this->solution_v_at_QPs_, bdry_index, oc_index);
        
        auto& boundary_flux_function = *this->boundary_flux_functions_[bdry_index];
        const auto& normals = this->set_of_normals_[bdry_index];
        const auto num_QPs = normals.size();

        auto& bdry_flux_QPs_m = this->set_of_boundary_flux_QPs_m_[bdry_index];
        for (ushort q = 0; q < num_QPs; ++q)
        {
            bdry_flux_QPs_m.change_column(q, boundary_flux_function.calculate(this->solution_v_at_QPs_[q], normals[q]));
        }

        const auto num_values = discrete_solution.num_values(oc_index);
        this->reset_residual_values(num_values);

        ms::gemm(bdry_flux_QPs_m, this->set_of_oc_side_QWs_basis_m_[bdry_index], this->residual_values_.data());
        residual.update_rhs(oc_index, this->residual_values_.data());
    }
}

void Boundaries_DG::reset_residual_values(const uint num_values) const
{
    std::fill(this->residual_values_.begin(), this->residual_values_.begin() + num_values, 0.0);
}