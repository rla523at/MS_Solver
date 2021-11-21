#include "../INC/Boundaries.h"

Boundaries_DG::Boundaries_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux)
{
    SET_TIME_POINT;

    this->num_boundaries_ = grid.num_boundaries();
    this->num_equations_ = discrete_solution.num_equations();

    this->oc_indexes_.resize(this->num_boundaries_);
    this->boundary_flux_functions_.reserve(this->num_boundaries_);
    this->set_of_oc_side_QWs_basis_m_.reserve(this->num_boundaries_);
    this->set_of_normals_.reserve(this->num_boundaries_);

    std::vector<Quadrature_Rule> quadrature_rules(this->num_boundaries_);

    for (int bdry_index = 0; bdry_index < this->num_boundaries_; ++bdry_index)
    {
        //set oc index
        const auto oc_index = grid.boundary_owner_cell_index(bdry_index);
        this->oc_indexes_[bdry_index] = oc_index;

        //set boundary flux function
        const auto boundary_type = grid.boundary_type(bdry_index);
        this->boundary_flux_functions_.push_back(Boundary_Flux_Function_Factory::make_unique(boundary_type, numerical_flux));

        //set oc_side_basis_QPs_m_ & oc_side_QWs_basis_m_
        const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
        const auto integrand_degree = 2 * oc_solution_degree + 1;
        auto quadrature_rule = grid.boundary_quadrature_rule(bdry_index, integrand_degree);

        const auto& QPs = quadrature_rule.points;
        const auto& QWs = quadrature_rule.weights;
        const auto num_QPs = QPs.size();

        const auto num_basis = discrete_solution.num_basis(oc_index);

        Matrix QWs_basis_m(num_QPs, num_basis);
        for (int q = 0; q < num_QPs; ++q)
        {
            QWs_basis_m.change_row(q, discrete_solution.calculate_basis_point_v(oc_index, QPs[q]) * QWs[q]);
        }
        QWs_basis_m *= -1.0; //owner side

        this->set_of_oc_side_QWs_basis_m_.push_back(std::move(QWs_basis_m));

        //set normals
        this->set_of_normals_.push_back(grid.boundary_normals(bdry_index, oc_index, QPs));

        quadrature_rules[bdry_index] = std::move(quadrature_rule);
    }

    discrete_solution.precalculate_basis_bdry_QPs(this->oc_indexes_, quadrature_rules);

    LOG << std::left << std::setw(50) << "@ Boundaries DG precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n" << LOG.print_;
}

void Boundaries_DG::calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_soltuion) const
{
    for (int bdry_index = 0; bdry_index < this->num_boundaries_; ++bdry_index)
    {
        const auto oc_index = this->oc_indexes_[bdry_index];
        const auto solution_at_QPs = discrete_soltuion.calculate_solution_at_bdry_QPs(bdry_index, oc_index);
        const auto num_QPs = solution_at_QPs.size();

        auto& boundary_flux_function = *this->boundary_flux_functions_[bdry_index];
        const auto& normals = this->set_of_normals_[bdry_index];

        Matrix boundary_flux_quadrature(this->num_equations_, num_QPs);

        for (ushort q = 0; q < num_QPs; ++q)
        {
            boundary_flux_quadrature.change_column(q, boundary_flux_function.calculate(solution_at_QPs[q], normals[q]));
        }

        const auto delta_rhs = boundary_flux_quadrature * this->set_of_oc_side_QWs_basis_m_[bdry_index];
        residual.update_rhs(oc_index, delta_rhs);
    }
}