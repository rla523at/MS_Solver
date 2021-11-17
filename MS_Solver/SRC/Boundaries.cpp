#include "../INC/Boundaries.h"

Boundaries_DG::Boundaries_DG(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux)
{
    SET_TIME_POINT;

    this->num_boundaries_ = grid.num_boundaries();
    this->num_equations_ = discrete_solution.num_equations();

    this->oc_indexes_.resize(this->num_boundaries_);
    this->boundary_flux_functions_.reserve(this->num_boundaries_);
    this->set_of_oc_side_basis_QPs_m_.reserve(this->num_boundaries_);
    this->set_of_oc_side_QWs_basis_m_.reserve(this->num_boundaries_);
    this->set_of_normals_.reserve(this->num_boundaries_);

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
        const auto quadrature_rule = grid.boundary_quadrature_rule(bdry_index, integrand_degree);

        const auto& QPs = quadrature_rule.points;
        const auto& QWs = quadrature_rule.weights;
        const auto num_QPs = QPs.size();

        auto basis_QPs_m = discrete_solution.calculate_basis_points_m(oc_index, QPs);
        const auto num_basis = basis_QPs_m.num_row();

        Matrix QWs_basis_m(num_QPs, num_basis);
        for (int q = 0; q < num_QPs; ++q)
        {
            QWs_basis_m.change_row(q, discrete_solution.calculate_basis_point_v(oc_index, QPs[q]) * QWs[q]);
        }
        QWs_basis_m *= -1.0; //owner side

        this->set_of_oc_side_basis_QPs_m_.push_back(std::move(basis_QPs_m));
        this->set_of_oc_side_QWs_basis_m_.push_back(std::move(QWs_basis_m));

        //set normals
        this->set_of_normals_.push_back(grid.boundary_normals(bdry_index, oc_index, QPs));
    }

    LOG << std::left << std::setw(50) << "@ Boundaries DG precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n" << LOG.print_;
}

void Boundaries_DG::calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_soltuion) const
{
    for (int i = 0; i < this->num_boundaries_; ++i)
    {
        const auto oc_index = this->oc_indexes_[i];
        const auto solution_at_QPs = discrete_soltuion.calculate_solution_at_points(oc_index, this->set_of_oc_side_basis_QPs_m_[i]);
        const auto num_QPs = solution_at_QPs.size();

        auto& boundary_flux_function = *this->boundary_flux_functions_[i];
        const auto& normals = this->set_of_normals_[i];

        Matrix boundary_flux_quadrature(this->num_equations_, num_QPs);

        for (ushort q = 0; q < num_QPs; ++q)
        {
            boundary_flux_quadrature.change_column(q, boundary_flux_function.calculate(solution_at_QPs[q], normals[q]));
        }

        const auto delta_rhs = boundary_flux_quadrature * this->set_of_oc_side_QWs_basis_m_[i];
        residual.update_rhs(oc_index, delta_rhs);
    }
}