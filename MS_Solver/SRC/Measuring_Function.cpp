#include "../INC/Measuring_Function.h"

Scaled_Average_Difference_Measurer::Scaled_Average_Difference_Measurer(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    : criterion_solution_index_(criterion_solution_index)
    , num_infcs_(grid.num_inner_faces())
    , infc_index_to_oc_nc_index_pair_table_(grid.inner_face_index_to_oc_nc_index_pair_table())
{
    //for precalculation
    discrete_solution.precalculate_cell_P0_basis_values();
};

std::vector<double> Scaled_Average_Difference_Measurer::measure_infc_index_to_scaled_average_difference_table(const Discrete_Solution_DG& discrete_solution) const
{
    std::vector<double> infc_index_to_scaled_avg_diff_table(this->num_infcs_);
    for (uint infc_index = 0; infc_index < this->num_infcs_; ++infc_index)
    {
        const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];

        const auto oc_avg_value = discrete_solution.calculate_P0_nth_solution(oc_index, this->criterion_solution_index_);
        const auto nc_avg_value = discrete_solution.calculate_P0_nth_solution(nc_index, this->criterion_solution_index_);

        const auto avg_diff = std::abs(oc_avg_value - nc_avg_value);
        const auto min_value = (std::min)(oc_avg_value, nc_avg_value);
        const auto scaled_diff = avg_diff / min_value;

        infc_index_to_scaled_avg_diff_table[infc_index] = scaled_diff;
    }

    return infc_index_to_scaled_avg_diff_table;
}


Extrapolation_Differences_Measurer::Extrapolation_Differences_Measurer(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    : criterion_solution_index_(criterion_solution_index)
    , num_cells_(grid.num_cells())
    , cell_index_to_face_share_cell_indexes_table_(grid.cell_index_to_face_share_cell_indexes_table_consider_pbdry())
{
    this->cell_index_to_volume_reciprocal_table_.resize(this->num_cells_);

    for (uint i = 0; i < this->num_cells_; ++i)
    {
        const auto cell_volume = grid.cell_volume(i);
        this->cell_index_to_volume_reciprocal_table_[i] = 1.0 / cell_volume;
    }

    this->cell_index_to_QW_v_table_.resize(this->num_cells_);

    const auto cell_index_to_face_share_cell_indexes_table_ignore_pbdry = grid.cell_index_to_face_share_cell_indexes_table_ignore_pbdry();

    std::vector<std::vector<Euclidean_Vector>> cell_index_to_order_n_QPs_table(this->num_cells_);

    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto solution_degree = discrete_solution.solution_degree(cell_index);
        const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, solution_degree);

        cell_index_to_order_n_QPs_table[cell_index] = quadrature_rule.points;
        this->cell_index_to_QW_v_table_[cell_index] = quadrature_rule.weights;
    }

    //for precalculation
    const auto pbdry_pair_index_to_oc_nc_index_pair_table = grid.pbdry_pair_index_to_oc_nc_index_pair_table();
    const auto pbdry_pair_index_to_ocs_to_ncs_v_table = grid.pbdry_pair_index_to_ocs_to_ncs_v_table();
    const auto num_pbdry_pair = pbdry_pair_index_to_oc_nc_index_pair_table.size();

    std::vector<std::pair<std::vector<Euclidean_Vector>, std::vector<Euclidean_Vector>>> pbdry_pair_index_to_translated_oc_nc_order_n_QPs_pair_table(num_pbdry_pair);

    for (uint pbdry_pair_index = 0; pbdry_pair_index < num_pbdry_pair; ++pbdry_pair_index)
    {
        const auto [oc_index, nc_index] = pbdry_pair_index_to_oc_nc_index_pair_table[pbdry_pair_index];
        const auto& ocs_to_ncs_v = pbdry_pair_index_to_ocs_to_ncs_v_table.at(pbdry_pair_index);

        const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
        const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);

        const auto& oc_quadrature_rule = grid.get_cell_quadrature_rule(oc_index, oc_solution_degree);
        const auto& nc_quadrature_rule = grid.get_cell_quadrature_rule(nc_index, nc_solution_degree);

        auto translated_oc_QPs = oc_quadrature_rule.points;
        for (auto& QP : translated_oc_QPs)
        {
            QP += ocs_to_ncs_v;
        }

        auto translated_nc_QPs = nc_quadrature_rule.points;
        for (auto& QP : translated_nc_QPs)
        {
            QP -= ocs_to_ncs_v;
        }

        pbdry_pair_index_to_translated_oc_nc_order_n_QPs_pair_table[pbdry_pair_index] = { std::move(translated_oc_QPs),std::move(translated_nc_QPs) };
    }

    discrete_solution.precalculate_set_of_cell_index_to_target_cell_basis_QPs_m_(
        cell_index_to_order_n_QPs_table,
        cell_index_to_face_share_cell_indexes_table_ignore_pbdry,
        pbdry_pair_index_to_oc_nc_index_pair_table,
        pbdry_pair_index_to_translated_oc_nc_order_n_QPs_pair_table);
    //
}


std::vector<std::vector<double>> Extrapolation_Differences_Measurer::measure_cell_index_to_extrapolation_differences(const Discrete_Solution_DG& discrete_solution) const
{
    std::vector<std::vector<double>> cell_index_to_extrapolation_differences_table(this->num_cells_);

    static constexpr ushort max_num_QPs = 100;
    std::array<double, max_num_QPs> nth_extrapolated_solution_at_cell_QPs = { 0 };

    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        auto& extrapolation_differences_table = cell_index_to_extrapolation_differences_table[cell_index];

        const auto one_over_volume = this->cell_index_to_volume_reciprocal_table_[cell_index];
        const auto P0_value = discrete_solution.calculate_P0_nth_solution(cell_index, this->criterion_solution_index_);

        const auto& QWs_v = this->cell_index_to_QW_v_table_[cell_index];
        const auto num_QPs = static_cast<int>(QWs_v.size());

        const auto& face_share_cell_indexes = this->cell_index_to_face_share_cell_indexes_table_[cell_index];
        const auto num_face_share_cells = face_share_cell_indexes.size();

        extrapolation_differences_table.resize(num_face_share_cells);

        for (ushort j = 0; j < num_face_share_cells; ++j)
        {
            const auto face_share_cell_index = face_share_cell_indexes[j];

            discrete_solution.calculate_nth_extrapolated_solution_at_cell_QPs(nth_extrapolated_solution_at_cell_QPs.data(), face_share_cell_index, cell_index, this->criterion_solution_index_);

            const auto P0_value_by_extrapolate = ms::BLAS::x_dot_y(num_QPs, nth_extrapolated_solution_at_cell_QPs.data(), QWs_v.data());

            extrapolation_differences_table[j] = std::abs(P0_value_by_extrapolate * one_over_volume - P0_value);
        }
    }

    return cell_index_to_extrapolation_differences_table;
}

Divergence_Velocity_Measurer::Divergence_Velocity_Measurer(const Grid& grid, Discrete_Solution_DG& discrete_solution)
{
    this->num_cells_ = grid.num_cells();
    this->cell_index_to_num_QPs_.resize(this->num_cells_);

    //construction optimization
    solution_at_cell_QPs.fill(Euclidean_Vector(discrete_solution.num_solutions()));
    ddx_GE_solution_at_cell_QPs.fill(Euclidean_Vector(discrete_solution.num_solutions()));
    ddy_GE_solution_at_cell_QPs.fill(Euclidean_Vector(discrete_solution.num_solutions()));

    //for precalculation
    std::vector<std::vector<Euclidean_Vector>> cell_index_to_QPs(this->num_cells_);

    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto solution_degree = discrete_solution.solution_degree(cell_index);
        const auto volume_integral_degree = 2 * solution_degree;

        const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, volume_integral_degree);
        cell_index_to_QPs[cell_index] = quadrature_rule.points;
        this->cell_index_to_num_QPs_[cell_index] = static_cast<ushort>(quadrature_rule.points.size());
    }    

    discrete_solution.precalcualte_cell_QPs_ddx_basis_values(cell_index_to_QPs);
    discrete_solution.precalcualte_cell_QPs_ddy_basis_values(cell_index_to_QPs);
}

std::vector<std::vector<double>> Divergence_Velocity_Measurer::measure_cell_index_to_divergence_velocities_table(const Discrete_Solution_DG& discrete_solution) const
{
    std::vector<std::vector<double>> cell_index_to_divergence_velocities_table(this->num_cells_);

    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        discrete_solution.calculate_solution_at_cell_QPs(solution_at_cell_QPs.data(), cell_index);
        discrete_solution.calculate_ddx_GE_solution_at_cell_QPs(ddx_GE_solution_at_cell_QPs.data(), cell_index);
        discrete_solution.calculate_ddy_GE_solution_at_cell_QPs(ddy_GE_solution_at_cell_QPs.data(), cell_index);
        
        const auto num_QPs = this->cell_index_to_num_QPs_[cell_index];        

        auto& divergence_velocities = cell_index_to_divergence_velocities_table[cell_index];
        divergence_velocities.resize(num_QPs);

        for (ushort j = 0; j < num_QPs; ++j)
        {
            const auto rho = solution_at_cell_QPs[j][0];
            const auto u = solution_at_cell_QPs[j][4];
            const auto v = solution_at_cell_QPs[j][5];

            const auto ddx_rho = ddx_GE_solution_at_cell_QPs[j][0];
            const auto ddx_rhou = ddx_GE_solution_at_cell_QPs[j][1];

            const auto ddy_rho = ddy_GE_solution_at_cell_QPs[j][0];
            const auto ddy_rhov = ddy_GE_solution_at_cell_QPs[j][2];

            const auto one_over_rho = 1.0 / rho;

            const auto ddx_u = one_over_rho * (ddx_rhou - ddx_rho * u);
            const auto ddy_v = one_over_rho * (ddy_rhov - ddy_rho * v);

            const auto div_velocity = ddx_u + ddy_v;

            divergence_velocities[j] = div_velocity;
        }
    }

    return cell_index_to_divergence_velocities_table;
}

Average_Solution_Jump_Measurer::Average_Solution_Jump_Measurer(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    :criterion_solution_index_(criterion_solution_index)
{
    const auto space_dimension = grid.space_dimension();
    this->num_infcs_ = grid.num_inner_faces();

    this->infc_index_to_reciprocal_volume_table_.resize(this->num_infcs_);
    this->infc_index_to_oc_nc_index_pair_table_.resize(this->num_infcs_);
    this->infc_index_to_jump_QWs_v_table_.resize(this->num_infcs_);

    std::vector<uint> oc_indexes(this->num_infcs_);
    std::vector<uint> nc_indexes(this->num_infcs_);
    std::vector<std::vector<Euclidean_Vector>> infc_index_to_ocs_QPs_table_(this->num_infcs_);
    std::vector<std::vector<Euclidean_Vector>> infc_index_to_ncs_QPs_table_(this->num_infcs_);

    for (uint i = 0; i < this->num_infcs_; ++i)
    {
        this->infc_index_to_reciprocal_volume_table_[i] = 1.0 / grid.inner_face_volume(i);
        this->infc_index_to_oc_nc_index_pair_table_[i] = grid.inner_face_oc_nc_index_pair(i);

        const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[i];
        const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
        const auto nc_solution_degree = discrete_solution.solution_degree(oc_index);
        const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

        const auto& [ocs_quadrature_rule, ncs_quadrature_rule] = grid.inner_face_quadrature_rule(i, max_solution_degree);
        this->infc_index_to_jump_QWs_v_table_[i] = ocs_quadrature_rule.weights;

        //for precalculation
        oc_indexes[i] = oc_index;
        nc_indexes[i] = nc_index;
        infc_index_to_ocs_QPs_table_[i] = ocs_quadrature_rule.points;
        infc_index_to_ncs_QPs_table_[i] = ncs_quadrature_rule.points;
    }

    //precalculation
    discrete_solution.precalculate_infs_ocs_jump_QPs_basis_values(oc_indexes, infc_index_to_ocs_QPs_table_);
    discrete_solution.precalculate_infs_ncs_jump_QPs_basis_values(nc_indexes, infc_index_to_ncs_QPs_table_);
}

std::vector<double> Average_Solution_Jump_Measurer::measure_infc_index_to_scaled_average_solution_jump_table(const Discrete_Solution_DG& discrete_solution)
{
    std::vector<double> infc_index_to_average_solution_jump_table(this->num_infcs_);

    for (uint infc_index = 0; infc_index < this->num_infcs_; ++infc_index)
    {
        const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];
        discrete_solution.calculate_nth_solution_at_infc_ocs_jump_QPs(this->value_at_ocs_jump_QPs_.data(), infc_index, oc_index, this->criterion_solution_index_);
        discrete_solution.calculate_nth_solution_at_infc_ncs_jump_QPs(this->value_at_ncs_jump_QPs_.data(), infc_index, nc_index, this->criterion_solution_index_);

        auto& jump_QWs_v = this->infc_index_to_jump_QWs_v_table_[infc_index];        
        const auto num_QPs = static_cast<int>(jump_QWs_v.size());

        ms::BLAS::x_minus_y(num_QPs, this->value_at_ocs_jump_QPs_.data(), this->value_at_ncs_jump_QPs_.data(), this->value_diff_at_jump_QPs_.data());
        ms::BLAS::abs_x(num_QPs, this->value_diff_at_jump_QPs_.data());
        const auto jump = ms::BLAS::x_dot_y(num_QPs, this->value_diff_at_jump_QPs_.data(), jump_QWs_v.data());

        const auto one_over_volume = this->infc_index_to_reciprocal_volume_table_[infc_index];        
        const auto avg_sol_jump = jump * one_over_volume;
        const auto scaled_avg_sol_jump = this->calculate_scail_factor(discrete_solution, infc_index) * avg_sol_jump;

        infc_index_to_average_solution_jump_table[infc_index] = scaled_avg_sol_jump;
    }

    return infc_index_to_average_solution_jump_table;
}

//
//std::vector<double> Difference_Of_Extrapolatiion_Difference::measure_infc_index_to_scaled_average_difference_table(const Discrete_Solution_DG& discrete_solution) const
//{
//    const auto cell_index_to_extrapolate_diff = this->extrapolate_difference_.measure_cell_index_to_extrapolation_differences(discrete_solution);
//
//    std::vector<double> cell_index_to_avg_diff_of_extrapolate_diff(this->num_cells_);
//
//    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
//    {
//        double diff_sum = 0.0;
//
//        const auto& face_share_cell_indexes = this->cell_index_to_face_share_cell_indexes_table_[cell_index];
//        const auto num_face_share_cells = face_share_cell_indexes.size();
//
//        const auto extrapolate_diff = cell_index_to_extrapolate_diff[cell_index];
//
//        for (ushort j = 0; j < num_face_share_cells; ++j)
//        {
//            const auto face_share_cell_index = face_share_cell_indexes[j];
//            const auto face_neighbor_extrapolate_diff = cell_index_to_extrapolate_diff[face_share_cell_index];
//
//            const auto diff = std::abs(extrapolate_diff - face_neighbor_extrapolate_diff);
//            diff_sum += diff;
//        }
//
//        cell_index_to_avg_diff_of_extrapolate_diff[cell_index] = diff_sum / num_face_share_cells;
//    }
//
//    return cell_index_to_avg_diff_of_extrapolate_diff;
//}