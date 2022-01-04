#include "../INC/Indicating_Function_Impl.h"

void Heuristic_Discontinuity_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution)
{
    std::fill(this->cell_index_to_has_discontinuity_table_.begin(), this->cell_index_to_has_discontinuity_table_.end(), false);

    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto my_value = discrete_solution.calculate_P0_nth_solution(cell_index, this->criterion_solution_index_);

        const auto& face_share_cell_indexes = this->cell_index_to_face_share_cell_indexes_table_[cell_index];

        for (const auto face_share_cell_index : face_share_cell_indexes)
        {
            const auto other_value = discrete_solution.calculate_P0_nth_solution(face_share_cell_index, this->criterion_solution_index_);

            if (std::abs(my_value - other_value) >= 0.1 * (std::min)(my_value, other_value))
            {
                this->cell_index_to_has_discontinuity_table_[cell_index] = true;
                break;
            }
        }
    }
};

Extrapolation_Discontinuity_Indicator::Extrapolation_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    : Discontinuity_Indicator(grid, criterion_solution_index)
{
    this->cell_index_to_has_discontinuity_table_.resize(this->num_cells_, false);

    this->cell_index_to_volume_reciprocal_table_.resize(this->num_cells_);
    this->cell_index_to_threshold_value_table_.resize(this->num_cells_);
    const auto space_dimension = grid.space_dimension();

    for (uint i = 0; i < this->num_cells_; ++i)
    {
        const auto cell_volume = grid.cell_volume(i);

        this->cell_index_to_threshold_value_table_[i] = std::pow(cell_volume, 1.0 / space_dimension);
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
};

void Extrapolation_Discontinuity_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution) 
{
    static constexpr ushort max_num_QPs = 100;
    std::array<double, max_num_QPs> nth_extrapolated_solution_at_cell_QPs = { 0 };

    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        double discontinuity_indicator_value = 0.0;

        const auto one_over_volume = this->cell_index_to_volume_reciprocal_table_[cell_index];
        const auto P0_value = discrete_solution.calculate_P0_nth_solution(cell_index, this->criterion_solution_index_);

        const auto& QWs_v = this->cell_index_to_QW_v_table_[cell_index];
        const auto num_QPs = static_cast<int>(QWs_v.size());

        const auto& face_share_cell_indexes = this->cell_index_to_face_share_cell_indexes_table_[cell_index];
        const auto num_face_share_cells = face_share_cell_indexes.size();

        for (ushort j = 0; j < num_face_share_cells; ++j)
        {
            const auto face_share_cell_index = face_share_cell_indexes[j];

            discrete_solution.calculate_nth_extrapolated_solution_at_cell_QPs(nth_extrapolated_solution_at_cell_QPs.data(), face_share_cell_index, cell_index, this->criterion_solution_index_);

            const auto P0_value_by_extrapolate = ms::BLAS::x_dot_y(num_QPs, nth_extrapolated_solution_at_cell_QPs.data(), QWs_v.data());

            discontinuity_indicator_value += std::abs(P0_value_by_extrapolate * one_over_volume - P0_value);
        }

        discontinuity_indicator_value /= num_face_share_cells;

        const auto threshold_value = this->cell_index_to_threshold_value_table_[cell_index];
        if (threshold_value < discontinuity_indicator_value)
        {
            this->cell_index_to_has_discontinuity_table_[cell_index] = true;
        }

        //this->discontinuity_factor_[cell_index] = std::log(discontinuity_indicator_value) / std::log(characteristic_length);
    }
};