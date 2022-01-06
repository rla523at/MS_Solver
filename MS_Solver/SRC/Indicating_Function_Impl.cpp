#include "../INC/Indicating_Function_Impl.h"

Discontinuity_Indicating_Function_Base::Discontinuity_Indicating_Function_Base(const std::vector<double>& cell_index_to_threshold_value_table, std::unique_ptr<Discontinuity_Measuring_Function>&& measuring_function)
    : cell_index_to_threshold_value_table_(cell_index_to_threshold_value_table)
    , measuring_function_(std::move(measuring_function))
{
    this->num_cells_ = static_cast<uint>(cell_index_to_threshold_value_table_.size());
    this->cell_index_to_has_discontinuity_table_.resize(this->num_cells_, false);
}

void Discontinuity_Indicating_Function_Base::precalculate(const Discrete_Solution_DG& discrete_solution)
{    
    std::fill(this->cell_index_to_has_discontinuity_table_.begin(), this->cell_index_to_has_discontinuity_table_.end(), false);

    const auto cell_index_to_measuring_value_table = measuring_function_->measure(discrete_solution);

    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto threshold_value = cell_index_to_threshold_value_table_[cell_index];
        const auto measuring_value = cell_index_to_measuring_value_table[cell_index];

        if (this->compare(measuring_value, threshold_value))
        {
            this->cell_index_to_has_discontinuity_table_[cell_index] = true;
        }
    }
}