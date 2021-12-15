#pragma once
#include "Discrete_Solution.h"

#include <unordered_map>

class MLP_u1_Limiting_Strategy
{
private:
    MLP_u1_Limiting_Strategy(void) = delete;

public:
    static double limiter_function(const std::vector<double>& P1_projected_solutions, const double P0_criterion_value, const std::vector<std::pair<double,double>>& set_of_allowable_min_max)
    {
        const auto num_vertices = P1_projected_solutions.size();

        double limiting_value = 1.0;
        for (ushort j = 0; j < num_vertices; ++j)
        {
            const auto [allowable_min, allowable_max] = set_of_allowable_min_max[j];

            const auto P1_mode_criterion_value = P1_projected_solutions[j] - P0_criterion_value;

            limiting_value = (std::min)(limiting_value, MLP_u1_Limiting_Strategy::limiter_function(P1_mode_criterion_value, P0_criterion_value, allowable_min, allowable_max));
        }

        return limiting_value;
    };

    static double limiter_function(const double P1_mode_solution, const double P0_solution, const double allowable_min, const double allowable_max) {
        if (P1_mode_solution == 0)
        {
            return 1.0;
        }

        if (P1_mode_solution < 0)
        {
            return (std::min)((allowable_min - P0_solution) / P1_mode_solution, 1.0);
        }
        else
        {
            return (std::min)((allowable_max - P0_solution) / P1_mode_solution, 1.0);
        }
    };
};

class P1_Projected_MLP_Condition
{
private:
    P1_Projected_MLP_Condition(void) = delete;

public:
    static bool is_satisfy(const double P1_projected_value, const double allowable_min, const double allowable_max) 
    {
        return allowable_min <= P1_projected_value && P1_projected_value <= allowable_max;
    }
};


class MLP_Smooth_Extrema_Detector
{
private:
    MLP_Smooth_Extrema_Detector(void) = delete;

public:
    static bool is_smooth_extrema(const double solution, const double higher_mode_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) 
    {
        if (P1_mode_solution > 0 && higher_mode_solution < 0 && solution > allowable_min)
            return true;
        else if (P1_mode_solution < 0 && higher_mode_solution > 0 && solution < allowable_max)
            return true;
        else
            return false;
    }
};


class Constant_Region_Detector
{
private:
    Constant_Region_Detector(void) = delete;

public:
    static bool is_constant(const double solution, const double P0_solution, const double volume) 
    {
        const auto constant_criterion = (std::max)(1.0E-3 * std::abs(P0_solution), volume);
        return std::abs(solution - P0_solution) <= constant_criterion;
    }
};

class HOM_Reconstruction
{
public://Query
	virtual void reconstruct(Discrete_Solution_DG& discrete_solution) const abstract;
};

class HOM_No_Reconstruction : public HOM_Reconstruction
{
public://Query
	void reconstruct(Discrete_Solution_DG& discrete_solution) const {};
};

class HOM_hMLP_Reconstruction : public HOM_Reconstruction
{
public:
	void reconstruct(Discrete_Solution_DG& discrete_solution)
	{
        this->precalculate(discrete_solution);

        const auto num_cell = this->P0_criterion_values_.size();
        for (uint cell_index = 0; cell_index < num_cell; ++cell_index) 
        {
            auto solution_degree = discrete_solution.solution_degree(cell_index);

            const auto P0_criterion_value = this->P0_criterion_values_[cell_index];
            const auto& allowable_min_max_criterion_values = this->set_of_allowable_min_max_criterion_values_[cell_index];

            while (true) 
            {
                if (this->is_trouble(discrete_solution, cell_index))
                {               
                    if (solution_degree <= 2)
                    {                        
                        const auto P1_projected_solution_at_vertices = discrete_solution.calculate_P1_projected_nth_solution_at_vertices(cell_index);
                        const auto num_vertices = P1_projected_solution_at_vertices.size();

                        for (ushort i = 0; i < num_vertices; ++i)
                        {
                            this->P1_projected_criterion_values_[i] = P1_projected_solution_at_vertices[i][this->criterion_variable_index_];
                        }

                        const auto limiting_value = MLP_u1_Limiting_Strategy::limiter_function(this->P1_projected_criterion_values_, P0_criterion_value, allowable_min_max_criterion_values);

                        solution_degree = 1;
                        discrete_solution.project_to_Pn_space(cell_index, solution_degree);
                        discrete_solution.limiting_slope(cell_index, limiting_value);
                        break;
                    }
                    else
                    {
                        discrete_solution.project_to_Pn_space(cell_index, --solution_degree);
                    }
                }
                else
                {
                    break;
                }
            }
        }
	};

private:
	void precalculate(const Discrete_Solution_DG& discrete_solution)
	{
        const auto P0_solutions = discrete_solution.calculate_P0_solutions();
        const auto num_solutions = P0_solutions.size();

        for (uint i = 0; i < num_solutions; ++i)
        {
            this->P0_criterion_values_[i] = P0_solutions[i][this->criterion_variable_index_];
        }

		//const auto num_vnode = this->vnode_index_to_share_cell_index_set_.size();

		//this->vnode_index_to_allowable_min_max_criterion_value_.reserve(num_vnode);

		for (const auto& [vnode_index, share_cell_index_set] : this->vnode_index_to_share_cell_index_set_) 
		{
			const auto num_share_cell = share_cell_index_set.size();
			std::vector<double> criterion_variables;
			criterion_variables.reserve(num_share_cell);

			for (const auto cell_index : share_cell_index_set)
			{
				criterion_variables.push_back(this->P0_criterion_values_[cell_index]);
			}

			const auto min_criterion_value = *std::min_element(criterion_variables.begin(), criterion_variables.end());
			const auto max_criterion_value = *std::max_element(criterion_variables.begin(), criterion_variables.end());

            this->vnode_index_to_allowable_min_max_criterion_value_[vnode_index] = std::make_pair(min_criterion_value, max_criterion_value);
		}

        for (uint i = 0; i < this->num_cell_; ++i)
        {
            const auto& vnode_indexes = this->set_of_vnode_indexes_[i];
            const auto num_vnode_indexes = vnode_indexes.size();

            for (ushort j = 0; j < num_vnode_indexes; ++j)
            {
                const auto vnode_index = vnode_indexes[j];
                this->set_of_allowable_min_max_criterion_values_[i][j] = this->vnode_index_to_allowable_min_max_criterion_value_[vnode_index];
            }
        }
	}

    bool is_trouble(const Discrete_Solution_DG& discrete_solution, const uint cell_index) const
    {
        const auto& vnode_indexes = this->set_of_vnode_indexes_[cell_index];
        const auto num_vnodes = vnode_indexes.size();

        const auto solution_at_vertices = discrete_solution.calculate_solution_at_vertices(cell_index);
        const auto P1_projected_solution_at_vertices = discrete_solution.calculate_P1_projected_nth_solution_at_vertices(cell_index);
        const auto P0_criterion_value = this->P0_criterion_values_[cell_index];
        const auto volume = this->volumes_[cell_index];

        for (ushort j = 0; j < num_vnodes; ++j)
        {
            const auto vnode_index = vnode_indexes[j];
            const auto [allowable_min, allowable_max] = this->vnode_index_to_allowable_min_max_criterion_value_.at(vnode_index);

            const auto& criterion_value = solution_at_vertices[j][this->criterion_variable_index_];
            const auto& P1_projected_criterion_value = P1_projected_solution_at_vertices[j][this->criterion_variable_index_];

            const auto higher_mode_criterion_value = criterion_value - P1_projected_criterion_value;
            const auto P1_mode_criterion_value = P1_projected_criterion_value - P0_criterion_value;

            if (!Constant_Region_Detector::is_constant(criterion_value, P0_criterion_value, volume) &&
                !P1_Projected_MLP_Condition::is_satisfy(P1_projected_criterion_value, allowable_min, allowable_max) &&
                !MLP_Smooth_Extrema_Detector::is_smooth_extrema(criterion_value, higher_mode_criterion_value, P1_mode_criterion_value, allowable_min, allowable_max))
            {
                return true;
            }
        }

        return false;
    }

    void calculate_P1_projected_criterion_values(const Discrete_Solution_DG& discrete_solution, const uint cell_index)
    {

    }

private:
	ushort criterion_variable_index_ = 0;	
	uint num_cell_=0;
	
	const std::unordered_map<uint, std::set<uint>>& vnode_index_to_share_cell_index_set_;
    std::unordered_map<uint, std::pair<double, double>> vnode_index_to_allowable_min_max_criterion_value_;
    std::vector<std::vector<std::pair<double, double>>> set_of_allowable_min_max_criterion_values_;

    std::vector<double> P0_criterion_values_;
    std::vector<double> P1_projected_criterion_values_;

	std::vector<std::vector<uint>> set_of_vnode_indexes_;
    std::vector<double> volumes_;
};