#pragma once
#include "Discrete_Solution.h"
#include "Indicator.h"
#include "Limiter_Function.h"

class Reconstruction_DG
{
public://Query
	virtual void reconstruct(Discrete_Solution_DG& discrete_solution) abstract;
};

class No_Reconstruction_DG : public Reconstruction_DG
{
public://Query
	void reconstruct(Discrete_Solution_DG& discrete_solution) override {};
};

class hMLP_Reconstruction_DG : public Reconstruction_DG
{
public:
    hMLP_Reconstruction_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
    void reconstruct(Discrete_Solution_DG& discrete_solution) override;

private:
    bool is_trouble(const Discrete_Solution_DG& discrete_solution, const uint cell_index) const;
    void precalculate(const Discrete_Solution_DG& discrete_solution);

private:
	ushort criterion_variable_index_ = 0;	
	uint num_cells_=0;
	
	const std::unordered_map<uint, std::set<uint>>& vnode_index_to_share_cell_index_set_;
    std::vector<std::vector<uint>> set_of_vnode_indexes_;
    std::vector<double> volumes_;

    //precalculate
    std::vector<double> P0_criterion_values_;
    std::vector<std::vector<double>> set_of_P1_projected_criterion_values_at_verticies_;
    std::vector<std::vector<std::pair<double, double>>> set_of_allowable_min_max_criterion_values_;

    //construction optimization
    std::unordered_map<uint, std::pair<double, double>> vnode_index_to_allowable_min_max_criterion_value_;
};

class Reconstruction_DG_Factory//static class
{
public://Query
    static std::unique_ptr<Reconstruction_DG> make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        const auto& reconstruction_scheme = configuration.get_reconstruction_scheme();

        if (ms::contains_icase(reconstruction_scheme, "none"))
        {
            return std::make_unique<No_Reconstruction_DG>();
        }
        else if (ms::contains_icase(reconstruction_scheme, "hMLP"))
        {
            return std::make_unique<hMLP_Reconstruction_DG>(grid, discrete_solution);
        }
        else
        {
            EXCEPTION("reconstruction shceme in configuration file is not supported");
            return nullptr;
        }
    }

private:
    Reconstruction_DG_Factory(void) = delete;
};

//