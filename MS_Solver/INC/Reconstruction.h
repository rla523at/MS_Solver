#pragma once
#include "Limiter.h"

class Reconstruction_DG
{
public://Query
	virtual void reconstruct(Discrete_Solution_DG& discrete_solution) const abstract;
};

class No_Reconstruction_DG : public Reconstruction_DG
{
public://Query
	void reconstruct(Discrete_Solution_DG& discrete_solution) const override {};
};

class Hierarchical_Limiting_DG : public Reconstruction_DG
{
public:
	Hierarchical_Limiting_DG(std::unique_ptr<MLP_Criterion_Base>&& stability_criterion, std::unique_ptr<Indicator>&& indicator, std::unique_ptr<Limiter>&& limiter)
		: stability_criterion_(std::move(stability_criterion))
		, indicator_(std::move(indicator))
		, limiter_(std::move(limiter)) {};

public://Query
    void reconstruct(Discrete_Solution_DG& discrete_solution) const override
    {
		const auto num_cells = discrete_solution.num_cells();

		this->stability_criterion_->precaclulate(discrete_solution);
		this->indicator_->precalculate(discrete_solution);
		this->limiter_->precalculate(discrete_solution);

		for (uint cell_index = 0; cell_index < num_cells; ++cell_index)
		{
			auto projection_degree = discrete_solution.solution_degree(cell_index);

			while (true)
			{
				const auto cell_type = this->indicator_->indicate(discrete_solution, cell_index, *this->stability_criterion_);
				this->limiter_->limit(cell_index, cell_type, discrete_solution, projection_degree, *this->stability_criterion_);

				if (this->limiter_->is_end())
				{
					break;
				}
			}
		}
    }


private:
	std::unique_ptr<MLP_Criterion_Base> stability_criterion_;
    std::unique_ptr<Indicator> indicator_;
    std::unique_ptr<Limiter> limiter_;
};


//
//class hMLP_Reconstruction_DG : public Reconstruction_DG
//{
//public:
//    hMLP_Reconstruction_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution);//멤버 변수의 구체클레스가 결정되어 있기 때문에 의존성 주입을 할 필요가 없다
//    //hMLP_Reconstruction_DG(MLP_Criterion&& stability_criterion, MLP_Indicator&& MLP_indicator, MLP_u1_Limiter&& MLP_u1_limiter, const Grid& grid);
//    
//public:
//    void reconstruct(Discrete_Solution_DG& discrete_solution) override;
//
//private:
//    void precalculate(const Discrete_Solution_DG& discrete_solution);
//
//private:
//    uint num_cells_ = 0;
//	
//    MLP_Criterion stability_criterion_;
//    MLP_Indicator MLP_indicator_;
//    MLP_u1_Limiter MLP_u1_limiter_;
//
//    //test
//    Discontinuity_Indicator discontinuity_indicator_;    
//    //test
//
//    //precalculate
//    std::vector<std::vector<double>> set_of_P1_projected_criterion_values_at_verticies_;
//};
//
//class hMLP_BD_Reconstruction_DG : public Reconstruction_DG
//{
//public:
//    hMLP_BD_Reconstruction_DG(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::string& governing_equation_name);
//
//    hMLP_BD_Reconstruction_DG(Simplex_Decomposed_MLP_Criterion&& stability_criterion, MLP_Indicator&& MLP_indicator, MLP_u1_Limiter&& MLP_u1_limiter,
//        Subcell_Oscillation_Indicator&& boundary_indicator, Shock_Indicator&& shock_indicator,  const Grid& grid);
//
//public:
//    void reconstruct(Discrete_Solution_DG& discrete_solution) override;
//
//private:
//    void precalculate(const Discrete_Solution_DG& discrete_solution);
//    void limit_solution(Discrete_Solution_DG& discrete_solution, const uint cell_index, ushort& solution_degree) const;
//
//private:
//    uint num_cells_ = 0;
//    mutable bool is_end_limting_ = false;
//
//    Simplex_Decomposed_MLP_Criterion stability_criterion_;
//    MLP_Indicator MLP_indicator_;
//    MLP_u1_Limiter MLP_u1_limiter_;
//    Subcell_Oscillation_Indicator subcell_oscillation_indicator;
//    Shock_Indicator shock_indicator_;
//
//    //test
//    Discontinuity_Indicator discontinuity_indicator_;
//    //
//
//    //precalculate
//    std::vector<std::vector<double>> set_of_P1_projected_criterion_values_at_verticies_;
//};
//
//
//#include "../INC/Post_Processor.h"
//class Test_Reconstuction_DG : public Reconstruction_DG
//{
//public:
//    Test_Reconstuction_DG(Discontinuity_Indicator&& discontinuity_indicator)
//        : discontinuity_indicator_(discontinuity_indicator) {};
//
//public:
//    void reconstruct(Discrete_Solution_DG& discrete_solution) override;
//
//private:
//    Discontinuity_Indicator discontinuity_indicator_;
//};
//
