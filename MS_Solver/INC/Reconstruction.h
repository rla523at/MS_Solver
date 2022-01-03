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
	void reconstruct(Discrete_Solution_DG& discrete_solution) const override;
private:
	std::unique_ptr<MLP_Criterion_Base> stability_criterion_;
    std::unique_ptr<Indicator> indicator_;
    std::unique_ptr<Limiter> limiter_;
};