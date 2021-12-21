#pragma once
#include "Configuration.h"
#include "Semi_Discrete_Equation.h"

class Time_Discrete_Scheme
{
public:
	virtual void update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const abstract;    
};

class SSPRK33 : public Time_Discrete_Scheme
{
public:
    void update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const override;
};

class SSPRK54 : public Time_Discrete_Scheme
{
public:
    void update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const override;
};

class Time_Discrete_Scheme_Factory//static class
{
public:
    static std::unique_ptr<Time_Discrete_Scheme> make_unique(const Configuration& configuration);

private:
    Time_Discrete_Scheme_Factory(void) = delete;
};