#pragma once
enum class ElementType
{
	cell, face,
	slip_wall,
	supersonic_inlet1, supersonic_inlet2,
	supersonic_outlet,
	initial_constant_BC,
	periodic,
	not_in_list
};

