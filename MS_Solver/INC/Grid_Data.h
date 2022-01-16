#pragma once
#include "Figure.h"
#include "Element_Type.h"
#include "Euclidean_Vector.h"

#include <vector>

using uint = unsigned int;

struct Element_Data
{
	ElementType element_type;
	Figure figure;
	ushort figure_order;
	std::vector<uint> node_indexes;
};

struct Node_Data
{
	Euclidean_Vector coordinate_v;
};