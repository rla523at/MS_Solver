#pragma once
#include "Exception.h"
#include "Figure.h"

using ushort = unsigned short;


namespace Gmsh
{
	enum class FigureType
	{
		POINT = 0,
		LINE_P1 = 1, LINE_P2 = 8, LINE_P3 = 26, LINE_P4 = 27, LINE_P5 = 28, LINE_P6 = 62,
		TRIS_P1 = 2, TRIS_P2 = 9, TRIS_P3 = 21, TRIS_P4 = 23, TRIS_P5 = 25,
		QUAD_P1 = 3, QUAD_P2 = 10, QUAD_P3 = 36, QUAD_P4 = 37, QUAD_P5 = 38, QUAD_P6 = 47,
		TETS_P1 = 4, TETS_P2 = 11, TETS_P3 = 29, TETS_P4 = 30, TETS_P5 = 31,
		HEXA_P1 = 5, HEXA_P2 = 12, HEXA_P3 = 92, HEXA_P4 = 93, HEXA_P5 = 94,
		PRIS_P1 = 6, PRIS_P2 = 13, PRIS_P3 = 90, PRIS_P4 = 91, PRIS_P5 = 106,
		PYRA_P1 = 7, PYRA_P2 = 14, PYRA_P3 = 118, PYRA_P4 = 119
	};

	short figure_type_index_to_figure_order(const ushort figure_type_index)
	{
		switch (static_cast<FigureType>(figure_type_index))
		{
		case FigureType::POINT:			return 0;
		case FigureType::LINE_P1:
		case FigureType::TRIS_P1:
		case FigureType::QUAD_P1:
		case FigureType::TETS_P1:
		case FigureType::HEXA_P1:
		case FigureType::PRIS_P1:
		case FigureType::PYRA_P1:		return 1;
		case FigureType::LINE_P2:
		case FigureType::TRIS_P2:
		case FigureType::QUAD_P2:
		case FigureType::TETS_P2:
		case FigureType::HEXA_P2:
		case FigureType::PRIS_P2:
		case FigureType::PYRA_P2:		return 2;
		case FigureType::LINE_P3:
		case FigureType::TRIS_P3:
		case FigureType::QUAD_P3:
		case FigureType::TETS_P3:
		case FigureType::HEXA_P3:
		case FigureType::PRIS_P3:
		case FigureType::PYRA_P3:		return 3;
		case FigureType::LINE_P4:
		case FigureType::TRIS_P4:
		case FigureType::QUAD_P4:
		case FigureType::TETS_P4:
		case FigureType::HEXA_P4:
		case FigureType::PRIS_P4:
		case FigureType::PYRA_P4:		return 4;
		case FigureType::LINE_P5:
		case FigureType::TRIS_P5:
		case FigureType::QUAD_P5:
		case FigureType::TETS_P5:
		case FigureType::HEXA_P5:
		case FigureType::PRIS_P5:		return 5;
		case FigureType::LINE_P6:
		case FigureType::QUAD_P6:		return 6;
		default:
			EXCEPTION("invalid element type index");
			return -1;
		}
	}
	Figure figure_type_index_to_figure(const ushort figure_type_index)
	{
		switch (static_cast<FigureType>(figure_type_index))
		{
		case FigureType::POINT:			return Figure::point;
		case FigureType::LINE_P1:
		case FigureType::LINE_P2:
		case FigureType::LINE_P3:
		case FigureType::LINE_P4:
		case FigureType::LINE_P5:
		case FigureType::LINE_P6:		return Figure::line;
		case FigureType::TRIS_P1:
		case FigureType::TRIS_P2:
		case FigureType::TRIS_P3:
		case FigureType::TRIS_P4:
		case FigureType::TRIS_P5:		return Figure::triangle;
		case FigureType::QUAD_P1:
		case FigureType::QUAD_P2:
		case FigureType::QUAD_P3:
		case FigureType::QUAD_P4:
		case FigureType::QUAD_P5:
		case FigureType::QUAD_P6:		return Figure::quadrilateral;
		case FigureType::TETS_P1:
		case FigureType::TETS_P2:
		case FigureType::TETS_P3:
		case FigureType::TETS_P4:
		case FigureType::TETS_P5:		return Figure::tetrahedral;
		case FigureType::HEXA_P1:
		case FigureType::HEXA_P2:
		case FigureType::HEXA_P3:
		case FigureType::HEXA_P4:
		case FigureType::HEXA_P5:		return Figure::hexahedral;
		case FigureType::PRIS_P1:
		case FigureType::PRIS_P2:
		case FigureType::PRIS_P3:
		case FigureType::PRIS_P4:
		case FigureType::PRIS_P5:		return Figure::prism;
		case FigureType::PYRA_P1:
		case FigureType::PYRA_P2:
		case FigureType::PYRA_P3:
		case FigureType::PYRA_P4:		return Figure::pyramid;
		default:
			EXCEPTION("invalid element type index");
			return Figure::not_in_list;
		}
	}
}