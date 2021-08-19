#pragma once

// SETTING CODE

// GRID FILE TYPE
#define __GMSH__						10
// GOVERNING EQUATION
#define __LINEAR_ADVECTION__			20
#define __BURGERS__						21
#define __EULER__						22
// INITIAL CONDITION
#define __SINE_WAVE__					30
#define __SQUARE_WAVE__					31
#define __MODIFIED_SOD__				32
#define __CONSTANT1__					33
// SPATIAL DISCRETE METHOD
#define __FVM__							40
#define __HOM__							41
// RECONSTRUCTION METHOD
#define __CONSTANT_RECONSTRUCTION__		50
#define __LINEAR_RECONSTRUCTION__		51
#define __MLP_u1_RECONSTRUCTION__		52
#define __ANN_RECONSTRUCTION__			53
#define __POLYNOMIAL_RECONSTRUCTION__	54
#define __hMLP_RECONSTRUCTION__			55
// NUMERICAL FLUX
#define __LLF__							60
// TIME INTEGRAL METHOD
#define __SSPRK33__						70
#define __SSPRK54__						71
// TIME STEP METHOD
#define __CFL__							80
#define __CONSTANT_DT__					81
// SOLVE END CONDITION
#define __END_BY_TIME__					90
#define __END_BY_ITER__					91
// SOLVE POST CONDITION
#define __POST_BY_TIME__				100
#define __POST_BY_ITER__				101
// GRADIENT METHOD
#define __VERTEX_LEAST_SQUARE__			110
#define __FACE_LEAST_SQUARE__			111


//FORMAT SETTER
#define FORMAT1(x,y) x ## _ ## y ## D
#define SET_FORMAT1(x,y) FORMAT1(x,y) 
#define FORMAT2(x,y) x ## _ ## y 
#define SET_FORMAT2(x,y) FORMAT2(x,y)
#define STRING(x) #x
#define TO_STRING(x) STRING(x) 