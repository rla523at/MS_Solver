#pragma once

// ##########################################SETTING##################################################################
// choose one of options and comment out other options

//--------------- PATH ---------------//
// USER DEFINE OPTION
#define PATH "E:/Code/Result/MS_Solver/" + GOVERNING_EQUATION::name() + "/" + INITIAL_CONDITION::name() + "/" + SPATIAL_DISCRETE_METHOD::name() + "_" + RECONSTRUCTION_METHOD::name()  + "/" + GRID_FILE_NAME + "/"

//--------------- DIMENSION ---------------//
// USER DEFINE OPTION
#define DIMENSION 2

//--------------- GRID FILE TYPE -----------------//
#define __GMSH__

// USER DEFINE OPTION
#define GRID_FILE Quad50

//--------------- GOVERNING EQUATION -----------------//
#define __LINEAR_ADVECTION__
//#define __BURGERS__
//#define __EULER__

//USER DEFINE OPTION
#ifdef __LINEAR_ADVECTION__

#if DIMENSION >= 2
#define X_ADVECTION_SPEED 1.0
#define Y_ADVECTION_SPEED 1.0
#elif DIMENSION >= 3
#define Z_ADVECTION_SPEED 0.5
#endif

#define ERROR_CALCULATION_MODE

#endif



//--------------- INITIAL CONDITION -----------------//
#ifndef __EULER__
#define __SINE_WAVE__
//#define SQUARE_WAVE
#endif

#ifdef __EULER__
#define MODIFIED_SOD
#endif

//USER DEFINE OPTION
#ifdef __SINE_WAVE__
#if DIMENSION >= 2
#define X_WAVE_LENGTH 1.0
#define Y_WAVE_LENGTH 1.0
#elif DIMENSION >= 3
#define Z_WAVE_LENGTH 1
#endif
#endif

//--------------- SPATIAL DISCRETE METHOD -----------------//
//#define __FVM__
#define __HOM__

//--------------- RECONSTRUCTION METHOD -----------------//
#ifdef __FVM__
//#define __CONSTANT_RECONSTRUCTION__
//#define __LINEAR_RECONSTRUCTION__
#define __MLP_U1_RECONSTRUCTION__
//#define __ANN_RECONSTRUCTION__
#endif

#ifdef __HOM__
#define __POLYNOMIAL_RECONSTRUCTION__
#endif

//USER DEFINE OPTION
#ifdef __HOM__
#define SOLUTION_ORDER 1
#endif

//--------------- GRADIENT METHOD -----------------//
#ifdef __FVM__
#ifndef __CONSTANT_RECONSTRUCTION__
#define __VERTEX_LEAST_SQUARE__
//#define __FACE_LEAST_SQUARE__
#endif
#endif

//--------------- NUMERICAL FLUX -----------------//
#define __LLF__

//--------------- TIME INTEGRAL METHOD -----------------//
#define __SSPRK33__

//--------------- TIME STEP METHOD -----------------//
#define __CFL__
//#define __CONSTANT_DT__

//USER DEFINE OPTION
#define TIME_STEP_CONSTANT 0.9

//--------------- SOLVE END CONDITION -----------------//
#define __END_BY_TIME__
//#define __END_BY_ITER__

//USER DEFINE OPTION
#define END_CONDITION_CONSTANT 1.0

//--------------- SOLVE POST CONDITION -----------------//
#define __POST_BY_TIME__
//#define __POST_BY_ITER__

//USER DEFINE OPTION
#define POST_CONDITION_CONSTANT 0.1
#define POST_ORDER 1

// ########################################## SETTING  END ##################################################################





// ########################################## MACRO SETTING ##################################################################

//FORMAT SETTER
#define FORMAT1(x,y) x ## _ ## y ## D
#define SET_FORMAT1(x,y) FORMAT1(x,y) 
#define FORMAT2(x,y) x ## _ ## y 
#define SET_FORMAT2(x,y) FORMAT2(x,y)
#define STRING(x) #x
#define TO_STRING(x) STRING(x) 


#ifdef __GMSH__
#define GRID_FILE_TYPE	Gmsh
#endif

#define GRID_FILE_NAME TO_STRING(GRID_FILE)

#ifdef __LINEAR_ADVECTION__
#if DIMENSION == 2
#define GOVERNING_EQUATION		SET_FORMAT1(Linear_Advection, DIMENSION)<X_ADVECTION_SPEED,Y_ADVECTION_SPEED>
#endif
#endif
#ifdef __BURGERS__
#define GOVERNING_EQUATION		SET_FORMAT1(Burgers, DIMENSION)
#endif
#ifdef __EULER__
#define GOVERNING_EQUATION		SET_FORMAT1(Euler, DIMENSION)
#endif

#ifdef __SINE_WAVE__
#if DIMENSION == 2
#define INITIAL_CONDITION	SET_FORMAT1(Sine_Wave, DIMENSION)<X_WAVE_LENGTH, Y_WAVE_LENGTH>
#endif
#endif
#ifdef __SQUARE_WAVE__
#define INITIAL_CONDITION	SET_FORMAT1(Square_Wave, DIMENSION)
#endif
#ifdef __MODIFIED_SOD__
#define INITIAL_CONDITION	SET_FORMAT1(Modified_SOD, DIMENSION)
#endif

#ifdef __FVM__
#define SPATIAL_DISCRETE_METHOD	FVM
#endif
#ifdef __HOM__
#define SPATIAL_DISCRETE_METHOD	HOM
#endif

#ifndef __CONSTANT_RECONSTRUCTION__
#ifdef __VERTEX_LEAST_SQUARE__ 
#define GRADIENT_METHOD		Vertex_Least_Square<GOVERNING_EQUATION::num_equation(), DIMENSION>
#endif
#ifdef __FACE_LEAST_SQUARE__ 
#define GRADIENT_METHOD		Face_Least_Square<GOVERNING_EQUATION::num_equation(), DIMENSION>
#endif
#endif

#ifdef __CONSTANT_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD Constant_Reconstruction 
#endif
#ifdef __LINEAR_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD Linear_Reconstruction<GRADIENT_METHOD>
#endif
#ifdef __MLP_U1_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD MLP_u1<GRADIENT_METHOD>
#endif
#ifdef __ANN_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD ANN_limiter<GRADIENT_METHOD>
#endif

#ifdef __POLYNOMIAL_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD Polynomial_Reconstruction<DIMENSION, SOLUTION_ORDER>
#endif

#ifdef __LLF__
#define NUMERICAL_FLUX_FUNCTION	LLF<GOVERNING_EQUATION>
#endif

#ifdef __SSPRK33__
#define TIME_INTEGRAL_METHOD SSPRK33
#endif

#ifdef __CFL__
#define TIME_STEP_METHOD CFL<TIME_STEP_CONSTANT>
#endif
#ifdef __CONSTANT_DT__
#define TIME_STEP_METHOD Constant_Dt<TIME_STEP_CONSTANT>
#endif

#ifdef __END_BY_TIME__
#define SOLVE_END_CONDITION End_By_Time<END_CONDITION_CONSTANT>
#endif
#ifdef __END_BY_ITER__
#define SOLVE_END_CONDITION End_By_Iter<END_CONDITION_CONSTANT>
#endif

#ifdef __POST_BY_TIME__
#define SOLVE_POST_CONDITION Post_By_Time<POST_CONDITION_CONSTANT>
#endif
#ifdef __POST_BY_ITER__
#define SOLVE_POST_CONDITION Post_By_Iter<POST_CONDITION_CONSTANT>
#endif

// ########################################## MACRO SETTING END ##################################################################