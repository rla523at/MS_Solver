#pragma once
#include "Setting_Base.h"

// ########################################## OPTION ##################################################################

#define __DEFAULT_PATH__				"E:/Code/Result/MS_Solver/" + GOVERNING_EQUATION::name() + "/" + INITIAL_CONDITION::name() + "/" + SPATIAL_DISCRETE_METHOD::name() + "_" + RECONSTRUCTION_METHOD::name()  + "/debug/"

#define __DIMENSION__							2
#define __GRID_FILE_TYPE__						__GMSH__
#define __GOVERNING_EQUATION__					__LINEAR_ADVECTION__
#define __INITIAL_CONDITION__					__SQUARE_WAVE__
#define __SPATIAL_DISCRETE_METHOD__				__HOM__
#define __RECONSTRUCTION_METHOD__				__hMLP_RECONSTRUCTION__

#if		__SPATIAL_DISCRETE_METHOD__ ==	__FVM__ 
#if		__RECONSTRUCTION_METHOD__	!=			__CONSTANT_RECONSTRUCTION__
#define __GRADIENT_METHOD__						__VERTEX_LEAST_SQUARE__
#endif
#endif

#if		__SPATIAL_DISCRETE_METHOD__ ==	__HOM__
#define __SOLUTION_ORDER__						3
#endif 

#define __NUMERICAL_FLUX__						__LLF__
#define __TIME_INTEGRAL_METHOD__				__SSPRK54__
#define __TIME_STEP_METHOD__					__CFL__
#define __TIME_STEP_CONSTANT__					0.9
#define __SOLVE_END_CONDITION__					__END_BY_TIME__
#define __SOLVE_END_CONDITION_CONSTANT__		1.0
#define __SOLVE_POST_CONDITION__				__POST_BY_ITER__
#define __SOLVE_POST_CONDITION_CONSTANT__		1
#define __POST_ORDER__							1

// AVAILABLE OPTIONS
// __GRID_FILE_TYPE__				__GMSH__
// __GOVERNING_EQUATION__			__LINEAR_ADVECTION__, __BURGERS__, __EULER__
// __INITIAL_CONDITION__			__SINE_WAVE__, __SQUARE_WAVE__, __MODIFIED_SOD__, __CONSTANT1__
// __SPATIAL_DISCRETE_METHOD__		__FVM__, __HOM__
// __RECONSTRUCTION_METHOD__		__CONSTANT_RECONSTRUCTION__, __LINEAR_RECONSTRUCTION__,  __MLP_u1_RECONSTRUCTION__, __ANN_RECONSTRUCTION__
//									__POLYNOMIAL_RECONSTRUCTION__, __hMLP_RECONSTRUCTION__, __hMLP_BD_RECONSTRUCTION__
// __GRADIENT_METHOD__				__VERTEX_LEAST_SQUARE__, __FACE_LEAST_SQUARE__ 
// __NUMERICAL_FLUX__				__LLF__
// __TIME_INTEGRAL_METHOD__			__SSPRK33__, __SSPRK54__
// __TIME_STEP_METHOD__				__CFL__, __CONSTANT_DT__
// __SOLVE_END_CONDITION__			__END_BY_TIME__, __END_BY_ITER__
// __SOLVE_POST_CONDITION__			__POST_BY_TIME__, __POST_BY_ITER__

// ######################################### OPTION END ################################################################

// #################################### USER DEFINE SETTING ############################################################
// Linear Advection
#define X_ADVECTION_SPEED				1.0
#define Y_ADVECTION_SPEED				2.0
#define Z_ADVECTION_SPEED				0.5

// Sine Wave
#define X_WAVE_LENGTH					1.0
#define Y_WAVE_LENGTH					1.0
#define Z_WAVE_LENGTH					1

// Mode (comment out == turn off)
//#define ERROR_CALCULATION_MODE
#if __GOVERNING_EQUATION__		== 		__EULER__
#if	__SPATIAL_DISCRETE_METHOD__ ==		__HOM__
#define PRESSURE_FIX_MODE
#endif
#endif



// ################################# USER DEFINE SETTING END #########################################################
 







// ########################################## MACRO SETTING ##################################################################

#if		__GRID_FILE_TYPE__ == __GMSH__
#define GRID_FILE_TYPE	Gmsh
#endif

#define GRID_FILE_NAME TO_STRING(__GRID_FILE__)


#if		__GOVERNING_EQUATION__ == __LINEAR_ADVECTION__
#if		__DIMENSION__ == 2
#define GOVERNING_EQUATION		SET_FORMAT1(Linear_Advection, __DIMENSION__)<X_ADVECTION_SPEED,Y_ADVECTION_SPEED>
#endif
#ifdef  ERROR_CALCULATION_MODE
#define ERROR_CALCULATION
#endif
#endif
#if		__GOVERNING_EQUATION__ ==__BURGERS__
#define GOVERNING_EQUATION		SET_FORMAT1(Burgers, __DIMENSION__)
#endif
#if		__GOVERNING_EQUATION__ ==__EULER__
#define GOVERNING_EQUATION		SET_FORMAT1(Euler, __DIMENSION__)
#endif

#if		__INITIAL_CONDITION__ == __SINE_WAVE__
#if		__DIMENSION__ == 2
#define INITIAL_CONDITION	SET_FORMAT1(Sine_Wave, __DIMENSION__)<X_WAVE_LENGTH, Y_WAVE_LENGTH>
#endif
#endif
#if		__INITIAL_CONDITION__ == __SQUARE_WAVE__
#define INITIAL_CONDITION	SET_FORMAT1(Square_Wave, __DIMENSION__)
#endif
#if		__INITIAL_CONDITION__ == __MODIFIED_SOD__
#define INITIAL_CONDITION	SET_FORMAT1(Modified_SOD, __DIMENSION__)
#endif
#if		__INITIAL_CONDITION__ == __CONSTANT1__
#define INITIAL_CONDITION	SET_FORMAT1(Constant1, __DIMENSION__)
#endif

#if		__SPATIAL_DISCRETE_METHOD__ == __FVM__
#define SPATIAL_DISCRETE_METHOD	FVM
#endif
#if		__SPATIAL_DISCRETE_METHOD__ == __HOM__
#define SPATIAL_DISCRETE_METHOD	HOM
#endif

#if		__SPATIAL_DISCRETE_METHOD__ == __FVM__
#if		__RECONSTRUCTION_METHOD__ != __CONSTANT_RECONSTRUCTION__
#if		__GRADIENT_METHOD__ == __VERTEX_LEAST_SQUARE__ 
#define GRADIENT_METHOD		Vertex_Least_Square<GOVERNING_EQUATION::num_equation(), __DIMENSION__>
#endif
#if		__GRADIENT_METHOD__ == __FACE_LEAST_SQUARE__ 
#define GRADIENT_METHOD		Face_Least_Square<GOVERNING_EQUATION::num_equation(), __DIMENSION__>
#endif
#endif
#endif

#if		__RECONSTRUCTION_METHOD__ == __CONSTANT_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD Constant_Reconstruction 
#endif
#if		__RECONSTRUCTION_METHOD__ == __LINEAR_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD Linear_Reconstruction<GRADIENT_METHOD>
#endif
#if		__RECONSTRUCTION_METHOD__ == __MLP_u1_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD MLP_u1<GRADIENT_METHOD>
#endif
#if		__RECONSTRUCTION_METHOD__ == __ANN_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD ANN_limiter<GRADIENT_METHOD>
#endif

#if		__RECONSTRUCTION_METHOD__ == __POLYNOMIAL_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD Polynomial_Reconstruction<__DIMENSION__, __SOLUTION_ORDER__>
#endif
#if		__RECONSTRUCTION_METHOD__ == __hMLP_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD hMLP_Reconstruction<__DIMENSION__, __SOLUTION_ORDER__>
#endif
#if		__RECONSTRUCTION_METHOD__ == __hMLP_BD_RECONSTRUCTION__
#define RECONSTRUCTION_METHOD hMLP_BD_Reconstruction<__DIMENSION__, __SOLUTION_ORDER__>
#endif


#if		__NUMERICAL_FLUX__ == __LLF__
#define NUMERICAL_FLUX_FUNCTION	LLF<GOVERNING_EQUATION>
#endif

#if		__TIME_INTEGRAL_METHOD__ == __SSPRK33__
#define TIME_INTEGRAL_METHOD SSPRK33
#endif
#if		__TIME_INTEGRAL_METHOD__ == __SSPRK54__
#define TIME_INTEGRAL_METHOD SSPRK54
#endif

#if		__TIME_STEP_METHOD__ == __CFL__
#define TIME_STEP_METHOD CFL<__TIME_STEP_CONSTANT__>
#endif
#if		__TIME_STEP_METHOD__ == __CONSTANT_DT__
#define TIME_STEP_METHOD Constant_Dt<__TIME_STEP_CONSTANT__>
#endif

#if		__SOLVE_END_CONDITION__ == __END_BY_TIME__
#define SOLVE_END_CONDITION		Controll_Condition::by_time, __SOLVE_END_CONDITION_CONSTANT__
#endif
#if		__SOLVE_END_CONDITION__ == __END_BY_ITER__
#define SOLVE_END_CONDITION		Controll_Condition::by_iter, __SOLVE_END_CONDITION_CONSTANT__
#endif

#if		__SOLVE_POST_CONDITION__ == __POST_BY_TIME__
#define SOLVE_POST_CONDITION	Controll_Condition::by_time, __SOLVE_POST_CONDITION_CONSTANT__
#endif
#if		__SOLVE_POST_CONDITION__ == __POST_BY_ITER__
#define SOLVE_POST_CONDITION	Controll_Condition::by_iter, __SOLVE_POST_CONDITION_CONSTANT__
#endif

// ########################################## MACRO SETTING END ##################################################################