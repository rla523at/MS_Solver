#pragma once
#include "Setting_Base.h"


// ########################################## OPTION ##################################################################

#define __DEFAULT_PATH__						"E:/CodeData/Result/MS_Solver/_Temp/" + GOVERNING_EQUATION::name() + "/" + INITIAL_CONDITION::name() + "/" + SPATIAL_DISCRETE_METHOD::name() + "_" + RECONSTRUCTION_METHOD::name() + "/"
#define __DIMENSION__							2
#define __GRID_FILE_TYPE__						__GMSH__
#define __GRID_FILE_NAMES__						Shocktube_Quad_400x4
#define __GOVERNING_EQUATION__					__EULER__
#define __INITIAL_CONDITION__					__BLAST_WAVE_INTERACTION__
#define __SPATIAL_DISCRETE_METHOD__				__HOM__
#define __RECONSTRUCTION_METHOD__				__hMLP_BD_RECONSTRUCTION__
#define __NUMERICAL_FLUX__						__LLF__
#define __TIME_INTEGRAL_METHOD__				__SSPRK54__
#define __TIME_STEP_METHOD__					__CFL__
#define __TIME_STEP_CONSTANT__					0.9
#define __SOLVE_END_CONDITION__					__BY_TIME__
#define __SOLVE_END_CONDITION_CONSTANT__		0.038
#define __SOLVE_POST_CONDITION__				__BY_ITER__
#define __SOLVE_POST_CONDITION_CONSTANT__		50
#define __POST_ORDER__							4
#define __POST_FILE_FORMAT__					__BINARY__

// CONDITIONAL OPTIONS
#if  __RECONSTRUCTION_METHOD__ == __ANN_RECONSTRUCTION__
#define __ANN_MODEL__							model5_1
#endif

#if		__SPATIAL_DISCRETE_METHOD__ ==	__HOM__
#define __SOLUTION_ORDER__						4
#endif 

//temp
#define __hMLP_BD_TYPE__				BD_Type::standard

// AVAILABLE OPTIONS
// __GRID_FILE_TYPE__				__GMSH__
// __GOVERNING_EQUATION__			__LINEAR_ADVECTION__, __BURGERS__, __EULER__
// __INITIAL_CONDITION__			__SINE_WAVE__, __SQUARE_WAVE__, __CIRCLE_WAVE__, __GAUSSIAN_WAVE__, __CONSTANT1__,
//									__SOD__, __MODIFIED_SOD__, __SUPERSONIC_EXPANSION__, __BLAST_WAVE_PROBLEM__, __DOUBLE_STRONG_SHOCK_PROBLEM__, __SLOWLY_MOVING_CONTACT__
//									__SHU_OSHER__, __HARTEN_LAX_PROBLEM__, __BLAST_WAVE_INTERACTION__, __EXPLOSION_PROBLEM__
// __SPATIAL_DISCRETE_METHOD__		__FVM__, __HOM__
// __RECONSTRUCTION_METHOD__		__CONSTANT_RECONSTRUCTION__, __LINEAR_RECONSTRUCTION__,  __MLP_u1_RECONSTRUCTION__, __ANN_RECONSTRUCTION__
//									__POLYNOMIAL_RECONSTRUCTION__, __hMLP_RECONSTRUCTION__, __hMLP_BD_RECONSTRUCTION__
// __NUMERICAL_FLUX__				__LLF__
// __TIME_INTEGRAL_METHOD__			__SSPRK33__, __SSPRK54__
// __TIME_STEP_METHOD__				__CFL__, __CONSTANT_DT__ 
// __SOLVE_END_CONDITION__			__END_BY_TIME__, __END_BY_ITER__
// __SOLVE_POST_CONDITION__			__POST_BY_TIME__, __POST_BY_ITER__
// __POST_MODE__					__ASCII__, __BINARY__

// Reference Constant
// CFL		: Modified SOD(0.9), Supersonic Expansion(0.5), Blast wave problem(0.6), Double strong shock problem(0.8), Slowly-moving contact(0.6)
//			   Harten Lax(0.5), Shu_Osher(0.5), Blast Wave Interaction(0.5)	
// END TIME : Modified SOD(0.2), Supersonic Expansion(0.15), Blast wave problem(0.012), Double strong shock problem(0.035), Slowly-moving contact(0.012)
//			  Harten Lax(0.15), Shu_Osher(0.178), Explosion(0.25), Blast Wave Interaction(0.038)


// ######################################### OPTION END ################################################################

// #################################### USER DEFINED SETTING ############################################################

// MODE
// comment out if you do not want to use
#define __USE_SCAILING_METHOD__

// Linear Advection
#define X_ADVECTION_SPEED				1.0
#define Y_ADVECTION_SPEED				1.0
#define Z_ADVECTION_SPEED				1.0

// Sine Wave
#define X_WAVE_LENGTH					1.0
#define Y_WAVE_LENGTH					1.0
#define Z_WAVE_LENGTH					1.0

// Supersonic Inlet inflow value
#define INFLOW1_VALUES					1.0, 0.0, 0.0, 2500
#define INFLOW2_VALUES					1.0, 0.0, 0.0, 250

// Reference Constants
// Modified SOD				1.0, 0.75, 0.0, 2.78125 
//							0.125, 0.0, 0.0, 0.25
// Supersonic Expansion		1.0, -2.0, 0.0, 3.0
//							1.0, 2.0, 0.0, 3.0
// Blast Wave				1.0, 0.0, 0.0, 2500
//							1.0, 0.0, 0.0, 0.025
// Double strong shock		5.99924, 117.5701059, 0.0, 2304.275075187625134
//							5.99242, -37.1310118186, 0.0, 230.2755012309728784
// Slowly moving contact	1.0, -19.59745, 0.0, 2692.03002325125
//							1.0, -19.59745, 0.0, 192.05502325125
// Harten Lax				0.445, 0.31061, 0.0, 8.92840289
//							0.5, 0.0, 0.0, 1.4275
// Shu Osher				3.857143, 10.1418522328, 0.0, 39.1666684317
//							1.0, 0.0, 0.0, 2.5
// Blast wave interaction	1.0, 0.0, 0.0, 2500
//							1.0, 0.0, 0.0, 250


// ################################# USER DEFINED SETTING END #########################################################


// ########################################## MACRO SETTING ##################################################################


#if		__GRID_FILE_TYPE__ == __GMSH__
#define GRID_FILE_TYPE	Gmsh
#endif

#define GRID_FILE_NAMES TO_STRING(__GRID_FILE_NAMES__)


#if		__GOVERNING_EQUATION__ == __LINEAR_ADVECTION__
#define GOVERNING_EQUATION		Linear_Advection<__DIMENSION__>
#endif
#if		__GOVERNING_EQUATION__ ==__BURGERS__
#define GOVERNING_EQUATION		Burgers<__DIMENSION__>
#endif
#if		__GOVERNING_EQUATION__ ==__EULER__
#define GOVERNING_EQUATION		Euler<__DIMENSION__>
#endif

#if		__INITIAL_CONDITION__ == __SINE_WAVE__
#define INITIAL_CONDITION	Sine_Wave<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __SQUARE_WAVE__
#define INITIAL_CONDITION	Square_Wave<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __CIRCLE_WAVE__
#define INITIAL_CONDITION	Circle_Wave<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __GAUSSIAN_WAVE__
#define INITIAL_CONDITION	Gaussian_Wave<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __CONSTANT1__
#define INITIAL_CONDITION	Constant1<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __SOD__
#define INITIAL_CONDITION	SOD<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __MODIFIED_SOD__
#define INITIAL_CONDITION	Modified_SOD<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __SHU_OSHER__
#define INITIAL_CONDITION	Shu_Osher<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __EXPLOSION_PROBLEM__
#define INITIAL_CONDITION	Explosion_Problem<__DIMENSION__>
#endif 
#if		__INITIAL_CONDITION__ == __DOUBLE_RAREFACTION_WAVE__
#define INITIAL_CONDITION	Double_Rarefaction_Wave<__DIMENSION__>
#endif 
#if		__INITIAL_CONDITION__ == __HARTEN_LAX_PROBLEM__
#define INITIAL_CONDITION	Harten_Lax_Problem<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __SUPERSONIC_EXPANSION__
#define INITIAL_CONDITION	Supersonic_Expansion<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __BLAST_WAVE_PROBLEM__
#define INITIAL_CONDITION	Blast_Wave_Problem<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __DOUBLE_STRONG_SHOCK_PROBLEM__
#define INITIAL_CONDITION	Double_Strong_Shock_Problem<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __SLOWLY_MOVING_CONTACT__
#define INITIAL_CONDITION	Slowly_Moving_Contact_Problem<__DIMENSION__>
#endif
#if		__INITIAL_CONDITION__ == __BLAST_WAVE_INTERACTION__
#define INITIAL_CONDITION	Blast_Wave_Interaction<__DIMENSION__>
#endif


#if		__SPATIAL_DISCRETE_METHOD__ == __FVM__
#define SPATIAL_DISCRETE_METHOD	FVM
#endif
#if		__SPATIAL_DISCRETE_METHOD__ == __HOM__
#define SPATIAL_DISCRETE_METHOD	HOM
#endif

#if		__SPATIAL_DISCRETE_METHOD__ == __FVM__
#if		__RECONSTRUCTION_METHOD__ != __CONSTANT_RECONSTRUCTION__

#define GRADIENT_METHOD		Face_Least_Square<GOVERNING_EQUATION::num_equation(), __DIMENSION__>

//#if		__GRADIENT_METHOD__ == __VERTEX_LEAST_SQUARE__ 
//#define GRADIENT_METHOD		Vertex_Least_Square<GOVERNING_EQUATION::num_equation(), __DIMENSION__>
//#endif
//#if		__GRADIENT_METHOD__ == __FACE_LEAST_SQUARE__ 
//#define GRADIENT_METHOD		Face_Least_Square<GOVERNING_EQUATION::num_equation(), __DIMENSION__>
//#endif
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

#if		__SOLVE_END_CONDITION__ == __BY_TIME__
#define SOLVE_END_CONDITION		Controll_Condition::by_time, __SOLVE_END_CONDITION_CONSTANT__
#endif
#if		__SOLVE_END_CONDITION__ == __BY_ITER__
#define SOLVE_END_CONDITION		Controll_Condition::by_iter, __SOLVE_END_CONDITION_CONSTANT__
#endif

#if		__SOLVE_POST_CONDITION__ == __BY_TIME__
#define SOLVE_POST_CONDITION	Controll_Condition::by_time, __SOLVE_POST_CONDITION_CONSTANT__
#endif
#if		__SOLVE_POST_CONDITION__ == __BY_ITER__
#define SOLVE_POST_CONDITION	Controll_Condition::by_iter, __SOLVE_POST_CONDITION_CONSTANT__
#endif

#if		__POST_FILE_FORMAT__ == __ASCII__
#define POST_FILE_FORMAT  Post_File_Format::ASCII
#endif
#if		__POST_FILE_FORMAT__ == __BINARY__
#define POST_FILE_FORMAT  Post_File_Format::binary
#endif

#ifdef __USE_SCAILING_METHOD__
#define SCAILING_METHOD_FLAG true
#else
#define SCAILING_METHOD_FLAG false
#endif

// ########################################## MACRO SETTING END ##################################################################

