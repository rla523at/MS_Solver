#pragma once
#include "Inner_Faces_FVM.h"
#include "Spatial_Discrete_Method.h"


template <typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
class Inner_Faces;


template <typename Numerical_Flux_Function>
class Inner_Faces<FVM, Constant_Reconstruction, Numerical_Flux_Function> : public Inner_Faces_FVM_Constant<Numerical_Flux_Function>
{
private:
    Inner_Faces(void) = delete;
};


template<typename Reconstruction_Method, typename Numerical_Flux_Function>
class Inner_Faces<FVM, Reconstruction_Method, Numerical_Flux_Function> : public Inner_Faces_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>
{
private:
    Inner_Faces(void) = delete;
};


//template<typename Reconstruction_Method, typename Numerical_Flux_Function>
//class Inner_Faces<HOM, Reconstruction_Method, Numerical_Flux_Function> : public Inner_Faces_HOM<Reconstruction_Method, Numerical_Flux_Function>
//{
//private:
//    Inner_Faces(void) = delete;
//};