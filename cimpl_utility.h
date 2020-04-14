//
//  cimpl_utility.h
//  cimpl
//
//  Created by Nicholas Dwork starting on 11/16/17.
//  Copyright © 2017 Nicholas Dwork.
//

#ifndef cimpl_utility_h
#define cimpl_utility_h

#ifdef __cplusplus
extern "C" {
#endif
  
#include "stdlib.h"

void cimpl_bubbleSortf( float* arr, size_t n );
void cimpl_bubbleSortst( size_t* arr, size_t n );

void cimpl_halfMSE( float const * in, float const * b, size_t nIn, float *out );
void cimpl_halfMSE_grad( float const * in, float const * b, size_t N, float *out );
void cimpl_halfNormL2Sq( float const * const in, size_t const N, float *out );

void cimpl_linspace( float const minValue, float const maxValue, size_t n, float * const out );

void cimpl_normL2( float const * const in, size_t const N, float *out );

int cimpl_randui( int const low, int const high );
float cimpl_randuf( float const low, float const high );

float cimpl_sigmoid( float x );
float cimpl_dSigmoid( float x );

#ifdef __cplusplus
}
#endif


#endif  /* cimpl_utility_h */
