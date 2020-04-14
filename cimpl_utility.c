//
//  cimpl_utility.c
//  cimpl
//
//  Created by Nicholas Dwork starting on 11/16/17.
//  Copyright © 2017 Nicholas Dwork.
//

#include <assert.h>
#include <math.h>
#include "cimpl_utility.h"


void cimpl_bubbleSortf( float* arr, size_t n ){
  // Highly inefficient in-place sorting algorithm.
  float tmp;
  for( size_t i=0; i<n-1; ++i ){
    for( size_t j=0; j < n-i-1; ++j ){
      if( arr[j] > arr[j+1] ){
        tmp = arr[j]; arr[j] = arr[j+1]; arr[j+1] = tmp;  // swap
  } } }
}

void cimpl_bubbleSortst( size_t* arr, size_t n ){
  // Highly inefficient in-place sorting algorithm.
  int i, j;
  size_t tmp;
  for( i=0; i<n-1; ++i ){
    for( j=0; j < n-i-1; ++j ){
      if( arr[j] > arr[j+1] ){
        tmp = arr[j]; arr[j] = arr[j+1]; arr[j+1] = tmp;  // swap
  } } } 
}

int cimpl_randui( int const low, int const high ){
  int tmp;

  tmp = rand() / ( 1.0 + RAND_MAX ) * ( high - low + 1 );
  return tmp + low;
}

void cimpl_halfMSE( float const * in, float const * b, size_t nIn, float *out ){
  // Computes half the Mean Square Error between in and b
  size_t i;
  float diff;

  *out = 0;
  for( i=0; i<nIn; ++i ){
    diff = in[i] - b[i];
    *out += diff * diff;
  }
  *out /= ( 2 * nIn );
}

void cimpl_halfMSE_grad( float const * in, float const * b, size_t N, float *out ){
  // Computes the gradient of cimpl_halfMSE evaluated at in
  size_t i;

  for( i=0; i<N; ++i )
    out[i] = ( in[i] - b[i] ) / (float) N;
}

void cimpl_halfNormL2Sq( float const * const in, size_t const N, float *out ){
  *out=0;
  size_t i;
  for( i=0; i<N; ++i ){
    *out += in[i] * in[i];
  }
  *out = 0.5 * (*out) ;
}

void cimpl_linspace( float const minValue, float const maxValue, size_t n, float * const out ){
  float dx = ( maxValue - minValue ) / (float)( n - 1 );
  for( size_t i=0; i<n; ++i )
    out[i] = minValue + i * dx;
}

void cimpl_normL2( float const * const in, size_t const N, float *out ){
  *out=0;
  size_t i;
  for( i=0; i<N; ++i ){
    *out += in[i] * in[i];
  }
  *out = sqrtf( *out );
}

float cimpl_randuf( float const low, float const high ){
  // Create a uniformly distributed random variable in [low,high]
  float out;

  out = (float) rand( ) / (float) RAND_MAX * ( high - low );
  out += low;
  return out;
}

float cimpl_dSigmoid( float const x ){
  float sx, out;
  sx = cimpl_sigmoid( x );
  out = sx * ( 1 - sx );
  return out;
}

float cimpl_sigmoid( float const x ){
  float out;
  out = 1 / ( 1 + expf(-x) );
  return out;
}
