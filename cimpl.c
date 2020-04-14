//
//  cimpl.c
//  cimpl
//
//  Created by Nicholas Dwork starting on 9/18/16.
//  Copyright © 2016 Nicholas Dwork.
//


#ifdef _MSC_VER
#define _USE_MATH_DEFINE
#endif

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#include <malloc.h>
#endif

#include "cimpl.h"
#include "cimpl_utility.h"

#define cimpl_min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define cimpl_max(X, Y) (((X) > (Y)) ? (X) : (Y))


// ----- Static Functions -----

// +,-,*,/ of complex numbers are not needed for gcc
static inline CIMPL_COMPLEX cimpl_addCmp(CIMPL_COMPLEX const x, CIMPL_COMPLEX const y ) {
#ifdef __GNUC__
  return x+y;
#elif defined( _MSC_VER )
  CIMPL_COMPLEX out;
  out._Val[0] = x._Val[0] + y._Val[0];
  out._Val[1] = x._Val[1] + y._Val[1];
  return out;
#else
  CIMPL_COMPLEX out;
  float* outPtr = (float*) &out;
  *outPtr = crealf( x ) + crealf( y );
  *(outPtr+1) = cimagf( x ) + cimagf( y );
  return out;
#endif
}

static inline CIMPL_COMPLEX cimpl_divideCmp(CIMPL_COMPLEX const x, CIMPL_COMPLEX const y ){
#ifdef __GNUC__
  return x / y;
#elif defined( _MSC_VER )
  CIMPL_COMPLEX out;
  float tmp = y._Val[0] * y._Val[0] + y._Val[1] * y._Val[1];
  out._Val[0] = (x._Val[0] * y._Val[0] + x._Val[1] * y._Val[1]) / tmp;
  out._Val[1] = (x._Val[1] * y._Val[0] - x._Val[0] * y._Val[1]) / tmp;
  return out;
#else
  CIMPL_COMPLEX out;
  float* outPtr = (float*) &out;
  float tmp = crealf(y) * crealf(y) + cimagf(y) * cimagf(y);
  *outPtr = ( crealf(x) * crealf(y) + cimagf(x) * cimagf(y) ) / tmp;
  *(outPtr+1) = ( cimagf(y) * crealf(y) - crealf(x) * cimagf(y) ) / tmp;
#endif
}

static inline float cimpl_imag( CIMPL_COMPLEX const x ){
#ifdef __GNUC__
  return __imag__ x;
#elif defined( _MSC_VER )
  return x._Val[1];
#else
  return cimagf( x );
#endif
}

static inline CIMPL_COMPLEX cimpl_multiplyCmp( CIMPL_COMPLEX const x, CIMPL_COMPLEX const y ){
#ifdef __GNUC__
  return x * y;
#elif defined( _MSC_VER )
  CIMPL_COMPLEX out;
  out._Val[0] = x._Val[0] * y._Val[0] - x._Val[1] * y._Val[1];
  out._Val[1] = x._Val[1] * y._Val[0] + x._Val[0] * y._Val[1];
  return out;
#else
  CIMPL_COMPLEX out;
  float* outPtr = (float*) &out;
  *outPtr = crealf(x) * crealf(y) - cimagf(x) * cimagf(y);
  *(outPtr+1) = cimagf(x) * crealf(y) + crealf(x) * cimagf(y);
#endif
}

static inline float cimpl_real( CIMPL_COMPLEX const x ){
#ifdef __GNUC__
  return __real__ x;
#elif defined( _MSC_VER )
  return x._Val[0];
#else
  return crealf( x );
#endif
}

static inline float cimpl_absCmp( CIMPL_COMPLEX const x ){
  return sqrtf( cimpl_real(x) * cimpl_real(x) + cimpl_imag(x) * cimpl_imag(x) );
}

static inline float cimpl_sign( float const in ){
  // Returns 1 if positive, -1 if negative, and 0 if 0.
  return ( in > 0 ) - ( in < 0 );
}

static inline CIMPL_COMPLEX cimpl_subtractCmp(CIMPL_COMPLEX const x, CIMPL_COMPLEX const y) {
#ifdef __GNUC__
  return x - y;
#elif defined( _MSC_VER )
  CIMPL_COMPLEX out;
  out._Val[0] = x._Val[0] - y._Val[0];
  out._Val[1] = x._Val[1] - y._Val[1];
  return out;
#else
  CIMPL_COMPLEX out;
  float* outPtr = (float*) &out;
  *outPtr = crealf( x ) - crealf( y );
  *(outPtr+1) = cimagf( x ) - cimagf( y );
  return out;
#endif
}


// ----- Public Functions -----


void cimpl_absCmpImg( cimpl_cmpImg const in, cimpl_img * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  for( size_t i=0; i<in.nPix; ++i )
    out->data[i] = cabsf( in.data[i] );
}

void cimpl_absCmpVol( cimpl_cmpVol const in, cimpl_vol * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );  assert( out->s == in.s );
  for( size_t i=0; i<in.nVox; ++i )
    out->data[i] = cabsf( in.data[i] );
}

void cimpl_absImg( cimpl_img const in, cimpl_img * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  int i;

  for( i=0; i<in.nPix; ++i )
    out->data[i] = fabsf( in.data[i] );
}

void cimpl_absVol( cimpl_vol const in, cimpl_vol * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  int i;

  for( i=0; i<in.nVox; ++i )
    out->data[i] = fabsf( in.data[i] );
}

void cimpl_addCmpImgs( cimpl_cmpImg const img1, cimpl_cmpImg const img2,
  cimpl_cmpImg * const out ){
  assert( out->h == img1.h );  assert( out->w == img1.w );
  assert( out->h == img2.h );  assert( out->w == img2.w );
  int i;
  for( i=0; i<img1.nPix; ++i )
    out->data[i] = cimpl_addCmp( img1.data[i], img2.data[i] );
}

void cimpl_addImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->w == img1.w );  assert( out->h == img1.h );
  assert( out->w == img2.w );  assert( out->h == img2.h );
  int i;
  float* outData = out->data;

  for( i=0; i<img1.nPix; ++i, ++outData )
    *outData = img1.data[i] + img2.data[i];
}

void cimpl_addVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->w == vol1.w );  assert( out->h == vol1.h );  assert( out->s == vol1.s );
  assert( out->w == vol2.w );  assert( out->h == vol2.h );  assert( out->s == vol2.s );
  int i;
  float* outData = out->data;
  for( i=0; i<vol1.nVox; ++i, ++outData )
    *outData = vol1.data[i] + vol2.data[i];
}

void cimpl_addScalar2Img( float const scalar, cimpl_img const in, cimpl_img * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  int i;
  float* outData = out->data;

  for( i=0; i<in.nPix; ++i, ++outData )
    *outData = in.data[i] + scalar;
}

void cimpl_addScalar2Vol( float const scalar, cimpl_vol const in, cimpl_vol * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );  assert( out->s == in.s );
  int i;

  for( i=0; i<in.nVox; ++i )
    out->data[i] = in.data[i] + scalar;
}

void cimpl_argCmpImg( cimpl_cmpImg const in, cimpl_img * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  for( size_t i=0; i<in.nPix; ++i )
    out->data[i] = cargf( in.data[i] );
}

void cimpl_argCmpVol( cimpl_cmpVol const in, cimpl_vol * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );  assert( out->s == in.s );
  for( size_t i=0; i<in.nVox; ++i )
    out->data[i] = cargf( in.data[i] );
}

float cimpl_besseli0( float x ){
  // Modified Bessel function of the first kind and order 0 at input x
  float ax,ans;
  float y;
  if ((ax=fabsf(x)) < 3.75) {
    y=x/3.75; y=y*y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
      +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    y=3.75/ax;
    ans=(expf(ax)/sqrtf(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+
      y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
      +y*0.392377e-2))))))));
  }
  return ans;
}

float cimpl_besselj0( float x ){
  // Bessel function of the first kind and order 0 at input x
  float ax,z;
  float xx,y,ans,ans1,ans2;

  if ((ax=fabsf(x)) < 8.0){
    y=x*x;
    ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
      +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
      +y*(59272.64853+y*(267.8532712+y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-0.785398164;
    ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4 +
      y*(-0.2073370639e-5+y*0.2093887211e-6) ) );
    ans2 = -0.1562499995e-1+y*(0.1430488765e-3 + y * ( -0.6911147651e-5 +
      y*( 0.7621095161e-6 - y*0.934935152e-7) ) );
    ans = sqrtf(0.636619772/ax) * (cosf(xx)*ans1-z*sinf(xx)*ans2);
  }
  return ans;
}

void cimpl_ceilImg( cimpl_img const img, cimpl_img * const out ){
  assert( out->h == img.h );  assert( out->w == img.w );
  size_t i;
  for( i=0; i<img.nPix; ++i )
    out->data[i] = ceilf( img.data[i] );
}

void cimpl_ceilVol( cimpl_vol const vol, cimpl_vol * const out ){
  assert( out->h == vol.h );  assert( out->w == vol.w );  assert( out->s == vol.s );
  size_t i;

  for( i=0; i<vol.nVox; ++i )
    out->data[i] = ceilf( vol.data[i] );
}

void cimpl_circShiftImg( cimpl_img const in, long hShift, long vShift, cimpl_img * const out ){
  assert( out->w == in.w );  assert( out->h == in.h );
  assert( out->data != in.data );
  long inX, inY;

  for( long x=0; x < (long) in.w; ++x ){
    inX = x - hShift;
    inX = inX % in.w;
    for( long y=0; y < (long) in.h; ++y ){
      inY = y - vShift;
      inY = inY % in.h;

      out->data[y+x*in.h] = in.data[inY+inX*in.h];
  } }
}

void cimpl_circShiftVol( cimpl_vol const in, int hShift, int vShift, int sShift,
  cimpl_vol * const out ){
  assert( out->w == in.w );  assert( out->h == in.h );  assert( out->s == in.s );
  assert( out->data != in.data );
  long inX, inY, inZ;

  for( long x=0; x < (long) in.w; ++x ){
    inX = x + hShift;
    inX = inX % in.w;
    for( long y=0; y < (long) in.h; ++y ){
      inY = y + vShift;
      inY = inY % in.h;
      for( long z=0; z < (long) in.s; ++z ){
        inZ = z + sShift;
        inZ = inZ % in.s;

        out->data[y+x*in.h+z*in.h*in.w] = in.data[inY+inX*in.h+inZ*in.h*in.w];
  } } }
}

void cimpl_cmpEqImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->h == img1.h );   assert( out->w == img1.w );
  assert( out->h == img2.h );   assert( out->w == img2.w );
  int i;

  for( i=0; i<img1.nPix; ++i )
    out->data[i] = ( img1.data[i] == img2.data[i] );
}

void cimpl_cmpEqVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->h == vol1.h );   assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );   assert( out->w == vol2.w );  assert( out->s == vol2.s );
  int i;

  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = ( vol1.data[i] == vol2.data[i] );
}

void cimpl_cmpGeImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->h == img1.h );   assert( out->w == img1.w );
  assert( out->h == img2.h );   assert( out->w == img2.w );
  int i;

  for( i=0; i<img1.nPix; ++i )
    out->data[i] = ( img1.data[i] >= img2.data[i] );
}

void cimpl_cmpGeVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->h == vol1.h );   assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );   assert( out->w == vol2.w );  assert( out->s == vol2.s );
  int i;

  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = ( vol1.data[i] >= vol2.data[i] );
}

void cimpl_cmpGtImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->h == img1.h );   assert( out->w == img1.w );
  assert( out->h == img2.h );   assert( out->w == img2.w );
  int i;

  for( i=0; i<img1.nPix; ++i )
    out->data[i] = ( img1.data[i] > img2.data[i] );
}

void cimpl_cmpGtVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->h == vol1.h );   assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );   assert( out->w == vol2.w );  assert( out->s == vol1.s );
  int i;

  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = ( vol1.data[i] > vol2.data[i] );
}

void cimpl_cmpLeImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->h == img1.h );   assert( out->w == img1.w );
  assert( out->h == img2.h );   assert( out->w == img2.w );
  int i;

  for( i=0; i<img1.nPix; ++i )
    out->data[i] = ( img1.data[i] <= img2.data[i] );
}

void cimpl_cmpLeVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->h == vol1.h );   assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );   assert( out->w == vol2.w );  assert( out->s == vol2.s );
  int i;

  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = ( vol1.data[i] <= vol2.data[i] );
}

void cimpl_cmpLtImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->h == img1.h );   assert( out->w == img1.w );
  assert( out->h == img2.h );   assert( out->w == img2.w );
  int i=0;

  for( i=0; i<img1.nPix; ++i )
    out->data[i] = ( img1.data[i] < img2.data[i] );
}

void cimpl_cmpLtVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_img * const out ){
  assert( out->h == vol1.h );   assert( out->w == vol1.w );
  assert( out->h == vol2.h );   assert( out->w == vol2.w );
  int i;

  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = ( vol1.data[i] < vol2.data[i] );
}

void cimpl_cmpNeqImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->h == img1.h );   assert( out->w == img1.w );
  assert( out->h == img2.h );   assert( out->w == img2.w );
  int i;

  for( i=0; i<img1.nPix; ++i )
    out->data[i] = ( img1.data[i] != img2.data[i] );
}

void cimpl_cmpNeqVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->h == vol1.h );   assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );   assert( out->w == vol2.w );  assert( out->s == vol2.s );
  size_t i;

  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = ( vol1.data[i] != vol2.data[i] );
}

void cimpl_concatCmpImgsH( cimpl_cmpImg const img1, cimpl_cmpImg const img2, cimpl_cmpImg * const out ){
  assert( img1.h == img2.h );  assert( out->h == img1.h );
  assert( out->w == img1.w + img2.w );

  for( size_t x=0; x<out->w; ++x ){
    memcpy( out->data+x*out->w, img1.data+x*img1.h, img1.h*sizeof(CIMPL_COMPLEX) );
    memcpy( out->data+x*out->w, img2.data+x*img2.h, img2.h*sizeof(CIMPL_COMPLEX) );
  }
}

void cimpl_concatCmpImgsW( cimpl_cmpImg const img1, cimpl_cmpImg const img2, cimpl_cmpImg * const out ){
  assert( img1.h == img2.h );  assert( out->h == img1.h );
  assert( out->w == img1.w + img2.w );

  memcpy( out->data, img1.data, img1.nPix*sizeof(CIMPL_COMPLEX) );
  memcpy( out->data, img2.data, img2.nPix*sizeof(CIMPL_COMPLEX) );
}

void cimpl_concatImgsH( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->w == img1.w );  assert( out->w == img2.w );
  assert( out->h == img1.h + img2.h );

  for( size_t x=0; x<out->w; ++x ){
    memcpy( out->data+x*out->w, img1.data+x*img1.h, img1.h*sizeof(float) );
    memcpy( out->data+x*out->w+img1.h, img2.data+x*img2.h, img2.h*sizeof(float) );
  }
}

void cimpl_concatImgsW( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( img1.h == img2.h );  assert( out->h == img1.h );
  assert( out->w == img1.w + img2.w );

  memcpy( out->data, img1.data, img1.nPix*sizeof(float) );
  memcpy( out->data + img1.nPix, img2.data, img2.nPix*sizeof(float) );
}

void cimpl_concatVolsH( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->w == vol2.w );  assert( out->s == vol2.s );
  assert( out->h == vol1.h + vol2.h );
  
  size_t sliceSize1 = vol1.h * vol1.w;
  size_t sliceSize2 = vol2.h * vol2.w;
  size_t sliceSizeOut = out->h * out->w;
  
  for( size_t x=0; x<out->w; ++x ){
    for( size_t s=0; s<out->s; ++s ){
      memcpy( out->data + x*out->w + s*sliceSizeOut, vol1.data + x*vol1.w + s*sliceSize1,
        vol1.h*sizeof(float) );
      memcpy( out->data + x*out->w + s*sliceSizeOut + vol2.h, vol2.data + x*vol2.w + s*sliceSize2,
        vol2.h*sizeof(float) );
  } }
}

void cimpl_concatVolsS( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( vol1.h == vol2.h );  assert( vol1.w == vol2.w );
  assert( out->h == vol1.h );  assert( out->w == vol1.w );
  assert( out->s == vol1.s + vol2.s );
  
  memcpy( out->data, vol1.data, vol1.nVox*sizeof(float) );
  memcpy( out->data+vol1.nVox, vol2.data, vol2.nVox*sizeof(float) );
}

void cimpl_concatCmpVolsW( cimpl_cmpVol const vol1, cimpl_cmpVol const vol2, cimpl_cmpVol * const out ){
  assert( out->h == vol1.h );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );  assert( out->s == vol2.s );
  assert( out->w == vol1.w + vol2.w );

  size_t sliceSize1 = vol1.h*vol1.w;
  size_t sliceSize2 = vol2.h*vol2.w;

  for( size_t s=0; s<out->h*out->w; ++s ){
    memcpy( out->data+s*out->h*out->w, vol1.data+s*sliceSize1, sizeof(CIMPL_COMPLEX)*sliceSize1 );
    memcpy( out->data+s*out->h*out->w + sliceSize1*sizeof(CIMPL_COMPLEX), vol2.data+s*sliceSize2,
      sizeof(CIMPL_COMPLEX)*sliceSize2 );
  }
}

void cimpl_concatCmpVolsH( cimpl_cmpVol const vol1, cimpl_cmpVol const vol2, cimpl_cmpVol * const out ){
  assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->w == vol2.w );  assert( out->s == vol2.s );
  assert( out->h == vol1.h + vol2.h );

  size_t sliceSize1 = vol1.h * vol1.w;
  size_t sliceSize2 = vol2.h * vol2.w;
  size_t sliceSizeOut = out->h * out->w;

  for( size_t x=0; x<out->w; ++x ){
    for( size_t s=0; s<out->s; ++s ){
      memcpy( out->data + x*out->w + s*sliceSizeOut, vol1.data + x*vol1.w + s*sliceSize1,
        vol1.h*sizeof(CIMPL_COMPLEX) );
      memcpy( out->data + x*out->w + s*sliceSizeOut + vol2.h, vol2.data + x*vol2.w + s*sliceSize2,
        vol2.h*sizeof(CIMPL_COMPLEX) );
  } }
}

void cimpl_concatCmpVolsS( cimpl_cmpVol const vol1, cimpl_cmpVol const vol2, cimpl_cmpVol * const out ){
  assert( out->h == vol1.h );  assert( out->w == vol1.w );
  assert( out->h == vol2.h );  assert( out->w == vol2.w );
  assert( out->s == vol1.s + vol2.s );

  size_t sizeVol1 = vol1.h*vol1.s*vol1.w;
  size_t sizeVol2 = vol2.h*vol2.s*vol2.w;

  memcpy( out->data, vol1.data, sizeVol1*sizeof(CIMPL_COMPLEX) );
  memcpy( out->data+sizeVol1, vol2.data, sizeVol2*sizeof(CIMPL_COMPLEX) );
}

void cimpl_concatVolsW( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->h == vol1.h );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );  assert( out->s == vol2.s );
  assert( out->w == vol1.w + vol2.w );

  size_t sliceSize1 = vol1.h*vol1.w;
  size_t sliceSize2 = vol2.h*vol2.w;

  for( size_t s=0; s<out->h*out->w; ++s ){
    memcpy( out->data+s*out->h*out->w, vol1.data+s*sliceSize1, sizeof(float)*sliceSize1 );
    memcpy( out->data+s*out->h*out->w + sliceSize1*sizeof(float), vol2.data+s*sliceSize2,
      sizeof(float)*sliceSize2 );
  }
}

void cimpl_conjCmpImg( cimpl_cmpImg const img, cimpl_cmpImg * const out ){
  assert( out->h = img.h );  assert( out->w = img.w );
  size_t i;

  for( i=0; i<img.nPix; ++i )
    out->data[i] = conjf( img.data[i] );
}

void cimpl_copyCmpImg( cimpl_cmpImg const in, cimpl_cmpImg * const copy ){
  assert( copy->h == in.h );  assert( copy->w == in.w );
  assert( copy->data != in.data );  assert( copy->data != NULL );
  memcpy( copy->data, in.data, in.nPix * sizeof(CIMPL_COMPLEX) );
}

void cimpl_copyImg( cimpl_img const in, cimpl_img * const copy ){
  assert( copy->h == in.h );  assert( copy->w == in.w );
  assert( copy->data != in.data );  assert( copy->data != NULL );
  memcpy( copy->data, in.data, in.nPix * sizeof(float) );
}

void cimpl_cropImg( cimpl_img const in, cimpl_img * const out ){
  // Crops data to size of out
  assert( out->h <= in.h );  assert( out->w <= in.w );
  size_t halfH, halfW, minH, minW;
  size_t colOffset, minColOffset;

  if( in.h % 2 == 0 )
    halfH = in.h/2;
  else
    halfH = (in.h-1)/2;
  if( in.w % 2 == 0 )
    halfW = in.w/2;
  else
    halfW = (in.w-1)/2;

  if( out->h % 2 == 0 )
    minH = halfH - out->h/2;
  else
    minH = halfH - (out->h-1)/2;
  if( out->w % 2 == 0 )
    minW = halfW - out->w/2;
  else
    minW = halfW - (out->w-1)/2;

  for( size_t x=0; x<out->w; ++x ){
    colOffset = x*out->h;
    minColOffset = (minW+x)*in.h;

    memcpy( out->data+colOffset, in.data+minH+minColOffset, out->h*sizeof(float) );
  }
}

void cimpl_cropVol( cimpl_vol const in, cimpl_vol * const out ){
  // Crops data to size of out
  assert( out->w <= in.w );  assert( out->h <= in.h );  assert( out->s <= in.s );
  assert( out->data != in.data );
  size_t halfH, halfS, halfW, minH, minS, minW;
  size_t minColOffset, minSliceOffset;
  size_t colOffset, sliceOffset;

  if( in.h % 2 == 0 )
    halfH = in.h/2;
  else
    halfH = (in.h-1)/2;
  if( in.w % 2 == 0 )
    halfW = in.w/2;
  else
    halfW = (in.w-1)/2;
  if( in.s % 2 == 0 )
    halfS = in.s/2;
  else
    halfS = (in.s-1)/2;

  if( out->h % 2 == 0 )
    minH = halfH - out->h/2;
  else
    minH = halfH - (out->h-1)/2;
  if( out->w % 2 == 0 )
    minW = halfW - out->w/2;
  else
    minW = halfW - (out->w-1)/2;
  if( out->s % 2 == 0 )
    minS = halfS - out->s/2;
  else
    minS = halfS - (out->s-1)/2;

  for( size_t z=0; z<out->s; ++z ){
    sliceOffset = z*out->h*out->w;
    minSliceOffset = (minS+z)*in.h*in.w;

    for( size_t x=0; x<out->w; ++x ){
      colOffset = x*out->h;
      minColOffset = (minW+x)*in.h;

      memcpy(out->data+colOffset+sliceOffset, in.data+minH+minColOffset+minSliceOffset,
        out->h*sizeof(float) );
  } }
}

void cimpl_d1hImg( cimpl_img const in, cimpl_img * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  size_t i, j;

  for( i=0; i<in.nPix-1; ++i )
    out->data[i] = out->data[i+1] - out->data[i];

  assert(1);  assert(0);  // NICK, THERE'S A BUG HERE; YOU REALLY NEED TO DECIDE WHICH DIMENSIONS IS WHICH
  for( j=0; j<in.w; ++j )
    out->data[in.w-1 + j*in.w] = 0;
}

void cimpl_d1hVol( cimpl_vol const in, cimpl_vol * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );  assert( out->s == in.s );
  assert(1);  assert(0);
  //NICK, THERE IS A BUG IN THIS FUNCTION.  WHAT'S sqrt DOING HERE

  float *next, *current;
  current = in.data;
  next = in.data + 1;
  for( size_t z=0; z<in.s; ++z ){
    for( size_t x=0; x<in.w; ++x, ++current, ++next ){
      for( size_t y=0; y<in.h-1; ++y, ++current, ++next )
        out->data[y+x*in.h] = *next - *current;

      out->data[(in.h-1)+x*in.h] = out->data[(in.h-2)+x*in.h];
  } }
}

void cimpl_d1wImg( cimpl_img const in, cimpl_img * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );

  float *next, *current;
  current = in.data;
  next = in.data + in.h;
  for( size_t x=0; x<in.w-1; ++x, ++current, ++next ){
    for( size_t y=0; y<in.h; ++y, ++current, ++next )
      out->data[y+x*in.h] = *next - *current;
  }
  
  memcpy( out->data+(in.w-1)*in.h, out->data+(in.w-2)*in.h, in.h*sizeof(float) );
}

void cimpl_d1wVol( cimpl_vol const in, cimpl_vol * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );  assert( out->s == in.s );
  float *next, *current;

  for( size_t z=0; z<in.s; ++z ){
    current = in.data + z*in.h*in.w;
    next = current + in.h;

    for( size_t x=0; x<in.w-1; ++x, ++current, ++next ){
      for( size_t y=0; y<in.h; ++y, ++current, ++next )
        out->data[y+x*in.h] = *next - *current;
    }

    memcpy( out->data+(in.w-1)*in.h+z*in.h*in.s, out->data+(in.w-2)*in.h+z*in.h*in.s,
      in.h*sizeof(float) );
  }
}

void cimpl_d1sVol( cimpl_vol const in, cimpl_vol * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );  assert( out->s == in.s );
  float *next, *current;

  for( size_t z=0; z<in.s-1; ++z ){
    current = in.data + z*in.h*in.w;
    next = current + in.h*in.w;

    for( size_t x=0; x<in.w; ++x, ++current, ++next ){
      for( size_t y=0; y<in.h; ++y, ++current, ++next )
        out->data[y+x*in.h+z*in.h*in.w] = *next - *current;
  } }

  memcpy( out->data+(in.s-1)*in.h*in.w, in.data+(in.s-1)*in.h*in.w, in.h*in.w*sizeof(float) );
}

float cimpl_dotProd( float const * const x, float const * const y, size_t const N ){
  // N is the number of elements in the arrays
  size_t i;
  float out=0;

  for( i=0; i < N; ++i ){
    out += x[i] * y[i];
  }

  return out;
}

void cimpl_divideCmpImgs( cimpl_cmpImg const img1, cimpl_cmpImg const img2,
  cimpl_cmpImg * const out ){
  assert( out->w == img1.w );  assert( out->h == img1.h );
  assert( out->w == img2.w );  assert( out->h == img2.h );
  size_t i;

  for (i=0; i < out->nPix; ++i)
    out->data[i] = cimpl_divideCmp( img1.data[i], img2.data[i] );
}

void cimpl_divideCmpVols( cimpl_cmpVol const vol1, cimpl_cmpVol const vol2,
  cimpl_cmpVol * const out ){
  assert( out->w == vol1.w );  assert( out->h == vol1.h );  assert( out->s == vol1.s );
  assert( out->w == vol2.w );  assert( out->h == vol2.h );  assert( out->s == vol2.s );
  size_t i;

  for (i=0; i < out->nVox; ++i)
    out->data[i] = cimpl_divideCmp( vol1.data[i], vol2.data[i] );
}

void cimpl_divideImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->w == img1.w );  assert( out->h == img1.h );
  assert( out->w == img2.w );  assert( out->h == img2.h );
  int i;

  for( i=0; i<img1.nPix; ++i )
    out->data[i] = img1.data[i] / img2.data[i];
}

void cimpl_divideVols( cimpl_vol const vol1, cimpl_vol const vol2,
  cimpl_vol * const out ){
  assert( out->h == vol1.h );  assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );  assert( out->w == vol2.w );  assert( out->s == vol2.s );
  int i;

  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = vol1.data[i] / vol2.data[i];
}

void cimpl_divideImgByScalar( cimpl_img const in, float const scalar,
  cimpl_img * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  int i;

  for( i=0; i<in.nPix; ++i )
    out->data[i] = in.data[i] / scalar;
}

void cimpl_divideVolByScalar( cimpl_vol const in, float const scalar,
  cimpl_vol * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  int i;

  for( i=0; i<in.nVox; ++i )
    out->data[i] = in.data[i] / scalar;
}

float cimpl_dotImgs( cimpl_img const img1, cimpl_img const img2 ){
  assert( img1.w == img2.w );  assert( img1.h == img2.h );
  size_t i;
  float out=0;
  for( i=0; i<img1.h*img1.w; ++i )
    out += img1.data[i] * img2.data[i];

  return out;
}

void cimpl_downsampleCmpImg( cimpl_cmpImg const img1, int downsample,
  cimpl_cmpImg * const out ){
  assert( floorf( (float)img1.w/downsample ) == out->w );
  assert( floorf( (float)img1.h/downsample ) == out->h );
  size_t i, j;
  for( i=0; i<out->w; ++i ){
    for( j=0; j<out->h; ++j ){
      out->data[j+i*out->h] = img1.data[(j*downsample)+(i*downsample*img1.h)];
  } }
}

void cimpl_downsampleImg( cimpl_img const img1, int downsample, cimpl_img * const out ){
  assert( floorf( (float)img1.w/downsample ) == out->w );
  assert( floorf( (float)img1.h/downsample ) == out->h );
  size_t i, j;
  for( i=0; i<out->w; ++i ){
    for( j=0; j<out->h; ++j ){
      out->data[j+i*out->h] = img1.data[(j*downsample)+(i*downsample*img1.h)];
  } }
}

void cimpl_downsampleVol( cimpl_vol const vol1, int downsample, cimpl_vol * const out ){
  assert( floorf( (float)vol1.w/downsample ) == out->w );
  assert( floorf( (float)vol1.h/downsample ) == out->h );
  assert( floorf( (float)vol1.s/downsample ) == out->s );
  size_t u, i, j;
  for( u=0; u<out->s; ++u ){
    for( i=0; i<out->w; ++i ){
      for( j=0; j<out->h; ++j ){
        out->data[j+i*out->h+u*(out->h*out->w)] = 
          vol1.data[(j*downsample)+(i*downsample*vol1.h)+(u*downsample*vol1.h*vol1.w)];
  } } }
}

int cimpl_equalCmpImgs( cimpl_cmpImg const img1, cimpl_cmpImg const img2 ){
  size_t i;
  if( img1.h != img2.h )
    return 0;
  if( img1.w != img2.w )
    return 0;
  for( i=0; i<img1.nPix; ++i){
    if( crealf(img1.data[i]) != crealf(img2.data[i]) || cimagf(img1.data[i]) != cimagf(img2.data[i]) ){
      return 0;
  } }
  return 1;
}

int cimpl_equalCmpVols( cimpl_cmpVol const vol1, cimpl_cmpVol const vol2 ){
  size_t i;
  if( vol1.h != vol2.h ) return 0;
  if( vol1.w != vol2.w ) return 0;
  if( vol1.s != vol2.s ) return 0;
  for( i=0; i<vol1.nVox; ++i){
    if( crealf(vol1.data[i]) != crealf(vol2.data[i]) ||
        cimagf(vol1.data[i]) != cimagf(vol2.data[i]) )
      return 0;
  }
  return 1;
}

int cimpl_equalImgs( cimpl_img const img1, cimpl_img const img2 ){
  size_t i;
  if( img1.h != img2.h )
    return 0;
  if( img1.w != img2.w )
    return 0;
  for( i=0; i<img1.nPix; ++i){
    if( img1.data[i] != img2.data[i] ) return 0;
  }
  return 1;
}

int cimpl_equalVols( cimpl_vol const vol1, cimpl_vol const vol2 ){
  size_t i;
  if( vol1.h != vol2.h ) return 0;
  if( vol1.w != vol2.w ) return 0;
  if( vol1.s != vol2.s ) return 0;
  for( i=0; i<vol1.nVox; ++i){
    if( vol1.data[i] != vol2.data[i] ) return 0;
  }
  return 1;
}

void cimpl_flipImgLR( cimpl_img const in, cimpl_img * const out ){
  // Flip image left/right
  size_t x, y;
  assert( out->h == in.h );  assert( out->w == in.w );
  assert( out->data != in.data );
  for( x=0; x<in.w; ++x ){
    for( y=0; y<in.h; ++y ){
      out->data[y+x*in.h] = in.data[y+(in.w-x-1)*in.h];
  } }
}

void cimpl_flipImgUD( cimpl_img const in, cimpl_img * const out ){
  // Flip image up/down
  assert( out->h == in.h );  assert( out->w == in.w );
  assert( out->data != in.data );
  for( size_t x=0; x<in.w; ++x ){
    for( size_t y=0; y<in.h; ++y ){
      out->data[y+x*in.h] = in.data[(in.h-y-1)+x*in.h];
  } }
}

void cimpl_floorImg( cimpl_img const img, cimpl_img * const out ){
  assert( out->h == img.h );  assert( out->w == img.w );
  size_t i;
  for( i=0; i<img.nPix; ++i )
    out->data[i] = floorf( img.data[i] );
}

void cimpl_floorVol( cimpl_vol const vol, cimpl_vol * const out ){
  assert( out->h == vol.h );  assert( out->w == vol.w );  assert( out->s == vol.s );
  size_t i;
  for( i=0; i<vol.nVox; ++i )
    out->data[i] = floorf( vol.data[i] );
}

void cimpl_free( void * in ){
#if defined( _MSC_VER ) || defined( __MINGW32__ )
  _aligned_free( in );
#else
  free( in );
#endif
}

void cimpl_freeCmpImg( cimpl_cmpImg * in ){
  in->h = in->w = 0;
  cimpl_free( in->data );
  in->data = NULL;
}

void cimpl_freeCmpVol( cimpl_cmpVol * in ){
  in->h = in->w = in->s = 0;
  cimpl_free( in->data );
  in->data = NULL;
}

void cimpl_freeImg( cimpl_img *in ){
  in->h = in->w = 0;
  cimpl_free( in->data );
  in->data = NULL;
}

void cimpl_freeVol( cimpl_vol * in ){
  in->h = in->w = in->s = 0;
  cimpl_free( in->data );
  in->data = NULL;
}

void cimpl_imagImg( cimpl_cmpImg const in, cimpl_img * const out ){
  assert( out->h == in.h ); assert( out->w == in.w );
  CIMPL_COMPLEX* inPtr = in.data;
  float* outPtr = out->data;
  for( size_t i=0; i<out->nPix; ++i, ++inPtr, ++outPtr )
    *outPtr = cimpl_imag( *inPtr );
}

void cimpl_imagVol( cimpl_cmpVol const in, cimpl_vol * const out ){
  assert( out->h == in.h ); assert( out->w == in.w );
  CIMPL_COMPLEX* inPtr = in.data;
  float* outPtr = out->data;
  for( size_t i=0; i<out->nVox; ++i, ++inPtr, ++outPtr )
    *outPtr = cimpl_imag( *inPtr );
}

float cimpl_linInterp( size_t const N, float const * const x, float const * const y,
  float const outOfBounds, float const q ){
  // N is the number of values in the x and y arrays
  // x is an array of inputs (or domain values) in ascending order
  // y is an array of outputs (or range values) where y[i] corresponds to x[i]
  // outOfBounds is the value to return if extrapolating
  // q is the query domain value.

  float out=outOfBounds;
  size_t lowIndx=0;
  size_t highIndx=N-1;
  size_t midIndx;

  if( q < x[0] || q > x[N-1] ) return outOfBounds;

  // Perform a binary search
  while( highIndx - lowIndx > 1 ){
    midIndx = lowIndx + (highIndx-lowIndx)/2;
    if( q < x[midIndx] ){
      highIndx = midIndx;
    } else {
      lowIndx = midIndx;
    }
  }
  out = y[lowIndx] + (y[highIndx]-y[lowIndx])/(x[highIndx]-x[lowIndx]) * (q-x[lowIndx]);

  return out;
}

void cimpl_linInterps( size_t const N, float const * const x, float const * const y,
                      float const outOfBounds, size_t const M, float const * const q,
                      float * const out ){
  // This function takes advantage of order of q for greater efficiency than calling
  //   cimpl_linInterp multiple times.
  // N is the number of values in the x and y arrays
  // x is an array of inputs (or domain values) in ascending order
  // y is an array of outputs (or range values) where y[i] corresponds to x[i]
  // outOfBounds is the value to return if extrapolating
  // M is the number of query values
  // q are the query domain values in ascending order.

  size_t qIndx=0;
  size_t xIndx;

  while( q[qIndx] < x[0] ){
    out[qIndx] = outOfBounds;
    qIndx++;
    if( qIndx >= M )
      return;
  }

  for( xIndx=0; xIndx<N-1; ++xIndx ){
    while( q[qIndx] > x[xIndx] && q[qIndx] < x[xIndx+1] ){
      out[qIndx] = (y[xIndx+1]-y[xIndx])/(x[xIndx+1]-x[xIndx]) * (q[qIndx]-x[xIndx]);
      ++qIndx;
  } }

  while( qIndx < M ){
    out[qIndx] = outOfBounds;
    ++qIndx;
  }

  return;
}

void cimpl_linInterpImg( cimpl_img const img, size_t const N, float const * const xq,
  float const * const yq, float const outOfBounds, float * const out ){
  size_t x1, x2, y1, y2;
  float tmpX1, tmpX2;
  for( size_t i=0; i<N; ++i ){
    if( xq[i]<0 || xq[i]>=img.w || yq[i]<0 || yq[i]>=img.h ){
      out[i] = outOfBounds;
    } else {
      x1 = (size_t) floorf( xq[i] );  x2 = (size_t) ceilf( xq[i] );
      y1 = (size_t) floorf( yq[i] );  y2 = (size_t) ceilf( yq[i] );
      tmpX1 = (y2-yq[i])*img.data[y1+x1*img.h] + (yq[i]-y1)*img.data[y2+x1*img.h];
      tmpX2 = (y2-yq[i])*img.data[y1+x2*img.h] + (yq[i]-y1)*img.data[y2+x2*img.h];
      out[i] = (x2-xq[i])*tmpX1 + (xq[i]-x1)*tmpX2;
  } }
}

void cimpl_malloc( size_t nBytes, void ** ptrPtr ){
  int boundary = 16; // word alignment for improved efficiency

#if defined( _MSC_VER ) || defined( __MINGW32__ )
  ptrPtr = _aligned_malloc( nBytes, boundary );
#else
  posix_memalign( ptrPtr, boundary, nBytes );
#endif
  // *ptrPtr = malloc( nBytes );  // unaligned allocation
}

cimpl_cmpImg cimpl_mallocCmpImg( size_t const h, size_t const w ){
  cimpl_cmpImg out;
  out.h = h;  out.w = w;  out.nPix = h*w;
  cimpl_malloc( sizeof(CIMPL_COMPLEX) * out.nPix, (void**) &(out.data) );
  return out;
}

cimpl_cmpVol cimpl_mallocCmpVol( size_t const h, size_t const w , size_t const s ){
  cimpl_cmpVol out;
  out.h = h;  out.w = w;  out.s = s;  out.nVox = h*w*s;
  cimpl_malloc( sizeof(CIMPL_COMPLEX) * out.nVox, (void**) &(out.data) );
  return out;
}

cimpl_img cimpl_mallocImg( size_t const h, size_t const w ){
  cimpl_img out;
  out.h = h;  out.w = w;  out.nPix = h*w;
  cimpl_malloc( sizeof(float) * out.nPix, (void**) &(out.data) );
  return out;
}

cimpl_pts2 cimpl_mallocPts2( size_t const nPts ){
  cimpl_pts2 out;
  out.nPts = nPts;
  cimpl_malloc( sizeof(float) * nPts*2, (void**) &(out.pts) );
  return out;
}

cimpl_pts3 cimpl_mallocPts3( size_t const nPts ){
  cimpl_pts3 out;
  out.nPts = nPts;
  cimpl_malloc( sizeof(float) * nPts*3, (void**) &(out.pts) );
  return out;
}

cimpl_ptsN cimpl_mallocPtsN( size_t const nPts, const unsigned int N ){
  cimpl_ptsN out;
  out.nPts = nPts;
  cimpl_malloc( sizeof(float) * nPts*N, (void**) &(out.pts) );
  return out;
}

cimpl_vol cimpl_mallocVol( size_t const h, size_t const w, size_t const s ){
  cimpl_vol out;
  out.h = h;  out.w = w;  out.s = s;  out.nVox = h*w*s;
  cimpl_malloc( sizeof(float) * out.nVox, (void**) &(out.data) );
  return out;
}

void cimpl_maxImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->h == img1.h );  assert( out->w == img1.w );
  assert( out->h == img2.h );  assert( out->w == img2.w );
  int i;
  for( i=0; i<out->nPix; ++i )
    out->data[i] = fmaxf( img1.data[i], img2.data[i] );
}

void cimpl_maxVols( cimpl_vol vol1, cimpl_vol vol2, cimpl_vol * out ){
  assert( out->h == vol1.h );  assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );  assert( out->w == vol2.w );  assert( out->s == vol2.s );
  int i;
  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = fmaxf( vol1.data[i], vol2.data[i] );
}

void cimpl_meanCmpImg( cimpl_cmpImg img, CIMPL_COMPLEX * out ){
  size_t i;

  *out = 0;
  for( i=0; i < img.nPix; ++i ){
    *out = cimpl_addCmp( *out, img.data[i] );
  }
  cimpl_divideCmp( *out, (CIMPL_COMPLEX) img.nPix );
}

void cimpl_meanCmpVol( cimpl_cmpImg vol, CIMPL_COMPLEX * out ){
  size_t i;

  *out = 0;
  for( i=0; i < vol.nPix; ++i ){
    *out = cimpl_addCmp( *out, vol.data[i] );
  }
  cimpl_divideCmp( *out, (CIMPL_COMPLEX) vol.nPix );
}

void cimpl_meanImg( cimpl_img img, float *mean ){
  size_t i;

  *mean = 0;
  for( i=0; i < img.nPix; ++i ){
    *mean += (img.data)[i];
  }
  *mean /= (float) img.nPix;
}

void cimpl_meanVol( cimpl_vol vol, float *mean ){
  size_t i;

  *mean = 0;
  for( i=0; i < vol.nVox; ++i ){
    *mean += (vol.data)[i];
  }
  *mean /= (float) vol.nVox;
}

void cimpl_minImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->h == img1.h );  assert( out->w == img1.w );
  assert( out->h == img2.h );  assert( out->w == img2.w );
  int i;
  for( i=0; i<out->nPix; ++i )
    out->data[i] = fminf( img1.data[i], img2.data[i] );
}

void cimpl_minVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->h == vol1.h );  assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );  assert( out->w == vol2.w );  assert( out->s == vol2.s );
  int i;
  for( i=0; i<out->nVox; ++i )
    out->data[i] = fminf( vol1.data[i], vol2.data[i] );
}

void cimpl_multiplyCmpImgByRealScalar( cimpl_cmpImg const in, float const scalar,
  cimpl_cmpImg * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  int i;
  float *inData, *outData;
  inData = (float*) in.data;
  outData = (float*) out->data;
  for( i=0; i<in.nPix*2; ++i, ++inData, ++outData )
    *outData = *inData * scalar;
}

void cimpl_multiplyCmpImgs( cimpl_cmpImg const img1, cimpl_cmpImg const img2, cimpl_cmpImg * const out ){
  assert( out->w == img1.w );  assert( out->h == img1.h );
  assert( out->w == img2.w );  assert( out->h == img2.h );

  for (size_t i = 0; i < out->nPix; ++i)
    out->data[i] = cimpl_multiplyCmp( img1.data[i], img2.data[i] );
}

void cimpl_multiplyImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  assert( out->w == img1.w );  assert( out->h == img1.h );
  assert( out->w == img2.w );  assert( out->h == img2.h );
  int i;
  for( i=0; i<img1.nPix; ++i )
    out->data[i] = img1.data[i] * img2.data[i];
}

void cimpl_multiplyVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  assert( out->h == vol1.h );  assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );  assert( out->w == vol2.w );  assert( out->s == vol2.s );
  int i;
  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = vol1.data[i] * vol2.data[i];
}

void cimpl_multiplyImgByScalar( cimpl_img const in, float const scalar,
  cimpl_img * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  int i;
  for( i=0; i<in.nPix; ++i )
    out->data[i] = in.data[i] * scalar;
}

void cimpl_multiplyVolByScalar( cimpl_vol const in, float const scalar,
  cimpl_vol * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );  assert( out->s == in.s );
  int i;
  for( i=0; i<in.nVox; ++i )
    out->data[i] = in.data[i] * scalar;
}

float cimpl_normCmpImgL2( cimpl_cmpImg const in ){
  size_t i;
  float out;

  out = 0;
  for( i=0; i<in.nPix; ++i ){
    out += cimpl_multiplyCmp( in.data[i], conjf( in.data[i] ) );
  }
  out = sqrtf( out );
  return out;
}

float cimpl_normImgL1( cimpl_img const in ){
  size_t i;
  float out;

  out=0;
  for( i=0; i<in.nPix; ++i ){
    out += fabs( in.data[i] );
  }
  return out;
}

float cimpl_normImgL2( cimpl_img const in ){
  size_t i;
  float out;

  out=0;
  for( i=0; i<in.nPix; ++i ){
    out += in.data[i] * in.data[i];
  }
  out = sqrtf( out );
  return out;
}


cimpl_pts2 cimpl_poissonDisc2( float const r, float const * const bounds ){
  // Written according to "Fast Poisson Disk Sampling in Arbitrary Dimensions" by Robert Bridson
  // and "Poisson Disk Sampling" written by Herman Tulleken located at
  // http://devmag.org.za/2009/05/03/poisson-disk-sampling/

  // Attempts to create the number of points in the pts array
  // Gets cutoff if the number of iterations to do so exceeds maxIter
  // We leave it to the user to choose a method and appropriate action for compensating for this
  // latter case.  For example, the user may intialize all points with values that are out of bounds.
  // This could then be checked after processing.  Alternatively, the user could initialize the points
  // array to random values in the space and retain these values when the maximum number of
  // iterations is reached.

  // Inputs:
  // r - the poisson Disc parameter
  // bounds - a four element array specifying [ xLower yLower xUpper yUpper ]
  //   If bounds is NULL, then the bounds are [ -0.5 -0.5 0.5 0.5 ]
  // maxIter - the maximum number of iterations to run
  //
  // Outputs:
  // pd - the realized poisson disc points

  bool validPtFound;
  cimpl_pts2 pd;
  float xL, yL, xU, yU, distX, distY, dist2NearbyPt;
  float cellSize, kDist, kAngle;
  float pt[2], kPt[2], otherPt[2];
  int nNearCells;
  size_t kPtGridX, kPtGridY, gcX, gcY, gxIndx, gyIndx, nPdPts;
  size_t i, pListIndx, opIndx, ptIndx, nGridX, nGridY, nPList, pListSize;
  size_t *pList, *tmpList;
  int *bGrid;
  const size_t incSize = 128;
  const size_t nKPts = 10;

  if( bounds == NULL ){
    xL = yL = -0.5;  xU = yU = 0.5;
  } else {
    xL = bounds[0];  xU = bounds[2];
    yL = bounds[1];  yU = bounds[3];
  }

  cellSize = r / sqrt(2);
  distX = xU - xL;  nGridX = ceilf( distX / cellSize );
  distY = yU - yL;  nGridY = ceilf( distY / cellSize );

  cimpl_malloc( nGridX * nGridY * sizeof(int), (void**) &bGrid );  // Background grid
  for( i=0; i < nGridX * nGridY; ++i ) bGrid[i] = -1;

  pListSize = sizeof(size_t) * incSize;
  cimpl_malloc( pListSize, (void**) &pList );

  pd = cimpl_mallocPts2( incSize );

  // Create the first point and initialize the processing list
  nPdPts = 0;
  pt[0] = cimpl_randuf( xL, xU );
  pt[1] = cimpl_randuf( yL, yU );
  cimpl_pts2_assignPt( &pd, nPdPts, pt );
  ++nPdPts;
  pList[0] = 0;
  nPList = 1;
  // Set the corresponding element of the background grid to 1
  gcX = floorf( ( pt[0] - xL ) / cellSize );
  gcY = floorf( ( pt[1] - yL ) / cellSize );
  bGrid[ gcX + gcY * nGridX ] = 1;

  while( nPList > 0 ){
    // Choose a point at random in the processing list
    pListIndx = cimpl_randui( 0, (int) nPList-1 );
    ptIndx = pList[ pListIndx ];
    cimpl_pts2_getPt( &pd, ptIndx, pt );

    //Make k random points in the annulus between r and 2r from the active point
    validPtFound = 0;
    for( i=0; i<nKPts; ++i ){
      kDist = cimpl_randuf( r, 2*r );
      kAngle = cimpl_randuf( 0, 2 * CIMPL_PI );
      kPt[0] = pt[0] + kDist * cos( kAngle );
      kPt[1] = pt[1] + kDist * sin( kAngle );

      // Check to see if the candidate point is valid
      if( kPt[0] < xL || kPt[0] > xU || kPt[1] < yL || kPt[1] > yU ) continue;

      // Find the minimum distance to other points
      kPtGridX = floorf( ( kPt[0] - xL ) / cellSize );
      kPtGridY = floorf( ( kPt[1] - yL ) / cellSize );
      nNearCells = ceil( r / cellSize );
      for( gyIndx = cimpl_max( kPtGridY - nNearCells, 0 );
           gyIndx < cimpl_min( kPtGridY + nNearCells, nGridY-1); ++gyIndx ){
        for( gxIndx = cimpl_max( kPtGridX - nNearCells, 0 );
             gxIndx < cimpl_min( kPtGridX + nNearCells, nGridX-1 ); ++gxIndx ){
          if( bGrid[ gxIndx + gyIndx * nGridX ] < 0 ) continue;
          if( gxIndx == kPtGridX && gyIndx == kPtGridY ) continue;

          // Find the distance from kPt to the point in the background grid
          opIndx = pList[ bGrid[ gxIndx + gyIndx * nGridX ] ];
          cimpl_pts2_getPt( &pd, opIndx, otherPt );  // Other point
          dist2NearbyPt = sqrtf( (kPt[0]-otherPt[0])*(kPt[0]-otherPt[0]) +
                                 (kPt[1]-otherPt[1])*(kPt[1]-otherPt[1]) );
          if( dist2NearbyPt < r ) continue;

          validPtFound = 1;

          // Add this point to the list of points
          if( nPdPts == pd.nPts ) cimpl_pts2_grow( &pd, pd.nPts + incSize );
          cimpl_pts2_assignPt( &pd, nPdPts, kPt );

          // Add this point to the processing list
          if( nPList == pListSize ){
            // Increase the size of the processing list
            tmpList = pList;
            pListSize += incSize;
            cimpl_malloc( sizeof(size_t) * pListSize, (void**) &pList );
            memcpy( pList, tmpList, sizeof(size_t) * (pListSize-incSize) );
            cimpl_free( tmpList );
          }
          pList[ nPList ] = nPdPts;
          ++nPList;

          // Add this point to the background list
          bGrid[ kPtGridX + kPtGridY * nGridX ] = (int) nPdPts;
          ++nPdPts;
      } }
    }

    if( validPtFound == 0 ){
      // Remove this point from the processing list
      pList[ ptIndx ] = pList[ nPList - 1 ];
      --nPList;
    }
  }

  cimpl_free( bGrid );
  cimpl_free( pList );
  return pd;
}

void cimpl_poissonDisc2_dwork( float const min_r, float (*r)(float const * ), float const * const bounds,
  cimpl_pts2 * pd){
  // Poisson Disc 2 with varying distance parameter

  // Attempts to create the number of points in the pts array
  // Gets cutoff if the number of iterations to do so exceeds maxIter
  // We leave it to the user to choose a method and appropriate action for compensating for this
  // latter case.  For example, the user may intialize all points with values that are out of bounds.
  // This could then be checked after processing.  Alternatively, the user could initialize the points
  // array to random values in the space and retain these values when the maximum number of
  // iterations is reached.
  
  // Inputs:
  // max_r - the maximum value of the poisson Disc parameter observable within bounds
  // bounds - a four element array specifying [ xLower yLower xUpper yUpper ]
  //   If bounds is NULL, then the bounds are [ -0.5 -0.5 0.5 0.5 ]
  // maxIter - the maximum number of iterations to run
  //
  // Outputs:
  // pd - the realized poisson disc points

  bool validPtFound, pointTooClose;
  float xL, yL, xU, yU, distX, distY, dist2NearbyPt;
  float cellSize, kDist, kAngle, rActive, kr;
  float pt[2], activePt[2], bPt[2], kPt[2];
  float *tmpFloats;
  size_t kIndx, gc, kc, gxIndx, gyIndx, nPdPts;
  int gcX, gcY, nNearCells;
  size_t nGridX, nGridY, kcX, kcY, incSize, nKPts;
  size_t pListIndx, ptIndx, nPList, pListSize;
  size_t bIndx, bPtIndx;
  int cL, cR, cB, cT;  // Left, Right, Bottom, Top indices
  size_t *pList, *pTmp, *bTmp, *nBGrid, *bGrid, szBGrid;
  float *rPds = NULL;
  incSize = 256;
  nKPts = 10;

  if( bounds == NULL ){
    xL = yL = -0.5;  xU = yU = 0.5;
  } else {
    xL = bounds[0];  xU = bounds[2];
    yL = bounds[1];  yU = bounds[3];
  }

  cellSize = min_r / sqrt(2);
  distX = xU - xL;  nGridX = ceilf( distX / cellSize );
  distY = yU - yL;  nGridY = ceilf( distY / cellSize );

  szBGrid = incSize;
  cimpl_malloc( nGridX * nGridY * szBGrid * sizeof(size_t), (void**) &bGrid );
  cimpl_malloc( nGridX * nGridY * sizeof(size_t), (void**) &nBGrid );
  memset( nBGrid, 0, nGridX * nGridY );
  cimpl_malloc( incSize * sizeof(size_t), (void**) &pList );

  pListSize = incSize;
  *pd = cimpl_mallocPts2( 100 * incSize );  // The poisson disc samples
  cimpl_malloc( sizeof(float) * 100 * incSize, (void**) &rPds );

  // Create the first point and initialize the processing list
  nPdPts = 0;
  pt[0] = cimpl_randuf( xL, xU );
  pt[1] = cimpl_randuf( yL, yU );
  cimpl_pts2_assignPt( pd, nPdPts, pt );
  ++nPdPts;
  rPds[0] = r( pt );
  pList[0] = 0;
  nPList = 1;

  // Include this point in the elements of the backgroung grid that are within threshold distance
  gcX = floorf( ( pt[0] - xL ) / cellSize );
  gcY = floorf( ( pt[1] - yL ) / cellSize );
  nNearCells = ceil( rPds[0] / cellSize );
  cL = gcX - nNearCells;  cB = gcY - nNearCells;
  cR = gcX + nNearCells;  cT = gcY + nNearCells;
  for( gyIndx = cimpl_max( cB, 0 ); gyIndx < cimpl_min( cT, nGridY ); ++gyIndx ){
    for( gxIndx = cimpl_max( cL, 0 ); gxIndx < cimpl_min( cR, nGridX ); ++gxIndx ){
      if( ( gxIndx == cL || gxIndx == cR ) &&
          ( gyIndx == cB || gyIndx == cT ) ) continue;  // corner

      gc = gxIndx + gyIndx * nGridX;
      bGrid[ gc ] = 0;  // Add 0th point to background grid
      nBGrid[ gc ] = 1;
  } }
  
  while( nPList > 0 ){
    // Choose a point at random in the processing list
    pListIndx = cimpl_randui( 0, (int) nPList - 1 );
    ptIndx = pList[ pListIndx ];
    cimpl_pts2_getPt( pd, ptIndx, activePt );
    rActive = rPds[ ptIndx ];

    //Make k random points in the annulus between r and 2r from the active point
    validPtFound = 0;
    for( kIndx=0; kIndx<nKPts; ++kIndx ){
      kDist = cimpl_randuf( rActive, 2 * rActive );
      kAngle = cimpl_randuf( 0, 2 * CIMPL_PI );
      kPt[0] = activePt[0] + kDist * cos( kAngle );
      kPt[1] = activePt[1] + kDist * sin( kAngle );
      
      // Check to see if the candidate point is valid
      if( kPt[0] < xL || kPt[0] > xU || kPt[1] < yL || kPt[1] > yU ) continue;
      kr = r( kPt );
      kcX = floorf( ( kPt[0] - xL ) / cellSize );
      kcY = floorf( ( kPt[1] - yL ) / cellSize );
      kc = kcX + kcY * nGridX;
      pointTooClose = 0;
      for( bIndx = 0; bIndx < nBGrid[ kc ]; ++bIndx ){
        bPtIndx = bGrid[ kc + bIndx * nGridX * nGridY ];
        cimpl_pts2_getPt( pd, bPtIndx, bPt );
        dist2NearbyPt = sqrtf( ( bPt[0] - kPt[0] ) * ( bPt[0] - kPt[0] ) +
                               ( bPt[1] - kPt[1] ) * ( bPt[1] - kPt[1] ) );
        if( dist2NearbyPt < kr ){
          pointTooClose = 1;
          break;
      } }
      if( pointTooClose == 1 ) continue;

      validPtFound = 1;
      
      // Add this point to the list of points
      // Add the corresponding distance threshold to rPds
      if( nPdPts == pd->nPts ){
        tmpFloats = rPds;
        cimpl_malloc( sizeof(float) * ( pd->nPts + 10 * incSize ), (void**) &rPds );
        memcpy( rPds, tmpFloats, sizeof(float) * pd->nPts );
        cimpl_free( tmpFloats );
        cimpl_pts2_grow( pd, pd->nPts + 10 * incSize );
      }
      rPds[ nPdPts ] = kr;
      cimpl_pts2_assignPt( pd, nPdPts, kPt );

      // Add this point to the processing list
      if( nPList == pListSize ){
        // Increase the size of the processing list
        pTmp = pList;
        pListSize += incSize;
        cimpl_malloc( sizeof(size_t) * pListSize, (void**) &pList );
        memcpy( pList, pTmp, sizeof(size_t) * (pListSize-incSize) );
        cimpl_free( pTmp );
      }
      pList[ nPList ] = nPdPts;
      ++nPList;

      // Add this point into the background grid
      nNearCells = ceil( rActive / cellSize );
      cL = (int) kcX - (int) nNearCells;  cB = (int) kcY - (int) nNearCells;
      cR = (int) kcX + (int) nNearCells;  cT = (int) kcY + (int) nNearCells;
      for( gyIndx = cimpl_max( cB, 0 ); gyIndx <= cimpl_min( cT, nGridY-1 ); ++gyIndx ){
        for( gxIndx = cimpl_max( cL, 0 ); gxIndx <= cimpl_min( cR, nGridX-1 ); ++gxIndx ){
          if( ( gxIndx == cL || gxIndx == cR ) && ( gyIndx == cB || gyIndx == cT ) ) continue;  // corner

          gc = gxIndx + gyIndx * nGridX;
          if( nBGrid[ gc ] == szBGrid ){
            bTmp = bGrid;
            cimpl_malloc( sizeof( size_t ) * nGridX * nGridY * ( szBGrid + incSize ), (void**) &bGrid );
            memcpy( bGrid, bTmp, sizeof( size_t ) * nGridX * nGridY * szBGrid );
            szBGrid += incSize;
            cimpl_free( bTmp );
          }
          bGrid[ gc + nGridX * nGridY * nBGrid[ gc ] ] = nPdPts;
          ++nBGrid[ gc ];
      } }

      ++nPdPts;
    }
    if( validPtFound == 0 ){
      // Remove this point from the processing list
      pList[ pListIndx ] = pList[ nPList - 1 ];
      --nPList;
    }
  }

  if( pd->nPts != nPdPts ) cimpl_pts2_shrink( pd, nPdPts );

  cimpl_free( rPds );
  cimpl_free( nBGrid );
  cimpl_free( bGrid );
  cimpl_free( pList );
}

void cimpl_poissonDisc2_tulleken( float const max_r, float (*r)(float const * ),
  float const * const bounds, cimpl_pts2 * pd ){
  // Poisson Disc 2 with varying distance parameter

  // Written according to "Fast Poisson Disk Sampling in Arbitrary Dimensions" by Robert Bridson
  // and "Poisson Disk Sampling" written by Herman Tulleken located at
  // http://devmag.org.za/2009/05/03/poisson-disk-sampling/

  // Attempts to create the number of points in the pts array
  // Gets cutoff if the number of iterations to do so exceeds maxIter
  // We leave it to the user to choose a method and appropriate action for compensating for this
  // latter case.  For example, the user may intialize all points with values that are out of bounds.
  // This could then be checked after processing.  Alternatively, the user could initialize the points
  // array to random values in the space and retain these values when the maximum number of
  // iterations is reached.

  // Inputs:
  // max_r - the maximum value of the poisson Disc parameter observable within bounds
  // bounds - a four element array specifying [ xLower yLower xUpper yUpper ]
  //   If bounds is NULL, then the bounds are [ -0.5 -0.5 0.5 0.5 ]
  // maxIter - the maximum number of iterations to run
  //
  // Outputs:
  // pd - the realized poisson disc points

  bool validPtFound, pointTooClose;
  float xL, yL, xU, yU, distX, distY, dist2NearbyPt;
  float cellSize, kDist, kAngle, rActive, kr;
  float pt[2], activePt[2], jPt[2], kPt[2];
  float *tmpFloats;
  size_t gcX, gcY, gc, gxIndx, gyIndx, nPdPts;
  size_t nGridX, nGridY;
  int kPtGridX, kPtGridY, nNearCells;
  size_t i, j, jPtIndx, pListIndx, ptIndx, nPList, pListSize;
  size_t *pList, *pTmp, *bTmp, *nBGrid, *szBGrid;
  size_t **bGrid;  // Background grid
  float *rPds = NULL;
  const size_t incSize = 1024;
  const size_t nKPts = 10;

  if( bounds == NULL ){
    xL = yL = -0.5;  xU = yU = 0.5;
  } else {
    xL = bounds[0];  xU = bounds[2];
    yL = bounds[1];  yU = bounds[3];
  }

  cellSize = max_r / sqrt(2);
  distX = xU - xL;  nGridX = ceilf( distX / cellSize );
  distY = yU - yL;  nGridY = ceilf( distY / cellSize );

  cimpl_malloc( nGridX * nGridY * sizeof(size_t*), (void**) &bGrid );
  cimpl_malloc( nGridX * nGridY * sizeof(size_t), (void**) &nBGrid );
  cimpl_malloc( nGridX * nGridY * sizeof(size_t), (void**) &szBGrid );
  for( i=0; i < nGridX * nGridY; ++i ){
    bGrid[i] = NULL;  nBGrid[i] = szBGrid[i] = 0;
  }

  cimpl_malloc( sizeof(size_t) * incSize, (void**) &pList );
  cimpl_malloc( sizeof(float) * incSize, (void**) &rPds );
  pListSize = incSize;
  *pd = cimpl_mallocPts2( incSize );

  // Create the first point and initialize the processing list
  pt[0] = cimpl_randuf( xL, xU );
  pt[1] = cimpl_randuf( yL, yU );
  cimpl_pts2_assignPt( pd, 0, pt );
  nPdPts = 1;
  rPds[0] = r( pt );
  pList[0] = 0;
  nPList = 1;
  // Include this point in the corresponding element of the background grid
  gcX = floorf( ( pt[0] - xL ) / cellSize );
  gcY = floorf( ( pt[1] - yL ) / cellSize );
  gc = gcX + gcY * nGridX;
  cimpl_malloc( sizeof(size_t) * incSize, (void**) &bGrid[ gc ] );
  bGrid[ gc ][0] = 0;
  nBGrid[ gc ] = 1;
  szBGrid[ gc ] = incSize;

  while( nPList > 0 ){
    // Choose a point at random in the processing list
    pListIndx = cimpl_randui( 0, (int) nPList-1 );
    ptIndx = pList[ pListIndx ];
    cimpl_pts2_getPt( pd, ptIndx, activePt );
    rActive = rPds[ ptIndx ];

    //Make k random points in the annulus between r and 2r from the active point
    validPtFound = 0;
    for( i=0; i<nKPts; ++i ){
      kDist = cimpl_randuf( rActive, 2 * rActive );
      kAngle = cimpl_randuf( 0, 2 * CIMPL_PI );
      kPt[0] = activePt[0] + kDist * cos( kAngle );
      kPt[1] = activePt[1] + kDist * sin( kAngle );

      // Check to see if the candidate point is valid
      if( kPt[0] < xL || kPt[0] > xU || kPt[1] < yL || kPt[1] > yU ) continue;
      kr = r( kPt );

      // Find the minimum distance to other points
      kPtGridX = floorf( ( kPt[0] - xL ) / cellSize );
      kPtGridY = floorf( ( kPt[1] - yL ) / cellSize );
      nNearCells = ceil( kr / cellSize );

      pointTooClose = 0;
      for( gyIndx = cimpl_max( kPtGridY - nNearCells, 0 );
        gyIndx <= cimpl_min( kPtGridY + nNearCells, nGridY-1); ++gyIndx ){

        for( gxIndx = cimpl_max( kPtGridX - nNearCells, 0 );
          gxIndx <= cimpl_min( kPtGridX + nNearCells, nGridX-1 ); ++gxIndx ){

          gc = gxIndx + gyIndx * nGridX;
          for( j = 0; j < nBGrid[ gc ]; ++j ){
            jPtIndx = bGrid[ gc ][j];
            cimpl_pts2_getPt( pd, jPtIndx, jPt );
            dist2NearbyPt = sqrtf( ( jPt[0] - kPt[0] ) * ( jPt[0] - kPt[0] ) +
                                   ( jPt[1] - kPt[1] ) * ( jPt[1] - kPt[1] ) );
            if( dist2NearbyPt < kr ){
              pointTooClose = 1;
              goto pointTooCloseBreak;
      } } } }
      pointTooCloseBreak:
      if( pointTooClose == 1 ) continue;

      validPtFound = 1;

      // Add this point to the processing list
      if( nPList == pListSize ){
        // Increase the size of the processing list
        pTmp = pList;
        cimpl_malloc( sizeof(size_t) * ( pListSize + incSize ), (void**) &pList );
        memcpy( pList, pTmp, sizeof(size_t) * pListSize );
        pListSize += incSize;
        cimpl_free( pTmp );
      }
      pList[ nPList ] = nPdPts;
      ++nPList;

      // Add this point to the background list
      gc = kPtGridX + kPtGridY * nGridX;
      if( nBGrid[ gc ] == szBGrid[ gc ] ){
        bTmp = bGrid[ gc ];
        cimpl_malloc( sizeof(size_t) * ( szBGrid[ gc ] + incSize ), (void**) &bGrid[ gc ] );
        if( bTmp != NULL ){
          memcpy( bGrid[gc], bTmp, sizeof(size_t) * szBGrid[ gc ] );
          cimpl_free( bTmp );
        }
        szBGrid[ gc ] += incSize;
      }
      bGrid[ gc ][ nBGrid[ gc ] ] = nPdPts;
      ++nBGrid[ gc ];

      // Add this point to the list of points
      // Add the corresponding distance threshold to rPds
      if( nPdPts == pd->nPts ){
        tmpFloats = rPds;
        cimpl_malloc( sizeof(float) * ( pd->nPts + incSize ), (void**) &rPds );
        memcpy( rPds, tmpFloats, sizeof(float) * pd->nPts );
        cimpl_free( tmpFloats );
        cimpl_pts2_grow( pd, pd->nPts + incSize );
      }
      rPds[ nPdPts ] = kr;
      cimpl_pts2_assignPt( pd, nPdPts, kPt );
      ++nPdPts;
    }

    if( validPtFound == 0 ){
      // Remove this point from the processing list
      pList[ pListIndx ] = pList[ nPList - 1 ];
      --nPList;
    }
  }

  if( pd->nPts != nPdPts ) cimpl_pts2_shrink( pd, nPdPts );

  cimpl_free( rPds );
  for( i=0; i < nGridX * nGridY; ++i ) cimpl_free( bGrid[i] );
  cimpl_free( nBGrid );
  cimpl_free( szBGrid );
  cimpl_free( bGrid );
  cimpl_free( pList );
}

void cimpl_powCmpImg( cimpl_cmpImg const in, CIMPL_COMPLEX const y, cimpl_cmpImg * const out ){
  CIMPL_COMPLEX *inPtr, *outPtr;
  inPtr = in.data;  outPtr = out->data;
  for( size_t i=0; i<in.nPix; ++i, ++inPtr, ++outPtr ){
    *inPtr = cpowf( *outPtr, y );
  }
}

void cimpl_powCmpVol( cimpl_cmpVol const in, CIMPL_COMPLEX const y, cimpl_cmpVol * const out ){
  CIMPL_COMPLEX *inPtr, *outPtr;
  inPtr = in.data;  outPtr = out->data;
  for( size_t i=0; i<in.nVox; ++i, ++inPtr, ++outPtr ){
    *inPtr = cpowf( *outPtr, y );
  }
}

void cimpl_powImg( cimpl_img const in, float const y, cimpl_img * const out ){
  float *inPtr, *outPtr;
  inPtr = in.data;  outPtr = out->data;
  for( size_t i=0; i<in.nPix; ++i, ++inPtr, ++outPtr ){
    *inPtr = powf( *outPtr, y );
  }
}

void cimpl_powVol( cimpl_vol const in, float const y, cimpl_vol * const out ){
  float *inPtr, *outPtr;
  inPtr = in.data;  outPtr = out->data;
  for( size_t i=0; i<in.nVox; ++i, ++inPtr, ++outPtr ){
    *inPtr = powf( *outPtr, y );
  }
}

void cimpl_printCmpImg( cimpl_cmpImg const in ){
  float realPart, imagPart;
  for( size_t y=0; y<in.h; ++y ){
    for( size_t x=0; x<in.w; ++x ){
      realPart = cimpl_real(in.data[y+x*in.h]);
      imagPart = cimpl_imag(in.data[y+x*in.h]);
      if( imagPart < 0 ){
        printf( "%f-i%f ", realPart, fabs(imagPart) );
      } else {
        printf( "%f+i%f ", realPart, imagPart );
      }
    }
    printf( "\n" );
  }
}

void cimpl_printCmpVol( cimpl_cmpVol const in ){
  float realPart, imagPart;
  size_t u, y, x;
  for( u=0; u<in.w; ++u ){
    printf( "Slice %lu:\n", u );
    for( y=0; y<in.h; ++y ){
      for( x=0; x<in.w; ++x ){
        realPart = cimpl_real( in.data[y+x*in.h+u*in.h*in.w] );
        imagPart = cimpl_imag( in.data[y+x*in.h+u*in.h*in.w] );
        if( imagPart < 0 ){
          printf( "%f-i%f ", realPart, fabs(imagPart) );
        } else {
          printf( "%f+i%f ", realPart, imagPart );
        }
      }
      printf( "\n" );
    }
    printf( "\n" );
  }
}

void cimpl_printImg( cimpl_img const in ){
  size_t x, y;
  for( y=0; y<in.h; ++y ){
    for( x=0; x<in.w; ++x ){
      printf( "%f ", in.data[y+x*in.h] );
    }
    printf( "\n" );
  }
}

void cimpl_printVol( cimpl_vol const in ){
  size_t u, x, y;
  for( u=0; u<in.w; ++u ){
    printf( "Slice %lu:\n", u );
    for( y=0; y<in.h; ++y ){
      for( x=0; x<in.w; ++x ){
        printf( "%f ", in.data[y+x*in.h+u*in.h*in.w] );
      }
      printf( "\n" );
    }
    printf( "\n" );
  }
}

void cimpl_pts2_assignPt( cimpl_pts2 * const out, size_t ptIndx, float const * const pt ){
  // pt is a two element array with first/second coordinate's value
  out->pts[ ptIndx ] = pt[0];
  out->pts[ ptIndx + out->nPts ] = pt[1];
}

void cimpl_pts2_free( cimpl_pts2 * const in ){
  in->nPts = 0;
  cimpl_free( in->pts );
  in->pts = NULL;
}

float cimpl_pts2_getMaxX( cimpl_pts2 const * pts ){
  float maxX;
  size_t i;

  maxX = pts->pts[0];
  for( i=1; i<pts->nPts; ++i ){
    maxX = maxX < pts->pts[i] ? maxX : pts->pts[i];
  }
  return maxX;
}

float cimpl_pts2_getMaxY( cimpl_pts2 const * pts ){
  float maxY;
  size_t i, nPts;

  nPts = pts->nPts;
  maxY = pts->pts[0];
  for( i=0; i<nPts; ++i ){
    maxY = maxY < pts->pts[i+nPts] ? maxY : pts->pts[i];
  }
  return maxY;
}

float cimpl_pts2_getMinX( cimpl_pts2 const * pts ){
  float minX;
  size_t i;

  minX = pts->pts[0];
  for( i=0; i<pts->nPts; ++i ){
    minX = minX < pts->pts[i] ? minX : pts->pts[i];
  }
  return minX;
}

float cimpl_pts2_getMinY( cimpl_pts2 const * pts ){
  float minY;
  size_t i, nPts;

  nPts = pts->nPts;
  minY = pts->pts[0];
  for( i=0; i<nPts; ++i ){
    minY = minY < pts->pts[i+nPts] ? minY : pts->pts[i];
  }
  return minY;
}


void cimpl_pts2_getPt( cimpl_pts2 const * const pts, size_t const ptIndx, float * const pt ){
  // pt is a two element array
  assert( ptIndx < pts->nPts );
  pt[0] = pts->pts[ ptIndx ];
  pt[1] = pts->pts[ ptIndx + pts->nPts ];
}

void cimpl_pts2_grow( cimpl_pts2 * const in, size_t const newNumPts ){
  float *newPts, *tmp;
  if( newNumPts < in->nPts ) return;
  cimpl_malloc( sizeof(float) * newNumPts * 2, (void**) &newPts );
  memcpy( newPts, in->pts, sizeof(float) * in->nPts );
  memcpy( newPts + newNumPts, in->pts + in->nPts, sizeof(float) * in->nPts );
  tmp = in->pts;
  in->pts = newPts;
  cimpl_free( tmp );
  in->nPts = newNumPts;
}

void cimpl_pts2_set( float const in, cimpl_pts2 * const out ){
  size_t i;
  for( i=0; i<out->nPts*2; ++i )
    out->pts[i] = in;
}

void cimpl_pts2_shrink( cimpl_pts2 * const in, size_t const newNumPts ){
  float *newPts, *tmp;
  if( newNumPts > in->nPts ) return;
  cimpl_malloc( sizeof(float) * newNumPts * 2, (void**) &newPts );
  memcpy( newPts, in->pts, sizeof(float) * newNumPts );
  memcpy( newPts + newNumPts, in->pts + in->nPts, sizeof(float) * newNumPts );
  tmp = in->pts;
  in->pts = newPts;
  cimpl_free( tmp );
  in->nPts = newNumPts;
}

void cimpl_pts3_set( float const in, cimpl_pts3 * const out ){
  size_t i;
  for( i=0; i<out->nPts*3; ++i )
    out->pts[i] = in;
}

void cimpl_ptsN_set( float const in, cimpl_ptsN * const out ){
  size_t i;
  for( i=0; i<out->nPts*out->N; ++i )
    out->pts[i] = in;
}

void cimpl_randCmpImg( cimpl_cmpImg * const out ){
  // Puts random values between 0 and 1 into all pixels of out
  CIMPL_COMPLEX* ptr;
  assert( out->data != NULL );

  ptr = out->data;
  float tmpR, tmpI;
  for( size_t i=0; i<out->nPix; ++i, ++ptr ){
    tmpR = (float) rand() / (float) RAND_MAX;
    tmpI = (float) rand() / (float) RAND_MAX;
    *ptr = tmpR + tmpI * I;
  }
}

void cimpl_randCmpVol( cimpl_cmpVol * const out ){
  // Puts random values between 0 and 1 into all pixels of out
  size_t i;
  CIMPL_COMPLEX* ptr;
  assert( out->data != NULL );

  ptr = out->data;
  float tmpR, tmpI;
  for( i=0; i<out->nVox; ++i, ++ptr ){
    tmpR = (float) rand() / (float) RAND_MAX;
    tmpI = (float) rand() / (float) RAND_MAX;
    *ptr = tmpR + tmpI * I;
  }
}

void cimpl_randImg( cimpl_img * const out ){
  // Puts random values between 0 and 1 into all pixels of out
  size_t i;
  assert( out->data != NULL );

  float* ptr = out->data;
  for( i=0; i<out->nPix; ++i, ++ptr ){
    *ptr = (float) rand() / (float) RAND_MAX;
  }
}

void cimpl_randVol( cimpl_vol * const out ){
  // Puts random values between 0 and 1 into all voxels of out
  size_t i;
  assert( out->data != NULL );

  float* ptr = out->data;
  for( i=0; i<out->nVox; ++i, ++ptr )
    *ptr = (float) rand() / (float) RAND_MAX;
}

void cimpl_rcpImg( cimpl_img const in, cimpl_img * const out ){
  // Computes the reciprocals of floating point numbers
  int i;
  assert( out->h == in.h );  assert( out->w == in.w );

  for( i=0; i<in.nPix; ++i )
    out->data[i] = 1 / in.data[i];
}

void cimpl_rcpVol( cimpl_vol const in, cimpl_vol * const out ){
  // Computes the reciprocals of floating point numbers
  size_t i;
  assert( out->h == in.h );  assert( out->w == in.w );  assert( out->s == in.s );

  for( i=0; i<in.nVox; ++i )
    out->data[i] = 1 / in.data[i];
}

void cimpl_realImg( cimpl_cmpImg const in, cimpl_img * const out ){
  CIMPL_COMPLEX* inPtr;
  assert( out->h == in.h ); assert( out->w == in.w );

  inPtr = in.data;
  float* outPtr = out->data;
  for( size_t i=0; i<out->nPix; ++i, ++inPtr, ++outPtr )
    *outPtr = cimpl_real( *inPtr );
}

void cimpl_realVol( cimpl_cmpVol const in, cimpl_vol * const out ){
  CIMPL_COMPLEX* inPtr;
  float* outPtr;
  assert( out->h == in.h ); assert( out->w == in.w ); assert( out->s == in.s );

  inPtr = in.data;
  outPtr = out->data;
  for( size_t i=0; i<out->nVox; ++i, ++inPtr, ++outPtr )
    *outPtr = cimpl_real( *inPtr );
}

void cimpl_reshapeCmpImg( size_t H, size_t W, cimpl_cmpImg * const img ){
  assert( H*W == img->h*img->w );
  img->h = H;  img->w = W;
}

void cimpl_reshapeCmpVol( size_t H, size_t W, size_t S, cimpl_vol * const out ){
  assert( H*W*S == out->h*out->w*out->s );
  out->h = H;  out->w = W;  out->s = S;
}

void cimpl_reshapeImg( size_t H, size_t W, cimpl_img * const out ){
  assert( H*W == out->h*out->w );
  out->h = H;  out->w = W;
}

void cimpl_reshapeVol( size_t H, size_t W, size_t S, cimpl_vol * const out ){
  assert( H*W*S == out->h*out->w*out->s );
  out->h = H;  out->w = W;  out->s = S;
}

//void cimpl_rot( cimpl_img const in, float const angle, cimpl_img * const out ){
//  // angle has units of radians
//  float sAngle, cAngle;
//  size_t cx, cy;
//
//  sAngle = sin( angle );
//  cAngle = cos( angle );
//  cx = ceil( (float) out->w / 2 );
//  cy = ceil( (float) out->h / 2 );
//
//  for( int x=0; x<out->w; ++x){
//    for( int y=0; y<out->h; ++y){
//
//  } }
//}

void cimpl_rgb2gray( cimpl_img const red, cimpl_img const green, cimpl_img const blue,
  cimpl_img * const out ){
  assert( out->h = red.h );  assert( out->h = green.h );  assert( out->h = blue.h );
  assert( out->w = red.w );  assert( out->w = green.w );  assert( out->w = blue.w );

  for( size_t x=0; x<red.nPix; ++x ){
      out->data[x] = 0.2989 * red.data[x] +
                     0.5870 * green.data[x] +
                     0.1140 * blue.data[x];
  }
}

void cimpl_rot90( cimpl_img const in, cimpl_img * const out ){
  assert( out->h == in.w );  assert( out->w == in.h );
  assert( out->data != in.data );
  for( size_t x=0; x<out->w; ++x ){
    for( size_t y=0; y<out->h; ++y ){
      out->data[y+x*out->h] = in.data[x+(out->h-y-1)*in.h];
  } }
}

void cimpl_rot180( cimpl_img const in, cimpl_img * const out ){
  assert( out->h == in.h );  assert( out->w == in.w );
  assert( out->data != in.data );
  for( size_t i=0; i<in.nPix; ++i )
    out->data[i] = in.data[in.nPix-i-1];
}

void cimpl_rot270( cimpl_img const in, cimpl_img * const out ){
  assert( out->h == in.w );  assert( out->w == in.h );
  assert( out->data != in.data );
  for( size_t x=0; x<out->w; ++x ){
    for( size_t y=0; y<out->h; ++y )
      out->data[y+x*out->h] = in.data[(out->w-x-1)+y*in.h];
  }
}

void cimpl_roundImg( cimpl_img const img, cimpl_img * const out ){
  assert( out->h == img.h );  assert( out->w == img.w );
  int i;
  for( i=0; i<img.nPix; ++i )
    out->data[i] = roundf( img.data[i] );
}

void cimpl_roundVol( cimpl_vol const vol, cimpl_vol * const out ){
  assert( out->h == vol.h );  assert( out->w == vol.w );  assert( out->s == vol.s );
  int i;
  for( i=0; i<vol.nVox; ++i )
    out->data[i] = roundf( vol.data[i] );
}

void cimpl_setImg( float const in, cimpl_img * const out ){
  size_t i;
  for( i=0; i<out->nPix; ++i )
    out->data[i] = in;
}

void cimpl_setVol( float const in, cimpl_vol * const out ){
  size_t i;
  for( i=0; i<out->nVox; ++i )
    out->data[i] = in;
}

float cimpl_sinc( float x ){
  float tmp;
  if( x==0 ) return 1;
  tmp = CIMPL_PI * x;
  return sinf( tmp ) / tmp;
}

void cimpl_sliceImgX( cimpl_img const in, size_t xIndx, float * const out ){
  assert( xIndx < in.w );
  memcpy( out, in.data+xIndx*in.h, in.h*sizeof(float) );
}

void cimpl_sliceImgY( cimpl_img const in, size_t yIndx, float * const out ){
  assert( yIndx < in.h );
  for( size_t x=0; x<in.h; ++x )
    out[x] = in.data[yIndx+x*in.h];
}

void cimpl_size2imgCoordinates( size_t const N, float * const coords ){
  size_t i;
  float tmp;
  assert( N > 0 );

  tmp = floorf( 0.5*N );
  for( i=0; i<N; ++i )
    coords[i] = i - tmp;
}

void cimpl_sliceX( cimpl_vol const in, size_t xIndx, cimpl_img * const out ){
  size_t x;
  assert( in.h == out->h );  assert( in.s == out->w );
  assert( xIndx < in.w );
  for( x=0; x<out->w; ++x )
    memcpy( out->data+x*out->h, in.data+xIndx*in.h+x*out->h*out->w,
            out->h*sizeof(float) );
}

void cimpl_sliceXZ( cimpl_vol const in, size_t xIndx, size_t zIndx,
  float * const out ){
  assert( xIndx < in.w );  assert( zIndx < in.s );
  memcpy( out, in.data+xIndx*in.h+zIndx*in.h*in.w, in.h*sizeof(float) );
}

void cimpl_sliceY( cimpl_vol const in, size_t yIndx, cimpl_img * const out ){
  assert( in.w == out->h );  assert( in.s == out->w );
  assert( yIndx < in.w );
  for( size_t x=0; x<out->w; ++x ){
    for( size_t y=0; y<out->h; ++y )
      out->data[y+x*out->h] = in.data[ yIndx+x*in.h+y*out->h*out->w ];
  }
}

void cimpl_sliceYX( cimpl_vol const in, size_t yIndx, size_t xIndx,
  float * const out ){
  assert( yIndx < in.h );
  assert( xIndx < in.w );
  for( size_t z=0; z<in.s; ++z )
    out[z] = in.data[yIndx+xIndx*in.h+z*in.h*in.w];
}

void cimpl_sliceYZ( cimpl_vol const in, size_t yIndx, size_t zIndx,
  float * const out ){
  assert( yIndx < in.h );
  assert( zIndx < in.s );
  for( size_t x=0; x<in.w; ++x )
    out[x] = in.data[yIndx+x*in.h+zIndx*in.h*in.w];
}

void cimpl_sliceZ( cimpl_vol const in, size_t zIndx, cimpl_img * const out ){
  assert( in.h == out->h );  assert( in.w == out->w );
  assert( zIndx < in.s );
  for( size_t x=0; x<out->w; ++x )
    memcpy(out->data+x*out->h, in.data+x*out->h+zIndx*in.h*in.w, out->h*sizeof(float) );
}

void cimpl_spaceConvImgTemplate( cimpl_img const img1, cimpl_img const t, cimpl_img * const out ){
  float tmp;
  assert( out->h == img1.h );  assert( out->w == img1.w );
  assert( t.h % 2 == 1 );
  assert( t.w % 2 == 1 );

  long hh = (long) floorf( (float) t.h * 0.5 );
  long hw = (long) floorf( (float) t.w * 0.5 );

  for( long x=0; x < (long) out->w; ++x ){
    for( long y=0; y < (long) out->h; ++y ){
      tmp = 0;

      for( long tx=-hw; tx <= hw; ++tx ){
        if( x+tx < 0 || x+tx > (long) out->w ) continue;

        for( long ty=-hh; ty <= hh; ++ty ){
          if( y+ty < 0 || y+ty > (long) out->h ) continue;

          tmp += img1.data[(y+ty)+(x+tx)*img1.w] * t.data[(ty+hh)+(tx+hw)*t.w];
      } }

      out->data[y+x*out->w] = tmp;
  } }
}

void cimpl_sqrtImg( cimpl_img const in, cimpl_img * const out ){
  int i;
  assert( out->h == in.h );  assert( out->w == in.w );

  for( i=0; i<in.nPix; ++i )
    out->data[i] = sqrtf( in.data[i] );
}

void cimpl_sqrtVol( cimpl_vol const in, cimpl_vol * const out ){
  int i;
  assert( out->h == in.h );  assert( out->w == in.w );  assert( out->s == in.s );
  for( i=0; i<in.nVox; ++i )
    out->data[i] = sqrt( in.data[i] );
}

void cimpl_softThresh( cimpl_img in, float const thresh, cimpl_img * const out ){
  int i;
  assert( out->h == in.h );  assert( out->w == in.w );
  for( i=0; i<in.nPix; ++i )
    out->data[i] = cimpl_sign( in.data[i] ) * cimpl_max( fabsf( in.data[i] ) - thresh, 0 );
}

void cimpl_softThreshCmpImg( cimpl_cmpImg in, float const thresh, cimpl_cmpImg * const out ){
  int i;
  float magIn, scalingFactor;
  assert( out->h == in.h );  assert( out->w == in.w );
  for( i=0; i<in.nPix; ++i ){
    magIn = cimpl_absCmp( in.data[i] );
    if( magIn <= thresh ){
      out->data[i] = 0;
    } else {
      scalingFactor = thresh / magIn;
      out->data[i] = cimpl_subtractCmp( in.data[i], in.data[i] * scalingFactor );
  } }
}

void cimpl_subAssignImg( cimpl_img * const img, cimpl_img const * const toAssign,
  size_t const left, size_t const low ){
  size_t N;
  for( size_t x=left; x<cimpl_min(left+toAssign->w,img->w-1); ++x ){
    N = cimpl_min( toAssign->h, img->h-low );
    memcpy(img->data+low+x*img->h, toAssign->data+(x-left)*toAssign->h,
      N*sizeof(float) );
  }
}

void cimpl_subImg( cimpl_img const in, size_t const h1, size_t const w1,
  cimpl_img * const out ){
  size_t x;
  assert( h1+out->h < in.h );  assert( w1+out->w < in.w );
  assert( out->data != in.data );
  for( x=0; x<out->w; ++x )
    memcpy( out->data+x*out->h, in.data+h1+(w1+x)*in.h, out->h*sizeof(float) );
}

void cimpl_subtractCmpImgs( cimpl_cmpImg const img1, cimpl_cmpImg const img2,
  cimpl_cmpImg * const out ){
  int i;
  assert( out->h == img1.h );  assert( out->w == img1.w );
  assert( out->h == img2.h );  assert( out->w == img2.w );
  for( i=0; i<img1.nPix; ++i )
    out->data[i] = cimpl_subtractCmp( img1.data[i], img2.data[i] );
}

void cimpl_subtractImgFromScalar( cimpl_img const in, float const scalar, cimpl_img * const out ){
  int i;
  assert( out->w == in.w );  assert( out->h == in.h );
  for( i=0; i<in.nPix; ++i )
    out->data[i] = scalar - in.data[i];
}

void cimpl_subtract( float const * const x, float const * const y, size_t const N, float * const diff ){
  int i;
  float * out = diff;
  for( i=0; i < N; ++i )
    out[i] = x[i] - y[i];
}

void cimpl_subtractImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out ){
  int i;
  assert( out->w == img1.w );  assert( out->h == img1.h );
  assert( out->w == img2.w );  assert( out->h == img2.h );
  for( i=0; i<img1.nPix; ++i )
    out->data[i] = img1.data[i] - img2.data[i];
}

void cimpl_subtractVolFromScalar( cimpl_vol const in, float const scalar, cimpl_vol * const out ){
  int i;
  assert( out->w == in.w );  assert( out->h == in.h );

  for( i=0; i<in.nVox; ++i )
    out->data[i] = scalar - in.data[i];
}

void cimpl_subtractVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out ){
  int i;
  assert( out->h == vol1.h );  assert( out->w == vol1.w );  assert( out->s == vol1.s );
  assert( out->h == vol2.h );  assert( out->w == vol2.w );  assert( out->s == vol2.s );

  for( i=0; i<vol1.nVox; ++i )
    out->data[i] = vol1.data[i] - vol2.data[i];
}

void cimpl_subtractScalarFromImg( cimpl_img const in, float const scalar, cimpl_img * const out ){
  cimpl_addScalar2Img( -scalar, in, out );
}

void cimpl_subtractScalarFromVol( cimpl_vol const in, float const scalar, cimpl_vol * const out ){
  cimpl_addScalar2Vol( -scalar, in, out );
}

void cimpl_sumImg( cimpl_img * const in, float * out ){
  size_t i;

  *out = 0;
  for( i=0; i<in->nPix; ++i )
    *out += in->data[i];
}

float cimpl_sumVol( cimpl_vol const in ){
  size_t i;
  float out;

  out = 0;
  for( i=0; i<in.nVox; ++i ){
    out += in.data[i];
  }
  return out;
}

void cimpl_transposeImg( cimpl_img const in, cimpl_img * const out ){
  size_t x, y;
  assert( out->h == in.w );  assert( out->w == in.h );

  for( x=0; x<out->w; ++x ){
    for( y=0; y<out->h; ++y )
      out->data[y+x*out->h] = in.data[x+y*in.h];
  }
}

void cimpl_upsampleCmpImg( cimpl_cmpImg const in, int const upsample, 
  cimpl_cmpImg * const out ){
  size_t x, y;
  assert( out->h == in.h * upsample ); assert( out->w == in.w * upsample );

  cimpl_zeroCmpImg( out );
  for( x=0; x<in.w; ++x ){
    for( y=0; y<in.h; ++y ){
      out->data[(y*upsample)+(x*upsample*out->h)] = in.data[y+(x*in.h)];
  } }
}

void cimpl_upsampleImg( cimpl_img const in, int const upsample, cimpl_img * const out ){
  size_t x, y;
  assert( out->h == in.h * upsample );
  assert( out->w == in.w * upsample );

  cimpl_zeroImg( out );
  for( x=0; x<in.w; ++x ){
    for( y=0; y<in.h; ++y ){
      out->data[(y*upsample)+(x*upsample*out->h)] = in.data[y+(x*in.h)];
  } }
}

void cimpl_upsampleCmpVol( cimpl_cmpVol const in, int const upsample, 
  cimpl_cmpVol * const out ){
  size_t u, x, y;
  assert( out->h == in.h * upsample );
  assert( out->w == in.w * upsample );
  assert( out->s == in.s * upsample );

  cimpl_zeroCmpVol( out );
  for( u=0; u<in.s; ++u ){
    for( x=0; x<in.w; ++x ){
      for( y=0; y<in.h; ++y ){
        out->data[(y*upsample)+(x*upsample*out->h)+(u*upsample*out->h*out->w)] =
            in.data[y+(x*in.h)];
  } } }
}

void cimpl_upsampleVol( cimpl_vol const in, int const upsample, cimpl_vol * const out ){
  size_t u, x, y;
  assert( out->h == in.h * upsample );
  assert( out->w == in.w * upsample );
  assert( out->s == in.s * upsample );

  cimpl_zeroVol( out );
  for( u=0; u<in.s; ++u ){
    for( x=0; x<in.w; ++x ){
      for( y=0; y<in.h; ++y ){
        out->data[(y*upsample)+(x*upsample*out->h)+(u*upsample*out->h*out->w)] =
            in.data[y+(x*in.h)];
  } } }
}

void cimpl_zeroCmpImg( cimpl_cmpImg * const img  ){
  // Sets all values of complex image to 0
  int i;
  float * inData;
  assert( img->data != NULL );

  inData = (float*) img->data;
  for( i=0; i<img->nPix*2; ++i ) inData[i] = 0;
}

void cimpl_zeroCmpVol( cimpl_cmpVol * const vol  ){
  // Sets all values of complex image to 0
  int i;
  assert( vol->data != NULL );

  for( i=0; i<vol->nVox*2; ++i ) vol->data[i] = 0;
}

void cimpl_zeroImg( cimpl_img * const img  ){
  // Sets all values of image to 0
  int i;
  assert( img->data != NULL );

  for( i=0; i<img->nPix; ++i ) img->data[i] = 0;
}

void cimpl_zeroVol( cimpl_vol * const vol  ){
  size_t i;
  assert( vol->data != NULL );

  for( i=0; i<vol->nVox; ++i ) vol->data[i] = 0;
}

