//
//  cimpl.h
//  cimpl
//
//  Created by Nicholas Dwork starting on 9/18/16.
//  Copyright © 2016 Nicholas Dwork.
//

#ifndef cimpl_h
#define cimpl_h

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include <stddef.h>

#ifdef _MSC_VER
  #define CIMPL_COMPLEX _C_float_complex
#elif __MINGW32__
  #define CIMPL_COMPLEX _Complex
#else
  #define CIMPL_COMPLEX float complex
#endif

#define CIMPL_PI 3.1415926535897932384626

typedef enum { cimpl_false=0, cimpl_true } cimpl_bool;

typedef struct {
  size_t h;  // height
  size_t w;  // width
  size_t nPix;  // number of pixels in image
  float* data;  // column major ordering
} cimpl_img;

typedef struct {
  size_t h;  // height
  size_t w;  // width
  size_t nPix;  // number of pixels in image
  CIMPL_COMPLEX* data;
} cimpl_cmpImg;

typedef struct {
  size_t h;  // height
  size_t w;  // width
  size_t s;  // number of slices
  size_t nVox;  // number of voxels in volume
  float* data;  // column major ordering
} cimpl_vol;

typedef struct {
  size_t h;  // height
  size_t w;  // width
  size_t s;  // number of slices
  size_t nVox;  // number of voxels in volume
  CIMPL_COMPLEX* data;
} cimpl_cmpVol;

typedef struct {
  size_t nPts; // number of points
  float *pts;  // Coordinates of points, first dimension coordinates followed by second dimension
} cimpl_pts2;  // Points in two dimensions

typedef struct {
  size_t nPts; // number of points
  float* pts;  // Coordinates of points, first dimension coordinates followed by second dimension
} cimpl_pts3;  // Points in three dimensions

typedef struct {
  size_t nPts; // number of points
  unsigned int N;  // Number of dimensions
  float* pts;  // Coordinates of points, first dimension coordinates followed by second dimension
} cimpl_ptsN;  // Points in N dimensions

typedef enum {
  SUCCESS = 0,
  ERROR,
  FILE_ERROR,
  FILE_OPEN_ERROR,
  FILE_TYPE_ERROR,
  IO_ERROR,
  READ_ERROR,
  WRITE_ERROR
} cimpl_status;

typedef enum {
  ZERO = 0,
  REPEATED,
  SYMMETRIC
} cimpl_boundaryCondition;


void cimpl_absImg( cimpl_img const in, cimpl_img * const out );
void cimpl_absCmpImg( cimpl_cmpImg const in, cimpl_img * const out );
void cimpl_absCmpVol( cimpl_cmpVol const in, cimpl_vol * const out );
void cimpl_absVol( cimpl_vol const in, cimpl_vol * const out );
void cimpl_addCmpImgs( cimpl_cmpImg const img1, cimpl_cmpImg const img2,
  cimpl_cmpImg * const out );
void cimpl_addImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_addVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_addScalar2Img( float const scalar, cimpl_img const in, cimpl_img * const out );
void cimpl_addScalar2Vol( float const scalar, cimpl_vol const in, cimpl_vol * const out );
void cimpl_argCmpImg( cimpl_cmpImg const in, cimpl_img * const out );
void cimpl_argCmpVol( cimpl_cmpVol const in, cimpl_vol * const out );
float cimpl_besseli0( float x );
float cimpl_besselj0( float x );
void cimpl_ceilImg( cimpl_img const img, cimpl_img * const out );
void cimpl_ceilVol( cimpl_vol const vol, cimpl_vol * const out );
void cimpl_circShiftImg( cimpl_img const in, long hShift, long vShift, cimpl_img * const out );
void cimpl_circShiftVol( cimpl_vol const in, int hShift, int vShift, int sShift,
  cimpl_vol * const out );
void cimpl_cmpEqImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_cmpEqVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_cmpGeImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_cmpGeVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_cmpGtImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_cmpGtVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_cmpLeImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_cmpLeVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_cmpLtImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_cmpLtVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_img * const out );
void cimpl_cmpLtEqImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_cmpLtEqVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_cmpNeqImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_cmpNeqVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_concatCmpImgsH( cimpl_cmpImg const img1, cimpl_cmpImg const img2, cimpl_cmpImg * const out );
void cimpl_concatCmpImgsW( cimpl_cmpImg const img1, cimpl_cmpImg const img2, cimpl_cmpImg * const out );
void cimpl_concatImgsH( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_concatImgsW( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_concatVolsH( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_concatVolsS( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_concatVolsW( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_concatCmpVolsH( cimpl_cmpVol const vol1, cimpl_cmpVol const vol2, cimpl_cmpVol * const out );
void cimpl_concatCmpVolsS( cimpl_cmpVol const vol1, cimpl_cmpVol const vol2, cimpl_cmpVol * const out );
void cimpl_concatCmpVolsW( cimpl_cmpVol const vol1, cimpl_cmpVol const vol2, cimpl_cmpVol * const out );
void cimpl_conjCmpImg( cimpl_cmpImg const in, cimpl_cmpImg * const out );
void cimpl_copyCmpImg( cimpl_cmpImg in, cimpl_cmpImg * copy );
void cimpl_copyImg( cimpl_img in, cimpl_img * copy );
void cimpl_cropImg( cimpl_img in, cimpl_img * out );
void cimpl_cropVol( cimpl_vol in, cimpl_vol * out );
void cimpl_d1hImg( cimpl_img const in, cimpl_img * const out );
void cimpl_d1hVol( cimpl_vol in, cimpl_vol * out );
void cimpl_d1sVol( cimpl_vol in, cimpl_vol * out );
void cimpl_d1wImg( cimpl_img in, cimpl_img * out );
void cimpl_d1wVol( cimpl_vol in, cimpl_vol * out );
void cimpl_divideCmpImgs( cimpl_cmpImg img1, cimpl_cmpImg img2, cimpl_cmpImg * out );
void cimpl_divideCmpVols( cimpl_cmpVol vol1, cimpl_cmpVol vol2, cimpl_cmpVol * out );
void cimpl_divideImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_divideVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_divideImgByScalar( cimpl_img const in, float const scalar, cimpl_img * const out );
void cimpl_divideVolByScalar( cimpl_vol const in, float const scalar, cimpl_vol * const out );
float cimpl_dotProd( float const * const x, float const * const y, size_t const N );
float cimpl_dotImgs( cimpl_img const img1, cimpl_img const img2 );
void cimpl_downsampleCmpImg( cimpl_cmpImg const img1, int downsample, cimpl_cmpImg * const out );
void cimpl_downsampleImg( cimpl_img const img1, int downsample, cimpl_img * const out );
void cimpl_downsampleVol( cimpl_vol const vol1, int downsample, cimpl_vol * const out );
int cimpl_equalCmpImgs( cimpl_cmpImg const img1, cimpl_cmpImg const img2 );
int cimpl_equalCmpVols( cimpl_cmpVol const vol1, cimpl_cmpVol const vol2 );
int cimpl_equalImgs( cimpl_img const img1, cimpl_img const img2 );
int cimpl_equalVols( cimpl_vol const vol1, cimpl_vol const vol2 );
void cimpl_flipImgLR( cimpl_img const in, cimpl_img * const out );
void cimpl_flipImgUD( cimpl_img const in, cimpl_img * const out );
void cimpl_floorImg( cimpl_img const img, cimpl_img * const out );
void cimpl_floorVol( cimpl_vol const vol, cimpl_vol * const out );
void cimpl_free( void * in );
void cimpl_freeCmpImg( cimpl_cmpImg * in );
void cimpl_freeCmpVol( cimpl_cmpVol * in );
void cimpl_freeImg( cimpl_img * const in );
void cimpl_freeVol( cimpl_vol * const in );
void cimpl_imagImg( cimpl_cmpImg const in, cimpl_img * const out );
void cimpl_imagVol( cimpl_cmpVol const in, cimpl_vol * const out );
float cimpl_linInterp( size_t const N, float const * const x, float const * const y,
  float const outOfBounds, float const q );
void cimpl_linInterps( size_t const N, float const * const x, float const * const y,
  float const outOfBounds, size_t const M, float const * const q, float * const out );
void cimpl_linInterpImg( cimpl_img const img, size_t const N, float const * const xq,
  float const * const yq, float const outOfBounds, float * const out );
void cimpl_malloc( size_t nBytes, void ** ptrPtr );
cimpl_cmpImg cimpl_mallocCmpImg( size_t const h, size_t const w );
cimpl_cmpVol cimpl_mallocCmpVol( size_t const h, size_t const w , size_t const s );
cimpl_img cimpl_mallocImg( size_t const h, size_t const w );
cimpl_vol cimpl_mallocVol( size_t const h, size_t const w, size_t const s );
void cimpl_maxImgs( cimpl_img img1, cimpl_img img2, cimpl_img * out );
void cimpl_maxVols( cimpl_vol vol1, cimpl_vol vol2, cimpl_vol * out );
void cimpl_meanCmpImg( cimpl_cmpImg img, CIMPL_COMPLEX * out );
void cimpl_meanImg( cimpl_img img, float * out );
void cimpl_meanVol( cimpl_vol vol, float * out );
void cimpl_minImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_minVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_movePts2ToGrid( float * const mins, float * const maxs, cimpl_pts2 * pts );
void cimpl_multiplyCmpImgByRealScalar( cimpl_cmpImg const in, float const scalar,
  cimpl_cmpImg * const out );
void cimpl_multiplyCmpImgs( cimpl_cmpImg const img1, cimpl_cmpImg const img2, cimpl_cmpImg * const out );
void cimpl_multiplyImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_multiplyVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_multiplyImgByScalar( cimpl_img const in, float const scalar, cimpl_img * const out );
void cimpl_multiplyVolByScalar( cimpl_vol const in, float const scalar, cimpl_vol * const out );
float cimpl_normCmpImgL2( cimpl_cmpImg const in );
float cimpl_normImgL1( cimpl_img const in );
float cimpl_normImgL2( cimpl_img const in );
cimpl_pts2 cimpl_poissonDisc2( float r, float const * bounds );
void cimpl_poissonDisc2_dwork( float const min_r, float (*r)(float const * ), float const * const bounds,
  cimpl_pts2 * pd );
void cimpl_poissonDisc2_tulleken( float const max_r, float (*r)(float const * ),
  float const * const bounds, cimpl_pts2 * pd );
void cimpl_powCmpImg( cimpl_cmpImg in, CIMPL_COMPLEX y, cimpl_cmpImg * out );
void cimpl_powCmpVol( cimpl_cmpVol in, CIMPL_COMPLEX y, cimpl_cmpVol * out );
void cimpl_powImg( cimpl_img in, float y, cimpl_img * out );
void cimpl_powVol( cimpl_vol in, float y, cimpl_vol * out );
void cimpl_printCmpImg( cimpl_cmpImg in );
void cimpl_printCmpVol( cimpl_cmpVol in );
void cimpl_printImg( cimpl_img in );
void cimpl_printVol( cimpl_vol in );
void cimpl_pts2_assignPt( cimpl_pts2 * out, size_t ptIndx, float const * pt );
void cimpl_pts2_free( cimpl_pts2 * const in );
void cimpl_pts2_getPt( cimpl_pts2 const * pts, size_t const ptIndx, float * pt );
float cimpl_pts2_getMaxX( cimpl_pts2 const * pts );
float cimpl_pts2_getMaxY( cimpl_pts2 const * pts );
float cimpl_pts2_getMinX( cimpl_pts2 const * pts );
float cimpl_pts2_getMinY( cimpl_pts2 const * pts );
void cimpl_pts2_grow( cimpl_pts2 * const in, size_t const newNumPts );
void cimpl_pts2_set( float const in, cimpl_pts2 * const out );
void cimpl_pts2_shrink( cimpl_pts2 * const in, size_t const newNumPts );
void cimpl_pts3_set( float const in, cimpl_pts3 * const out );
void cimpl_ptsN_set( float const in, cimpl_ptsN * const out );
void cimpl_randCmpImg( cimpl_cmpImg * const out );
void cimpl_randCmpVol( cimpl_cmpVol * const out );
void cimpl_randImg( cimpl_img * const out );
void cimpl_rcpImg( cimpl_img const in, cimpl_img * const out );
void cimpl_rcpVol( cimpl_vol const in, cimpl_vol * const out );
void cimpl_realImg( cimpl_cmpImg const in, cimpl_img * const out );
void cimpl_realVol( cimpl_cmpVol const in, cimpl_vol * const out );
void cimpl_reshapeCmpImg( size_t H, size_t W, cimpl_cmpImg * const img );
void cimpl_reshapeImg( size_t H, size_t W, cimpl_img * const out );
void cimpl_reshapeVol( size_t H, size_t W, size_t S, cimpl_vol * const out );
void cimpl_rgb2gray( cimpl_img const red, cimpl_img const green, cimpl_img const blue,
  cimpl_img * const out );
//void cimpl_rot( cimpl_img const in, float const angle, cimpl_img * const out );
void cimpl_rot90( cimpl_img const in, cimpl_img * const out );
void cimpl_rot180( cimpl_img const in, cimpl_img * const out );
void cimpl_rot270( cimpl_img const in, cimpl_img * const out );
void cimpl_roundImg( cimpl_img const img, cimpl_img * const out );
void cimpl_roundVol( cimpl_vol const vol, cimpl_vol * const out );
void cimpl_setImg( float const in, cimpl_img * const out );
void cimpl_setVol( float const in, cimpl_vol * const out );
float cimpl_sinc( float x );
void cimpl_size2imgCoordinates( size_t const N, float * const coords );
void cimpl_sliceX( cimpl_vol const in, size_t xIndx, cimpl_img * const out );
void cimpl_sliceXZ( cimpl_vol const in, size_t xIndx, size_t zIndx, float * const out );
void cimpl_sliceY( cimpl_vol const in, size_t yIndx, cimpl_img * const out );
void cimpl_sliceYX( cimpl_vol const in, size_t yIndx, size_t xIndx, float * const out );
void cimpl_sliceYZ( cimpl_vol const in, size_t yIndx, size_t zIndx, float * const out );
void cimpl_sliceZ( cimpl_vol const in, size_t zIndx, cimpl_img * const out );
void cimpl_spaceConvImgTemplate( cimpl_img const img1, cimpl_img const t,
  cimpl_img * const out );
void cimpl_sqrtImg( cimpl_img const in, cimpl_img * const out );
void cimpl_sqrtVol( cimpl_vol const in, cimpl_vol * const out );
void cimpl_softThresh( cimpl_img in, float const thresh, cimpl_img * const out );
void cimpl_subAssignImg( cimpl_img * const img, cimpl_img const * const toAssign,
  size_t const left, size_t const low );
void cimpl_subtract( float const * const x, float const * const y, size_t const N, float * const diff );
void cimpl_subImg( cimpl_img const in, size_t const h1, size_t const v1,
  cimpl_img * const out );
void cimpl_subtractCmpImgs( cimpl_cmpImg const img1, cimpl_cmpImg const img2,
  cimpl_cmpImg * const out );
void cimpl_subtractImgFromScalar( cimpl_img const in, float const scalar, cimpl_img * const out );
void cimpl_subtractImgs( cimpl_img const img1, cimpl_img const img2, cimpl_img * const out );
void cimpl_subtractVolFromScalar( cimpl_vol const in, float const scalar, cimpl_vol * const out );
void cimpl_subtractVols( cimpl_vol const vol1, cimpl_vol const vol2, cimpl_vol * const out );
void cimpl_subtractScalarFromImg( cimpl_img const in, float const scalar, cimpl_img * const out );
void cimpl_subtractScalarFromVol( cimpl_vol const in, float const scalar, cimpl_vol * const out );
void cimpl_sumImg( cimpl_img * const in, float * out );
float cimpl_sumVol( cimpl_vol const in );
void cimpl_transposeImg( cimpl_img const in, cimpl_img * const out );
void cimpl_upsampleCmpImg( cimpl_cmpImg const in, int const upsample, cimpl_cmpImg * const out );
void cimpl_upsampleCmpVol( cimpl_cmpVol const in, int const upsample, cimpl_cmpVol * const out );
void cimpl_upsampleImg( cimpl_img const in, int const upsample, cimpl_img * const out );
void cimpl_upsampleVol( cimpl_vol const in, int const upsample, cimpl_vol * const out );
void cimpl_zeroCmpImg( cimpl_cmpImg * const img  );
void cimpl_zeroCmpVol( cimpl_cmpVol * const vol  );
void cimpl_zeroImg( cimpl_img * const in  );
void cimpl_zeroVol( cimpl_vol * const in  );


#ifdef __cplusplus
}
#endif

#endif /* cimpl_h */
