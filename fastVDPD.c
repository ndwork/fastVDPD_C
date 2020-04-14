//
//  fastVDPD.c
//  cimpl
//
//  Created by Nicholas Dwork on 3/20/20.
//  Copyright © 2020 ndwork. All rights reserved.
//

#include "fastVDPD.h"

#include <math.h>
#include <time.h>

#include "cimpl_opt.h"
#include "cimpl_utility.h"
#include "cimpl_io.h"


#include <stdio.h>



float runFastVDPD_gamma;  // Global variables.  Ugly, I know.  It works.  So there.
float fastVDPD_accelerationRate;
size_t fastVDPD_sizeX, fastVDPD_sizeY;

float r( float const * in ){
  // In is two element float array
  float delta = 0.15;
  float gamma = runFastVDPD_gamma;
  float dist, out;

  dist = sqrtf( in[0]*in[0] + in[1]*in[1] );
  out = ( dist + delta ) / gamma;
  return out;
}


void foDwork( float const * in, float * out ){
  float corner[2];
  float bounds[4] = { -0.5, -0.5, 0.5, 0.5 };
  float dx, dy, nSamples;
  float accelerationRate;
  size_t nPts, ptIndx, imgIndx;
  cimpl_img mask;
  corner[0] = bounds[2];
  corner[1] = bounds[3];

  //float max_r = r( corner );
  float smallest[2] = { 0, 0 };
  float min_r = r( smallest );

  cimpl_pts2 pts;
  cimpl_poissonDisc2_dwork( min_r, r, bounds, &pts );
  dx = ( bounds[2] - bounds[0] ) / fastVDPD_sizeX;
  dy = ( bounds[3] - bounds[1] ) / fastVDPD_sizeY;
  nPts = pts.nPts;
  mask = cimpl_mallocImg( fastVDPD_sizeX, fastVDPD_sizeY );
  cimpl_zeroImg( &mask );
  for( ptIndx=0; ptIndx<nPts; ++ptIndx ){
    pts.pts[ptIndx] = floorf( ( pts.pts[ptIndx] - bounds[0] ) / dx );
    pts.pts[ptIndx+nPts] = floorf( ( pts.pts[ptIndx+nPts] - bounds[1] ) / dy );
    imgIndx = pts.pts[ptIndx] + fastVDPD_sizeX * pts.pts[ptIndx+nPts];
    mask.data[ imgIndx ] = 1;
  }
  cimpl_sumImg( &mask, &nSamples );
  accelerationRate = (float) fastVDPD_sizeX * (float) fastVDPD_sizeY / nSamples;
  printf( "Dwork Gamma: %f, Acceleration Rate: %f\n", runFastVDPD_gamma, accelerationRate );
  *out = accelerationRate - fastVDPD_accelerationRate;

  cimpl_pts2_free( &pts );
  cimpl_freeImg( &mask );
}


void foTulleken( float const * in, float * out ){
  float corner[2];
  float bounds[4] = { -0.5, -0.5, 0.5, 0.5 };
  float dx, dy, nSamples;
  float accelerationRate;
  cimpl_img mask;
  corner[0] = bounds[2];
  corner[1] = bounds[3];
  float max_r = r( corner );
  size_t nPts, ptIndx, imgIndx;

  //float smallest[2] = { 0, 0 };
  //float min_r = r( smallest );
  cimpl_pts2 pts;
  cimpl_poissonDisc2_tulleken( max_r, r, bounds, &pts );
  dx = ( bounds[2] - bounds[0] ) / fastVDPD_sizeX;
  dy = ( bounds[3] - bounds[1] ) / fastVDPD_sizeY;
  nPts = pts.nPts;
  mask = cimpl_mallocImg( fastVDPD_sizeX, fastVDPD_sizeY );
  cimpl_zeroImg( &mask );
  for( ptIndx=0; ptIndx<nPts; ++ptIndx ){
    pts.pts[ptIndx]      = floorf( ( pts.pts[ptIndx]      - bounds[0] ) / dx );
    pts.pts[ptIndx+nPts] = floorf( ( pts.pts[ptIndx+nPts] - bounds[1] ) / dy );
    imgIndx = pts.pts[ptIndx] + fastVDPD_sizeX * pts.pts[ptIndx+nPts];
    mask.data[ imgIndx ] = 1;
  }
  cimpl_sumImg( &mask, &nSamples );
  accelerationRate = (float) fastVDPD_sizeX * (float) fastVDPD_sizeY / nSamples;
  printf( "Tulleken Gamma: %f, Acceleration Rate: %f\n", runFastVDPD_gamma, accelerationRate );
  *out = accelerationRate - fastVDPD_accelerationRate;

  cimpl_pts2_free( &pts );
  cimpl_freeImg( &mask );
}


void runFastVDPD( void ){
  size_t ptIndx, imgIndx, nPts;
  cimpl_img mask_tulleken, mask_dwork;
  cimpl_pts2 pd_tulleken, pd_dwork;
  float bounds[4];
  //void* bounds = NULL;
  float corner[2];
  float smallest[2] = { 0, 0 };
  float dx, dy;
  float undersampleX, undersampleY;
  float sumTulleken, samplingPercentageTulleken, timeTulleken;
  float sumDwork, samplingPercentageDwork, timeDwork;
  clock_t startTulleken, startDwork, endTulleken, endDwork;


  fastVDPD_sizeX = fastVDPD_sizeY = 256;
  undersampleX = 1;
  undersampleY = 1;
  runFastVDPD_gamma = 150;

  bounds[0] = -0.5 / undersampleX;  bounds[2] = 0.5 / undersampleX;
  bounds[1] = -0.5 / undersampleY;  bounds[3] = 0.5 / undersampleY;
  dx = ( bounds[2] - bounds[0] ) / fastVDPD_sizeX;
  dy = ( bounds[3] - bounds[1] ) / fastVDPD_sizeY;

  mask_dwork.data = NULL;    mask_tulleken.data = NULL;
  pd_dwork.pts = NULL;       pd_tulleken.pts = NULL;

  char* tullekenFile = "/Users/nicholasdwork/Dropbox/school/Stanford/ee391/cimpl/Output/tullekenMask.png";
  char* dworkFile = "/Users/nicholasdwork/Dropbox/school/Stanford/ee391/cimpl/Output/dworkMask.png";
  char* tullekenAccelFile =
  "/Users/nicholasdwork/Dropbox/school/Stanford/ee391/cimpl/Output/tullekenAcceleration.png";
  char* dworkAccelFile =
  "/Users/nicholasdwork/Dropbox/school/Stanford/ee391/cimpl/Output/dworkAcceleration.png";

  corner[0] = bounds[2];
  corner[1] = bounds[3];
  float max_r = r( corner );
  float min_r = r( smallest );

  startTulleken = clock();
  cimpl_poissonDisc2_tulleken( max_r, r, bounds, &pd_tulleken );
  endTulleken = clock();
  nPts = pd_tulleken.nPts;
  mask_tulleken = cimpl_mallocImg( fastVDPD_sizeX, fastVDPD_sizeY );
  cimpl_zeroImg( &mask_tulleken );
  for( ptIndx=0; ptIndx<nPts; ++ptIndx ){
    pd_tulleken.pts[ptIndx]      = floorf( ( pd_tulleken.pts[ptIndx] - bounds[0] ) / dx );
    pd_tulleken.pts[ptIndx+nPts] = floorf( ( pd_tulleken.pts[ptIndx+nPts] - bounds[1] ) / dy );

    imgIndx = pd_tulleken.pts[ptIndx] + fastVDPD_sizeX * pd_tulleken.pts[ptIndx+nPts];
    mask_tulleken.data[ imgIndx ] = 255;
  }
  cimpl_writePng( &mask_tulleken, NULL, NULL, NULL, tullekenFile );
  cimpl_sumImg( &mask_tulleken, &sumTulleken );
  samplingPercentageTulleken = sumTulleken / 255 / mask_tulleken.nPix;
  timeTulleken = ( (float) (endTulleken - startTulleken) ) / CLOCKS_PER_SEC;
  printf( "Sampling Fraction Tulleken: %f\n", samplingPercentageTulleken );
  printf( "Processing Time Tulleken: %f\n", timeTulleken );
  cimpl_pts2_free( &pd_tulleken );

  startDwork = clock();
  cimpl_poissonDisc2_dwork( min_r, r, bounds, &pd_dwork );
  endDwork = clock();
  nPts = pd_dwork.nPts;
  mask_dwork = cimpl_mallocImg( fastVDPD_sizeX, fastVDPD_sizeY );
  cimpl_zeroImg( &mask_dwork );
  for( ptIndx=0; ptIndx<nPts; ++ptIndx ){
    pd_dwork.pts[ptIndx]      = floorf( ( pd_dwork.pts[ptIndx]      - bounds[0] ) / dx );
    pd_dwork.pts[ptIndx+nPts] = floorf( ( pd_dwork.pts[ptIndx+nPts] - bounds[1] ) / dy );

    imgIndx = pd_dwork.pts[ptIndx] + fastVDPD_sizeX * pd_dwork.pts[ptIndx+nPts];
    mask_dwork.data[ imgIndx ] = 255;
  }
  cimpl_writePng( &mask_dwork, NULL, NULL, NULL, dworkFile );
  cimpl_sumImg( &mask_dwork, &sumDwork );
  samplingPercentageDwork = sumDwork / 255 / mask_tulleken.nPix;
  timeDwork = ( (float) (endDwork - startDwork) ) / CLOCKS_PER_SEC;
  printf( "Sampling Fraction Dwork: %f\n", samplingPercentageDwork );
  printf( "Processing Time Dwork: %f\n", timeDwork );
  printf( "Time Dwork / Time Tulleken: %f\n", timeDwork / timeTulleken );
  cimpl_pts2_free( &pd_dwork );
  

  // Now generate samples for specific acceleration rates
  float minGamma, maxGamma;
  float tolerance;
  unsigned int nMax;
  fastVDPD_accelerationRate = 4.5;
  minGamma = 0.15 * (float) fastVDPD_sizeX;
  maxGamma = 500;
  nMax = 100;
  tolerance = 0.01;

  runFastVDPD_gamma = 0;
  startDwork = clock();
  cimpl_binarySearch( foDwork, minGamma, maxGamma, tolerance, nMax, &runFastVDPD_gamma );
  endDwork = clock();
  timeDwork = ( (float) (endDwork - startDwork) ) / CLOCKS_PER_SEC;
  printf( "Dwork gamma, time(s): %f, %f\n", runFastVDPD_gamma, timeDwork );

  cimpl_poissonDisc2_dwork( min_r, r, bounds, &pd_dwork );
  nPts = pd_dwork.nPts;
  cimpl_zeroImg( &mask_dwork );
  for( ptIndx=0; ptIndx<nPts; ++ptIndx ){
    pd_dwork.pts[ptIndx]      = floorf( ( pd_dwork.pts[ptIndx]      - bounds[0] ) / dx );
    pd_dwork.pts[ptIndx+nPts] = floorf( ( pd_dwork.pts[ptIndx+nPts] - bounds[1] ) / dy );

    imgIndx = pd_dwork.pts[ptIndx] + fastVDPD_sizeX * pd_dwork.pts[ptIndx+nPts];
    mask_dwork.data[ imgIndx ] = 255;
  }
  float foDworkResult;
  foDwork( NULL, &foDworkResult );
  cimpl_writePng( &mask_dwork, NULL, NULL, NULL, dworkAccelFile );
  cimpl_pts2_free( &pd_dwork );
  cimpl_freeImg( &mask_dwork );

  
  runFastVDPD_gamma = 0;
  startTulleken = clock();
  cimpl_binarySearch( foTulleken, minGamma, maxGamma, tolerance, nMax, &runFastVDPD_gamma );
  endTulleken = clock();
  timeTulleken = ( (float) (endTulleken - startTulleken) ) / CLOCKS_PER_SEC;
  printf( "Tulleken gamma, time(s): %f, %f\n", runFastVDPD_gamma, timeTulleken );

  cimpl_poissonDisc2_tulleken( max_r, r, bounds, &pd_tulleken );
  nPts = pd_tulleken.nPts;
  cimpl_zeroImg( &mask_tulleken );
  for( ptIndx=0; ptIndx<nPts; ++ptIndx ){
    pd_tulleken.pts[ptIndx]      = floorf( ( pd_tulleken.pts[ptIndx]      - bounds[0] ) / dx );
    pd_tulleken.pts[ptIndx+nPts] = floorf( ( pd_tulleken.pts[ptIndx+nPts] - bounds[1] ) / dy );

    imgIndx = pd_tulleken.pts[ptIndx] + fastVDPD_sizeX * pd_tulleken.pts[ptIndx+nPts];
    mask_tulleken.data[ imgIndx ] = 255;
  }
  float foTullekenResult, tullekenSum;
  cimpl_sumImg( &mask_tulleken, &tullekenSum );
  foTulleken( NULL, &foTullekenResult );
  cimpl_writePng( &mask_tulleken, NULL, NULL, NULL, tullekenAccelFile );
  cimpl_pts2_free( &pd_tulleken );
  cimpl_freeImg( &mask_tulleken );
}
