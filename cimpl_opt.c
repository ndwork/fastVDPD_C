//
//  cimpl_opt.c
//  cimpl
//
//  Created by Nicholas Dwork on 4/2/20.
//  Copyright © 2020 ndwork. All rights reserved.
//

#include "cimpl_opt.h"

#include <assert.h>

void cimpl_binarySearch( void (*f)(float const *, float* ), float boundL, float boundU,
  float tol, unsigned int const nMax, float *param ){
  float L, U;  // lower, mid, and upper domain values
  float fL, fParam;  // lower, mid, and upper function values
  float err;
  unsigned int i;
  assert( boundU >= boundL );

  *param = L = boundL;
  U = boundU;
  f( &L, &fL );

  for( i = 0; i < nMax; ++i ){
    err = 0.5 * ( U - L );
    *param = L + err;
    if( err < tol ) break;

    f( param, &fParam );
    if( fParam == 0 ) break;

    if( fL * fParam > 0 ){
      L = *param;
    } else {
      U = *param;
    }
  }
}


