//
//  cimpl_opt.h
//  cimpl
//
//  Created by Nicholas Dwork on 4/2/20.
//  Copyright © 2020 ndwork. All rights reserved.
//

#ifndef cimpl_opt_h
#define cimpl_opt_h



void cimpl_binarySearch( void (*r)(float const *, float* ), float boundL, float boundU, float tol,
  unsigned int const nMax, float* param );

#endif /* cimpl_opt_h */
