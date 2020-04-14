//
//  cimpl_io.h
//  cimpl
//
//  Created by Nicholas Dwork on 8/29/17.
//  Copyright © 2017 Nicholas Dwork. All rights reserved.
//

#ifndef cimpl_io_h
#define cimpl_io_h

#ifdef __cplusplus
extern "C" {
#endif

#include "cimpl.h"

cimpl_status cimpl_readJpeg( char const * const inFilename,
  cimpl_img * outRed, cimpl_img * outGreen, cimpl_img * outBlue );

cimpl_status cimpl_readPng( char const * const inFilename,
  cimpl_img * outRed, cimpl_img * outGreen,
  cimpl_img * outBlue, cimpl_img * outAlpha );

cimpl_status cimpl_writeJpeg( cimpl_img const * const redImg,
  cimpl_img const * const greenImg, cimpl_img const * const blueImg,
  int quality, char const * const outFilename );

cimpl_status cimpl_writePng( cimpl_img const * const redImg,
  cimpl_img const * const greenImg, cimpl_img const * const blueImg,
  cimpl_img const * const alphaImg, char const * const outFilename );

#ifdef __cplusplus
}
#endif

#endif /* cimpl_png_h */
