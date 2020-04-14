
//
//  cimpl_png.c
//  cimpl
//
//  Created by Nicholas Dwork on 8/29/17.
//  Copyright © 2017 Nicholas Dwork. All rights reserved.
//
//  Based on code written by Guillaume Cottenceau distributed under the X11 license 
//  and located at http://zarb.org/~gc/html/libpng.html.
//  Also based on example.c located at https://dev.w3.org/Amaya/libpng/example.c
//
//

#include "cimpl_io.h"

#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <png.h>
#include <jpeglib.h>


cimpl_status cimpl_readPng( char const * const filename,
  cimpl_img * const redImg, cimpl_img * const greenImg,
  cimpl_img * const blueImg, cimpl_img * const alphaImg ){

  FILE* fp;
  cimpl_bool is_alpha;
  cimpl_bool is_gray;
  float *ptr_r, *ptr_g, *ptr_b, *ptr_a;
  int bit_depth, color_type, interlace_type;
  png_bytep* imgData;
  unsigned char header[8];    // 8 is the maximum size that can be checked
  unsigned int * bits_per_pixel;

  assert( redImg->data == NULL );   assert( greenImg->data == NULL );
  assert( blueImg->data == NULL );  assert( alphaImg->data == NULL );

  fp = fopen( filename, "rb" );
  if( fp == NULL ) return FILE_OPEN_ERROR;

  fread( header, 1, 8, fp );
  if( png_sig_cmp(header, 0, 8) ) return FILE_TYPE_ERROR;


  // Setup PNG structures for read
  png_structp png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
  if( !png_ptr ){
    fclose(fp);
    return READ_ERROR;
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if( info_ptr == NULL ){
    fclose(fp);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    return READ_ERROR;
  }

  if( setjmp(png_jmpbuf(png_ptr)) ){
    png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
    fclose(fp);
    return READ_ERROR;
  }

  png_init_io( png_ptr, fp );
  png_set_sig_bytes( png_ptr, 8 );

  // Get PNG Header Info up to data block
  png_read_info(png_ptr,info_ptr);
  png_uint_32 W, H;
  is_gray = cimpl_false;
  //png_get_IHDR(png_ptr,info_ptr,&W,&H,&bit_depth,&color_type,&interlace_type,(int*)0,(int*)0);
  W = png_get_image_width( png_ptr, info_ptr );
  H = png_get_image_height( png_ptr, info_ptr );
  bit_depth = png_get_bit_depth( png_ptr, info_ptr );
  color_type = png_get_color_type( png_ptr, info_ptr );
  interlace_type = png_get_interlace_type( png_ptr, info_ptr );
  if (bits_per_pixel) *bits_per_pixel = (unsigned int) bit_depth;
  
  // Transforms to unify image data
  if (color_type==PNG_COLOR_TYPE_PALETTE) {
    png_set_palette_to_rgb(png_ptr);
    color_type = PNG_COLOR_TYPE_RGB;
    bit_depth = 8;
  }
  if (color_type==PNG_COLOR_TYPE_GRAY && bit_depth<8) {
    png_set_expand_gray_1_2_4_to_8(png_ptr);
    is_gray = cimpl_true;
    bit_depth = 8;
  }
  if (png_get_valid(png_ptr,info_ptr,PNG_INFO_tRNS)) {
    png_set_tRNS_to_alpha(png_ptr);
    color_type |= PNG_COLOR_MASK_ALPHA;
  }
  if (color_type==PNG_COLOR_TYPE_GRAY || color_type==PNG_COLOR_TYPE_GRAY_ALPHA) {
    png_set_gray_to_rgb(png_ptr);
    color_type |= PNG_COLOR_MASK_COLOR;
    is_gray = cimpl_true;
  }
  
  png_read_update_info(png_ptr,info_ptr);

  if (bit_depth!=8 && bit_depth!=16) {
    png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
    fclose(fp);
    return READ_ERROR;
  }
  const int byte_depth = bit_depth >> 3;
  
  // Allocate Memory for Image Read
  imgData = (png_bytep*) malloc( sizeof(png_bytep) * H );
  for( unsigned int row = 0; row<H; ++row ){
    imgData[row] = (png_byte*) malloc( (size_t)byte_depth*4*W );
  }
  png_read_image(png_ptr,imgData);
  png_read_end(png_ptr,info_ptr);
    

  // Read pixel data
  if (color_type!=PNG_COLOR_TYPE_RGB && color_type!=PNG_COLOR_TYPE_RGB_ALPHA) {
    png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
    fclose(fp);
    return READ_ERROR;
  }
  is_alpha = (color_type==PNG_COLOR_TYPE_RGBA);

  redImg->h = H; redImg->w = W; redImg->nPix = H*W;
  cimpl_malloc( sizeof(float)*redImg->nPix,   (void**) &( redImg->data ) );

  if( !is_gray ){
    greenImg->h = H; greenImg->w = W; greenImg->nPix = H*W;
    blueImg->h  = H; blueImg->w  = W; blueImg->nPix  = H*W;
    cimpl_malloc( sizeof(float)*greenImg->nPix, (void**) &( greenImg->data ) );
    cimpl_malloc( sizeof(float)*blueImg->nPix,  (void**) &( blueImg->data ) );
  }

  if( is_alpha ){
    alphaImg->h = H; alphaImg->w = W; alphaImg->nPix = H*W;
    cimpl_malloc( sizeof(float)*alphaImg->nPix,  (void**) &( alphaImg->data ) );
  }

  ptr_r = redImg->data;
  ptr_g = greenImg->data;
  ptr_b = blueImg->data;
  ptr_a = alphaImg->data;
  if (bit_depth == 8) {

    for( int y=0; y<H; y++ ){
      const unsigned char *ptrs = (unsigned char*)( imgData[y] );

      for( int x=0; x<W; x++ ){
        ptr_r[y+x*H] = (float) *(ptrs++);
        if (ptr_g != NULL) ptr_g[y+x*H] = (float) *(ptrs++);
        if (ptr_b != NULL) ptr_b[y+x*H] = (float) *(ptrs++);
        if (ptr_a != NULL) ptr_a[y+x*H] = (float) *(ptrs++);
    } }

  } else {

    for( int y=0; y<H; y++ ){
      const unsigned short *ptrs = (unsigned short*)( imgData[y] );

      for( int x=0; x<W; x++ ){
        *(ptr_r++) = (float)*(ptrs++);
        if (ptr_g) *(ptr_g++) = (float) *(ptrs++);
        if (ptr_b) *(ptr_b++) = (float) *(ptrs++);
        if (ptr_a) *(ptr_a++) = (float) *(ptrs++);
    } }

  }
    
  /* Cleanup */
  for( int y=0; y<H; y++ ) free( imgData[y] );
  free( imgData );
  png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
  fclose(fp);

  return SUCCESS;
}


cimpl_status cimpl_writePng( cimpl_img const * const redImg,
  cimpl_img const * const greenImg, cimpl_img const * const blueImg,
  cimpl_img const * const alphaImg, char const * const filename ){
  // For channels that don't exist, pass in NULL

  FILE* fp;
  int color_type, height, width, byte_depth, bit_depth, numChan;
  unsigned int bytes_per_pixel = 1;
  png_bytep* imgData;
  png_infop info_ptr;
  png_structp png_ptr;

  width = (int) redImg->w; height = (int) redImg->h;

#ifndef NDEBUG
  if( blueImg && blueImg->data ){
    assert( greenImg && greenImg->data );
    assert( greenImg->w == width ); assert( greenImg->h == height );
    assert( blueImg->w == width ); assert( blueImg->h == height );
  }
#endif // NDEBUG

  if( alphaImg && alphaImg->data ){
    assert( alphaImg->w == width ); assert( alphaImg->h == height );
  }

  fp = fopen(filename, "wb");
  if( fp == NULL ) return FILE_OPEN_ERROR;

  png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
  if( png_ptr == NULL ){
    fclose(fp);
    return WRITE_ERROR;
  }

  info_ptr = png_create_info_struct( png_ptr );
  if( info_ptr == NULL ){
    fclose(fp);
    png_destroy_write_struct( &png_ptr, NULL );
    return WRITE_ERROR;
  }

  png_init_io(png_ptr, fp);
  
  bit_depth = bytes_per_pixel ? ( bytes_per_pixel * 8 ) : 8;
  
  if( alphaImg && alphaImg->data ){
    // An alpha channel exists
    if( blueImg->data ){
      color_type = PNG_COLOR_TYPE_RGB_ALPHA;
      numChan = 4;
    } else {
      color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
      numChan = 2;
    }
  } else {
    // No alpha channel exists
    if( blueImg && blueImg->data ){
      color_type = PNG_COLOR_TYPE_RGB;
      numChan = 3;
    } else {
      color_type = PNG_COLOR_TYPE_GRAY;
      numChan = 1;
    }
  }

  int interlace_type = PNG_INTERLACE_NONE;
  int compression_type = PNG_COMPRESSION_TYPE_DEFAULT;
  int filter_method = PNG_FILTER_TYPE_DEFAULT;
  png_set_IHDR( png_ptr, info_ptr, width, height, bit_depth, color_type,
      interlace_type, compression_type, filter_method );
  png_write_info( png_ptr, info_ptr );

  byte_depth = bit_depth >> 3;
  imgData = (png_bytep*) malloc( sizeof(png_bytep*) * height );
  for( unsigned int row = 0; row<height; ++row ){
    imgData[row] = (png_byte*) malloc( (size_t)byte_depth * numChan * width );
  }

  if( bit_depth == 8 ){
    if( alphaImg && alphaImg->data ){
      if( blueImg && blueImg->data ){
        // 8 bit color image with alpha channel
        float* pC0 = redImg->data;
        float* pC1 = greenImg->data;
        float* pC2 = blueImg->data;
        float* pC3 = alphaImg->data;
        for( int row=0; row<height; ++row ){
          unsigned char* ptrd = imgData[row];
          for( int col=0; col<width; ++col ){
            *(ptrd++) = (unsigned char) pC0[ row + col*height ];
            *(ptrd++) = (unsigned char) pC1[ row + col*height ];
            *(ptrd++) = (unsigned char) pC2[ row + col*height ];
            *(ptrd++) = (unsigned char) pC3[ row + col*height ];
          }
        }
      } else {
        // 8 bit grayscale image with alpha channel
        float* pC0 = redImg->data;
        float* pC1 = alphaImg->data;
        for( int row=0; row<height; ++row ){
          unsigned char* ptrd = imgData[row];
          for( int col=0; col<width; ++col ){
            *(ptrd++) = (unsigned char) pC0[ row + col*height ];
            *(ptrd++) = (unsigned char) pC1[ row + col*height ];
          }
        }
      }
    } else {
      if( blueImg && blueImg->data ){
        // 8 bit color image without alpha channel
        float* pC0 = redImg->data;
        float* pC1 = greenImg->data;
        float* pC2 = blueImg->data;
        for( int row=0; row<height; ++row ){
          unsigned char* ptrd = imgData[row];
          for( int col=0; col<width; ++col ){
            *(ptrd++) = (unsigned char) pC0[ row + col*height ];
            *(ptrd++) = (unsigned char) pC1[ row + col*height ];
            *(ptrd++) = (unsigned char) pC2[ row + col*height ];
          }
        }
      } else {
        // 8 bit grayscale image without alpha channel
        float* pC0 = redImg->data;
        for( int row=0; row<height; ++row ){
          unsigned char* ptrd = imgData[row];
          for( int col=0; col<width; ++col ){
            *(ptrd++) = (unsigned char) pC0[ row + col*height ];
          }
        }
      }
    }
  } else if( bit_depth == 16 ) {
    if( alphaImg && alphaImg->data ){
      if( blueImg && blueImg->data ){
        // 8 bit color image with alpha channel
        float* pC0 = redImg->data;
        float* pC1 = greenImg->data;
        float* pC2 = blueImg->data;
        float* pC3 = alphaImg->data;
        for( int row=0; row<height; ++row ){
          unsigned short* ptrd = (unsigned short*) imgData[row];
          for( int col=0; col<width; ++col ){
            *(ptrd++) = (unsigned short)*(pC0++);
            *(ptrd++) = (unsigned short)*(pC1++);
            *(ptrd++) = (unsigned short)*(pC2++);
            *(ptrd++) = (unsigned short)*(pC3++);
          }
        }
      } else {
        // 8 bit grayscale image with alpha channel
        float* pC0 = redImg->data;
        float* pC1 = alphaImg->data;
        for( int row=0; row<height; ++row ){
          unsigned short* ptrd = (unsigned short*) imgData[row];
          for( int col=0; col<width; ++col ){
            *(ptrd++) = (unsigned short)*(pC0++);
            *(ptrd++) = (unsigned short)*(pC1++);
          }
        }
      }
    } else {
      if( blueImg && blueImg->data ){
        // 8 bit color image without alpha channel
        float* pC0 = redImg->data;
        float* pC1 = greenImg->data;
        float* pC2 = blueImg->data;
        for( int row=0; row<height; ++row ){
          unsigned short* ptrd = (unsigned short*) imgData[row];
          for( int col=0; col<width; ++col ){
            *(ptrd++) = (unsigned short)*(pC0++);
            *(ptrd++) = (unsigned short)*(pC1++);
            *(ptrd++) = (unsigned short)*(pC2++);
          }
        }
      } else {
        // 8 bit grayscale image without alpha channel
        float* pC0 = redImg->data;
        for( int row=0; row<height; ++row ){
          unsigned short* ptrd = (unsigned short*) imgData[row];
          for( int col=0; col<width; ++col )
            *(ptrd++) = (unsigned short)*(pC0++);
        }
      }
    }
  }

  png_write_image( png_ptr, imgData );
  png_write_end( png_ptr, info_ptr );
  png_destroy_write_struct( &png_ptr, &info_ptr );
  
  // Deallocate Image Write Memory
  for( int row=0; row<height; row++ ) free( imgData[row] );
  free( imgData );
  fclose(fp);

  return SUCCESS;
}




cimpl_status cimpl_writeJpeg( cimpl_img const * const redImg,
  cimpl_img const * const greenImg, cimpl_img const * const blueImg,
  int quality, char const * const outFilename ){
  // set quality to 0 for default value

  FILE *fp;
  int image_height, image_width, nChannels;

  JSAMPLE *image_buffer;	/* Points to large array of R,G,B-order data */

  image_width = (int) redImg->w;
  image_height = (int) redImg->h;

#ifndef NDEBUG
  if( blueImg && blueImg->data ){
    assert( greenImg && greenImg->data );
    assert( greenImg->w == image_width );  assert( blueImg->w == image_width );
    assert( greenImg->h == image_height ); assert( blueImg->h == image_height );
  }
#endif // NDEBUG

  if( blueImg && blueImg->data ){
    nChannels = 3;
    image_buffer = malloc( sizeof(JSAMPLE) * image_width * image_height * nChannels );
    for( size_t i=0; i<redImg->w; ++i ){
      for( size_t j=0; j<redImg->h; ++j ){
        image_buffer[ 0 + i*nChannels + j*redImg->w*nChannels ] =
          (JSAMPLE) redImg->data[ j+i*redImg->h ];
        image_buffer[ 1 + i*nChannels + j*redImg->w*nChannels ] =
          (JSAMPLE) greenImg->data[ j+i*redImg->h ];
        image_buffer[ 2 + i*nChannels + j*redImg->w*nChannels ] =
          (JSAMPLE) blueImg->data[ j+i*redImg->h ];
    } }
  } else {
    nChannels = 1;
    image_buffer = malloc( sizeof(JSAMPLE) * image_width * image_height );
    for( size_t i=0; i<redImg->w; ++i ){
      for( size_t j=0; j<redImg->h; ++j ){
        image_buffer[i+j*redImg->w] = (JSAMPLE) redImg->data[j+i*redImg->h];
    } }
  }

  fp = fopen( outFilename, "wb" );
  if( fp == NULL ){
    free( image_buffer );
    return FILE_OPEN_ERROR;
  }

  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;

  JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
  int row_stride;		/* physical row width in image buffer */

  /* Step 1: allocate and initialize JPEG compression object */

  cinfo.err = jpeg_std_error( &jerr );
  jpeg_create_compress( &cinfo );
  jpeg_stdio_dest( &cinfo, fp );

  cinfo.image_width = image_width;  /* image width and height, in pixels */
  cinfo.image_height = image_height;
  cinfo.input_components = nChannels;  /* # of color components per pixel */
  if( nChannels == 3 ){
    cinfo.in_color_space = JCS_RGB;
  } else {
    cinfo.in_color_space = JCS_GRAYSCALE;
  }
  jpeg_set_defaults( &cinfo );
  if( quality != 0 ) jpeg_set_quality( &cinfo, quality, TRUE );
  jpeg_start_compress( &cinfo, TRUE );

  row_stride = image_width * cinfo.input_components;
  while( cinfo.next_scanline < cinfo.image_height ){
    row_pointer[0] = &image_buffer[ cinfo.next_scanline * row_stride ];
    (void) jpeg_write_scanlines( &cinfo, row_pointer, 1 );
  }

  jpeg_finish_compress(&cinfo);
  fclose(fp);
  jpeg_destroy_compress(&cinfo);

  free( image_buffer );
  return SUCCESS;
}


cimpl_status cimpl_readJpeg( char const * const filename,
  cimpl_img * outRed, cimpl_img * outGreen, cimpl_img * outBlue ){

  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  FILE* fp;
  JSAMPARRAY buffer;  /* Output row buffer */
  int row_stride;  /* physical row width in output buffer */
  int nChannels;
  size_t W, H, rowIndx;

  assert( outRed->data == NULL );
  assert( outGreen->data == NULL );
  assert( outBlue->data == NULL );
  fp = fopen( filename, "rb" );
  if( fp == NULL ) return FILE_OPEN_ERROR;

  cinfo.err = jpeg_std_error( &jerr );
  jpeg_create_decompress( &cinfo );
  jpeg_stdio_src( &cinfo, fp );
  (void) jpeg_read_header( &cinfo, TRUE );
  (void) jpeg_start_decompress( &cinfo );

  H = (size_t) cinfo.output_height;
  W = (size_t) cinfo.output_width;

  outRed->h = H; outRed->w = W; outRed->nPix = H*W;
  cimpl_malloc( sizeof(float)*outRed->nPix,   (void**) &( outRed->data ) );
  nChannels = cinfo.output_components;
  if( nChannels == 3 ){
    outGreen->h = H; outGreen->w = W; outGreen->nPix = H*W;
    outBlue->h  = H; outBlue->w  = W; outBlue->nPix  = H*W;
    cimpl_malloc( sizeof(float)*outGreen->nPix, (void**) &( outGreen->data ) );
    cimpl_malloc( sizeof(float)*outBlue->nPix,  (void**) &( outBlue->data ) );
  } else if( nChannels == 1 ){
  } else {
    fclose( fp );
    jpeg_destroy_decompress( &cinfo );
    return FILE_TYPE_ERROR;
  }

  row_stride = cinfo.output_width * cinfo.output_components;
  rowIndx = 0;
  while (cinfo.output_scanline < cinfo.output_height) {
    buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    for( size_t i=0; i<W; ++i ){
      if( nChannels == 3){
        outRed->data[rowIndx+i*H]   = (float) (*buffer)[3*i+0];
        outGreen->data[rowIndx+i*H] = (float) (*buffer)[3*i+1];
        outBlue->data[rowIndx+i*H]  = (float) (*buffer)[3*i+2];
      } else {
        outRed->data[rowIndx+i*H] = (float) (*buffer)[i];
      }
    }
    ++rowIndx;
  }

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress( &cinfo );
  fclose( fp );
  return SUCCESS;
}


