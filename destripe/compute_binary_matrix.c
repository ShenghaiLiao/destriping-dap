//!;+
//!; NAME:
//!;       COMPUTE_BINARY_MATRICE
//!;
//!; PURPOSE:
//!;       This function returns the binary matrice that determines the "homogeneous" spatial domain to use for destriping 
//!;       
//!; CALLING SEQUENCE:
//!;      binary_M = COMPUTE_BINARY_MATRICE(input_image, binary_M, thresh_x, thresh_y, sx, sy)
//!;      
//!; INPUTS:
//!;   input_image:  input image corresponding to brightness temperature
//!;   thresh_x: threshold for the horizontal gradient norm
//!;   thresh_y: threshold for the vertical gradient norm
//!;   Tmin: threshold for the minimum value 
//!;   Tmax: threshold for the maximum value 
//!;   sx:  size of input image and output matrix (horizontal direction)
//!;   sy:  size of input image and output matrix ( vertical  direction)
//!;                                 
//!; OUPUTS:
//!;   binary_M: binary matrice that determines the location of image discontinuities 
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!    2014-02-04 - Karlis Mikelsons (translation to C)
// !

#include "math.h"

void compute_binary_matrix(float ** input_image, int ** binary_M, float thresh_x, float thresh_y, float Tmin, float Tmax, int sx, int sy){

  int ix, iy;
  float gx, gy;
#pragma omp parallel for private(ix, iy, gx, gy)
  for(iy=0; iy<sy; iy++){
    for(ix=0; ix<sx; ix++){
      if(ix<(sx-1)) gx = input_image[iy][ix+1] - input_image[iy][ix]; else gx = 0.0;
      if(iy<(sy-1)) gy = input_image[iy+1][ix] - input_image[iy][ix]; else gy = 0.0;
      if( (fabs(gx)>thresh_x) || (fabs(gy)>thresh_y) || (input_image[iy][ix]<Tmin) ||  (input_image[iy][ix]>Tmax) ) {
        binary_M[iy][ix] = 1;
      } 
      else {
        binary_M[iy][ix] = 0;
      }
    }
  }

  return;
}
