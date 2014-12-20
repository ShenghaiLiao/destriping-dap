//!;+
//!; NAME:
//!;       COMPUTE_MESHGRID_DCT
//!;
//!; PURPOSE:
//!;       This function generates the meshgrid to be used for the DCT and IDCT in the Poisson reconstruction.  
//!;       
//!; CALLING SEQUENCE:
//!;      COMPUTE_MESHGRID_DCT(grid_dct, sx, sy)
//!;      
//!; INPUTS:
//!;   sx: number of lines of the image
//!;   sy: number of columns of the image
//!;                                 
//!; OUPUTS:
//!;   grid_dct: meshgrid for the DCT 
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!    2014-02-04 - Karlis Mikelsons (translated to C)

#include "destripe.h"
#include "math.h"


void compute_meshgrid_dct(float ** grid_dct, int sx, int sy)
{

  int ix, iy;
  double w_x[sx], w_y[sy];

  for(ix=0; ix<sx; ix++) { w_x[ix] = 2.0*cos((PI*ix)/sx) - 2.0; }
  for(iy=0; iy<sy; iy++) { w_y[iy] = 2.0*cos((PI*iy)/sy) - 2.0; }

  grid_dct[0][0] = 0.0;
//        ! we need to avoid (ix=0,iy=0) term, since then we may have division by zero
  for(iy=0; iy<1; iy++) {
    for(ix=1; ix<sx; ix++) {
      grid_dct[iy][ix] = (float) (1.0/((w_x[ix] + w_y[iy])*(4.0*sx*sy)));
//            ! 4*sx*sy is normalization parameter that is incorporated 
//            ! in grid_dct for optimization purposes
    }
  }
#pragma omp parallel for
  for(iy=1; iy<sy; iy++) {
    for(ix=0; ix<sx; ix++) {
      grid_dct[iy][ix] = (float) (1.0/((w_x[ix] + w_y[iy])*(4.0*sx*sy)));
//            ! 4*sx*sy is normalization parameter that is incorporated 
//            ! in grid_dct for optimization purposes
    }
  }

  return;
};
