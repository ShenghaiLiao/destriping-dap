//!;+
//!; NAME:
//!;       COMPUTE_POISSON_REC_CPU
//!;
//!; PURPOSE:
//!;       This function integrates a gradient field (dx, dy) using a poisson solver (projection in a cosine basis)  
//!;       
//!; CALLING SEQUENCE:
//!;      COMPUTE_POISSON_REC_GPU(dx, dy, grid_dct, outut_image, sx, sy)
//!;      
//!; INPUTS:
//!;       grid_dct = meshgrid to be used for the projection in the cosine basis  
//!;       sx, sy = size of matrices dx, dy, grid_dct, output_image                                
//!;                                 
//!; OUPUTS:
//!;   output_image: image obtained from the integration of the input gradient field 
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!;   2014-02-13 - Karlis Mikelsons (translated to C)

#include "destripe.h"
#include "fftw3.h"

void compute_poisson_rec(float ** input_image, float ** output_image, int ** binary_M, float ** grid_dct, int sx, int sy, 
  float ** Lap, float ** Lap_dct, fftwf_plan plan_forward, fftwf_plan plan_backward){

  int ix, iy;

//          ! calculate Laplacian
//          ! HEre, we calculate Laplacian in one sweep from the input image
//          ! and matrix M of binary values that is used to restrict the gradient along y direction
//          ! The boundary conditions are Neumann kind, meaning the first derivative at the boundary is assumed to be 0
//          ! Thus, the first and last rows and columns have to be treated separately

  Lap[0][0] = -input_image[0][0] + input_image[0][1];
  if(binary_M[0][0]!=0) Lap[0][0] += input_image[1][0] - input_image[0][0];

  for(ix=1; ix<(sx-1); ix++){
    Lap[0][ix] = input_image[0][ix-1] - 2.0*input_image[0][ix] + input_image[0][ix+1];
    if(binary_M[0][ix]!=0) Lap[0][ix] += input_image[1][ix] - input_image[0][ix];
  }

  Lap[0][sx-1] = input_image[0][sx-2] - input_image[0][sx-1];
  if(binary_M[0][sx-1]!=0) Lap[0][sx-1] += input_image[1][sx-1] - input_image[0][sx-1];

#pragma omp parallel for
  for(iy=1; iy<(sy-1); iy++){

    Lap[iy][0] = -input_image[iy][0] + input_image[iy][1];
    if(binary_M[iy  ][0]!=0) Lap[iy][0] += input_image[iy+1][0] - input_image[iy][0];
    if(binary_M[iy-1][0]!=0) Lap[iy][0] += input_image[iy-1][0] - input_image[iy][0];

    for(ix=1; ix<(sx-1); ix++){
      Lap[iy][ix] = input_image[iy][ix-1] - 2.0*input_image[iy][ix] + input_image[iy][ix+1]; 
      if(binary_M[iy  ][ix]!=0) Lap[iy][ix] += input_image[iy+1][ix] - input_image[iy][ix];
      if(binary_M[iy-1][ix]!=0) Lap[iy][ix] += input_image[iy-1][ix] - input_image[iy][ix];

    }

    Lap[iy][sx-1] = input_image[iy][sx-2] - input_image[iy][sx-1];
    if(binary_M[iy  ][sx-1]!=0) Lap[iy][sx-1] += input_image[iy+1][sx-1] - input_image[iy][sx-1];
    if(binary_M[iy-1][sx-1]!=0) Lap[iy][sx-1] += input_image[iy-1][sx-1] - input_image[iy][sx-1];
  }


  Lap[sy-1][0] = -input_image[sy-1][0] + input_image[sy-1][1];
  if(binary_M[sy-2][0]!=0) Lap[sy-1][0] += input_image[sy-2][0] - input_image[sy-1][0];

  for(ix=1; ix<(sx-1); ix++){
    Lap[sy-1][ix] = input_image[sy-1][ix-1] - 2.0*input_image[sy-1][ix] + input_image[sy-1][ix+1];
    if(binary_M[sy-2][ix]!=0) Lap[sy-1][ix] += input_image[sy-2][ix] - input_image[sy-1][ix];
  }

  Lap[sy-1][sx-1] = input_image[sy-1][sx-2] - input_image[sy-1][sx-1];
  if(binary_M[sy-2][sx-1]!=0) Lap[sy-1][sx-1] += input_image[sy-2][sx-1] - input_image[sy-1][sx-1];


          
  // use FFT to do a forward DCT 
  fftwf_execute_r2r(plan_forward, Lap[0], Lap_dct[0]);  /* repeat as needed */

#pragma omp parallel for
  for(ix=0; ix<sx*sy; ix++) { Lap_dct[0][ix] *= grid_dct[0][ix]; }

  // use FFT to do a backward DCT 
  fftwf_execute_r2r(plan_backward, Lap_dct[0], Lap[0]);  /* repeat as needed */

//         ! copy data to output array 
#pragma omp parallel for
  for(ix=0; ix<sx*sy; ix++) { output_image[0][ix] = Lap[0][ix]; }

  return;
};
