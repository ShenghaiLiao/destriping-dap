
/////////////////////////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------------------------- //
// This code was developed in NOAA/NESDIS/STAR SST and Ocean Color groups by Karlis Mikelsons, //
// and is partially based on an earlier version of destriping algorithm for use with the       //
// Sea Surface Temperature (SST) data written by Marouan Bouali at NOAA/NESDIS/STAR SST group. //
//                                                                                             //
// Please acknowledge any use of these codes in your presentations                             //
// and publications by citing the following publications:                                      //
//                                                                                             //
// 1. M. Bouali and A. Ignatov, "Adaptive Reduction of Striping                                //
//    for Improved Sea Surface Temperature Imagery                                             //
//    from Suomi National Polar-Orbiting Partnership (S-NPP)                                   //
//    Visible Infrared Imaging Radiometer Suite (VIIRS)",                                      //
//    J. Atmos. Oceanic Technol., 31, 150–163 (2014).                                          //
//                                                                                             //
// 2. K. Mikelsons, M. Wang, L. Jiang, and M. Bouali,                                          //
//    "Destriping algorithm for improved satellite-derived                                     //
//    ocean color product imagery," Opt. Express  22, 28058-28070 (2014).                      //
//                                                                                             //
// Please report any bugs to Karlis.Mikelsons@noaa.gov                                         //
// ------------------------------------------------------------------------------------------- //
/////////////////////////////////////////////////////////////////////////////////////////////////

//!; NAME:
//!;       compute_poisson_rec
//!;
//!; PURPOSE:
//!;       This function integrates a gradient field (dx, dy) using a poisson solver 
//!;       This amounts to 
//!;            1) calculation of Laplacian
//!;            2) forawrd Discrete Cosine Transform
//!;            3) scaling by inverse Laplacian
//!;            4) reverse Discrete Cosine Transform
//!;       
//!; INPUTS:
//!;       input_image[sy][sx]             input image
//!;       binary_M[sy][sx]                matrix that tells which y gradients to restrict
//!;       grid_dct[sy][sx]                meshgrid to be used for the projection in the cosine basis  
//!;       sx, sy                          size of matrices dx, dy, grid_dct, output_image                                
//!;       Lap[sy][sx], Lap_dct[sy][sx]    temporary work arrays
//!;       plan_forward, plan_backward     fftw plans for fast FFT execution
//!;       
//!;                                 
//!; OUPUTS:
//!;       output_image[sy][sx]            image obtained from the integration of the input gradient field 
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!;   2014-02-13 - Karlis Mikelsons (translated to C)

#include "destripe.h"
#include "fftw3.h"

void compute_poisson_rec(float ** input_image, float ** output_image, int ** binary_M, float ** grid_dct, int sx, int sy, 
                         float ** Lap, float ** Lap_dct, fftwf_plan plan_forward, fftwf_plan plan_backward)
{

    int ix, iy;

    //  calculate Laplacian [eq. 3 in ref. 2]
    //  Here, we calculate Laplacian in one sweep from the input image
    //  and matrix M of binary values that is used to restrict the gradient along y direction
    //  The boundary conditions are Neumann kind, meaning the first derivative at the boundary is assumed to be 0
    //  Thus, the first and last rows and columns have to be treated separately

    // first row
    Lap[0][0] = -input_image[0][0] + input_image[0][1];
    if(binary_M[0][0]!=0) Lap[0][0] += input_image[1][0] - input_image[0][0];

    for(ix=1; ix<(sx-1); ix++){
        Lap[0][ix] = input_image[0][ix-1] - 2.0*input_image[0][ix] + input_image[0][ix+1];
        if(binary_M[0][ix]!=0) Lap[0][ix] += input_image[1][ix] - input_image[0][ix];
    }

    Lap[0][sx-1] = input_image[0][sx-2] - input_image[0][sx-1];
    if(binary_M[0][sx-1]!=0) Lap[0][sx-1] += input_image[1][sx-1] - input_image[0][sx-1];

    // middle rows
#pragma omp parallel for
    for(iy=1; iy<(sy-1); iy++){

        Lap[iy][0] = -input_image[iy][0] + input_image[iy][1];
        if(binary_M[iy  ][0]!=0) Lap[iy][0] += input_image[iy+1][0] - input_image[iy][0];
        if(binary_M[iy-1][0]!=0) Lap[iy][0] += input_image[iy-1][0] - input_image[iy][0];

        for(ix=1; ix<(sx-1); ix++){
            Lap[iy][ix] = input_image[iy][ix-1] - 2.0*input_image[iy][ix] + input_image[iy][ix+1]; 
            if(binary_M[iy  ][ix]!=0) Lap[iy][ix] += input_image[iy+1][ix] - input_image[iy][ix];
            if(binary_M[iy-1][ix]!=0) Lap[iy][ix] += input_image[iy-1][ix] - input_image[iy][ix];
        } // ix

        Lap[iy][sx-1] = input_image[iy][sx-2] - input_image[iy][sx-1];
        if(binary_M[iy  ][sx-1]!=0) Lap[iy][sx-1] += input_image[iy+1][sx-1] - input_image[iy][sx-1];
        if(binary_M[iy-1][sx-1]!=0) Lap[iy][sx-1] += input_image[iy-1][sx-1] - input_image[iy][sx-1];
    } // iy

    // last row
    Lap[sy-1][0] = -input_image[sy-1][0] + input_image[sy-1][1];
    if(binary_M[sy-2][0]!=0) Lap[sy-1][0] += input_image[sy-2][0] - input_image[sy-1][0];

    for(ix=1; ix<(sx-1); ix++){
        Lap[sy-1][ix] = input_image[sy-1][ix-1] - 2.0*input_image[sy-1][ix] + input_image[sy-1][ix+1];
        if(binary_M[sy-2][ix]!=0) Lap[sy-1][ix] += input_image[sy-2][ix] - input_image[sy-1][ix];
    }

    Lap[sy-1][sx-1] = input_image[sy-1][sx-2] - input_image[sy-1][sx-1];
    if(binary_M[sy-2][sx-1]!=0) Lap[sy-1][sx-1] += input_image[sy-2][sx-1] - input_image[sy-1][sx-1];

    ////////   done calculating Laplacian  /////////////////////////////////////////////////////////
  
    // now we need to solve Discrete Poisson Equation  [Eq. 4 in ref. 2] 
    // use FFT to do a forward DCT 
    fftwf_execute_r2r(plan_forward, Lap[0], Lap_dct[0]);

    // solve discrete Poisson equation (algabraic equations in Fourier space, Eq. 5 in ref. 2)
    // by scaling Fourier components
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++) { Lap_dct[0][ix] *= grid_dct[0][ix]; }

    // use FFT to do a backward DCT, and get the solution 
    fftwf_execute_r2r(plan_backward, Lap_dct[0], Lap[0]);

    // copy data to output array 
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++) { output_image[0][ix] = Lap[0][ix]; }

    return;
};
