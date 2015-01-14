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
//!;       compute_meshgrid_dct
//!;
//!; PURPOSE:
//!;       This function generates the two dimensional matrix that 
//!;       contains scaled inverse Laplacian in Fourier space. 
//!;       The scaling takes care of factors introduced by forward and
//!;       backward Discrete Cosine Transforms.
//!;       It is used in the Poisson reconstruction step of the algorithm.  
//!;       See Eq. 26 in ref. 1 and Eq. 5 in ref. 2
//!;       
//!; INPUTS:
//!;   sx: number of lines of the image
//!;   sy: number of columns of the image
//!;                                 
//!; OUPUTS:
//!;   grid_dct[sy][sx]: meshgrid for the DCT 
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!    2014-02-04 - Karlis Mikelsons (translated to C, optimized and modified)

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
