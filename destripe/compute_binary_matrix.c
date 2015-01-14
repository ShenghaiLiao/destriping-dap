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
//!;       compute_binary_matrix
//!;
//!; PURPOSE:
//!;       This function returns the binary matrix that determines the
//!;       "homogeneous" spatial domain to use for destriping 
//!;       See Eq. 3 in ref. 1 and Eq. 2 in ref. 2
//!;       
//!; INPUTS:
//!;   input_image[sy][sx]   input image corresponding to brightness temperature
//!;   thresh_x              threshold for the horizontal gradient norm
//!;   thresh_y              threshold for the vertical gradient norm
//!;   Tmin                  threshold for the minimum value 
//!;   Tmax                  threshold for the maximum value 
//!;   sx                    size of input image and output matrix (horizontal direction)
//!;   sy                    size of input image and output matrix ( vertical  direction)
//!;                                 
//!; OUPUTS:
//!;   binary_M[sy][sx]      binary matrice that determines the location of image discontinuities 
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!    2014-02-04 - Karlis Mikelsons (translation to C)

#include "math.h"

void compute_binary_matrix(float ** input_image, int ** binary_M, float thresh_x, float thresh_y, float Tmin, float Tmax, int sx, int sy){

    int ix, iy;
    float gx, gy;
#pragma omp parallel for private(ix, iy, gx, gy)
    for(iy=0; iy<sy; iy++){
        for(ix=0; ix<sx; ix++){
            if(ix<(sx-1)) gx = input_image[iy][ix+1] - input_image[iy][ix]; else gx = 0.0;
            if(iy<(sy-1)) gy = input_image[iy+1][ix] - input_image[iy][ix]; else gy = 0.0;
            if( (fabs(gx)>thresh_x) || (fabs(gy)>thresh_y) ||
                (input_image[iy][ix]<Tmin) ||  (input_image[iy][ix]>Tmax) ) 
            {
                binary_M[iy][ix] = 1;
            } 
            else {
                binary_M[iy][ix] = 0;
            }
        } // ix
    } // iy 

    return;
}
