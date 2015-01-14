
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
//!;       compute_nly_filter
//!;
//!; PURPOSE:
//!;       This function applies a unidirectional non local filter (Yaroslavsky) to the input image  
//!;       
//!; INPUTS:
//!;       input_image[sy][sx]    input image 
//!;       binary_M[sy][sx]       binary matrix
//!;       NEdT                   nonlinear filter parameter
//!;                              if NEdT < 0.0, then |NEdT| is used as a upper limit 
//!;                              for adaptively determined filter parameter (used for Ocean Color bands)
//!;       number_detectors       number of detectors used to define the size of the NLY filter
//!;                              normally, it should be 10 for MODIS and 16 for VIIRS
//!;                              but can be doubled for bands that have different 
//!;                              sensor mirror reflectances (MODIS)
//!;       sx, sy                 size of input_image, output_image and binary_M
//!;                                 
//!; OUPUTS:
//!;       output_image[sy][sx]   image processed with the unidirectional non local filter (Yaroslavsky)
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!;   2014-02-13 - Karlis Mikelsons (translation into C, optimized and modified)

#include "math.h"

void compute_nly_filter(float ** input_image, float ** output_image, int ** binary_M, float NEdT, int number_detectors, int sx, int sy )
{

    int filter_size, ix, iy, iymin, iymax, iybuf, num_diff;
    double central_pixel, sum_weights, sum_weightximg, weight;
    double ave_diff, std_diff, NEdT_min = 0.0001;

    // Derive the size of the NLY filter from the number of detectors
    filter_size = (number_detectors/2);

    if(NEdT < 0.0){

        // if negative NEdT is supplied, the minus sign acts as a switch to 
        // perform an automatic evaluation of optimal value of NEdT to be used 
        // in the nonlinear filter. The absolute value of NEdT then acts as 
        // an upper bound to ensure reasonable results if some abnormal cases 
        // when automatic evaluation breaks down.
        // This is described by Eqs. 7-8 in ref. 2

        // calculate average difference from central pixel
        ave_diff = 0.0; num_diff = 0;
#pragma omp parallel for private(ix, iy, iymin, iymax, central_pixel) reduction(+:ave_diff, num_diff)
        for(iy=0; iy<sy; iy++){
            for(ix=0;ix<sx; ix++){
                if(binary_M[iy][ix]!=0) continue; 
                central_pixel = input_image[iy][ix];
                iymin = iy-filter_size; if(iymin<0)      iymin = 0;
                iymax = iy+filter_size; if(iymax>(sy-1)) iymax = (sy-1);

                //    ! perform sum of weights and weights*image amplitudes 
                for(iybuf=iymin; iybuf<iymax; iybuf++){
                    ave_diff += (double) (input_image[iybuf][ix]-central_pixel);
                    num_diff++;
                }
            } // ix
        } // iy

        if(num_diff>0) ave_diff /= num_diff;

        // calculate standard deviation of all differences
        std_diff = 0.0; num_diff = 0;
#pragma omp parallel for private(ix, iy, iymin, iymax, central_pixel) reduction(+:std_diff, num_diff)
        for(iy=0; iy<sy; iy++){
            for(ix=0;ix<sx; ix++){
                if(binary_M[iy][ix]!=0) continue; 
                central_pixel = input_image[iy][ix];
                iymin = iy-filter_size; if(iymin<0)      iymin = 0;
                iymax = iy+filter_size; if(iymax>(sy-1)) iymax = (sy-1);

                //    ! perform sum of weights and weights*image amplitudes 
                for(iybuf=iymin; iybuf<iymax; iybuf++){
                    std_diff += (double) (input_image[iybuf][ix]-central_pixel-ave_diff)*
                                         (input_image[iybuf][ix]-central_pixel-ave_diff);
                    num_diff++;
                }
            } // ix
        } // iy

        if(num_diff>0) std_diff = sqrt(fabs(std_diff/num_diff));
  
        NEdT = fabs(NEdT);
        if(NEdT>2.0*(4.0*std_diff)*(4.0*std_diff)) NEdT = 2.0*(4.0*std_diff)*(4.0*std_diff);
        printf("nly_filter: ave_diff = %e  std_diff = %e  NEdT = %e\n", ave_diff, std_diff, NEdT);

    } //  if(NEdT < 0.0)

    // make sure NEdT is not too small (since we need to divide by it)
    if(NEdT<NEdT_min) NEdT = NEdT_min;
  

    // apply nonlinear filter [Eq. 9 in ref. 1 or Eq. 6 in ref. 2]
#pragma omp parallel for private(ix, iy, iymin, iymax, central_pixel, sum_weights, sum_weightximg, weight)
    for(iy=0; iy<sy; iy++){
        for(ix=0;ix<sx; ix++){
            output_image[iy][ix] = input_image[iy][ix]; 
            if(binary_M[iy][ix]!=0) continue; 
            central_pixel = input_image[iy][ix];
            iymin = iy-filter_size; if(iymin<0)      iymin = 0;
            iymax = iy+filter_size; if(iymax>(sy-1)) iymax = (sy-1);

            //    ! perform sum of weights and weights*image amplitudes 
            sum_weights    = 0.0;
            sum_weightximg = 0.0;
            for(iybuf=iymin; iybuf<iymax; iybuf++){
                weight = exp((-(input_image[iybuf][ix]-central_pixel)*
                               (input_image[iybuf][ix]-central_pixel))/NEdT);

                sum_weightximg = sum_weightximg + weight*input_image[iybuf][ix];
                sum_weights    = sum_weights    + weight;
            }
            output_image[iy][ix] = sum_weightximg / sum_weights;
        } // ix
    } // iy

    return;
};
