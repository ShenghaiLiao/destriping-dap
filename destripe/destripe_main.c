
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
//!;       destripe_main
//!;
//!; PURPOSE:
//!;       This function is the main destriping routine. 
//!;       It applies the Directional Hierarchical Decomposition (DHD), 
//!;       the non local filtering (NLY) and returns the destriped image 
//!;       See Fig. 5 in ref. 2 for a flowchart of this function
//!;       
//!; INPUTS:
//!;   input_image[sy][sx]  Input image to be destriped. The image should represent brightness temperatures
//!;   number_iterations    number of iterations to be used in the Directional Hierarchical Decomposition (DHD)
//!;   NEdT                 nonlinear filter parameter used in non local filtering
//!;   thresh_x, thresh_y   gradient threshold parameters for x, y gradients
//!;   Tmin, Tmax           treshold paraeters for min, max values of image subject to destriping
//!;   sx, sy               size of image in x, y directions
//!;   number_detectors     number of detectors used in the instrument (10 for MODIS and 16 for VIIRS)
//!;                                 
//!; OUPUTS:
//!;   output_image[sy][sx] destriped image 
//!;   binary_M[sy][sx]     binary matrice that defines the spatial domain to be destriped
//!;
//!; RETURN VALUE:
//!;   on success, returns 0; other values indicate failure
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!;   2014-02-13 - Karlis Mikelsons (translated to FORTRAN, updated and optimized)

#include "fftw3.h"
#include "omp.h"
#include "destripe.h"
#include "allocate_2d.c"
#include "compute_meshgrid_dct.c"
#include "compute_binary_matrix.c"
#include "compute_poisson_rec.c"
#include "compute_nly_filter.c"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int destripe_main(float ** input_image, float ** output_image, int ** binary_M,  int number_detectors,
                   int number_iterations, float NEdT, float thresh_x, float thresh_y, float Tmin, float Tmax, int sx, int sy)
{

    int ix, iy, k;
    float max_input_image;
    double sumimg;
    float ** grid_dct, ** rec_tmp, **rx, **image_tmp, **Lap, **Lap_dct;
    fftwf_plan plan_forward, plan_backward;

    // allocate memory
    grid_dct  = allocate_2d_f(sy, sx);
    rec_tmp   = allocate_2d_f(sy, sx);
    rx        = allocate_2d_f(sy, sx);
    image_tmp = allocate_2d_f(sy, sx);
    Lap       = allocate_2d_f(sy, sx);
    Lap_dct   = allocate_2d_f(sy, sx);

    if( (grid_dct==NULL) || (rec_tmp==NULL) || (rx==NULL) || (image_tmp==NULL) || (Lap==NULL) || (Lap_dct==NULL) ){
        printf("ERROR: destripe_main cannot allocate memory\n");
        return -1;
    }

    // calculate inverse Laplacian in Fourier space 
    compute_meshgrid_dct(grid_dct, sx, sy);

    // find destriping domain [see Eq. 3 in ref. 1 and Eq. 2 in ref. 2]
    compute_binary_matrix(input_image, binary_M, thresh_x, thresh_y, Tmin, Tmax, sx, sy);

    // initialize threads for FFTW3 library
    int i = fftwf_init_threads();
    if(i==0) { printf("ERROR: Cannot init fftw threads\n"); return -1; }
    int nthreads = omp_get_max_threads();
    fftwf_plan_with_nthreads(nthreads);
    printf("fftw3 initialized on %i threads\n", nthreads);

    // create plans for fftw3 dct forward and backward transformations
    plan_forward  = fftwf_plan_r2r_2d(sy, sx, Lap[0], Lap_dct[0], FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    plan_backward = fftwf_plan_r2r_2d(sy, sx, Lap_dct[0], Lap[0], FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);

    // set the temporary buffer equal to input image
    // set rx equal to zero
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++) { 
        image_tmp[0][ix] = input_image[0][ix]; 
        rx[0][ix] = 0.0; 
    }


    // perform iterative image decomposition
    for(k=0; k<number_iterations; k++){
  
        printf("Iteration number: %i\n", k);

        compute_poisson_rec(image_tmp, rec_tmp, binary_M, grid_dct, sx, sy, Lap, Lap_dct, plan_forward, plan_backward);

        //!; Estimate vertical reconstruction using horizontal integration
#pragma omp parallel for
        for(ix=0; ix<sx*sy; ix++){  image_tmp[0][ix] -= rec_tmp[0][ix]; }
        //         !; update horizontal reconstruction 
        //         !; Add result to the previous reconstruction to increase ID
#pragma omp parallel for
        for(ix=0; ix<sx*sy; ix++){         rx[0][ix] += rec_tmp[0][ix]; }

    }

    // Start non local filtering only after the last DHD iteration
    compute_nly_filter(image_tmp, output_image, binary_M, NEdT, number_detectors, sx, sy);

    // Adjust mean value (unknown constant from gradient field integration)
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++) {  output_image[0][ix] += rx[0][ix]; }
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++) {  rec_tmp[0][ix]       = output_image[0][ix] - input_image[0][ix]; }
  
    sumimg = 0.0;
    for(ix=0; ix<sx*sy; ix++) {  sumimg += rec_tmp[0][ix]; }
    sumimg /= (sx*sy);
    printf("Average correction = %e\n", sumimg);
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++) {  output_image[0][ix] -= sumimg; } 

  
    // Account for missing lines due to bowtie deletion module (fill_in values)
    max_input_image = input_image[0][0]; 
    for(ix=0; ix<sx*sy; ix++) {
         if(input_image[0][ix]>max_input_image)  
              max_input_image = input_image[0][ix]; 
    }

#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++) { 
        if(fabs(max_input_image-input_image[0][ix])<0.0001) 
            output_image[0][ix] = input_image[0][ix]; 
    } 
 
    // clean up memory
    free(grid_dct[0]);  free(grid_dct);
    free(image_tmp[0]); free(image_tmp);
    free(rec_tmp[0]);   free(rec_tmp);
    free(rx[0]);        free(rx);
    free(Lap[0]);       free(Lap);
    free(Lap_dct[0]);   free(Lap_dct);

    fftwf_destroy_plan(plan_forward );
    fftwf_destroy_plan(plan_backward);
    fftwf_cleanup_threads();

    return 0;
};

