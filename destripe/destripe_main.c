//!;+
//!; NAME:
//!;       DESTRIPE_MAIN_ACSPO
//!;
//!; PURPOSE:
//!;       This function is the main destriping routine. 
//!;       It applies the Directional Hierarchical Decomposition (DHD), the non local filtering (NLY) and returns the destriped image 
//!;       
//!; CALLING SEQUENCE:
//!;      output = DESTRIPE_MAIN_ACSPO(input_image, number_iterations, binary_M, NEdT, number_detectors )
//!;      
//!; INPUTS:
//!;   input_image:  Input image to be destriped. The image should represent brightness temperatures
//!;   number_iterations: number of iterations to be used in the Directional Hierarchical Decomposition (DHD)
//!;   binary_M: binary matrice that defines the spatial domain to be destriped
//!;   NEdT: Noise equivalent difference temperature to be used in the non local filtering NLY
//!;   number_detectors: number of detectors used in the instrument (10 for MODIS and 16 for VIIRS)
//!;                                 
//!; OUPUTS:
//!;   output_image: destriped image 
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!;   2014-02-13 - Karlis Mikelsons (translated to FORTRAN)

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
int destripe_main(float ** input_image, float ** output_image, int ** binary_M,  int number_detectors, int number_iterations, float NEdT, float thresh_x, float thresh_y, float Tmin, float Tmax, int sx, int sy)
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

//          !;;;;;;;;;;; Generate meshgrid for Poisson reconstruction
  compute_meshgrid_dct(grid_dct, sx, sy);
  compute_binary_matrix(input_image, binary_M, thresh_x, thresh_y, Tmin, Tmax, sx, sy);

  int i = fftwf_init_threads();
  if(i==0) { printf("Cannot init threads\n"); return -1; }
  int nthreads = omp_get_max_threads();
  fftwf_plan_with_nthreads(nthreads);
  printf("fftw3 initialized on %i threads\n", nthreads);

//          ! create plans for fftw3 dct forward and backward transformations
  plan_forward  = fftwf_plan_r2r_2d(sy, sx, Lap[0], Lap_dct[0], FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
  plan_backward = fftwf_plan_r2r_2d(sy, sx, Lap_dct[0], Lap[0], FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);

//          ! set the temporary buffer equal to input image
//          ! set rx equal to zero
#pragma omp parallel for
  for(ix=0; ix<sx*sy; ix++) { image_tmp[0][ix] = input_image[0][ix]; rx[0][ix] = 0.0; }


//          !;;;;;;;;;;;;;;;;;;; START ITERATIVE SCHEME (TNV) ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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

  // !; Start non local filtering only after the last DHD iteration
  compute_nly_filter(image_tmp, output_image, binary_M, NEdT, number_detectors, sx, sy);

  //!; Adjust mean value (unknown constant from gradient field integration
#pragma omp parallel for
  for(ix=0; ix<sx*sy; ix++) {   output_image[0][ix] += rx[0][ix]; }
#pragma omp parallel for
  for(ix=0; ix<sx*sy; ix++) {   rec_tmp[0][ix]      = output_image[0][ix] - input_image[0][ix]; }
  
  sumimg = 0.0;
  for(ix=0; ix<sx*sy; ix++) {  sumimg += rec_tmp[0][ix]; }
  sumimg /= (sx*sy);
  printf("Average correction = %e\n", sumimg);
#pragma omp parallel for
  for(ix=0; ix<sx*sy; ix++) {  output_image[0][ix] -= sumimg; } 

  
  //!; Account for missing lines due to bowtie deletion module (fill_in values)
  max_input_image = input_image[0][0]; 
  for(ix=0; ix<sx*sy; ix++) { if(input_image[0][ix]>max_input_image) max_input_image = input_image[0][ix]; }
#pragma omp parallel for
  for(ix=0; ix<sx*sy; ix++) { if(fabs(max_input_image-input_image[0][ix])<0.0001) output_image[0][ix] = input_image[0][ix]; } 
 
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



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double fsmooth(double x){  return 0.5*(1.0 - cos(PI*x)); }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// parse the interval of integers from n0 to n1 with increment inc
// and find the number whose largest prime factor is the smallest
// of all the numbers parsed
int get_best_composite_number(int n0, int n1, int inc){
  int i, j, i1, j1, best_num = n0;
  int largest_prime = n0;
  
  for(i=n0; i<n1; i+=inc){
    // factorize i
    i1 = i; 
    for(j=2; j*j<i; j++){
      while( (i1%j)==0 ){ i1 /= j; }
      if(i1==1) break;
    }
    if(i1>1) j1 = i1; else j1 = j;
    // now j1 is the largest prime factor of i
    // update largest prime and "best number" if j1 is smaller
    if(largest_prime>j1) { largest_prime = j1; best_num = i; }
  }
  return best_num;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int destripe_main_frame(float ** input_image, float ** output_image, int ** binary_M, int number_detectors, int number_iterations, float NEdT, float thresh_x, float thresh_y, float Tmin, float Tmax, int sx, int sy)
{

  int sx1 = sx, sy1 = sy;
  int ix, iy, dx, dy, ierr;

  float ** input_image1, ** output_image1, aveimg, r1, r2;
  int  ** binary_M1;
/*
  if(sx1==1354) sx1 = 1440;
  if(sx1==3200) sx1 = 3200;
  
  if(sy1==2030) sy1 = 2048;
  if(sy1==5392) sy1 = 5488;
  if(sy1==5408) sy1 = 5488;
*/
  sx1 = get_best_composite_number( sx + 16, (int) (sx*1.1), 2 );
  sy1 = get_best_composite_number( sy + 16, (int) (sy*1.1), 2 );

  dx = (sx1 - sx)/2;
  dy = (sy1 - sy)/2;

  binary_M1     = allocate_2d_i(sy1, sx1);
  input_image1  = allocate_2d_f(sy1, sx1);
  output_image1 = allocate_2d_f(sy1, sx1);
  if( (binary_M1==NULL) || (input_image1==NULL) || (output_image1==NULL) ){
    printf("ERROR: cannot allocate memory\n");
    return -1; 
  }

  aveimg = 0.0;
  for(ix=0; ix<sx; ix++) { aveimg += input_image[0][ix] + input_image[sy-1][ix]; }
  for(iy=0; iy<sy; iy++) { aveimg += input_image[iy][0] + input_image[iy][sx-1]; }
  aveimg /= 2.0*(sx+sy);

  // top, bottom
  for(iy=0; iy<dy; iy++){
    r1 = fsmooth(iy/(1.0*dy));
    for(ix=0; ix<sx; ix++){
      input_image1[    iy  ][ix+dx] = r1*input_image[0   ][ix] + (1.0-r1)*aveimg;
      input_image1[sy1-iy-1][ix+dx] = r1*input_image[sy-1][ix] + (1.0-r1)*aveimg;
    }
  }

  // sides
  for(ix=0; ix<dx; ix++){
    r1 = fsmooth(ix/(1.0*dx));
    for(iy=0; iy<sy; iy++){
      input_image1[iy+dy][    ix  ] = r1*input_image[iy][0   ] + (1.0-r1)*aveimg;
      input_image1[iy+dy][sx1-ix-1] = r1*input_image[iy][sx-1] + (1.0-r1)*aveimg;
    }
  }

  // corners
  for(iy=0; iy<dy; iy++){
    r1 = fsmooth(iy/(1.0*dy));
    for(ix=0; ix<dx; ix++){
      r2 = fsmooth(ix/(1.0*dx));
      input_image1[    iy  ][    ix  ] = r1*r2*input_image[   0][   0] + (1.0-r1*r2)*aveimg;
      input_image1[sy1-iy-1][    ix  ] = r1*r2*input_image[sy-1][   0] + (1.0-r1*r2)*aveimg;
      input_image1[    iy  ][sx1-ix-1] = r1*r2*input_image[   0][sx-1] + (1.0-r1*r2)*aveimg;
      input_image1[sy1-iy-1][sx1-ix-1] = r1*r2*input_image[sy-1][sx-1] + (1.0-r1*r2)*aveimg;
    }
  }

  // center
#pragma omp parallel for
  for(iy=0; iy<sy; iy++) {
    for(ix=0; ix<sx; ix++) {
      input_image1[iy+dy][ix+dx] = input_image[iy][ix];
    }
  }

  ierr = destripe_main(input_image1, output_image1, binary_M1, number_detectors, number_iterations, NEdT, thresh_x, thresh_y, Tmin, Tmax, sx1, sy1);

  // restore values to the original arrays
#pragma omp parallel for
  for(iy=0; iy<sy; iy++) {
    for(ix=0; ix<sx; ix++) {
      output_image[iy][ix] = output_image1[iy+dy][ix+dx];
      binary_M[iy][ix] = binary_M1[iy+dy][ix+dx];
    }
  }

  free(input_image1[0]);  free(input_image1);
  free(output_image1[0]); free(output_image1);
  free(binary_M1[0]);     free(binary_M1);

  return ierr;
}
