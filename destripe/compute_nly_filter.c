//!;+
//!; NAME:
//!;       COMPUTE_NLY_FILTER
//!;
//!; PURPOSE:
//!;       This function applies a unidirectional non local filter (Yaroslavsky) to the input image  
//!;       
//!; CALLING SEQUENCE:
//!;      COMPUTE_NLY_FILTER(input_image, output_image, binary_M, NEdT, number_detectors, sx, sy)
//!;      
//!; INPUTS:
//!;       input_image:       input image 
//!;       binary_M:          binary matrix
//!;       NEdT:              Noise equivalent difference temperature to be used in the non local filtering NLY
//!;       number_detectors:  number of detectors (10 for MODIS and 16 for VIIRS) used here to define the size of the NLY filter   
//!;       sx, sy :           size of input_image, output_image and binary_M
//!;                                 
//!; OUPUTS:
//!;       output_image: image processed with the unidirectional non local filter (Yaroslavsky)
//!;    
//!; HISTORY:
//!;   2012-03-01 - Marouan Bouali
//!;   2014-02-13 - Karlis Mikelsons (translation into C)
#include "math.h"

void compute_nly_filter(float ** input_image, float ** output_image, int ** binary_M, float NEdT, int number_detectors, int sx, int sy )
{

  int filter_size, ix, iy, iymin, iymax, iybuf;
  double central_pixel, sum_weights, sum_weightximg, weight;
  double ave_diff, std_diff, NEdT_min = 0.0001;
  int num_diff;

  // Derive the size of the NLY filter from the number of detectors
  filter_size = (number_detectors/2);

  if(NEdT < 0.0){
  // if negative NEdT is supplied, the minus sign acts as a switch to 
  // perform an automatic evaluation of optimal value of NEdT to be used 
  // in the nonlinear filter. The absolute value of NEdT then acts as 
  // an upper bound to ensure reasonable results if some abnormal cases 
  // when automatic evaluation breaks down.

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
    }
  }
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
        std_diff += (double) (input_image[iybuf][ix]-central_pixel-ave_diff)*(input_image[iybuf][ix]-central_pixel-ave_diff);
        num_diff++;
      }
    }
  }
  if(num_diff>0) std_diff = sqrt(fabs(std_diff/num_diff));
  
  NEdT = fabs(NEdT);
  if(NEdT>2.0*(4.0*std_diff)*(4.0*std_diff)) NEdT = 2.0*(4.0*std_diff)*(4.0*std_diff);
  printf("nly_filter: ave_diff = %e  std_diff = %e  NEdT = %e\n", ave_diff, std_diff, NEdT);

  } //  if(NEdT < 0.0)

  // make sure NEdT is not too small (since we need to divide by it)
  if(NEdT<NEdT_min) NEdT = NEdT_min;
  

//         ! loop over all image pixels
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
        weight = exp((-(input_image[iybuf][ix]-central_pixel)*(input_image[iybuf][ix]-central_pixel))/NEdT);
        sum_weightximg = sum_weightximg + weight*input_image[iybuf][ix];
        sum_weights    = sum_weights    + weight;
      }
      output_image[iy][ix] = sum_weightximg / sum_weights;
    }
  }

  return;
}
