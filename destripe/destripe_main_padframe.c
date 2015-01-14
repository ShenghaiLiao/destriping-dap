
/////////////////////////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------------------------- //
// This code was developed in NOAA/NESDIS/STAR SST and Ocean Color groups by Karlis Mikelsons. //
//                                                                                             //
// Please report any bugs to Karlis.Mikelsons@noaa.gov                                         //
// ------------------------------------------------------------------------------------------- //
/////////////////////////////////////////////////////////////////////////////////////////////////

#include "destripe_main.c"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  a smooth function in [0,1] interval and range [0,1] ///////////////////////////////////////////////////////////////////////////////////
double fsmooth(double x){  
    return 0.5*(1.0 - cos(PI*x)); 
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int get_best_composite_number(int n0, int n1, int inc){
    ///////////////////////////////////////////////////////////////////////////////////////////
    // parse the interval of integers from n0 to n1 with increment inc
    // and find the number whose largest prime factor is the smallest
    // of all the numbers parsed
    ///////////////////////////////////////////////////////////////////////////////////////////

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
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int destripe_main_frame(float ** input_image, float ** output_image, int ** binary_M, int number_detectors, 
                        int number_iterations, float NEdT, float thresh_x, float thresh_y, float Tmin, float Tmax, int sx, int sy)
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // this subroutine acts as a wrapper for  destripe_main(...) subroutine, 
    // it passes all arguments intact except for the original image,
    // which gets padded on all four sides with extra layers of interpolated pixels.
    // The purpose is to speed up Fourier transform by using good composite numbers as 
    // the size of the two - dimensional arrays, as well as to minimize any boundary effects.
    // After getting destriped image of the padded image, this subroutine then
    // extracts the image of the original dimentsions (by discarding the padded frame),
    // as well as corresponding binary matrix 
    // All input/output parameters are same as for destripe_main
    ///////////////////////////////////////////////////////////////////////////////////////////

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

    // destripe the padded image
    ierr = destripe_main(input_image1, output_image1, binary_M1, number_detectors,
                         number_iterations, NEdT, thresh_x, thresh_y, Tmin, Tmax, sx1, sy1);

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
};
