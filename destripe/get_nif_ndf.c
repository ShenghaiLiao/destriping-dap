
/////////////////////////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------------------------- //
// This code was developed in NOAA/NESDIS/STAR SST and Ocean Color groups by Karlis Mikelsons, //
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
//!;       get_nif_ndf
//!;
//!; PURPOSE:
//!;       This function calculates Normalized Improvement Factor (NIF)
//!;       and Normalized Distortion  Factor (NDF) and 
//!;       the size of destriping domain
//!;       
//!; INPUTS:
//!;   oldimg[sy][sx]       original  input  image
//!;   newimg[sy][sx]       destriped output image 
//!;   binary_M[sy][sx]     matrix specifying destriping domain
//!;   nx                   size of input image (horizontal direction)
//!;   ny                   size of input image ( vertical  direction)
//!;                                 
//!; OUPUTS:
//!;   nif                  Normalized Improvement Factor [Eq. 27 in ref. 1]
//!;   ndf                  Normalized Distortion  Factor [Eq. 28 in ref. 1]
//!;   nif                  size of destriping domain in pixels
//!;    
//!; HISTORY:
//!;   2014-09-30 - Karlis Mikelsons
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_nif_ndf(float ** oldimg, float ** newimg, int ** binary_M, int nx, int ny, float * nif, float * ndf, int * nwf){

    int ix, iy;
    float dxold, dxnew, dyold, dynew;
    float sumdiffy = 0.0, sumdiffx = 0.0, sumx = 0.0, sumy = 0.0;

    // assign default return values
    *nif = 0.0;
    *ndf = 1.0;
    *nwf = 0;

    // loop over all image
    for(iy=0; iy<(ny-1); iy++){
        for(ix=0; ix<(nx-1); ix++){

            // exclude points outside destriping domain
            if( (binary_M[iy][ix]!=0) || (binary_M[iy+1][ix]!=0) || (binary_M[iy][ix+1]!=0) )continue;

            // exclude points with zero values both before and after destriping
            if( (fabs(oldimg[iy][ix])<0.0001) &&  (fabs(newimg[iy][ix])<0.0001) ) continue; 

            dxold = oldimg[iy][ix+1] - oldimg[iy][ix];
            dxnew = newimg[iy][ix+1] - newimg[iy][ix];
            dyold = oldimg[iy+1][ix] - oldimg[iy][ix];
            dynew = newimg[iy+1][ix] - newimg[iy][ix];

            sumdiffy += fabs(dyold) - fabs(dynew);
            sumdiffx += fabs(dxold - dxnew);
            sumx += fabs(dxold);
            sumy += fabs(dyold);
            (*nwf)++;
        }
    }

    // printf("%10.5e %10.5e     %10.5e %10.5e\n", sumdiffy, sumy, sumdiffx, sumx);

    if(sumy>0.0) { *nif =       sumdiffy/sumy; }
    if(sumx>0.0) { *ndf = 1.0 - sumdiffx/sumx; }

    return;
};
