
/////////////////////////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------------------------- //
// This code was developed in NOAA/NESDIS/STAR SST and Ocean Color groups by Karlis Mikelsons. //
//                                                                                             //
// Please report any bugs to Karlis.Mikelsons@noaa.gov                                         //
// ------------------------------------------------------------------------------------------- //
/////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int fill_viirs_bowtie(float ** img, float ** mask, int nx, int ny, float minval, float maxval, float fillval){
    int ix, iy, ntot = 0, ndet = 16;

// #pragma omp parallel for 
    for(ix=0; ix<nx*ny; ix++) mask[0][ix] = fillval;
  
    // this routine is custom made for VIIRS with 3200 pixels per line
    if(nx!=3200) { return -1; }
 
    // fill rows 1 and (ny-2)
    for(ix=0; ix<nx; ix++){
        if( (ix >= 640) && (ix < 2560) ) continue;
        if( (img[2][ix] > minval) && (img[2][ix] < maxval) ) {
            mask[1][ix] = img[1][ix];
            img[1][ix] = img[2][ix];
            ntot++;
        } 
        if( (img[ny-3][ix] > minval) && (img[ny-3][ix] < maxval) ) {
            mask[ny-2][ix] = img[ny-2][ix];
            img[ny-2][ix] = img[ny-3][ix];
            ntot++;
        }
    }

    // fill rows 0 and (ny-1)
    for(ix=0; ix<nx; ix++){
        if( (ix >= 1008) && (ix < 2192) ) continue;
        if( (img[2][ix] > minval) && (img[2][ix] < maxval)  ) {
            mask[0][ix] = img[0][ix];
            img[0][ix] = img[2][ix];
            ntot++;
        }
        if( (img[ny-3][ix] > minval) && (img[ny-3][ix] < maxval) ) {
            mask[ny-1][ix] = img[ny-1][ix];
            img[ny-1][ix] = img[ny-3][ix];
            ntot++;
        }
    }
  
//#pragma omp parallel for 
    for(iy=16; iy<(ny-3); iy+=ndet) {
        for(ix=0; ix<nx; ix++){
            if( (ix>=1008) && (ix<2192) ) continue;
            if( (img[iy-3][ix] <= minval) || (img[iy-3][ix] >= maxval) ||
                (img[iy+2][ix] <= minval) || (img[iy+2][ix] >= maxval) ) continue;
            mask[iy][ix] = img[iy][ix];
            img[iy][ix] = 0.4*img[iy-3][ix] + 0.6*img[iy+2][ix];
            ntot++;
        }
    }

// #pragma omp parallel for 
    for(iy=17; iy<(ny-2); iy+=ndet) {
        for(ix=0; ix<nx; ix++){
            if( (ix>= 640) && (ix<2560) ) continue;
            if( (img[iy-4][ix] <= minval) || (img[iy-4][ix] >= maxval) ||
                (img[iy+1][ix] <= minval) || (img[iy+1][ix] >= maxval) ) continue;
            mask[iy][ix] = img[iy][ix];
            img[iy][ix] = 0.2*img[iy-4][ix] + 0.8*img[iy+1][ix];
            ntot++;
        }
    }
  
// #pragma omp parallel for 
    for(iy=14; iy<(ny-5); iy+=ndet) {
        for(ix=0; ix<nx; ix++){
            if( (ix>= 640) && (ix<2560) ) continue;
            if( (img[iy-1][ix] <= minval) || (img[iy-1][ix] >= maxval) ||
                (img[iy+4][ix] <= minval) || (img[iy+4][ix] >= maxval) ) continue;
            mask[iy][ix] = img[iy][ix];
            img[iy][ix] = 0.8*img[iy-1][ix] + 0.2*img[iy+4][ix];
            ntot++;
        }
    }
  
// #pragma omp parallel for 
    for(iy=15; iy<(ny-4); iy+=ndet) {
        for(ix=0; ix<nx; ix++){
            if( (ix>=1008) && (ix<2192) ) continue;
            if( (img[iy-2][ix] <= minval) || (img[iy-2][ix] >= maxval) ||
                (img[iy+3][ix] <= minval) || (img[iy+3][ix] >= maxval) ) continue;
            mask[iy][ix] = img[iy][ix];
            img[iy][ix] = 0.6*img[iy-2][ix] + 0.4*img[iy+3][ix];
            ntot++;
        }
    }

    // print how many pixels were changed
    // printf("fill_viirs_bowtie done ntot = %i\n", ntot);
 
/*
    // how many pixels have not been changed
    ntot = 0;
    for(ix=0; ix<nx*ny; ix++) { 
        if(fabs((mask[0][ix]) + 100.0) < 0.01) { 
            ntot++; 
        } 
    }
    printf("fill_viirs_bowtie done ntot = %i\n", ntot);

*/

    return 0;
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int restore_viirs_bowtie(float ** img, float ** mask, int nx, int ny, float fillval){
    int ix, ntot = 0;

    // this routine is custom made for VIIRS with 3200 pixels per line
    if(nx!=3200) { return -1; }

// #pragma omp parallel for 
    for(ix=0; ix<nx*ny; ix++) { 
        if(fabs((mask[0][ix])-fillval) > 0.001) { 
            img[0][ix] = mask[0][ix]; 
            ntot++; 
        } 
    }

    // printf("restore_viirs_bowtie done ntot = %i\n", ntot);
    return 0;
};


