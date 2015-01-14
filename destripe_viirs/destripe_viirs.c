
/////////////////////////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------------------------- //
// This code was developed in NOAA/NESDIS/STAR SST group by Karlis Mikelsons,                  //
// and is partially based on an earlier version of destriping algorithm for use with the       //
// Sea Surface Temperature (SST) data written by Marouan Bouali at NOAA/NESDIS/STAR SST group. //
//                                                                                             //
// Please acknowledge any use of these codes in your presentations                             //
// and publications by citing the following publication:                                       //
//                                                                                             //
// 1. M. Bouali and A. Ignatov, "Adaptive Reduction of Striping                                //
//    for Improved Sea Surface Temperature Imagery                                             //
//    from Suomi National Polar-Orbiting Partnership (S-NPP)                                   //
//    Visible Infrared Imaging Radiometer Suite (VIIRS)",                                      //
//    J. Atmos. Oceanic Technol., 31, 150–163 (2014).                                          //
//                                                                                             //
// Please report any bugs to Karlis.Mikelsons@noaa.gov                                         //
// ------------------------------------------------------------------------------------------- //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include "omp.h"
#include "readwrite_viirs.c"
#include "destripe.h"
#include "allocate_2d.c"
#include "get_nif_ndf.c"

////////////////////////////////////
// #include "resample_viirs.h"
////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {

    unsigned short *buffer1  = NULL;
    float          *bufferf1 = NULL;
    float          *bufferf2 = NULL;
    float          *bufferf3 = NULL;
    unsigned long long dims1[32];
    int status, i, j, is;
    float scale1, offset1;
    int ix, iy, sx, sy;
    double r1, scale, offset, lambda;

    int Ndet, Niter, nthreads;
    int Ndet_arr[40], Niter_arr[40], isband[40];
    float Qmin, Qmax, Thresh_x, Thresh_y, NEdQ, scalefact = 1000.0;
    float Qplotmin_chlor, Qplotmax_chlor, Qplotmin_kd490, Qplotmax_kd490, Qplotmin[40], Qplotmax[40];
    float Qmin_arr[40], Qmax_arr[40], Tx_arr[40], Ty_arr[40], NEdQ_arr[40];

    float ** img_in, **img_out, **lat, **lon;
    int   ** binary_M;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // make sure we have at least two command line arguments
    //////////////////////////////////////////////////////////////////////////////////////////////////
    if(argc<3) { 
        printf("Not enough arguments!\nUsage:\n");
        printf(" %s viirs_h5_file destriping_parameter_file_viirs.txt\n", argv[0]);
        return -9;
    }

    // echo command line
    printf("destripe_viirs %s %s\n", argv[1], argv[2]);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // set the number of OpenMP threads  
    //////////////////////////////////////////////////////////////////////////////////////////////////
    nthreads = omp_get_max_threads();
    omp_set_num_threads(nthreads); 
    printf("numthreads = %i maxthreads =  %i\n", omp_get_num_threads(), omp_get_max_threads());

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // open parameter file
    //////////////////////////////////////////////////////////////////////////////////////////////////
    FILE *fp = fopen(argv[2],"r");
    if(fp==NULL) { 
        printf("ERROR: Cannot open param. file\n");
        return -8;
    }
 
    // read parameters
    for(i=0; i<16; i++){ isband[i] = 0; } // initially set all band parameters as "absent"
    for(i=0; i<16; i++){                  // read at most 16 lines of parameters
        j = fscanf(fp,"%i ", &is);        // read band number
        if( (j!=1) || (is<1) || (is>16) ) break; // if band number not read, or if outside range, break

        // read destriping parameters for and is
        j = fscanf(fp,"%i %i %f %f %f %f %f\n", &(Ndet_arr[is]), &(Niter_arr[is]), &(NEdQ_arr[is]),
                                      &(Tx_arr[is]), &(Ty_arr[is]), &(Qmin_arr[is]), &(Qmax_arr[is]));

        if(j!=7) break;                   // if did not read all parameters, break

        // echo read parameters
        printf("%i %i %i %f %f %f %f %f\n", is, (Ndet_arr[is]), (Niter_arr[is]), (NEdQ_arr[is]),
                                       (Tx_arr[is]), (Ty_arr[is]), (Qmin_arr[is]), (Qmax_arr[is]));

        isband[is] = 1;                   // set band parameters as "present" for this band
    }

    // close param. file
    fclose(fp);
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // done reading parameters
    //////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////
    // extract the name of band from the file
    //////////////////////////////////////////////////////////////////////////////////////////////////
    is = 0; 
    for(i=(strlen(argv[1])-6); i>=0; i--) {

        // look for sequence "SVM", then following two chars give band number
        if( (argv[1][i]=='S') && (argv[1][i+1]=='V') && (argv[1][i+2]=='M') ){
            is = (argv[1][i+3]-'0')*10 + (argv[1][i+4]-'0');
            break;
        }
    }
    printf("Band = %i\n", is);

    // check that the band number extracted from file name is within the valid range
    if((is<1)||(is>16)) {
        printf("ERROR: Invalid band\n");   
        return -7;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // make sure we have parameters for destriping this band
    //////////////////////////////////////////////////////////////////////////////////////////////////
    if(isband[is]!=1) {
        printf("No destriping parameters for this band\n");
        return 0;
    }

    // print parameters for the band to be destriped
    printf("%i %i %i %f %f %f %f %f\n", is, (Ndet_arr[is]), (Niter_arr[is]), (NEdQ_arr[is]), 
                                       (Tx_arr[is]), (Ty_arr[is]), (Qmin_arr[is]), (Qmax_arr[is]));

  
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // prepare all data field strings for later use
    //////////////////////////////////////////////////////////////////////////////////////////////////
    char attrfieldstr[128], attrnamestr[128], btstr[128], latstr[128], lonstr[128], geofile[1024];

    // destriping attribute field
    sprintf(attrfieldstr,"Data_Products/VIIRS-M%i-SDR/VIIRS-M%i-SDR_Aggr\0", is, is);
    printf("Destriping atribute location = %s\n", attrfieldstr);

    // generate names of main data fields to be destriped
    // and the names of the corresponding destriping attributes
    if(is<12) { 

        // for M11 and below, destripe Reflectance
        sprintf(btstr, "All_Data/VIIRS-M%i-SDR_All/Reflectance\0", is);
        sprintf(attrnamestr, "DestripingReflectance\0");
    } 
    else {

        // for M12 and above, destripe Brightness Temperature
        sprintf(btstr, "All_Data/VIIRS-M%i-SDR_All/BrightnessTemperature\0", is);
        sprintf(attrnamestr, "DestripingBrightnessTemperature\0");
    }
    printf("Destriping atribute name = %s\n", attrnamestr);
    printf("Data location = %s\n", btstr);          // name of main data field to be destriped


    //////////////////////////////////////////////////////////////////////////////////////////////////
    // lat, lon data fields in geofile - normally not needed
    sprintf(latstr, "All_Data/VIIRS-MOD-GEO_All/Latitude\0");
    sprintf(lonstr, "All_Data/VIIRS-MOD-GEO_All/Longitude\0");

    // generate the name of the corresponding geofile
    // start with provided SVM file name
    sprintf(geofile, "%s\0", argv[1]);
  
    for(j=(strlen(argv[1])-6); j>=0; j--){
        i = j;
        if(argv[1][j]=='/') { 
            i = j + 1; 
            break; 
        }
    }
    

    // construct geolocation file name  - replace "SVMXY" with "GMODO"
    geofile[i+0] = 'G';
    geofile[i+1] = 'M';
    geofile[i+2] = 'O';
    geofile[i+3] = 'D';
    geofile[i+4] = 'O';

    printf("Corresponding geofile = %s\n", geofile);
    //////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////
    // read data
    //////////////////////////////////////////////////////////////////////////////////////////////////
    if(is!=13){  status = readwrite_viirs(&buffer1, dims1, &scale1, &offset1, argv[1], btstr, 0); } 
    else {       status = readwrite_viirs_float( &bufferf1, dims1, argv[1], btstr, 0); }
    if(status!=0) { 
        printf("ERROR: Cannot read VIIRS data!\n");  
        return 10*status; 
    }


  /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // read geolocation data
    ///////////////////////////////////////////////////////////////////////////////////////////////
    int isgeo = 1;
    status = readwrite_viirs_float( &bufferf2, dims1, geofile, latstr, 0);
    if(status!=0) { printf("Cannot read VIIRS (lat) geolocation data!\n"); isgeo = 0; }
    status = readwrite_viirs_float( &bufferf3, dims1, geofile, lonstr, 0);
    if(status!=0) { printf("Cannot read VIIRS (lon) geolocation data!\n"); isgeo = 0; }
    ///////////////////////////////////////////////////////////////////////////////////////////////
  */


    // extract scale, offset and dimensions info
    sy = dims1[0]; // height, along the track
    sx = dims1[1]; // width, across track, along scan line
    printf("nx = %i ny = %i\n", sx, sy);
    if(is!=13) {
        scale  = ((double) scale1); 
        offset = ((double) offset1);
        printf("scale = %f offset = %f\n", scale, offset);
    }

    // allocate temporary data arrays
    img_in   = allocate_2d_f(sy, sx);
    img_out  = allocate_2d_f(sy, sx);
    binary_M = allocate_2d_i(sy, sx);
    if( (img_in==NULL) || (img_out==NULL) || (binary_M==NULL) ) {
        printf("ERROR: Cannot allocate memory\n");
        return -1; 
    }

    // if needed, apply scale and offset to get physical data
    if(is!=13) {

        // apply scale and offset to get real physical units  
#pragma omp parallel for
        for(ix=0; ix<sx*sy; ix++){ img_in[0][ix] = scale*buffer1[ix] + offset; }
    }
    else {

        // no scaling for band 13
#pragma omp parallel for
        for(ix=0; ix<sx*sy; ix++){ img_in[0][ix] = bufferf1[ix]; }
    }


  /*
    /////////////////////////////////////////////////////////////
    if(isgeo==1) {

        // allocate geolocation arrays
        lat      = allocate_2d_f(sy, sx);
        lon      = allocate_2d_f(sy, sx);
        if( (lat==NULL) || (lon==NULL) ) {
            printf("ERROR: Cannot allocate memory\n"); return -1; 
        }

#pragma omp parallel for
        for(ix=0; ix<sx*sy; ix++){ lat[0][ix] = bufferf2[ix]; }
#pragma omp parallel for
        for(ix=0; ix<sx*sy; ix++){ lon[0][ix] = bufferf3[ix]; }

        free(bufferf2);
        free(bufferf3);

        resampling of image on sorted lon, lat grid
        resample_viirs(img_in, lat, lon, sx, sy, Qmin_arr[is], Qmax_arr[is]);
    }
    /////////////////////////////////////////////////////////////
  */

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // destripe data
    //////////////////////////////////////////////////////////////////////////////////////////////////
    status = destripe_main_frame(img_in, img_out, binary_M, Ndet_arr[is], Niter_arr[is], NEdQ_arr[is], 
                                 Tx_arr[is], Ty_arr[is], Qmin_arr[is], Qmax_arr[is], sx, sy);

    if(status!=0) {
        printf("ERROR: destriping failure\n"); 
        return 100*status; 
    }

    // calculate NIF and NDF - see Ref. 1 for definitions
    float nif, ndf;
    int nwf;
    get_nif_ndf( img_in, img_out, binary_M, sx, sy, &nif, &ndf, &nwf);
    printf("VIIRS_band    NIF        NDF      domain_size\n");
    printf("    %2i      %2.6f  %2.6f      %i\n", is, nif, ndf, nwf);
  
    // Restore all values outside of the destriping domain
    // and scale destriped data back to integers if band != M13
    if(is!=13) {
#pragma omp parallel for private(ix, j)
        for(ix=0; ix<sx*sy; ix++){ 

            // scale destriped data back to integer value
            j = (int) round((img_out[0][ix] - offset)/scale); 

            // check if integer is in the valid range
            if(j<0) { 
                printf("Output data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, j);
                j = 0;
            }
            if(j>65535) { 
                printf("Output data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, j);
                j = 65535;
            }

            // scale original data back to integer value
            i = (int) round(( img_in[0][ix] - offset)/scale); 

            // check if integer is in the valid range
            if(j<0) { 
                printf("Input  data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, i);
                i = 0;
            }
            if(j>65535) { 
                printf("Input  data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, i);
                i = 65535;
            }

            if(binary_M[0][ix]==0) { 

                // for data in destriping domain, return destriped data
                buffer1[ix] = (unsigned short) j; 
            }
            else                   { 

                // for data outside destriping domain, return original data
                buffer1[ix] = (unsigned short) i; 
            }

        } // for ix = all data in 2d array

        // write destriped data back to file as short int
        status = readwrite_viirs(&buffer1, dims1, &scale1, &offset1, argv[1], btstr, 1);
        free(buffer1);
    }
    else { 
        // no conversion for band M13
#pragma omp parallel for
        for(ix=0; ix<sx*sy; ix++) { 

            if(binary_M[0][ix]==0) { 

                // for data in destriping domain, return destriped data
                bufferf1[ix] = img_out[0][ix]; 
            } 
            else                   {

                // for data outside destriping domain, return original data
                bufferf1[ix] =  img_in[0][ix]; 
            }

        } // for ix = all data in 2d array

        // write destriped band M13 data back to file as float
        status = readwrite_viirs_float(&bufferf1, dims1, argv[1], btstr, 1);
        free(bufferf1);
    }

    if(status!=0) { 
        printf("ERROR: Cannot write VIIRS data!\n");
        return 10*status; 
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // write a destriping attribute
    //////////////////////////////////////////////////////////////////////////////////////////////////
    status = write_viirs_destriping_attribute(argv[1], attrfieldstr, attrnamestr, 1.0);
    if(status<0) { 
        printf("ERROR: Cannot write VIIRS attribute!\n");
        return 40*status;
    }

    // free memory
    free(img_in[0]);   free(img_in); 
    free(img_out[0]);  free(img_out); 
    free(binary_M[0]); free(binary_M); 

    //////////////////////////////////////////
    // free(lat[0]);      free(lat); 
    // free(lon[0]);      free(lon); 
    //////////////////////////////////////////

    return 0;
}


