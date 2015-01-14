
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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mfhdf.h"
#include "readwrite_modis.c"
#include "destripe.h"
#include "get_nif_ndf.c"

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {

    // there are four data fields in MODIS HDF files that contain all bands this code can process
    // here are the data field names for the four data fields
    char EV_250_Ref[] = "EV_250_Aggr1km_RefSB";
    char EV_500_Ref[] = "EV_500_Aggr1km_RefSB";
    char EV_1KM_Ref[] = "EV_1KM_RefSB";
    char EV_1KM_Emi[] = "EV_1KM_Emissive";
    char * dataFieldNames[4] = { EV_250_Ref, EV_500_Ref, EV_1KM_Ref, EV_1KM_Emi };

    // attribute names for the four data fields
    char AttrRef[] = "reflectance";
    char AttrRad[] = "radiance";
    char * attrBaseNames[4] = { AttrRef, AttrRef, AttrRef, AttrRad };

    // define band names for all bands
    char bandNames[38][4] = {"1\0\0\0", "2\0\0\0",                                      /* in "EV_250_Aggr1km_RefSB" */
                             "3\0\0\0", "4\0\0\0", "5\0\0\0", "6\0\0\0", "7\0\0\0",     /* in "EV_500_Aggr1km_RefSB" */
                             "8\0\0\0", "9\0\0\0", "10\0\0",  "11\0\0",  "12\0\0",  "13lo", "13hi", "14lo", "14hi", "15\0\0", "16\0\0", "17\0\0", "18\0\0", "19\0\0", "26\0\0",                       /* in "EV_1KM_RefSB" */
                             "20\0\0",  "21\0\0",  "22\0\0",  "23\0\0",  "24\0\0",  "25\0\0", "27\0\0", "28\0\0", "29\0\0", "30\0\0", "31\0\0", "32\0\0", "33\0\0", "34\0\0", "35\0\0", "36\0\0" };   /* in "EV_1KM_Emissive" */

    // indices in bandNames array of first band in a data field
    int bandIndex[5] = { 0, 2, 7, 22, 38};

    // define array of wavelengths and assign correct values for emissive bands
    double lambda[38];
    lambda[22] = 0.5*( 3.660 +  3.840)*1.0E-6; // band 20
    lambda[23] = 0.5*( 3.929 +  3.989)*1.0E-6; // band 21
    lambda[24] = 0.5*( 3.929 +  3.989)*1.0E-6; // band 22
    lambda[25] = 0.5*( 4.020 +  4.080)*1.0E-6; // band 23
    lambda[26] = 0.5*( 4.433 +  4.498)*1.0E-6; // band 24
    lambda[27] = 0.5*( 4.482 +  4.549)*1.0E-6; // band 25
    lambda[28] = 0.5*( 6.535 +  6.895)*1.0E-6; // band 27
    lambda[29] = 0.5*( 7.175 +  7.475)*1.0E-6; // band 28
    lambda[30] = 0.5*( 8.400 +  8.700)*1.0E-6; // band 29
    lambda[31] = 0.5*( 9.580 +  9.880)*1.0E-6; // band 30
    lambda[32] = 0.5*(10.780 + 11.280)*1.0E-6; // band 31
    lambda[33] = 0.5*(11.770 + 12.270)*1.0E-6; // band 32
    lambda[34] = 0.5*(13.185 + 13.485)*1.0E-6; // band 33
    lambda[35] = 0.5*(13.485 + 13.785)*1.0E-6; // band 34
    lambda[36] = 0.5*(13.785 + 14.085)*1.0E-6; // band 35
    lambda[37] = 0.5*(14.085 + 14.385)*1.0E-6; // band 36
  
    // k = Boltzmann gas constant (joules/Kelvin)
    float k_Boltz = 1.3806488e-23;
    // h = Planck’s constant (joule * second)
    float h_Planck = 6.62606957e-34;
    // c = speed of light in vacuum (m/s)
    float c_light = 299792458.0;
  

    unsigned short *buffer1 = NULL, *buff1, usi;
    float **inp_img  = NULL;
    float **outp_img = NULL;
    int   **binary_M = NULL;
    int    *maskNaN  = NULL;

    int nx1, ny1, is, ns, ix, iy, status, i, j, n, nthreads;

    int Ndet, Niter, ib, nb, nx, ny, iband,  iDataField, nmask;
    int Ndet_arr[40], Niter_arr[40], isBand[40];
    float Qmin, Qmax, Thresh_x, Thresh_y, NEdQ, r1, r2;
    float Qmin_arr[40], Qmax_arr[40], Tx_arr[40], Ty_arr[40], NEdQ_arr[40], Scale_arr[40], Offset_arr[40];


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // check if at least two command line arguments
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(argc<3) { 
        printf("Not enough arguments!\nUsage:\n %s MODIS_hdf_file  destriping_param_file.txt\n", argv[0]);
        return -9;
    }
    printf("destripe_modis %s %s\n", argv[1], argv[2]);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // set the number of OpenMP threads
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    nthreads = omp_get_max_threads();
    omp_set_num_threads(nthreads);
    printf("numthreads = %i maxthreads =  %i\n", omp_get_num_threads(), omp_get_max_threads());


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // read input parameters from parameter file
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("Input parameters\n");

    // open parameter file
    FILE *fp = fopen(argv[2],"r");
    if(fp==NULL) {  printf("ERROR: Cannot open parameter file\n");  return -8; }

    // read parameters
    char tmpbandname[128];
    for(i=0; i<40; i++){ isBand[i] = 0; }    // initially, set all band parameters as "absent"
    for(i=0; i<40; i++){                     // read at most 40 lines in parameter file

        j = fscanf(fp,"%s ", tmpbandname);   // read the first item in a line = band name
        if(j!=1) break;                      // if not read, then break

        // check if the band name is sensible
        is = -1;
        for(j=0;j<38;j++){ 

            // compare the band name to all possible band names in an array
            if(strcmp(tmpbandname, bandNames[j]) == 0) {
                is = j;  // if found a match, assign band index
                break;   // and stop comparing
            }
        }
        if(is==-1) break; // if band name not valid, break

        // read all destriping parameters for this band
        j = fscanf(fp,"%i %i %f %f %f %f %f\n", &(Ndet_arr[is]), &(Niter_arr[is]), &(NEdQ_arr[is]),
                                                &(Tx_arr[is]), &(Ty_arr[is]), &(Qmin_arr[is]), &(Qmax_arr[is]));

        if(j==7) { isBand[is]=1; } else break;  // if destriping parameters not read correctly, break

        // echo read destriping parameters
        printf(" %s \t %i %i %f %f %f %f %f\n", bandNames[is], (Ndet_arr[is]), (Niter_arr[is]), (NEdQ_arr[is]), 
                                               (Tx_arr[is]), (Ty_arr[is]), (Qmin_arr[is]), (Qmax_arr[is]));
    }

    fclose(fp); // close parameter file
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // done reading parameters
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // loop over all 4 data fields
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(iDataField=0; iDataField<4; iDataField++){
        printf("========================================================================\n");
        printf("Data_field number = %i   name = %s\n", iDataField, dataFieldNames[iDataField]);
        printf("------------------------------------------------------------------------\n");
        ib = bandIndex[iDataField];                              // index of first band of this data field in bandNames[] array
        nb = bandIndex[iDataField+1] - bandIndex[iDataField];    // number of bands in this data field

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // make sure we have at least one band to destripe in this data field
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int nreadwrite = 0;
        for(iband=0; iband<nb; iband++){ if(isBand[ib+iband]>0) nreadwrite++; }
        if(nreadwrite==0) continue;   // if no bands to destripe in this data field, continue

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // read data to destripe in this data field    
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        status = readwrite_modis( &buffer1, &nx, &ny, nb, &(Scale_arr[ib]), &(Offset_arr[ib]), &(isBand[ib]), 
                                  dataFieldNames[iDataField], attrBaseNames[iDataField], argv[1], 0);

        if(status<0) { 
            printf("ERROR: Cannot read data field %s\n", dataFieldNames[iDataField]); 
            return 10*status;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // allocate temporary arrays
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        inp_img  = allocate_2d_f(ny,nx);
        outp_img = allocate_2d_f(ny,nx);
        binary_M = allocate_2d_i(ny,nx);
        maskNaN  = (int *) malloc(ny*nx*sizeof(int)); 
        if( (maskNaN==NULL) || (inp_img==NULL) || (outp_img==NULL) || (binary_M==NULL) ) { 
            printf("ERROR: Cannot allocate memory\n");
            return -1;
        }
    

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // go through all the bands in the current data field
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int iBandIndx = 0;
        for(iband=0; iband<nb; iband++){

            // index of the current band in bandNames[] array is a sum of 
            //    index of first band of this data field in bandNames[] array and 
            //    index of band within the current data field
            is = ib + iband; 

            // printf("iband = %i isband = %i\n", iband, isBand[is]);
            if(isBand[is]==0) continue; // if no parameters for this band, then pass

            buff1 = &(buffer1[iBandIndx*nx*ny]); // location of the current band data to destripe
            iBandIndx++;                         // increment the index of next band data to destripe

            printf("Band = %i  MODIS_band_number = %s   scale = %e  offset = %e\n", iband, bandNames[is], Scale_arr[is], Offset_arr[is]);
 
            if(iDataField==3) {

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // emissive bands  - transfrom to brightness temperature
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                r1 = h_Planck*c_light/(k_Boltz*lambda[is]);
                r2 = lambda[is];
                r2 = 1.0e-6*(2.0*h_Planck*c_light*c_light)/(r2*r2*r2*r2*r2);

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // find the minimum valid radiance > 0, mask all pixels with radiance <= 0
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                unsigned short jmin = 65535;
                nmask = 0;
                for(ix=0; ix<nx*ny; ix++) { 

                    maskNaN[ix] = 0;      // originally assume data are physically valid

                    if(buff1[ix]<=Offset_arr[is]) { 

                        // found a pixel with negative radiance
                        maskNaN[ix] = 1;  // set the mask for unphysical data
                        nmask++;          // count such pixels
                        continue; 
                    }

                    // update minimum physical radiance - used later as fill in value for unphysical values
                    if(buff1[ix]<jmin) { jmin = buff1[ix]; }
                }
                // printf("nmask = %i nx*ny = %i jmin = %i\n", nmask, nx*ny, jmin);
 
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // convert radiance to brightness temperature
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                nmask = 0;
                double avebt = 0.0;

#pragma omp parallel for private(ix, j) reduction(+:nmask, avebt)
                for(ix=0; ix<nx*ny; ix++) { 
                    j = buff1[ix];

                    // if radiance less than smallest physical value, set it to fill in value
                    if(j<jmin) { j = jmin; nmask++; }

                    // scale integers to physical radiance values
                    inp_img[0][ix] = Scale_arr[is]*(j - Offset_arr[is]);

                    // calculate Brightness Temperature from Radiance
                    inp_img[0][ix] = r1/log(1.0 + r2/inp_img[0][ix]);

                    // counter for average BT - just to check sanity of data
                    avebt +=  inp_img[0][ix];
                }

                printf("Number of pixels with negative radiances on input = %i\n", nmask); 
                // printf("Average Brightness Temperature = %e \n", avebt/(nx*ny));
                printf("Parameters:\n");
                printf("%i %i %f   %f %f   %f %f\n",  Ndet_arr[is], Niter_arr[is], NEdQ_arr[is], Tx_arr[is], Ty_arr[is], Qmin_arr[is], Qmax_arr[is]);


                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // destripe here
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                status = destripe_main_frame(inp_img, outp_img, binary_M, Ndet_arr[is], Niter_arr[is], NEdQ_arr[is], 
                                             Tx_arr[is], Ty_arr[is], Qmin_arr[is], Qmax_arr[is], nx, ny);

                if(status!=0) { 
                    printf("ERROR: destriping failure\n"); 
                    return 100*status; 
                }
                printf("Destriping done\n");


                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // calculate NIF and NDF - see Ref. 1 for details
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                float nif, ndf;
                int nwf;
                get_nif_ndf( inp_img, outp_img, binary_M, nx, ny, &nif, &ndf, &nwf);
                printf("MODIS_band_number     NIF        NDF      domain_size\n");
                printf("         %s         %2.6f  %2.6f      %i\n", bandNames[is], nif, ndf, nwf);
   
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // convert output brightness temperature back to radiance
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                nmask = 0;
#pragma omp parallel for private(ix, j) reduction(+:nmask)
                for(ix=0; ix<nx*ny; ix++) { 

                    if(binary_M[0][ix]==1) nmask++;  // count the number of pixels outside destriping domain

                    // preserve the original data outside destriping domain
                    // or in case of unphysical negative radiance in original data (produces NaN in Brightness Temperature)
                    if( (binary_M[0][ix]==1) || (maskNaN[ix] == 1) ) continue;            

                    // get the radiance from brightness temperature
                    outp_img[0][ix] = r2/(exp(r1/outp_img[0][ix]) - 1.0);

                    // scale the radiance back to integer
                    j = (int) round( outp_img[0][ix]/Scale_arr[is] + Offset_arr[is] );

                    // check that integer is within valid bounds
                    if((j<0) || (j>65535)) {
                        printf("Value outside range %i\n", j);
                        j = 65535;
                    }

                    // copy value to output buffer
                    buff1[ix] = (unsigned short) j;
                }
                printf("binary_M==1 for %i out of %i pixels (%5.2f%)\n", nmask, nx*ny, (100.0*nmask)/(nx*ny));

            }  //  if(iDataField==3)
            else {

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // reflective band
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // convert scaled integers to physical reflectance
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
                for(ix=0; ix<nx*ny; ix++) { inp_img[0][ix] = Scale_arr[is]*( ((float) (buff1[ix])) - Offset_arr[is]); }

                printf("Parameters:\n");
                printf("%i %i %f   %f %f   %f %f\n",  Ndet_arr[is], Niter_arr[is], NEdQ_arr[is], Tx_arr[is], Ty_arr[is], Qmin_arr[is], Qmax_arr[is]);
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // destripe here
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                status = destripe_main_frame(inp_img, outp_img, binary_M, Ndet_arr[is], Niter_arr[is], NEdQ_arr[is], 
                                             Tx_arr[is], Ty_arr[is], Qmin_arr[is], Qmax_arr[is], nx, ny);

                if(status!=0) { 
                    printf("ERROR: destriping failure\n");
                    return 100*status;
                }
                printf("Destriping done\n");


                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // calculate NIF and NDF  - see Ref. 1 for details
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                float nif, ndf;
                int nwf;
                get_nif_ndf( inp_img, outp_img, binary_M, nx, ny, &nif, &ndf, &nwf);
                printf("MODIS_band_number     NIF        NDF      domain_size\n");
                printf("         %s         %2.6f  %2.6f      %i\n", bandNames[is], nif, ndf, nwf);

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // convert reflectance back to scaled integers
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                nmask = 0;
#pragma omp parallel for private(ix, j) reduction(+:nmask)
                for(ix=0; ix<nx*ny; ix++) { 

                    if(binary_M[0][ix]==1) { 
                        nmask++;   // count the number of pixels outside destriping domain
                        continue;  // but otherwise preserve the original data in these pixels
                    }

                    // scale the reflectance back to integer
                    j = (int) round( outp_img[0][ix]/Scale_arr[is] + Offset_arr[is] );

                    // check that integer is within valid bounds
                    if((j<0) || (j>65535)) {
                        printf("Value outside range %i\n", j);
                        j = 65535;
                    }
         
                    // copy value to output buffer
                    buff1[ix] = (unsigned short) j;
                } 
                printf("binary_M==1 for %i out of %i pixels (%5.2f%)\n", nmask, nx*ny, (100.0*nmask)/(nx*ny));

            }  //  if(iDataField!=3) meaning reflective band
    
            printf("------------------------------------------------------------------------\n");
     
        } // for iband

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // write destriped data back to hdf file, as well as set destriping attribute
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        status = readwrite_modis( &buffer1, &nx, &ny, nb, &(Scale_arr[ib]), &(Offset_arr[ib]), &(isBand[ib]), 
                                  dataFieldNames[iDataField], attrBaseNames[iDataField], argv[1], 1);

        if(status<0) { 
            printf("ERROR: Failed to write data\n"); 
            return 10*status; 
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // clean up
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(buffer1!=NULL) {  free(buffer1); buffer1 = NULL; }
        free(inp_img[0]);  free(inp_img);
        free(outp_img[0]); free(outp_img);
        free(binary_M[0]); free(binary_M);
        free(maskNaN);

    } // for(iDataField = ...

    return 0;
}
