#include <stdio.h>
#include "omp.h"
#include "readwrite_viirs.c"
#include "destripe.h"
#include "allocate_2d.c"
#include "resample_viirs.h"

void get_nif_ndf(float ** oldimg, float ** newimg, int ** binary_M, int nx, int ny, float * nif, float * ndf, int * nwf);

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

  // make sure we have at least two command line arguments
  if(argc<3) { 
    printf("Not enough arguments!\nUsage:\n %s viirs_h5_file destriping_parameter_file_viirs.txt\n", argv[0]);
    return 0;
  }
  printf("destripe_viirs %s %s\n", argv[1], argv[2]);

  // set the number of threads  
  nthreads = omp_get_max_threads();
  omp_set_num_threads(nthreads); 
  printf("numthreads = %i maxthreads =  %i\n", omp_get_num_threads(), omp_get_max_threads());

  // open parameter file
  FILE *fp = fopen(argv[2],"r");
  if(fp==NULL) { 
    printf("Cannot open param. file\n");
    return -1;
  }
 
  // read parameters
  is = 0;
  for(i=0; i<16; i++){ isband[is] = 0; }
  for(i=0; i<16; i++){
    j = fscanf(fp,"%i ", &is);
    if( (j!=1) || (is<1) || (is>16) ) break;
    j = fscanf(fp,"%i %i %f %f %f %f %f\n", &(Ndet_arr[is]), &(Niter_arr[is]), &(NEdQ_arr[is]), &(Tx_arr[is]), &(Ty_arr[is]), &(Qmin_arr[is]), &(Qmax_arr[is]));
    if(j!=7) break;
    isband[is] = 1;
  }

  // close param. file
  fclose(fp);
  // done reading parameters

  // extract the name of band from the file
  is = 0;
  for(i=(strlen(argv[1])-6); i>=0; i--) {
    if( (argv[1][i]=='S') && (argv[1][i+1]=='V') && (argv[1][i+2]=='M') ){
      is = (argv[1][i+3]-'0')*10 + (argv[1][i+4]-'0');
      break;
    }
  }
  printf("Band = %i\n", is);
  if((is<1)||(is>16)) {
    printf("Invalid band\n");   
    return -1;
  }
  // make sure we have parameters for destriping this band
  if(isband[is]!=1) {
    printf("No destriping parameters for this band\n");
    return 0;
  }

  // print parameters for this band
  printf("%i %i %i %f %f %f %f %f\n", is, (Ndet_arr[is]), (Niter_arr[is]), (NEdQ_arr[is]), (Tx_arr[is]), (Ty_arr[is]), (Qmin_arr[is]), (Qmax_arr[is]));

  
  // prepare all data field strings for later use
  char attrfieldstr[128], attrnamestr[128], btstr[128], latstr[128], lonstr[128], geofile[1024];
  sprintf(attrfieldstr,"Data_Products/VIIRS-M%i-SDR/VIIRS-M%i-SDR_Aggr\0", is, is);
  printf("attrfieldstr = %s\n", attrfieldstr);
  if(is<12) { 
    sprintf(btstr, "All_Data/VIIRS-M%i-SDR_All/Reflectance\0", is);
    sprintf(attrnamestr, "DestripingReflectance\0");
  } else {
    sprintf(btstr, "All_Data/VIIRS-M%i-SDR_All/BrightnessTemperature\0", is);
    sprintf(attrnamestr, "DestripingBrightnessTemperature\0");
  }
  sprintf(latstr, "All_Data/VIIRS-MOD-GEO_All/Latitude\0");
  sprintf(lonstr, "All_Data/VIIRS-MOD-GEO_All/Longitude\0");
  // sprintf(geofile, "GMODO%s\0", argv[1][5]);
  sprintf(geofile, "%s\0", argv[1]);
  
  // find if this is VIIRS or MODIS file
  // search for the first letter of filename
  for(i=(strlen(argv[1])-2); i>0; ){
    if(argv[1][i-1]=='/') break;
    i--;
  }
  geofile[i+0] = 'G';
  geofile[i+1] = 'M';
  geofile[i+2] = 'O';
  geofile[i+3] = 'D';
  geofile[i+4] = 'O';

  printf("btstr = %s\n", btstr);
  printf("geofile = %s\n", geofile);


  // read data
  if(is!=13){  status = readwrite_viirs(&buffer1, dims1, &scale1, &offset1, argv[1], btstr, 0); } 
  else {       status = readwrite_viirs_float( &bufferf1, dims1, argv[1], btstr, 0); }
  if(status!=0) { printf("Cannot read VIIRS data!\n");  return 1; }

  // read geolocation data
  int isgeo = 1;
  status = readwrite_viirs_float( &bufferf2, dims1, geofile, latstr, 0);
  if(status!=0) { printf("Cannot read VIIRS (lat) geolocation data!\n"); isgeo = 0; }
  status = readwrite_viirs_float( &bufferf3, dims1, geofile, lonstr, 0);
  if(status!=0) { printf("Cannot read VIIRS (lon) geolocation data!\n"); isgeo = 0; }


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
  lat      = allocate_2d_f(sy, sx);
  lon      = allocate_2d_f(sy, sx);
  if( (img_in==NULL) || (img_out==NULL) || (binary_M==NULL) || (lat==NULL) || (lon==NULL) ) {
    printf("ERROR: Cannot allocate memory\n"); return -1; 
  }
 
  if(is!=13) {
    // apply scale and offset to get real physical units  
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++){ img_in[0][ix] = scale*buffer1[ix] + offset; }
  }
  else { // no scaling for band 13
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++){ img_in[0][ix] = bufferf1[ix]; }
  }

  if(isgeo==1) {
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++){ lat[0][ix] = bufferf2[ix]; }
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++){ lon[0][ix] = bufferf3[ix]; }

    free(bufferf2);
    free(bufferf3);

    // resampling of image on sorted lon, lat grid
    resample_viirs(img_in, lat, lon, sx, sy, Qmin_arr[is], Qmax_arr[is]);
  }

  // destripe data
  destripe_main_frame(img_in, img_out, binary_M, Ndet_arr[is], Niter_arr[is], NEdQ_arr[is], Tx_arr[is], Ty_arr[is], Qmin_arr[is], Qmax_arr[is], sx, sy);

  // calculate NIF and NDF
  float nif, ndf;
  int nwf;
  get_nif_ndf( img_in, img_out, binary_M, sx, sy, &nif, &ndf, &nwf);
  printf("  NIF        NDF      domain_size\n %2.6f  %2.6f   %i\n", nif, ndf, nwf);

  // Convert destriped data back to integers
  // and Restore all values outside of the destriping domain
  if(is!=13) {
#pragma omp parallel for private(ix, j)
    for(ix=0; ix<sx*sy; ix++){ 
      j = (int) round((img_out[0][ix] - offset)/scale); 
      if(j<0) { 
        printf("Data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, j);
        j = 0;
      }
      if(j>65535) { 
        printf("Data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, j);
        j = 65535;
      }

      i = (int) round(( img_in[0][ix] - offset)/scale); 
      if(j<0) { 
        printf("Data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, i);
        i = 0;
      }
      if(j>65535) { 
        printf("Data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, i);
        i = 65535;
      }
      if(binary_M[0][ix]==0) { buffer1[ix] = (unsigned short) j; }
      else                   { buffer1[ix] = (unsigned short) i; }
    }
    status = readwrite_viirs(&buffer1, dims1, &scale1, &offset1, argv[1], btstr, 1);
    free(buffer1);
  }
  else { // special treatment for band 13
#pragma omp parallel for
    for(ix=0; ix<sx*sy; ix++){ 
      if(binary_M[0][ix]==0) { bufferf1[ix] = img_out[0][ix]; } 
      else                   { bufferf1[ix] =  img_in[0][ix]; }
    }
    status = readwrite_viirs_float(&bufferf1, dims1, argv[1], btstr, 1);
    free(bufferf1);
  }
  if(status!=0) { printf("Cannot write VIIRS data!\n");  return 1; }

  // write a destriping attribute
  status = write_viirs_destriping_attribute(argv[1], attrfieldstr, attrnamestr, 1.0);
  if(status==-1) { printf("Cannot write VIIRS attribute!\n");  return 1; }

  free(img_in[0]);   free(img_in); 
  free(img_out[0]);  free(img_out); 
  free(binary_M[0]); free(binary_M); 
  free(lat[0]);      free(lat); 
  free(lon[0]);      free(lon); 

  return 0;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_nif_ndf(float ** oldimg, float ** newimg, int ** binary_M, int nx, int ny, float * nif, float * ndf, int * nwf){
  // assign default return values
  *nif = 0.0;
  *ndf = 1.0;
  *nwf = 0;

  int ix, iy;
  float dxold, dxnew, dyold, dynew;
  float sumdiffy = 0.0, sumdiffx = 0.0, sumx = 0.0, sumy = 0.0;

  // loop over all image
  for(iy=0; iy<(ny-1); iy++){
    for(ix=0; ix<(nx-1); ix++){
      // exclude points outside destriping domain
      if( (binary_M[iy][ix]!=0) || (binary_M[iy+1][ix]!=0) || (binary_M[iy][ix+1]!=0) )continue;
      // exclude points with zero values both before and after destriping
      // if( (fabs(oldimg[iy][ix])<0.0001) &&  (fabs(newimg[iy][ix])<0.0001) ) continue; 
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
}
