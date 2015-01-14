
/////////////////////////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------------------------- //
// This code was developed in NOAA/NESDIS/STAR SST group by Karlis Mikelsons.                  //
//                                                                                             //
// Please report any bugs to Karlis.Mikelsons@noaa.gov                                         //
// ------------------------------------------------------------------------------------------- //
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "mfhdf.h"
#define MAX_STR_LEN 256

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// This subroutine reads/writes MODIS data from/to HDF4 file as unsigned short (2 byte unsigned integer)
// type array, to be scaled with scale and offset values.
//
// Arguments:
// unsigned short **    buffer      IN/OUT    If readwrite == 0, on output contains the pointer to an
//                                            newly allocated array with read data; can be undefined on input; 
//                                            but caller is responsible for release of this memory later.
//                                            If readwrite != 0, on input must contain the pointer to 
//                                            array with data to be written; unchanged on return.
//
// int *                nx          OUT       On output contains the size of dimension x of 2d array
//                                            [swath width in pixels of granule (should be 3200 for VIIRS)]
//
// int *                ny          OUT       On output contains the size of dimension y of 2d array
//                                            [height in pixels of granule (depends on length of granule)]
//
// int                  nband       IN        Number of bands in record, also size of isband array  
//
// float *              scales      OUT       Returns scale  for each band in data record
// float *              offsets     OUT       Returns offset for each band in data record
//                                            Physical value = scale*(buffer[0][i] - offset) for MODIS
//
// int *                isband      IN        if isband[i]==0, skip reading/writing band i
//                                            in data record
//
// char *               sds_name    IN        Name of data record in hdf4 file
//                                            (such as "EV_1KM_Emissive")
//
// char *               attr_name   IN        Scales and offsets attribute base name
//                                            (such as "radiance" for "radiance_scales" and "radiance_offsets")
//
// char *               filename    IN        Name of the hdf4 file to read of write to
//
// int                  readwrite   IN        if readwrite == 0, read data
//                                            if readwrite != 0, write data
//
// Return value: 
// Upon sucessful completion, returns a non-negative number equal to the number of bands read/written;
// negative return value indicates error;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int readwrite_modis(unsigned short ** buffer, int * nx, int * ny, int nband, float *scales, float *offsets,
                    int *isband, char * sds_name, char * attr_name, char * filename, int readwrite) 
{

    int32 sd_id, sds_index, sds_id; /* SD interface and data set identifiers */
    intn status; /* status returned by some routines; has value SUCCEED or FAIL */
    int i,j,k, iprint = 0;


    // find how many bands need to be read/written
    int nreadwrite = 0;
    for(i=0; i<nband; i++) { if(isband[i]>0) nreadwrite++; }
    if(iprint > 0) printf("nreadwrite = %i\n", nreadwrite);
    if(nreadwrite==0) return 0;


    // open file for read or write
    if (readwrite==0) sd_id = SDstart(filename, DFACC_READ);
    else              sd_id = SDstart(filename, DFACC_WRITE);
    if(sd_id==FAIL) {  
        if(iprint > 0) printf("Cannot open file %s with SDstart\n", filename);
        return -1;  
    }


    // find the index of data record
    sds_index = SDnametoindex(sd_id, sds_name);
    if(sds_index==FAIL) { 
        if(iprint > 0) printf("Cannot get index of %s with SDnametoindex\n", sds_name);
        return -1; 
    }


    // open the data record
    sds_id = SDselect(sd_id, sds_index);
    if(sds_id==FAIL) {
        if(iprint > 0) printf("Cannot select data set with SDselect\n");
        return -1;
    }

    char sds_name1[MAX_STR_LEN];
    int32 rank = 0, data_type = 0, num_attrs = 0, n_values = 0;


    // read the scales   for all bands in data record
    char full_attr_name[MAX_STR_LEN], full_attr_name1[MAX_STR_LEN];
    sprintf(full_attr_name, "%s_scales\0", attr_name);
    intn attr_index = SDfindattr (sds_id, full_attr_name);
    if(status==FAIL) {
        printf("Cannot find attr %s with SDfindattr\n", full_attr_name);
        return -1;  
    }

    status = SDattrinfo (sds_id, attr_index, full_attr_name1, &data_type, &n_values);
    if(status==FAIL) {
        printf("Cannot get attr info about %s with SDattrinfo\n", full_attr_name);
        return -1;
    }
    if(iprint > 0) printf("attribute  %s   data type = %i, nvalues = %i \n", full_attr_name, data_type, n_values);

    status = SDreadattr (sds_id, attr_index, scales);
    if(status==FAIL) { 
        printf("Cannot read attr %s with SDreadattr\n", full_attr_name);
        return -1;
    }
    if(iprint > 0) for(i=0;i<n_values;i++) { printf("%i %e\n", i, scales[i]); }


    // read the offsets for all bands in data record
    sprintf(full_attr_name, "%s_offsets\0", attr_name);
    attr_index = SDfindattr (sds_id, full_attr_name);
    if(status==FAIL) { printf("Cannot find attr %s with SDfindattr\n", full_attr_name);  return -1;  }
    status = SDattrinfo (sds_id, attr_index, full_attr_name1, &data_type, &n_values);
    if(status==FAIL) { printf("Cannot get attr info about %s with SDattrinfo\n", full_attr_name);  return -1;  }
    if(iprint > 0) printf("attribute  %s   data type = %i, nvalues = %i \n", full_attr_name, data_type, n_values);
    status = SDreadattr (sds_id, attr_index, offsets);
    if(status==FAIL) { printf("Cannot read attr %s with SDreadattr\n", full_attr_name);  return -1;  }
    if(iprint > 0) for(i=0;i<n_values;i++) { printf("%i %e\n", i, offsets[i]); }


    // get the dimension size info
    int dimsizes[32];
    for(i=0; i<32; i++) dimsizes[i] = 0;
    status = SDgetinfo(sds_id, sds_name1, &rank, dimsizes, &data_type, &num_attrs);
    if(status==FAIL) {
        if(iprint > 0) printf("Cannot get info about %s with SDgetinfo\n", sds_name1);
        return -1;
    }
    if(iprint > 0) printf("rank = %i\n", rank);
    if(rank!=3) {
        printf("ERROR: %s rank = %i != 3\n", sds_name, rank);
        return -1;
    }
    if(data_type!=DFNT_UINT16) {
        printf("ERROR: %s data_type = %i != DFNT_UINT16\n", sds_name, data_type);
         return -1;
    }
    if(iprint > 0) for(i=0; i<rank; i++) printf("%i %i\n",i,dimsizes[i]);
    if(iprint > 0) printf("datatype = %i\n", data_type);

    // dimsizes[0] is number of bands in data record
    *ny = dimsizes[1];
    *nx = dimsizes[2];
  

    // if reading, allocate buffer for data
    int ntot = dimsizes[1] * dimsizes[2] * nreadwrite; 
    if(readwrite==0){
        buffer[0] = NULL;
        buffer[0] = (unsigned short *) malloc(ntot*sizeof(unsigned short));
        if(buffer[0]==NULL) {
            if(iprint > 0) printf("Cannot allocate memory %i bytes\n", ntot*sizeof(unsigned short));
            return -1;
        } 
    }

    int32 start[3]  = { 0, 0, 0 };
    int32 stride[3] = { 1, 1, 1 };
    int32 edge[3]   = { 1, dimsizes[1], dimsizes[2] };


    // read / write bands 
    int nb = 0;
    for(i=0; i<nband; i++){
        if(isband[i]==0) continue;   // skip bands that are not needed
        start[0] = i;                // start index of band in data record in file
                                 // nb is here used as an index of band in memory
        if(readwrite == 0){
            status =  SDreaddata(sds_id, start, stride, edge, (VOIDP) &(buffer[0][nb*dimsizes[1]*dimsizes[2]]) );
            if(status==FAIL) {
                if(iprint > 0) printf("Cannot  read data with SDreaddata\n");
                return -1; 
            } 
        }
        else {
            status = SDwritedata(sds_id, start, stride, edge, (VOIDP) &(buffer[0][nb*dimsizes[1]*dimsizes[2]]) );
            if(status==FAIL) {
                 if(iprint > 0) printf("Cannot write data with SDwritedata\n");
                 return -1; 
            } 
        }
        nb++;
    }


    // now deal with destriping attribute
    if(readwrite==1) {

        // first, try to find an existing destriping attribute
        sprintf(full_attr_name, "Destriping\0");
        float attrbuff[32];
        for(i=0; i<nband; i++){ attrbuff[i] = 0; }
        attr_index = SDfindattr (sds_id, full_attr_name);
        if(attr_index!=FAIL) {

            // we need to read the values of this attribute and modify those of destriped bands
            status = SDreadattr (sds_id, attr_index, attrbuff);
            if(status==FAIL) {
                printf("Cannot read attr %s with SDreadattr\n", full_attr_name);
                return -2;
            }
        }

        // modify attribute for destriped bands  
        for(i=0; i<nband; i++){ 
            if(isband[i]==0) continue;
            if(attrbuff[i]>0) { printf("WARNING: This band %i in data field %s was already destriped\n", i, sds_name); }
            attrbuff[i] = (float) isband[i];
        }

        // write the modified attributes back to the data record
        status = SDsetattr(sds_id, full_attr_name, DFNT_FLOAT32, nband, attrbuff);
        if(status==FAIL) {
            if(iprint > 0) printf("Cannot write attribute with SDsetattr\n");
            return -2;
        }

    } // if(readwrite==1) 
    // done with destriping attribute 


    // close data record
    status = SDendaccess(sds_id);
    if(status==FAIL) {     if(iprint > 0) printf("Cannot end access with SDendaccess\n");     return -1;  }

    // close hdf4 file access
    status = SDend(sd_id);
    if(status==FAIL) {     if(iprint > 0) printf("Cannot end with SDend\n");     return -1;  }

    if(iprint > 0) printf("modis_readwrite done \n");
    return nb;
};

