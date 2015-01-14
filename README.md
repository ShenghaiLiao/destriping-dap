
## Destriping algorithm for VIIRS/MODIS Ocean Products.

This software was developed by the SST group at NOAA STAR (Lead:
A. Ignatov). It is designed to destripe VIIRS Level 1b (SDR) data.
See the following publication for more details on the algorithm:

    Bouali, Marouan, and Alexander Ignatov, 2014: Adaptive Reduction
    of Striping for Improved Sea Surface Temperature Imagery from Suomi
    National Polar-Orbiting Partnership (S-NPP) Visible Infrared Imaging
    Radiometer Suite (VIIRS).  J. Atmos. Oceanic Technol., 31, 150–163. doi:
    http://dx.doi.org/10.1175/JTECH-D-13-00035.1

### Dependencies

Besides standard C library functions, this code relies on the following
libraries:
* [fftw3 library](http://www.fftw.org/) (float version),
* [HDF5 library](http://www.hdfgroup.org/HDF5/) (and dependencies, such as libsz, libz),
* libpthread, and
* OpenMP library (such as libiomp5 or [libgomp](https://gcc.gnu.org/projects/gomp/))

The code is designed to run multiple threads using OpenMP. 
The number of threads can be set before execution by
setting environment variable `OMP_NUM_THREADS`.

### Installation

To install, run `make install` inside destripe_viirs directory. You may
need to adapt the `Makefile` for your system. The installed binary is named
`destripe_viirs`.

### Usage

On input, it requires two command line arguments: First argument should
be the name of hdf5 file containing VIIRS SDR data, Second argument is the
parameter file (see sample at
[destripe_viirs/destriping_viirs_param.txt](destripe_viirs/destriping_viirs_param.txt)).

On output, the hdf5 file with VIIRS SDR data is modified -
the original data are replaced by the destriped data.
In addition, it includes a destriping attribute under data field
`Data_Products/VIIRS-M??-SDR/VIIRS-M??-SDR_Aggr`.
The name of destriping attribute is `DestripingBrightnessTemperature` for bands 12-16,
and it is `DestripingReflectance` for bands 1-11.

Below is an example of destriping band 12 of a VIIRS granule. First, download the example granule we've made available:
```
$ curl -O ftp://www.star.nesdis.noaa.gov/pub/sod/osb/aignatov/ACSPO/destriping/GMODO_npp_d20140410_t2130000_e2140001_b12705_c20140414174333267547_star_dev.h5
$ curl -O ftp://www.star.nesdis.noaa.gov/pub/sod/osb/aignatov/ACSPO/destriping/SVM12_npp_d20140410_t2130000_e2140001_b12705_c20140414174333267547_star_dev.h5
```
This is the dataset that we want to destripe:
```
$ h5dump -H SVM12_npp_d20140410_t2130000_e2140001_b12705_c20140414174333267547_star_dev.h5
...
   GROUP "Data_Products" {
      GROUP "VIIRS-M12-SDR" {
         DATASET "VIIRS-M12-SDR_Aggr" {
            DATATYPE  H5T_REFERENCE { H5T_STD_REF_OBJECT }
            DATASPACE  SIMPLE { ( 5 ) / ( 5 ) }
         }
...
```
Next, run the VIIRS destriping program. One OpenMP thread is fast enough for us:
```
$ export OMP_NUM_THREADS=1
$ ./destripe_viirs SVM12_npp_d20140410_t2130000_e2140001_b12705_c20140414174333267547_star_dev.h5 SVM12_npp_d20140410_t2130000_e2140001_b12705_c20140414174333267547_star_dev.h5 destriping_viirs_param.txt
destripe_viirs SVM12_npp_d20140410_t2130000_e2140001_b12705_c20140414174333267547_star_dev.h5 destriping_viirs_param.txt
numthreads = 1 maxthreads =  1
5 16 8 0.001000 0.020000 0.050000 0.000000 0.250000
7 16 8 0.001000 0.020000 0.050000 0.000000 0.250000
10 16 8 0.001000 0.020000 0.050000 0.000000 0.250000
12 16 8 0.050000 0.300000 0.500000 270.000000 320.000000
13 16 8 0.200000 0.300000 0.500000 270.000000 320.000000
14 16 8 0.050000 0.300000 0.500000 270.000000 320.000000
15 16 8 0.050000 0.300000 0.500000 270.000000 320.000000
16 16 8 0.050000 0.300000 0.500000 270.000000 320.000000
Band = 12
12 16 8 0.050000 0.300000 0.500000 270.000000 320.000000
Destriping atribute location = Data_Products/VIIRS-M12-SDR/VIIRS-M12-SDR_Aggr
Destriping atribute name = DestripingBrightnessTemperature
Data location = All_Data/VIIRS-M12-SDR_All/BrightnessTemperature
Corresponding geofile = GMODO_npp_d20140410_t2130000_e2140001_b12705_c20140414174333267547_star_dev.h5
nx = 3200 ny = 5408
scale = 0.002518 offset = 203.000000
fftw3 initialized on 1 threads
Iteration number: 0
Iteration number: 1
Iteration number: 2
Iteration number: 3
Iteration number: 4
Iteration number: 5
Iteration number: 6
Iteration number: 7
Average correction = -9.832491e-05
VIIRS_band    NIF        NDF      domain_size
    12      0.071900  0.956226      3654949
```
Now, if we look at the dataset, we see that the attribute `DestripingBrightnessTemperature` was added to it, which means it was successfully destriped:
```
$ h5dump -H SVM12_npp_d20140410_t2130000_e2140001_b12705_c20140414174333267547_star_dev.h5
...
   GROUP "Data_Products" {
      GROUP "VIIRS-M12-SDR" {
         DATASET "VIIRS-M12-SDR_Aggr" {
            DATATYPE  H5T_REFERENCE { H5T_STD_REF_OBJECT }
            DATASPACE  SIMPLE { ( 5 ) / ( 5 ) }
            ATTRIBUTE "DestripingBrightnessTemperature" {
               DATATYPE  H5T_IEEE_F32LE
               DATASPACE  SCALAR
            }
         }
...
```

### License

Copyright 2014, NOAA STAR SST Team. All Rights Reserved.
See accompanying [LICENSE.txt](LICENSE.txt) file for use and distribution.
