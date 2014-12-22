
## Destriping algorithm for VIIRS/MODIS Ocean Products.

This software was developed by the SST group at NOAA STAR (Lead:
A. Ignatov). It is designed to destripe VIIRS Level 1b (SDR) data.
See the following publication for more details on the algorithm:

    Bouali, Marouan, and Alexander Ignatov, 2014: Adaptive Reduction
    of Striping for Improved Sea Surface Temperature Imagery from Suomi
    National Polar-Orbiting Partnership (S-NPP) Visible Infrared Imaging
    Radiometer Suite (VIIRS).  J. Atmos. Oceanic Technol., 31, 150–163. doi:
    http://dx.doi.org/10.1175/JTECH-D-13-00035.1

## Dependencies

Besides standard C library functions, this code relies on the following
libraries:
* [fftw3 library](http://www.fftw.org/) (float version),
* [HDF5 library](http://www.hdfgroup.org/HDF5/) (and dependencies, such as libsz, libz),
* libpthread, and
* OpenMP library (such as libiomp5 or [libgomp](https://gcc.gnu.org/projects/gomp/))

The code is designed to run multiple threads using OpenMP. 
The number of threads can be set before execution by
setting environment variable `OMP_NUM_THREADS`.

## Installation

To install, run `make install` inside destripe_viirs directory. You may
need to adapt the `Makefile` for your system. The installed binary is named
`destripe_viirs`.

## Usage

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

## License

Copyright 2014, NOAA STAR SST Team. All Rights Reserved.
See accompanying [LICENSE.txt](LICENSE.txt) file for use and distribution.
