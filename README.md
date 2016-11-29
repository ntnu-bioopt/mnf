Line-by-line noise removal of hyperspectral images using an algorithm based on the MNF
transform.

If you use this software, we would be grateful if you could cite the following paper:

[A. Bjorgan, L. L. Randeberg,
"Real-time noise removal for line-scanning hyperspectral devices using a
Minimum Noise Fraction-based approach", Sensors 15(2), pp. 3362-3378 (2015).
doi:10.3390/s150203362](http://www.mdpi.com/1424-8220/15/2/3362). 

Underlying theory can be found in the same reference. 


Running
-------

`mnf` is run from the command line. `./mnf filename --option1 --option2...`, or `mnf.exe filename --option1 --option2 ...` .  See `./mnf --help` for
options. BIL-interleaved ENVI images are assumed. 

The line-by-line algorithm is run when the `--line-by-line option` is set. When
it is not set, the ordinary MNF transform is run. 

The conventional MNF transform is based on Green, A. A., Berman, M., Switzer,
P., and Craig, M. D., "A transformation for ordering multispectral data in
terms of image quality with implications for noise removal", IEEE Transactions
on Geoscience and Remote Sensing 26(1), pp. 65-74 (1988).

Building
--------

This application is built using [CMake](http://cmake.org), producing `mnf`, the main application executable, and `libmnf`, for using the algorithms in other applications. 

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`
5. `make install`

The output of `cmake ..` should give indication of missing libraries or the need to
add libraries to the `CMakeLists.txt` file. 

Requirements
------------

This program requires a LAPACKE and BLAS implementation.

MKL can be enabled over the standard BLAS implementations available on the system by running `cmake -DUSE_MKL_LIBRARIES=True ..` and rebuilding.

GNU Regex is required for src/readimage.cpp. 
