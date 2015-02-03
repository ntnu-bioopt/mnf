Line-by-line noise removal of hyperspectral images using an algorithm based on the MNF
transform.

For use of this source code, please cite A. Bjorgan, L. L. Randeberg,
"Real-time noise removal for line-scanning hyperspectral devices using a
Minimum Noise Fraction-based approach", Sensors 15(2), pp. 3362-3378 (2015).
doi:10.3390/s150203362. URL: http://www.mdpi.com/1424-8220/15/2/3362. 
Theory can be found in the same paper. 

RUNNING 

mnf is run from the command line. ./mnf filename --option1 --option2
..., or mnf.exe filename --option1 --option2 ... .  See ./mnf --help for
options. BIL-interleaved ENVI images are assumed. 

The line-by-line algorithm is run when the --line-by-line option is set. When
it is not set, the ordinary MNF transform is run. 

The conventional MNF transform is based on Green, A. A., Berman, M., Switzer,
P., and Craig, M. D., "A transformation for ordering multispectral data in
terms of image quality with implications for noise removal", IEEE Transactions
on Geoscience and Remote Sensing 26(1), pp. 65-74 (1988).

BUILDING

This application is built using CMake (http://cmake.org). Create a new
subdirectory build/, navigate to it and run cmake .. to generate makefiles. Run
make or nmake.exe or equivalent. 

REQUIREMENTS

This program requires some LAPACKE and BLAS implementation. This is currently
configured to Intel MKL in the CMakeLists.txt file. 

GNU Regex is required for src/readimage.cpp. 
