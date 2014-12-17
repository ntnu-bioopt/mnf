#ifndef MNF_C_H_DEFINED
#define MNF_C_H_DEFINED

/* Author: Asgeir Bjørgan
asgeir.bjorgan@iet.ntnu.no
NTNU */

#include <string>
#include "hyperspectral.h"
#include <vector>
#include <gsl/gsl_matrix.h>
#include <vector>

//struct containing variables and arrays neccessary for performing the MNF transform
typedef struct{
	Hyperspectral *img;
	std::string basefilename; //base filename for output
	float *imgCov; //image convariance
	float *noiseCov; //noise covariance
	float *forwardTransf; //forward MNF matrix. NB: Stored as transposed. 
	float *inverseTransf; //inverse MNF matrix. NB: Stored as transposed. 
	float *means; //band means, removed from data, to be added again
	std::vector<float> wlens; //wavelengths, for writing to file
	float *R; //matrix with which to pre-multiply transformed image before inverse transform (first m elements bla bla)
	float *ones_samples; //vector of length samples consisting only of ones, for convenience
} MnfWorkspace;

void initializeMnfWorkspace(MnfWorkspace *workspace, int numBandsInInverse, int samples, int bands, int lines = -1, std::string basefilename = std::string(), float *hyData = NULL, std::vector<float> wlens = std::vector<float>());
void deinitializeMnfWorkspace(MnfWorkspace *workspace);

//enum specifying which direction to run the MNF transform
enum TransformDirection{RUN_BOTH, RUN_FORWARD, RUN_INVERSE};

//used in covariance estimation. ESTIMATE_NOISE for calculating noise instead of using the values directly, USUAL_VAL for using values from image array directly without no processing
enum WhatValue{ESTIMATE_NOISE, USUAL_VAL};

//data: assumed bil-interleaved (mostly abstracted away into the Hyperspectral class, but be aware that the forward and inverse transformations are explicitly using BIL-interleavedness)
//perform forward and mnf transform (whether both or one direction depends on the dir variable)
void run_mnf(TransformDirection dir, int numBands, float *data, std::vector<float> wlens, int samples, int bands, int lines, std::string basefilename);

//calculate forwardTransf in container based on imgCov and noiseCov
void calculateForwardTransfMatrix(MnfWorkspace *container);
void calculateInverseTransfMatrix(MnfWorkspace *container);

//returns matrix of means for later mean addition
//removes band mean from each band image
void removeMean(MnfWorkspace *container);

//add supplied band mean to each band image
void addMean(Hyperspectral img, float* means);

//noise estimation, shift difference
//not used, since the noise is implicitly calculated in the covariance estimation
Hyperspectral* estimateNoise(Hyperspectral img, int startSamp, int endSamp, int startLine, int endLine);

//calculate covariance matrix. WhatValue determines whether the noise covariance is to be found, and the noise implicitly calculated. 
//noise is estimated by shift differencing neighboring _columns_
void findCov(MnfWorkspace *workspace, Hyperspectral img, float *cov, WhatValue whichVal = USUAL_VAL);

#endif
