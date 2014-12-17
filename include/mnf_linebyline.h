#ifndef MNF_LINEBYLINE_H_DEFINED
#define MNF_LINEBYLINE_H_DEFINED

#include "mnf.h"

void mnf_linebyline(float *data, int lines, int samples, int bands, int numBands, std::string mnfOutFilename);



typedef struct{
	int n; //number of pixels so far summed over

	//float *y;
	float *yiyj;

	float *C; //not used?
	float *means_prev;
	double *means;
	int numSamplesInMeans;
} Covariance;

//run mnf on one line, incrementally updating the covariances by this line and denoising line in place
void mnf_oneline(float *line, MnfWorkspace *workspace, Covariance *imageCov, Covariance *noiseCov);

//create covariance arrays
void initializeCov(Covariance *cov, int bands, int samples);

//delete covariance arrays
void deinitializeCov(Covariance *cov);

//update covariance matrices with a new line
void updateCovariances(float *data, int samples, int bands, Covariance *imCov);

//update means with a new line
void updateMeans(MnfWorkspace *workspace, float *data, int samples, int bands, double *means, int *numMeans);

//update covariances and means using a numerically stable algorithm
void updateStatistics(MnfWorkspace *workspace, float *bilData, int samples, int bands, Covariance *cov);

//calculate total covariance from the partial sums and means
void calculateTotCovariance(Covariance *cov, float *means, int bands, float *outputCovMat);

//explicitly estimate noise by shift difference
//since we are doing this line by line, there is no point to doing it implicitly in order to save RAM, except for jumbling more shit inside the same for loop
float *estimateNoise(float *line, int samples, int bands);

//estimate the ENVI way using neighboring values in current line and the sample above in the previous line
float *estimateENVINoise(float *line, float *prevLine, int samples, int bands);

#endif
