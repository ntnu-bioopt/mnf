#include <stdio.h>
#include <string.h> 
#include <sys/time.h>
#include "mnf.h"
#include <iostream>
#include <sys/time.h>
#include <cblas.h>
#include "mnf_linebyline.h"
using namespace std;


void mnf_linebyline_estimate_noise(int bands, int samples, float *line, float **noise_est, int *noise_samples){
	*noise_samples = samples-1;
	*noise_est = new float[*noise_samples*bands];
	for (int i=0; i < bands; i++){
		for (int j=0; j < *noise_samples; j++){
			*noise_est[i*(*noise_samples) + j] = line[i*samples + j] - line[i*samples + (j + 1)];
		}
	}
}


void imagestatistics_get_means(ImageStatistics *stats, int numBands, float *means){
	for (int i=0; i < bands; i++){
		means[i] = stats->means[i];
	}
}

void imagestatistics_get_covariance(ImageStatistics *stats, int numBands, float *cov){
	for (int i=0; i < numBands*numBands; i++){
		cov[i] = stats->C[i]/(stats->n*1.0f);
	}
}

void mnf_linebyline_remove_mean(const MnfWorkspace *workspace, float *means, int bands, int samples, float *line){
	cblas_sger(CblasRowMajor, bands, samples, -1.0f, means, 1, workspace->ones_samples, 1, line, samples);
}

void mnf_linebyline_add_mean(const MnfWorkspace *workspace, float *means, int bands, int samples, float *line){
	cblas_sger(CblasRowMajor, bands, samples, 1.0f, means, 1, workspace->ones_samples, 1, line, samples);
}


void mnf_linebyline_run_oneline(MnfWorkspace *workspace, int bands, int samples, float *line, ImageStatistics *imageStats, ImageStatistics *noiseStats){
	float *noise;
	int noise_samples = 0;

	//image statistics
	imagestatistics_update_statistics_with_line(workspace, bands, samples, line, imageStats);
	
	//noise statistics
	mnf_linebyline_estimate_noise(bands, samples, line, &noise, &noise_samples);
	imagestatistics_update_statistics_with_line(workspace, bands, noise_samples, noise, noiseStats);

	//get means	
	float *means_float = new float[bands];
	float *noisemeans_float = new float[bands];
	imagestatistics_get_means(noiseStats, bands, noisemeans_float);
	imagestatistics_get_means(imageStats, bands, means_float);

	//get true covariance
	imagestatistics_get_covariance(imageStats, bands, workspace->imgCov);
	imagestatistics_get_covariance(noiseStats, bands, workspace->noiseCov);

	//find forward and inverse transformation matrices	
	mnf_calculate_forward_transf_matrix(bands, workspace->imgCov, workspace->noiseCov, workspace->forwardTransf);
	mnf_calculate_inverse_transf_matrix(bands, workspace->forwardTransf, workspace->inverseTransf);
	float *forwardTransfArr = workspace->forwardTransf;
	float *inverseTransfArr = workspace->inverseTransf;
	
	//remove mean from data
	mnf_linebyline_remove_mean(workspace, means_float, bands, samples, line);

	//pre-calculate matrix with which to multiply the image for noise removal
	float *invTrMultR = new float[bands*bands];
	cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, bands, bands, bands, 1.0f, inverseTransfArr, bands, workspace->R, bands, 0.0f, invTrMultR, bands);

	float *invTrMultRMultForTr = new float[bands*bands];
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, bands, bands, bands, 1.0f, invTrMultR, bands, forwardTransfArr, bands, 0.0f, invTrMultRMultForTr, bands);

	//actual noise removal	
	float *submatrNew = new float[bands*samples];
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, bands, samples, bands, 1.0f, invTrMultRMultForTr, bands, line, samples, 0.0f, submatrNew, samples);
	//cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, bands, samples, bands, 1.0f, forwardTransfArr, bands, line, samples, 0.0f, submatrNew, samples); //if we would only want the transform, which we don't anyway
	
	//add mean to data
	mnf_linebyline_add_mean(workspace, means_float, bands, samples, submatrNew);

	memcpy(line, submatrNew, sizeof(float)*samples*bands);
	delete [] submatrNew;

	delete [] invTrMultR;
	delete [] invTrMultRMultForTr;
	delete [] noise;
	delete [] means_float;
}

void mnf_linebyline_run_image(int bands, int samples, int lines, float *data, std::string mnfOutFilename){
	Covariance noiseCov;
	Covariance imageCov;

	initializeCov(&noiseCov, bands, samples);
	initializeCov(&imageCov, bands, samples);

	MnfWorkspace workspace;
	initializeMnfWorkspace(&workspace, numBands, samples, bands);

	for (int i=0; i < lines; i++){
		float *line = data + i*samples*bands;
		fprintf(stderr, "Line %d/%d\n", i, lines);
		mnf_oneline(line, &workspace, &imageCov, &noiseCov);
	}
	deinitializeMnfWorkspace(&workspace);	
	deinitializeCov(&noiseCov);
	deinitializeCov(&imageCov);

	Hyperspectral img(data, lines, samples, bands, BIL);
	img.writeToFile(mnfOutfilename + string("_inversetransformed"));
}

void imagestatistics_initialize(ImageStatistics *stats, int bands){
	stats->n = 0;
	stats->C = new float[bands*bands];
	for (int i=0; i < bands*bands; i++){
		stats->C[i] = 0.0f;
	}

	stats->means = new double[bands];
	for (int i=0; i < bands; i++){
		stats->means[i] = 0.0f;
	}
}

void imagestatistics_deinitialize(ImageStatistics *stats){
	delete [] stats->C;
	delete [] stats->means;
}


void imagestatistics_update_with_line(const MnfWorkspace *workspace, int numBands, int numSamples, float *bilData, ImageStatistics *stats){
	stats->n += samples;
	float *C = stats->C;
	double *means = stats->means;

	//copy line to temporary variable
	float *tempLine = new float[samples*bands];
	memcpy(tempLine, bilData, sizeof(float)*samples*bands);

	//estimate mean of single line
	float *meanTemp = new float[bands];
	cblas_sgemv(CblasRowMajor, CblasNoTrans, bands, samples, 1.0f/(1.0f*samples), tempLine, samples, workspace->ones_samples, 1, 0.0f, meanTemp, 1);

	//subtract mean from line
	removeMean(workspace, tempLine, meanTemp, samples, bands);

	//calculate covariance for current line and add to accumulated covariance
	cblas_ssyrk(CblasRowMajor, CblasLower, CblasNoTrans, bands, samples, 1.0f, tempLine, samples, 1.0f, C, bands);

	//find the difference between the mean of the current line and the old mean
	//update total means
	float *meanDiff = new float[bands];
	for (int i=0; i < bands; i++){
		meanDiff[i] = meanTemp[i] - means[i];
		means[i] = means[i] + 1.0*samples*(meanTemp[i] - means[i])/(1.0*stats->n);
	}

	//update to actual covariance
	cblas_ssyr(CblasRowMajor, CblasLower, bands, samples*(stats->n - samples)/(stats->n), meanDiff, 1, C, bands);



	delete [] meanDiff;
	delete [] meanTemp;
	delete [] tempLine;
}
