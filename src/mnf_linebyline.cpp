#include <stdio.h>
#include <string.h> 
#include <sys/time.h>
#include "mnf.h"
#include <iostream>
#include <sys/time.h>
#include <cblas.h>
#include "mnf_linebyline.h"
using namespace std;


float *estimateNoise(float *line, int samples, int bands){
	int newSamples = samples-1;
	float *noise = new float[newSamples*bands];
	for (int i=0; i < bands; i++){
		for (int j=0; j < newSamples; j++){
			noise[i*newSamples + j] = line[i*samples + j] - line[i*samples + (j + 1)];
		}
	}
	return noise;
}

float *estimateENVINoise(float *line, float *prevLine, int samples, int bands){
	int newSamples = samples-1;
	float *noise = new float[newSamples*bands];
	for (int i=0; i < bands; i++){
		for (int j=0; j < newSamples; j++){
			noise[i*newSamples + j] = line[i*samples + j] - 0.5*(line[i*samples + (j + 1)] + prevLine[i*samples + j]);
		}
	}
	return noise;
}



void removeMean(MnfWorkspace *workspace, float *line, float *means, int samples, int bands){
	float *ones = workspace->ones_samples;
	cblas_sger(CblasRowMajor, bands, samples, -1.0f, means, 1, ones, 1, line, samples);
}

void addMean(MnfWorkspace *workspace, float *line, float *means, int samples, int bands){
	float *ones = workspace->ones_samples;
	cblas_sger(CblasRowMajor, bands, samples, 1.0f, means, 1, ones, 1, line, samples);
}

void mnf_oneline(float *line, MnfWorkspace *workspace, Covariance *imageCov, Covariance *noiseCov){
	int bands = workspace->img->getBands();
	int samples = workspace->img->getPixels();

	float *noise;

	//covariances

	//image
	updateStatistics(workspace, line, samples, bands, imageCov);
	
	//noise
	noise = estimateNoise(line, samples, bands);
	updateStatistics(workspace, noise, samples-1, bands, noiseCov);
	

	
	float *means_float = new float[bands];
	float *noisemeans_float = new float[bands];
	for (int i=0; i < bands; i++){
		means_float[i] = imageCov->means[i];
		noisemeans_float[i] = noiseCov->means[i];
	}

	for (int i=0; i < bands*bands; i++){
		workspace->imgCov[i] = imageCov->C[i]/(imageCov->n*1.0f);
		workspace->noiseCov[i] = noiseCov->C[i]/(noiseCov->n*1.0f);
	}
	
	//find forward and inverse transformation matrices	
	calculateForwardTransfMatrix(workspace);
	float *forwardTransfArr = workspace->forwardTransf;
	
	calculateInverseTransfMatrix(workspace);
	float *inverseTransfArr = workspace->inverseTransf;
	
	//remove mean from data
	removeMean(workspace, line, means_float, samples, bands);

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
	addMean(workspace, submatrNew, means_float, samples, bands);

	memcpy(line, submatrNew, sizeof(float)*samples*bands);
	delete [] submatrNew;

	delete [] invTrMultR;
	delete [] invTrMultRMultForTr;
	delete [] noise;
	delete [] means_float;
}

void mnf_linebyline(float *data, int lines, int samples, int bands, int numBands, string mnfOutfilename){

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

void initializeCov(Covariance *cov, int bands, int samples){
	cov->n = 0;
	cov->C = new float[bands*bands];

	cov->yiyj = new float[bands*bands];

	for (int i=0; i < bands*bands; i++){
		cov->yiyj[i] = 0.0f;
		cov->C[i] = 0.0f;
	}

	cov->numSamplesInMeans = 0;
	cov->means = new double[samples];
	
	for (int i=0; i < samples; i++){
		cov->means[i] = 0.0f;
	}
}

void deinitializeCov(Covariance *cov){
	delete [] cov->yiyj;
	delete [] cov->means;
}


void updateStatistics(MnfWorkspace *workspace, float *bilData, int samples, int bands, Covariance *cov){
	cov->n += samples;
	float *C = cov->C;
	double *means = cov->means;

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
		means[i] = means[i] + 1.0*samples*(meanTemp[i] - means[i])/(1.0*cov->n);
	}

	//update to actual covariance
	cblas_ssyr(CblasRowMajor, CblasLower, bands, samples*(cov->n - samples)/(cov->n), meanDiff, 1, C, bands);



	delete [] meanDiff;
	delete [] meanTemp;
	delete [] tempLine;
}


//only used by conventional MNF. conventional MNF should use the mnf-line-by-line update algorithm as they are more numerically stable, instead of these.  
void updateCovariances(float *bilData, int samples, int bands, Covariance *imCov){
	imCov->n += samples;
	//yi yj
	//cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, bands, bands, samples, 1.0f, bilData, samples, bilData, samples, 1.0f, imCov->yiyj, bands);

	//gsl_matrix_float_view A = gsl_matrix_float_view_array(bilData, bands, samples);
	//gsl_matrix_float_view C = gsl_matrix_float_view_array(imCov->yiyj, bands, bands);

	//old version without dividing by N-tjafs, might be it is uneccessary anyway since N is so large that it will only be a shortening of the potens
	//gsl_blas_ssyrk(CblasLower, CblasNoTrans, 1.0f, &A.matrix, 1.0f, &C.matrix);

	//old version of update: cannot use since yiyj must be double due to precision issues
	//FIXME: This is not really numerically stable, and might be cause for future errors
	//alternativ måte å gjøre det på er å ta Cx = Ca + Cb + (bar{x_A} - bar{x_B})(bar{y_A} - bar{y_B})*nA*nB/nX

	//first an optimized matrix multiplication
	float *C = new float[bands*bands];

	//C = bilData * bilData^T
	cblas_ssyrk(CblasRowMajor, CblasLower, CblasNoTrans, bands, samples, 1.0f, bilData, samples, 0.0f, C, bands);
	
	//rather unoptimized update of the actual yiyj array
	for (int i=0; i < bands*bands; i++){
		imCov->yiyj[i] = (imCov->n - samples)*1.0f/(imCov->n*1.0f)*imCov->yiyj[i] + 1.0f/(imCov->n*1.0f)*C[i];
	}
	delete [] C;
}

void updateMeans(MnfWorkspace *workspace, float *data, int samples, int bands, double *means, int *numMeans){
	//find mean across current line, hope that float is enough precision
	
	float *meanTemp = new float[bands];
	cblas_sgemv(CblasRowMajor, CblasNoTrans, bands, samples, 1.0f, data, samples, workspace->ones_samples, 1, 0.0f, meanTemp, 1);

	*numMeans += samples;

	//update actual mean array (doing this explicitly, I don't trust sgemv update)
	//also some shit to partially ensure numerical stability
	for (int i=0; i < bands; i++){
		//means[i] = means[i]*1.0*(*numMeans - samples)/(1.0*(*numMeans)) + meanTemp[i]*1.0/(1.0*(*numMeans));
		means[i] = means[i] + 1.0*(meanTemp[i] - samples*means[i])/(1.0*(*numMeans));
	}
	delete [] meanTemp;
	/*for (int i=0; i < bands; i++){
		for (int j=0; j < samples; j++){
			if (i==0){
				(*numMeans)++;
			}
			means[i] = means[i] + (data[i*samples + j] - means[i])/(*numMeans*1.0);
		}
	}*/
}


void calculateTotCovariance(Covariance *imageCov, float *means, int bands, float *covMat){
	size_t size = bands*bands*sizeof(float);
	float *tempCov = (float*) malloc(size);
	memcpy(tempCov, imageCov->yiyj, size);

	cblas_ssyr(CblasRowMajor, CblasLower, bands, -1.0f, means, 1, imageCov->yiyj, bands);
	memcpy(covMat, imageCov->yiyj, size);
	memcpy(imageCov->yiyj, tempCov, size);
	free(tempCov);
}
