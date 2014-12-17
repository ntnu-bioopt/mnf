
/* Author: Asgeir Bjørgan
asgeir.bjorgan@iet.ntnu.no
NTNU */

#include <string>
#include <iostream>
#include <fstream>
#include "mnf.h"
#include "hyperspectral.h"
#include <string.h>
#include <sys/time.h>
#include <stdio.h>
#include <pthread.h>
#include <sstream>
extern "C"
{
#include <lapacke.h>
#include <cblas.h>
}
#include "mnf_linebyline.h"
using namespace std;


//find mean across band image for specified band
float findBandMean(Hyperspectral img, int band);

//print supplied covariance to file
void printCov(float *cov, int bands, string filename);
void printMeans(float *means, int bands, string filename);

//read covariances and band means from basefilename_imgcov.dat and basefilename_noisecov.dat		
void readStatistics(MnfWorkspace *container);

//for running forward and inverse MNF
void runForwardMNF(MnfWorkspace *container);
void runInverseMNF(MnfWorkspace *container);

//container functions in order to run covariance estimation inside pthreads
void *covImage(void *param){
	fprintf(stderr, "Image covariance\n");
	MnfWorkspace *input = (MnfWorkspace*)param; //throw all type safety out of the window!

	//find covariance matrix of the noisy image
	findCov(input, *(input->img), input->imgCov);
}

void *covNoise(void *param){
	fprintf(stderr, "Noise covariance\n");
	MnfWorkspace *input = (MnfWorkspace*)param;

	//find covariance matrix of the noise
	findCov(input, *(input->img), input->noiseCov, ESTIMATE_NOISE);
}

void run_mnf(TransformDirection dir, int numBands, float *data, vector<float> wlens, int samples, int bands, int lines, string basefilename){
	//initialization
	timeval t1,t2;
	gettimeofday(&t1, NULL);

	MnfWorkspace container;
	initializeMnfWorkspace(&container, numBands, samples, bands, lines, basefilename, data, wlens);

	
	switch (dir){
		case RUN_FORWARD:
			//estimate covariances, run forward transformation
			runForwardMNF(&container);
		break;
		case RUN_INVERSE:
			//read in covariances from file
			readStatistics(&container);
			calculateForwardTransfMatrix(&container);

			//run inverse transformation
			runInverseMNF(&container);
		break;
		case RUN_BOTH:
			//estimate covariances, run forward transformation
			runForwardMNF(&container);
			//run inverse transformation
			runInverseMNF(&container);
		break;	
	}

	//cleanup
	//not deleting the data array, main thread calling this function should handle it
	deinitializeMnfWorkspace(&container);
	gettimeofday(&t2, NULL);
	cout << (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec)*1.0e-06 << endl;
}	

void initializeMnfWorkspace(MnfWorkspace *workspace, int numBandsInInverse, int samples, int bands, int lines, string basefilename, float *hyData, std::vector<float> wlens){
	workspace->img = new Hyperspectral(hyData, lines, samples, bands, BIL);
	workspace->imgCov = new float[bands*bands];
	workspace->noiseCov = new float[bands*bands];
	workspace->forwardTransf = new float[bands*bands];
	workspace->inverseTransf = new float[bands*bands];
	workspace->wlens = wlens;
	workspace->basefilename = basefilename;
	workspace->means = new float[bands];

	workspace->ones_samples = new float[samples];
	for (int i=0; i < samples; i++){
		workspace->ones_samples[i] = 1.0f;
	}
	
	workspace->R = new float[bands*bands]();
	for (int i=0; i < numBandsInInverse; i++){
		workspace->R[i*bands + i] = 1.0f;
	}
}

void deinitializeMnfWorkspace(MnfWorkspace *workspace){
	delete workspace->img;
	delete [] workspace->imgCov;
	delete [] workspace->noiseCov;
	delete [] workspace->forwardTransf;
	delete [] workspace->inverseTransf;
	delete [] workspace->means;
	delete [] workspace->R;
	
}

void runForwardMNF(MnfWorkspace *container){
	//FIXME: if you want to use a smaller subset for noise estimation, you should call noiseEstimation with start and end samples and lines, and call findCov on the noise image without the NOISE_ESTIMATION-option

	float *imgCov = container->imgCov;
	float *noiseCov = container->noiseCov;
	Hyperspectral *img = container->img;
	float *forwardTransf = container->forwardTransf;
	float *inverseTransf = container->inverseTransf;
	int bands = img->getBands();
	
	
	fprintf(stderr, "Removing band means from image\n");
	removeMean(container); //remove mean from data
	printMeans(container->means, container->img->getBands(), container->basefilename + "_bandmeans.dat");

	//create threads for covariance estimation
	fprintf(stderr, "Estimate noise and image covariances\n");

	#pragma omp parallel
	{	
		#pragma omp sections
		{
			#pragma omp section
			{ 
				covImage((void*)(container));
			}
			#pragma omp section
			{
				covNoise((void*)(container));
			}
		}
	}
	
	/*pthread_t thImg, thNoise;
	pthread_create(&thImg, NULL, covImage, (void*)(container));
	pthread_create(&thNoise, NULL, covNoise, (void*)(container));

	pthread_join(thImg, NULL);
	pthread_join(thNoise, NULL);*/

	//print to file
	printCov(imgCov, img->getBands(), container->basefilename + "_imgcov.dat");
	printCov(noiseCov, img->getBands(), container->basefilename + "_noisecov.dat");

	fprintf(stderr, "Solve eigenvalue problem\n");
	calculateForwardTransfMatrix(container);
	
	

	//run forward transform
	float *data = img->getAllData();
	fprintf(stderr, "Running forward transformation\n");

	#pragma omp parallel for
	for (int i=0; i < img->getLines(); i++){
		//import submatrix into a gsl_view
		float *submatr = img->getAllData() + i*img->getPixels()*img->getBands();
		float *submatrNew = new float[img->getPixels()*img->getBands()];
		cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, img->getBands(), img->getPixels(), img->getBands(), 1.0f, forwardTransf, img->getBands(), submatr, img->getPixels(), 0.0f, submatrNew, img->getPixels());
		memcpy(data + i*img->getPixels()*img->getBands(), submatrNew, sizeof(float)*img->getPixels()*img->getBands());
		delete [] submatrNew;
	}

	//write forward transform to file
	fprintf(stderr, "Writing to file\n");
	img->writeToFile(container->basefilename + string("_transformed"), container->wlens);

}	

void calculateForwardTransfMatrix(MnfWorkspace *container){
	float *imgCov = container->imgCov;
	float *noiseCov = container->noiseCov;
	Hyperspectral *img = container->img;
	
	float *forwardTransf = container->forwardTransf;
	int bands = img->getBands();

	//copy covariances to temp variables, since dsygv will destroy them 
	float *noiseCovTemp = new float[bands*bands];
	float *imgCovTemp = new float[bands*bands];
	memcpy(noiseCovTemp, noiseCov, sizeof(float)*bands*bands);
	memcpy(imgCovTemp, imgCov, sizeof(float)*bands*bands);
	

	//solve eigenvalue problem
	int itype = 1; //solves Ax = lam * Bx
	char jobz = 'V'; //compute both eigenvalues and eigenvectors
	char uplo = 'L'; //using lower triangles of matrices
	float *eigvals = new float[bands];
	int err = LAPACKE_ssygv(LAPACK_ROW_MAJOR, itype, jobz, uplo, bands, noiseCov, bands, imgCov, bands, eigvals);
	//int err = LAPACKE_dsygv(LAPACK_ROW_MAJOR, itype, jobz, uplo, bands, imgCov, bands, noiseCov, bands, eigvals);
	if (err != 0){
		fprintf(stderr, "LAPACKE_dsygv failed: calculation of forward MNF transformation matrix, err code %d\n", err);
	}

	//print eigenvalues to file
	if (container->basefilename != string()){
		ofstream eigvalsFile;
		eigvalsFile.open(string(container->basefilename + "_eigvals.dat").c_str());
		for (int i=0; i < bands; i++){
			eigvalsFile << eigvals[i] << endl;
		}
		eigvalsFile.close();
	}

	//copy back covariances, move transfer matrix to transformation array
	memcpy(forwardTransf, noiseCov, sizeof(float)*bands*bands);
	memcpy(noiseCov, noiseCovTemp, sizeof(float)*bands*bands);
	memcpy(imgCov, imgCovTemp, sizeof(float)*bands*bands);

	delete [] noiseCovTemp;
	delete [] imgCovTemp;
	delete [] eigvals;
}

void calculateInverseTransfMatrix(MnfWorkspace *container){
	size_t size = container->img->getBands()*container->img->getBands()*sizeof(float);
	
	//keep forward transf matrix
	float *forwardTransfTemp = (float*)malloc(size);
	memcpy(forwardTransfTemp, container->forwardTransf, size);

	//find inverse of forward eigenvector transfer matrix
	int *ipiv = new int[container->img->getBands()];
	int err = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, container->img->getBands(), container->img->getBands(), container->forwardTransf, container->img->getBands(), ipiv);
	if (err != 0){
		fprintf(stderr, "LU decomposition failed: calculation of inverse transformation matrix\n");
	}
	err = LAPACKE_sgetri(LAPACK_ROW_MAJOR, container->img->getBands(), container->forwardTransf, container->img->getBands(), ipiv);
	if (err != 0){
		fprintf(stderr, "Inversion failed: calculation of inverse transformation matrix\n");
	}
	delete [] ipiv;

	//copy back
	memcpy(container->inverseTransf, container->forwardTransf, size);
	memcpy(container->forwardTransf, forwardTransfTemp, size);
	
	free(forwardTransfTemp);
}

void runInverseMNF(MnfWorkspace *container){
	//prepare inverse transformation matrix
	calculateInverseTransfMatrix(container);

	Hyperspectral *img = container->img;

	float *inverseTransfArr = container->inverseTransf;
	
	//post-multiply R matrix to get k < m first bands in inverse transformation
	float *R = container->R;
	float *inverseTransfArrShort = new float[img->getBands()*img->getBands()];
	cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, img->getBands(), img->getBands(), img->getBands(), 1.0f, inverseTransfArr, img->getBands(), R, img->getBands(), 0.0f, inverseTransfArrShort, img->getBands());

	//run inverse transformation
	float *data = img->getAllData();
	fprintf(stderr, "Running inverse transformation\n");
	#pragma omp parallel for
	for (int i=0; i < img->getLines(); i++){
		float *submatr = data + i*img->getPixels()*img->getBands();
		float *submatrNew = new float[img->getBands()*img->getPixels()];
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, img->getBands(), img->getPixels(), img->getBands(), 1.0f, inverseTransfArrShort, img->getBands(), submatr, img->getPixels(), 0.0f, submatrNew, img->getPixels());

		memcpy(data + i*img->getPixels()*img->getBands(), submatrNew, sizeof(float)*img->getPixels()*img->getBands());
		delete [] submatrNew;
	}
	
	//adding back band means
	fprintf(stderr, "Adding back band means\n");
	addMean(*(container->img), container->means);

	//writing inverse transform to file
	fprintf(stderr, "Writing to file\n");
	img->writeToFile(container->basefilename + string("_inversetransformed"), container->wlens);

	delete [] inverseTransfArrShort;
}

void removeMean(MnfWorkspace *container){
	Hyperspectral *img = container->img;
	int lines = img->getLines();
	int samples = img->getPixels();
	int bands = img->getBands();
	float *means = container->means;

	for (int j=0; j < bands; j++){
		//find band mean
		float bandMean = findBandMean(*img, j);

		//remove band mean
		for (int i=0; i < lines; i++){
			for (int k=0; k < samples; k++){
				img->write(i, k, j, img->read(i,k,j) - bandMean);
			}
		}
		means[j] = bandMean;
	}
}

void addMean(Hyperspectral img, float* means){
	int lines = img.getLines();
	int samples = img.getPixels();
	int bands = img.getBands();
	for (int i=0; i < lines; i++){
		for (int j=0; j < bands; j++){
			for (int k=0; k < samples; k++){
				float mean = means[j];
				img.write(i, k, j, img.read(i, k, j) + mean);
			}
			
		}
	}
}

float findBandMean(Hyperspectral img, int band){
	double mean = 0;
	int num=0;
	for (int line=0; line < img.getLines(); line++){
		for (int sample = 0; sample < img.getPixels(); sample++){
			num++;
			mean = mean + (img.read(line, sample, band) - mean)/(1.0*num);
		}
	}
	return mean*1.0f;
}

void readStatistics(MnfWorkspace *container){
	const char *imgCovStr = string(container->basefilename + "_imgcov.dat").c_str();
	const char *noiseCovStr = string(container->basefilename + "_noisecov.dat").c_str();

	ifstream imgCovFile;
	imgCovFile.open(imgCovStr);

	ifstream noiseCovFile;
	noiseCovFile.open(noiseCovStr);

	if (imgCovFile.fail() || imgCovFile.fail()){
		fprintf(stderr, "Could not find image and noise covariances. Ensure that %s and %s exist. Exiting.\n", imgCovStr, noiseCovStr);
		exit(1);
	}

	//copy to container
	for (int i=0; i < container->img->getBands(); i++){
		string imgLine, noiseLine;
		getline(imgCovFile, imgLine);
		stringstream sI(imgLine);
		stringstream sN(noiseLine);
		for (int j=0; j < container->img->getBands(); j++){
			float val;
			sI >> val;
			container->imgCov[i*container->img->getBands() + j] = val;
			sN >> val;
			container->noiseCov[i*container->img->getBands() + j] = val;
		}
	}
	noiseCovFile.close();
	imgCovFile.close();


	//read band means from file
	const char *meanStr = string(container->basefilename + "_bandmeans.dat").c_str();
	ifstream meanFile;
	meanFile.open(meanStr);

	if (meanFile.fail()){
		fprintf(stderr, "Could not find file containing band means, %s. Exiting\n", meanStr);
		exit(1);
	}

	//copy to container
	for (int i=0; i < container->img->getBands(); i++){
		meanFile >> container->means[i];
	}

	meanFile.close();
}

void findCov(MnfWorkspace *workspace, Hyperspectral img, float *covMat, WhatValue whichVal){
	float *data = img.getAllData();
	int bands = img.getBands();

	double *means = new double[bands];
	int numMeanSamples = 0; //used if mean is incrementally calculated for each line (i.e. for the noise)

	for (int i=0; i < img.getBands(); i++){
		if (whichVal == USUAL_VAL){
			means[i] = findBandMean(img, i);
		} else {
			means[i] = 0.0;
		}
	}


	//initialize to zero
	for (int i=0; i < bands; i++){
		for (int j=0; j < bands; j++){
			covMat[i*img.getBands() + j] = 0;
		}
	}

	Covariance cov;
	initializeCov(&cov, bands, img.getPixels());
	

	//run over lines
	int samples = img.getPixels();

	float *prevLine = NULL;
	if (whichVal == ESTIMATE_NOISE){
		prevLine = new float[samples*bands];
	}
	for (int i=0; i < img.getLines(); i++){
		float *line = data + i*img.getPixels()*bands;
		switch(whichVal){
			case USUAL_VAL:
				updateCovariances(line, samples, bands, &cov);
				//updateMeans(line, samples, bands, means, &numMeanSamples);
			break;
			case ESTIMATE_NOISE:
				#if 1
				//my usual way of calculating noise, shift difference neighboring columns
				float *noise = estimateNoise(line, samples, bands);
				updateCovariances(noise, samples-1, bands, &cov);
				updateMeans(workspace, noise, samples-1, bands, means, &numMeanSamples);
				delete [] noise;
				#else
				//envi's way of doing it
				//skip first line
				if (i >= 1){
					float *noise = estimateENVINoise(line, prevLine, samples, bands);
					updateCovariances(noise, samples-1, bands, &cov);
					updateMeans(workspace, noise, samples-1, bands, means, &numMeanSamples);
					delete [] noise;
				} 
				memcpy(prevLine, line, sizeof(float)*samples*bands);
				#endif
			break;
		}
	}

	delete [] prevLine;
	
	float *meansTemp = new float[bands];
	for (int i=0; i < bands; i++){
		meansTemp[i] = means[i];
	}

	calculateTotCovariance(&cov, meansTemp, bands, covMat);

	deinitializeCov(&cov);
	delete [] means;
	delete [] meansTemp;
}

Hyperspectral* estimateNoise(Hyperspectral img, int startSamp, int endSamp, int startLine, int endLine){
	int newSamples = endSamp - startSamp - 1;
	int newLines = endLine - startLine;
	float *noise = new float[img.getBands()*newLines*newSamples];
	for (int i=0; i < img.getBands(); i++){
		for (int j=startLine; j < endLine; j++){
			for (int k=startSamp; k < endSamp-1; k++){
				//noise[j*img.getBands()*newSamples + i*newSamples + k] = img.read(j, k+1, i) - 0.5*(img.read(j, k+1, i) + img.read(j, k, i)); //the envi way
				noise[(j-startLine)*img.getBands()*newSamples + i*newSamples + (k-startSamp)] = img.read(j, k, i) - img.read(j, k+1, i);
			}
		}
	}
	Hyperspectral *ret = new Hyperspectral(noise, newLines, newSamples, img.getBands(), BIL);
	return ret;
}

void printCov(float *cov, int bands, string filename){
	ofstream file;
	file.open(filename.c_str());
	for (int i=0; i < bands; i++){
		for (int j=0; j < bands; j++){
			file << cov[i*bands + j] << " ";
		}
		file << endl;
	}
	file.close();
	
}

void printMeans(float *means, int bands, string filename){
	ofstream file;
	file.open(filename.c_str());
	for (int i=0; i < bands; i++){
		file << means[i] << " ";
	}
	file << endl;
	file.close();
	
}

