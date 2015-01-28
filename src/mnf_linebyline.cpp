//=======================================================================================================
// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//=======================================================================================================

#include <stdio.h>
#include <string.h> 
#include <sys/time.h>
#include "mnf.h"
#include <iostream>
#include <sys/time.h>
#include <cblas.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "readimage.h"
#include "mnf_linebyline.h"
using namespace std;

///////////////////////////////
// MNF-LBL RELEVANT ROUTINES //
///////////////////////////////

void mnf_linebyline_estimate_noise(int bands, int samples, float *line, float **noise_est, int *noise_samples){
	*noise_samples = samples-1;
	*noise_est = new float[*noise_samples*bands];
	for (int i=0; i < bands; i++){
		for (int j=0; j < *noise_samples; j++){
			(*noise_est)[i*(*noise_samples) + j] = line[i*samples + j] - line[i*samples + (j + 1)];
		}
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
	imagestatistics_update_with_line(workspace, bands, samples, line, imageStats);
	
	//noise statistics
	mnf_linebyline_estimate_noise(bands, samples, line, &noise, &noise_samples);
	imagestatistics_update_with_line(workspace, bands, noise_samples, noise, noiseStats);

	//get means	
	float *means_float = new float[bands];
	imagestatistics_get_means(imageStats, bands, means_float);

	//find forward and inverse transformation matrices
	float *forwardTransfArr = new float[bands*bands];
	float *inverseTransfArr = new float[bands*bands];
	float *eigvals = new float[bands];
	mnf_get_transf_matrix(bands, imageStats, noiseStats, forwardTransfArr, inverseTransfArr, eigvals);
	delete [] eigvals;
	
	//remove mean from data
	mnf_linebyline_remove_mean(workspace, means_float, bands, samples, line);

	//pre-calculate matrix with which to multiply the image for noise removal
	float *invTrMultR = new float[bands*bands];
	cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, bands, bands, bands, 1.0f, inverseTransfArr, bands, workspace->R, bands, 0.0f, invTrMultR, bands);

	float *invTrMultRMultForTr = new float[bands*bands];
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, bands, bands, bands, 1.0f, invTrMultR, bands, forwardTransfArr, bands, 0.0f, invTrMultRMultForTr, bands);

	//actual noise removal	
	float *submatrNew = new float[bands*samples];

	switch(workspace->direction){
		case RUN_FORWARD:
			cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, bands, samples, bands, 1.0f, forwardTransfArr, bands, line, samples, 0.0f, submatrNew, samples);
			break;
		default:
			//inverse only doesn't make sense for the line by line method, so do both directions
			cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, bands, samples, bands, 1.0f, invTrMultRMultForTr, bands, line, samples, 0.0f, submatrNew, samples);
			
			//add mean to data
			mnf_linebyline_add_mean(workspace, means_float, bands, samples, submatrNew);
			break;
	}		
	

	memcpy(line, submatrNew, sizeof(float)*samples*bands);
	delete [] submatrNew;

	delete [] invTrMultR;
	delete [] invTrMultRMultForTr;
	delete [] noise;
	delete [] means_float;
	delete [] forwardTransfArr;
	delete [] inverseTransfArr;
}

void mnf_linebyline_run_image(MnfWorkspace *workspace, int bands, int samples, int lines, float *data, vector<float> wlens){
	ImageStatistics noiseStats;
	ImageStatistics imageStats;

	imagestatistics_initialize(&noiseStats, bands);
	imagestatistics_initialize(&imageStats, bands);


	for (int i=0; i < lines; i++){
		float *line = data + i*samples*bands;
		fprintf(stderr, "Line %d/%d\n", i, lines);
		
		timeval time1;
		gettimeofday(&time1, NULL);

		mnf_linebyline_run_oneline(workspace, bands, samples, line, &imageStats, &noiseStats);
		
		timeval time2;
		gettimeofday(&time2, NULL);
		cout << (time2.tv_sec + time2.tv_usec*1.0e-06) - (time1.tv_sec + time1.tv_usec*1.0e-06) << endl;
	}

	imagestatistics_deinitialize(&noiseStats);
	imagestatistics_deinitialize(&imageStats);

	//write image to file
	hyperspectral_write_header(string(workspace->basefilename + "_linebyline").c_str(), bands, samples, lines, wlens);
	hyperspectral_write_image(string(workspace->basefilename + "_linebyline").c_str(), bands, samples, lines, data);
}


///////////////////////////
// STATISTICS ESTIMATION //
///////////////////////////


void imagestatistics_get_means(ImageStatistics *stats, int numBands, float *means){
	for (int i=0; i < numBands; i++){
		means[i] = stats->means[i];
	}
}

void imagestatistics_get_cov(ImageStatistics *stats, int numBands, float *cov){
	for (int i=0; i < numBands*numBands; i++){
		cov[i] = stats->C[i]/(stats->n*1.0f);
	}
}

void imagestatistics_read_from_file(MnfWorkspace *workspace, int bands, ImageStatistics *imgStats, ImageStatistics *noiseStats){
	const char *imgCovStr = string(workspace->basefilename + "_imgcov.dat").c_str();
	const char *noiseCovStr = string(workspace->basefilename + "_noisecov.dat").c_str();

	ifstream imgCovFile;
	imgCovFile.open(imgCovStr);

	ifstream noiseCovFile;
	noiseCovFile.open(noiseCovStr);

	if (imgCovFile.fail() || imgCovFile.fail()){
		fprintf(stderr, "Could not find image and noise covariances. Ensure that %s and %s exist. Exiting.\n", imgCovStr, noiseCovStr);
		exit(1);
	}

	//copy to container
	for (int i=0; i < bands; i++){
		string imgLine, noiseLine;
		getline(imgCovFile, imgLine);
		stringstream sI(imgLine);
		stringstream sN(noiseLine);
		for (int j=0; j < bands; j++){
			float val;
			sI >> val;
			imgStats->C[i*bands + j] = val;
			sN >> val;
			noiseStats->C[i*bands + j] = val;
		}
	}

	imgStats->n = 1;
	noiseStats->n = 1;

	noiseCovFile.close();
	imgCovFile.close();


	//read band means from file
	const char *meanStr = string(workspace->basefilename + "_bandmeans.dat").c_str();
	ifstream meanFile;
	meanFile.open(meanStr);

	if (meanFile.fail()){
		fprintf(stderr, "Could not find file containing band means, %s. Exiting\n", meanStr);
		exit(1);
	}

	//copy to container
	for (int i=0; i < bands; i++){
		float val = 0;
		meanFile >> val;
		imgStats->means[i] = val;
	}

	meanFile.close();
}

void imagestatistics_write_cov_to_file(ImageStatistics *imgStats, int bands, std::string filename){
	ofstream file;
	file.open(filename.c_str());
	float *cov = new float[bands*bands];
	imagestatistics_get_cov(imgStats, bands, cov);
	for (int i=0; i < bands; i++){
		for (int j=0; j < bands; j++){
			file << cov[i*bands + j] << " ";
		}
		file << endl;
	}
	file.close();
	delete [] cov;
}

void imagestatistics_write_mean_to_file(ImageStatistics *imgStats, int bands, std::string filename){
	float *means = new float[bands];
	imagestatistics_get_means(imgStats, bands, means);
	ofstream file;
	file.open(filename.c_str());
	for (int i=0; i < bands; i++){
		file << means[i] << " ";
	}
	file << endl;
	file.close();
	delete [] means;
}

void imagestatistics_write_to_file(MnfWorkspace *workspace, int bands, ImageStatistics *imgStats, ImageStatistics *noiseStats){
	imagestatistics_write_cov_to_file(imgStats, bands, workspace->basefilename + "_imgcov.dat");
	imagestatistics_write_cov_to_file(noiseStats, bands, workspace->basefilename + "_noisecov.dat");
	imagestatistics_write_mean_to_file(imgStats, bands, workspace->basefilename + "_bandmeans.dat");
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


void imagestatistics_update_with_line(const MnfWorkspace *workspace, int bands, int samples, float *bilData, ImageStatistics *stats){
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
	mnf_linebyline_remove_mean(workspace, meanTemp, bands, samples, tempLine);

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
