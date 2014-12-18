
/* Author: Asgeir Bjørgan
asgeir.bjorgan@iet.ntnu.no
NTNU */

#include <string>
#include <iostream>
#include <fstream>
#include "mnf.h"
#include <string.h>
#include <sys/time.h>
#include <stdio.h>
#include <pthread.h>
#include <sstream>
#include "readimage.h"
extern "C"
{
#include <lapacke.h>
#include <cblas.h>
}
#include "mnf_linebyline.h"
using namespace std;


void mnf_run(MnfWorkspace *workspace, int bands, int samples, int lines, float *data, std::vector<float> wlens){
	ImageStatistics imgStats;
	ImageStatistics noiseStats;
	imagestatistics_initialize(&imgStats, bands);
	imagestatistics_initialize(&noiseStats, bands);

	if (workspace->direction == RUN_FORWARD || workspace->direction == RUN_BOTH){
		//estimate image statistics
		mnf_estimate_statistics(workspace, bands, samples, lines, data, &imgStats, &noiseStats);

		//run forward transform
		cout << "Run forward transform." << endl;
		mnf_run_forward(workspace, &imgStats, &noiseStats, bands, samples, lines, data);

		//write transformed data to file
		cout << "Write transformed data to file." << endl;
		hyperspectral_write_header(string(workspace->basefilename + "_transformed").c_str(), bands, samples, lines, wlens);
		hyperspectral_write_image(string(workspace->basefilename + "_transformed").c_str(), bands, samples, lines, data);

		//write image statistics to file
		imagestatistics_write_to_file(workspace, bands, &imgStats, &noiseStats);
	} else if (workspace->direction == RUN_INVERSE){
		//get statistics from file
		imagestatistics_read_from_file(workspace, bands, &imgStats, &noiseStats);
	} else if (workspace->direction == RUN_BOTH){
		//run inverse transform
		cout << "Run inverse transform." << endl;
		mnf_run_inverse(workspace, &imgStats, &noiseStats, bands, samples, lines, data);
		cout << "Write transformed data to file." << endl;
		hyperspectral_write_header(string(workspace->basefilename + "_inversetransformed").c_str(), bands, samples, lines, wlens);
		hyperspectral_write_image(string(workspace->basefilename + "_inversetransformed").c_str(), bands, samples, lines, data);
	}

	imagestatistics_deinitialize(&imgStats);
	imagestatistics_deinitialize(&noiseStats);
}	

void mnf_initialize(TransformDirection direction, int bands, int samples, int numBandsInInverse, MnfWorkspace *workspace, std::string basefilename){
	workspace->direction = direction;

	workspace->ones_samples = new float[samples];
	for (int i=0; i < samples; i++){
		workspace->ones_samples[i] = 1.0f;
	}
	
	workspace->R = new float[bands*bands]();
	for (int i=0; i < numBandsInInverse; i++){
		workspace->R[i*bands + i] = 1.0f;
	}

	workspace->basefilename = basefilename;
}

void mnf_deinitialize(MnfWorkspace *workspace){
	delete [] workspace->ones_samples;
	delete [] workspace->R;
}

void mnf_estimate_statistics(MnfWorkspace *workspace, int bands, int samples, int lines, float *img, ImageStatistics *imgStats, ImageStatistics *noiseStats){
	for (int i=0; i < lines; i++){
		//image statistics
		float *line = img + i*samples*bands;
		imagestatistics_update_with_line(workspace, bands, samples, line, imgStats);

		//noise statistics
		float *noise;
		int noise_samples = 0;
		mnf_linebyline_estimate_noise(bands, samples, line, &noise, &noise_samples);
		imagestatistics_update_with_line(workspace, bands, noise_samples, noise, noiseStats);
		delete [] noise;
	}
}


void mnf_run_forward(MnfWorkspace *workspace, ImageStatistics *imgStats, ImageStatistics *noiseStats, int bands, int samples, int lines, float *img){
	//get means for removing
	float *means = new float[bands];
	imagestatistics_get_means(imgStats, bands, means);

	//get transfer matrices
	float *forwardTransf = new float[bands*bands];
	float *inverseTransf = new float[bands*bands];
	float *eigvals = new float[bands*bands];
	mnf_get_transf_matrix(bands, imgStats, noiseStats, forwardTransf, inverseTransf, eigvals);

	//write eigenvalues to file
	ofstream eigvalsFile;
	eigvalsFile.open(string(workspace->basefilename + "_eigvals.dat").c_str());
	for (int i=0; i < bands; i++){
		eigvalsFile << eigvals[i] << endl;
	}
	eigvalsFile.close();
	
	//run transform
	for (int i=0; i < lines; i++){
		float *submatr = img + i*samples*bands;
		
		//remove means
		mnf_linebyline_remove_mean(workspace, means, bands, samples, submatr);

		//perform transform of single line
		float *submatrNew = new float[samples*bands];
		cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, bands, samples, bands, 1.0f, forwardTransf, bands, submatr, samples, 0.0f, submatrNew, samples);
		memcpy(img + i*samples*bands, submatrNew, sizeof(float)*samples*bands);
		delete [] submatrNew;
	}

	delete [] means;
	delete [] forwardTransf;
	delete [] inverseTransf;
	delete [] eigvals;
}

void mnf_run_inverse(MnfWorkspace *workspace, ImageStatistics *imgStats, ImageStatistics *noiseStats, int bands, int samples, int lines, float *img){
	//get means for adding
	float *means = new float[bands];
	imagestatistics_get_means(imgStats, bands, means);

	//get transfer matrices
	float *forwardTransf = new float[bands*bands];
	float *inverseTransf = new float[bands*bands];
	float *eigvals = new float[bands*bands];
	mnf_get_transf_matrix(bands, imgStats, noiseStats, forwardTransf, inverseTransf, eigvals);
	
	//post-multiply R matrix to get k < m first bands in inverse transformation
	float *R = workspace->R;
	float *inverseTransfArrShort = new float[bands*bands];
	cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, bands, bands, bands, 1.0f, inverseTransf, bands, R, bands, 0.0f, inverseTransfArrShort, bands);
	
	//run transform
	for (int i=0; i < lines; i++){
		float *submatr = img + i*samples*bands;

		//perform inverse transform of single line
		float *submatrNew = new float[samples*bands];
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, bands, samples, bands, 1.0f, inverseTransfArrShort, bands, submatr, samples, 0.0f, submatrNew, samples);
		
		//add means
		mnf_linebyline_add_mean(workspace, means, bands, samples, submatrNew);
		
		memcpy(img + i*samples*bands, submatrNew, sizeof(float)*samples*bands);
		delete [] submatrNew;
		
	}

	delete [] means;
	delete [] forwardTransf;
	delete [] inverseTransf;
	delete [] inverseTransfArrShort;
}


void mnf_get_transf_matrix(int bands, ImageStatistics *imgStats, ImageStatistics *noiseStats, float *forwardTransf, float *inverseTransf, float *eigvals){
	//get true covariances from supplied statistics
	float *imgCov = new float[bands*bands];
	imagestatistics_get_cov(imgStats, bands, imgCov);

	float *noiseCov = new float[bands*bands];
	imagestatistics_get_cov(noiseStats, bands, noiseCov);
	
	//estimate forward transformation matrix
	mnf_calculate_forward_transf_matrix(bands, imgCov, noiseCov, forwardTransf, eigvals);
	
	//estimate inverse transformation matrix
	mnf_calculate_inverse_transf_matrix(bands, forwardTransf, inverseTransf);
	
	delete [] imgCov;
	delete [] noiseCov;
}

void mnf_calculate_forward_transf_matrix(int bands, const float *imgCov, const float *noiseCov, float *forwardTransf, float *eigvals){
	//copy covariances to temp variables, since dsygv will destroy them 
	float *noiseCovTemp = new float[bands*bands];
	float *imgCovTemp = new float[bands*bands];
	memcpy(noiseCovTemp, noiseCov, sizeof(float)*bands*bands);
	memcpy(imgCovTemp, imgCov, sizeof(float)*bands*bands);

	//solve eigenvalue problem
	int itype = 1; //solves Ax = lam * Bx
	char jobz = 'V'; //compute both eigenvalues and eigenvectors
	char uplo = 'L'; //using lower triangles of matrices
	int err = LAPACKE_ssygv(LAPACK_ROW_MAJOR, itype, jobz, uplo, bands, noiseCovTemp, bands, imgCovTemp, bands, eigvals);
	if (err != 0){
		fprintf(stderr, "LAPACKE_dsygv failed: calculation of forward MNF transformation matrix, err code %d\n", err);
	}

	//move transfer matrix (stored in noiseCovTemp) to output variable
	memcpy(forwardTransf, noiseCovTemp, sizeof(float)*bands*bands);

	delete [] noiseCovTemp;
	delete [] imgCovTemp;
}

void mnf_calculate_inverse_transf_matrix(int bands, const float *forwardTransf, float *inverseTransf){
	size_t size = bands*bands*sizeof(float);
	
	//keep forward transf matrix
	float *forwardTransfTemp = (float*)malloc(size);
	memcpy(forwardTransfTemp, forwardTransf, size);

	//find inverse of forward eigenvector transfer matrix
	int *ipiv = new int[bands];
	int err = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, bands, bands, forwardTransfTemp, bands, ipiv);
	if (err != 0){
		fprintf(stderr, "LU decomposition failed: calculation of inverse transformation matrix\n");
	}
	err = LAPACKE_sgetri(LAPACK_ROW_MAJOR, bands, forwardTransfTemp, bands, ipiv);
	if (err != 0){
		fprintf(stderr, "Inversion failed: calculation of inverse transformation matrix\n");
	}
	delete [] ipiv;

	//copy back
	memcpy(inverseTransf, forwardTransfTemp, size);
	
	free(forwardTransfTemp);
}


