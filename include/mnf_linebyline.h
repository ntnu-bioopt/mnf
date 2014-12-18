#ifndef MNF_LINEBYLINE_H_DEFINED
#define MNF_LINEBYLINE_H_DEFINED

#include "mnf.h"

void mnf_linebyline_run_image(MnfWorkspace *workspace, int bands, int numSamples, int numLines, float *data, std::string mnfOutFilename);


typedef struct{
	int n; //number of pixels so far summed over
	float *C; 
	double *means; 

} ImageStatistics;

//(de)initialize covariance arrays
void imagestatistics_initialize(ImageStatistics *stats, int bands);
void imagestatistics_deinitialize(ImageStatistics *stats);

//run mnf on one line, incrementally updating the covariances by this line and denoising line in place
void mnf_linebyline_run_oneline(MnfWorkspace *workspace, int numBands, int numSamples, float *line, ImageStatistics *imageStats, ImageStatistics *noiseStats);

//update covariances and means using a numerically stable algorithm
void imagestatistics_update_with_line(const MnfWorkspace *workspace, int numBands, int numSamples, float *bilData, ImageStatistics *stats);

//explicitly estimate noise by shift difference
//allocates the noise estimate array
void mnf_linebyline_estimate_noise(int bands, int samples, float *line, float **noise_est, int *noise_samples);

void imagestatistics_get_means(ImageStatistics *stats, int numBands, float *means);
void imagestatistics_get_cov(ImageStatistics *stats, int numBands, float *cov);



void mnf_linebyline_remove_mean(const MnfWorkspace *workspace, float *means, int bands, int samples, float *line);
void mnf_linebyline_add_mean(const MnfWorkspace *workspace, float *means, int bands, int samples, float *line);



#endif
