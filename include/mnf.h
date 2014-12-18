#ifndef MNF_C_H_DEFINED
#define MNF_C_H_DEFINED

/* Author: Asgeir Bjørgan
asgeir.bjorgan@iet.ntnu.no
NTNU */

#include <string>
#include <vector>

//enum specifying which direction to run the MNF transform
enum TransformDirection{RUN_BOTH, RUN_FORWARD, RUN_INVERSE};

//struct containing variables and arrays neccessary for performing the MNF transform
typedef struct{
	std::string basefilename; //base filename for output
	float *R; //matrix with which to pre-multiply transformed image before inverse transform (first m elements bla bla)
	float *ones_samples; //vector of length samples consisting only of ones, for convenience
	
	TransformDirection direction;
} MnfWorkspace;

#include "mnf_linebyline.h"

//(de)initialize the mnf workspace
void mnf_initialize(TransformDirection dir, int bands, int samples, int numBandsInInverse, MnfWorkspace *workspace, std::string basefilename);
void mnf_deinitialize(MnfWorkspace *workSpace);

//run mnf, save results to file
void mnf_run(MnfWorkspace *workspace, int bands, int samples, int lines, float *data, std::vector<float> wlens);

//estimate mnf statistics from the image
void mnf_estimate_statistics(MnfWorkspace *workspace, int bands, int samples, int lines, float *img, ImageStatistics *imgStats, ImageStatistics *noiseStats);

//run forward mnf transform in place
void mnf_run_forward(MnfWorkspace *workspace, ImageStatistics *imgStats, ImageStatistics *noiseStats, int bands, int samples, int lines, float *img);

//run inverse mnf transform in place
void mnf_run_inverse(MnfWorkspace *workspace, ImageStatistics *imgStats, ImageStatistics *noiseStats, int bands, int samples, int lines, float *img);

//forward and inverse transformations based on covariances
void mnf_get_transf_matrix(int bands, ImageStatistics *imgStats, ImageStatistics *noiseStats, float *forwardTransf, float *inverseTransf, float *eigvals);
void mnf_calculate_forward_transf_matrix(int bands, const float *imgCov, const float *noiseCov, float *forwardTransf, float *eigvals);
void mnf_calculate_inverse_transf_matrix(int bands, const float *forwardTransf, float *inverseTransf);



#endif
