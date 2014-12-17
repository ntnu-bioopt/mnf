/* Author: Asgeir Bjørgan
asgeir.bjorgan@iet.ntnu.no
NTNU */


#include <fstream>
#include <iostream>
#include <math.h>
#include <cstring>
#include <sstream>
using namespace std;
#include "hyperspectral.h"
void Hyperspectral::writeToFile(string outfilename, vector<float> wlens){
	//prepare image file
	ostringstream imgFname;
	imgFname << outfilename << ".img";
	ofstream *hyspexOut = new ofstream(imgFname.str().c_str(),ios::out | ios::binary);
	
	//write image header
	ostringstream hdrFname;
	hdrFname << outfilename << ".hdr";
	ofstream hdrOut(hdrFname.str().c_str());
	hdrOut << "ENVI" << endl;
	hdrOut << "samples = " << numPixels << endl;
	hdrOut << "lines = " << numLines << endl;
	hdrOut << "bands = " << numBands << endl;
	hdrOut << "header offset = 0" << endl;
	hdrOut << "file type = ENVI Standard" << endl;
	hdrOut << "data type = 4" << endl;
	hdrOut << "interleave = bil" << endl;
	hdrOut << "default bands = {55,41,12}" << endl;
	hdrOut << "byte order = 0" << endl;
	hdrOut << "wavelength = {";
	for (int i=0; i < wlens.size(); i++){
		hdrOut << wlens[i] << " ";
	}
	hdrOut << "}" << endl;
	hdrOut.close();

	//write image
	for (int i=0; i < numLines; i++){
		float *data = getData(i);
		hyspexOut->write((char*)(data), sizeof(float)*numBands*numPixels);
		delete [] data;
	}

	hyspexOut->close();
	delete hyspexOut;
}

void Hyperspectral::rescale(int startLine, int endLine, int startPixel, int endPixel){
	int newLines = endLine - startLine;
	int newSamples = endPixel - startPixel;

	float *newData = new float[newLines*newSamples*numBands];
	for (int i=startLine; i < endLine; i++){
		for (int j=startPixel; j < endPixel; j++){
			for (int k=0; k < numBands; k++){
				newData[(i-startLine)*newSamples*numBands + k*newSamples + (j-startPixel)] = data[i*numPixels*numBands + k*numPixels + j];
			}
		}
	}
	delete [] data;
	data = newData;
	numLines = newLines;
	numPixels = newSamples;
}


void Hyperspectral::writeToFile(string file){
	vector<float> wlens;
	for (int i=0; i < numBands; i++){
		wlens.push_back(i);
	}
	writeToFile(file, wlens);
}
		
float Hyperspectral::read(int line, int pixel, int band){
	if ((line >= numLines) || (pixel >= numPixels) || (band >= numBands)){
		return 0.0f;
	}
	
	float val;
	switch (interleave){
		case BIL:
			val = data[numPixels*numBands*line + band*numPixels + pixel];
		break;
		case BSQ:
			val = data[band*numPixels*numLines + line*numPixels + pixel];
		break;
	}
	if (isnan(val)){
		val = 0.0f;
	}
	return val;
}
		
void Hyperspectral::write(int line, int pixel, int band, float val){
	switch (interleave){
		case BIL:
			data[numPixels*numBands*line + band*numPixels + pixel] = val;
		break;
		case BSQ:
			data[band*numPixels*numLines + line*numPixels + pixel] = val;
		break;

	}
}
		
float *Hyperspectral::getData(int line){
	float *retarr = new float[numBands*numPixels];
	switch (interleave){
		case BIL:
			memcpy(retarr, data + line*numPixels*numBands, sizeof(float)*numBands*numPixels);
		break;
		case BSQ:
			for (int i=0; i < numBands; i++){
				for (int j=0; j < numPixels; j++){
					retarr[i*numPixels + j] = read(line, j, i);
				}
			}
		break;
	}
	return retarr;
}
