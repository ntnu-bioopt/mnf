
/* Author: Asgeir Bjørgan
asgeir.bjorgan@iet.ntnu.no
NTNU */

#ifndef READIMAGE_H_DEFINED
#define READIMAGE_H_DEFINED
#include <vector>

typedef struct {
	int samples;
	int bands;
	int lines;
	int offset;
	std::vector<float> wlens;
	int datatype;
} HyspexHeader;

typedef struct {
	int startSamp;
	int endSamp;
	int startLine;
	int endLine;
} ImageSubset;
	
void readHeader(char *filename, HyspexHeader *header);
void readImage(char *filename, HyspexHeader *header, ImageSubset subset, float *data);


#endif
