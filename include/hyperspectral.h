
/* Author: Asgeir Bjørgan
asgeir.bjorgan@iet.ntnu.no
NTNU */

#ifndef HYPERSPECTRAL_H_DEFINED
#define HYPERSPECTRAL_H_DEFINED

#include <string>
#include <vector>

enum Interleave{BSQ, BIL, BIP};

class Hyperspectral{
	public:
		Hyperspectral(float *data, int lines, int pixels, int bands, Interleave interleave) : numLines(lines), numPixels(pixels), numBands(bands), data(data), interleave(interleave){}; 
		float read(int line, int pixel, int band);
		void write(int line, int pixel, int band, float val); //will write as BIL-interleaved
		int getBands(){return numBands;};
		int getPixels(){return numPixels;};
		int getLines(){return numLines;};
		float *getData(int line);
		float *getAllData(){return data;};

		void writeToFile(std::string file);
		void writeToFile(std::string file, std::vector<float> wlens);

		void rescale(int startLine, int endLine, int startPixel, int endPixel);
	private:
		int numPixels;
		int numLines;
		int numBands;
		float *data;
		Interleave interleave;

	
};


#endif
