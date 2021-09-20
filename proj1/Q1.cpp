// CS 474 Project 1, Question 1
// Team #8: Alexander Mathew and Tal Shahar Zemach
#include <iostream>
#include <fstream>
#include <sstream>
#include "image.h"
using namespace std;


// Function declarations
//void subSampleImage(ImageType, ImageType&, int);
int readImage(char fname[], ImageType& image);
int writeImage(char fname[], ImageType& image);
int readImageHeader(char fname[], int& N, int& M, int& Q, bool& type);


// function implementations
void subSampleImage(ImageType originalImage, ImageType& newImage, int factor){
	int orig_rows, orig_cols, orig_qlvl;
	int sampled_rows = 0;
	int sampled_cols = 0;
	int val;
	int pixelVal;

	// get image info
	originalImage.getImageInfo(orig_rows, orig_cols, orig_qlvl);
	
	// iterate over image
	for(int i = 0; i < orig_cols; i += factor){
		sampled_rows = 0;
		for(int j = 0; j < orig_rows; j += factor){
			originalImage.getPixelVal(i, j, val); // get pixel values from the original image
			newImage.setPixelVal(sampled_cols, sampled_rows, val); // set pixel values
			newImage.assignPixelVal(sampled_cols, sampled_rows, pixelVal); // assign pixel values
		
			for(int x = 0; x < factor; x++){
				for(int y = 0; y < factor; y++){
					newImage.resizePixel(i+x, j+y, pixelVal); //newImage[i+x][j+y] = newImage[][];
				}				
			}
		sampled_rows++;

		}
	sampled_cols++;
	}

}

int main() {
	int m, n, q;
	bool type;
	int ss_factor = 2;

	// read image
	char* readPath = "./Images/lenna.pgm";
	char *writePath = "./Images/lenna_sampled.pgm"; 	

	readImageHeader(readPath, n, m, q, type);
	ImageType ImageToSample(n, m, q);
	ImageToSample.getImageInfo(n, m, q);

	int x = m;
	int y = q;
	int z = n;

	ImageType resizedImage(z, x, y);

	readImage(readPath, ImageToSample);


	// sample image
	subSampleImage(ImageToSample, resizedImage, ss_factor);


	// write image
	writeImage(writePath, resizedImage);

	return 0;
}

