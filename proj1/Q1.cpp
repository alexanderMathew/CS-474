// CS 474 Project 1, Question 1
// Team #8: Alexander Mathew and Tal Shahar Zemach
#include <iostream>
#include <fstream>
#include <sstream>
#include "image.h"
using namespace std;


// Function declarations
void subSampleImage(ImageType originalImage, int factor);
int readImage(char fname[], ImageType& image);
int writeImage(char fname[], ImageType& image);
int readImageHeader(char fname[], int& N, int& M, int& Q, bool& type);

int main(int argc, char** argv) {


}

// function implementations
void subSampleImage(ImageType originalImage, ImageType newImage, int factor){
	int orig_rows, orig_cols, orig_qlvl;
	int sampled_rows = 0;
	int sampled_cols = 0;
	int val;

	// get image info
	originalImage.getImageInfo(orig_rows, orig_cols, orig_qlvl);
	
	// iterate over image
	for(int i = 0; i < orig_cols; i += factor){
		for(int j = 0; j < orig_rows; j += factor){
			originalImage.getPixelVal(i, j, val);
			newImage.setPixelVal(sampled_rows, sampled_cols, val);
			
		sampled_rows++;

		}
	sampled_cols++;
	}

}
