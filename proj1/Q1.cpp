// CS 474 Project 1, Question 1
// Team #8: Alexander Mathew and Tal Shahar Zemach
#include <iostream>
#include <fstream>
#include <sstream>
#include "image.h"
using namespace std;


// Function declarations
void subSampleImage(ImageType originalImage, ImageType& newImage, int factor);
int readImage(char fname[], ImageType& image);
int writeImage(char fname[], ImageType& image);

int main(int argc, char** argv) {
	

}

// function implementations
void subSampleImage(ImageType originalImage,ImageType& newImage, int factor){

	int orig_rows, orig_cols, orig_qlvls;
	int sampled_rows, sampled_cols;
	int val;

	// get image info
	originalImage.getImageInfo(orig_rows, orig_cols, orig_qlvls);

	// iterate over image
	for(int i = 0; i < orig_cols; i += factor){
		for(int j = 0; j < orig_rows; j+= factor){
			originalImage.getPixelVal(i, j, val);
			newImage.setPixelVal(sampled_rows, sampled_cols, val);
		sampled_rows++;
		}
	sampled_cols++;
	}
}
