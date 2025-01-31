//CS 474 Project 2
//Team #8: Alexander Mathew & Tal Shahar Zemach

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "image.h"

using namespace std;

// Function declarations
int readImage(char fname[], ImageType &image);
int writeImage(char fname[], ImageType &image);
int readImageHeader(char fname[], int &N, int &M, int &Q, bool &type);


//this function rescales values to [0-255]
void rescale(int n, int m, ImageType OriginalImage, char* writePath)
{

	int max = 0;
	int min = 0;

	int temp;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			OriginalImage.getPixelVal(i, j, temp);
			if (temp > max)
			{
				max = temp;
				OriginalImage.setPixelVal(i, j, max); //max=OriginalImage[i][j];
			}
			if (temp < min)
			{
				min = temp;
				OriginalImage.setPixelVal(i, j, min); //min=OriginalImage[i][j];
			}
		}
	}

	// map to [0,255] equation
	int x;
	int y;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
			OriginalImage.getPixelVal(i, j, x);
			y = 255.0 * (x - min) / (max - min); // The Equation: 255.0 is a DOUBLE (prevent overflow)
			OriginalImage.setPixelVal(i, j, y);
		}
	}

	writeImage(writePath, OriginalImage);
}


// Function to perform correlation on a given image
void correlation(ImageType& image, ImageType& mask, int mask_x, int mask_y, char * writePath){

	//Image originalImage = Image(image);
	int rows, cols, qlvl;
	int m_rows, m_cols, m_qlvl;

	image.getImageInfo(rows, cols, qlvl);
	mask.getImageInfo(m_rows, m_cols, m_qlvl);

	// Vector to store correlation values
	vector<int> corr_values;
	int temp1;
	int temp2;

	ImageType resizedImage(rows, cols, qlvl);

	// iterate through image pixels
    for(int i = 0; i < rows; i++){
   		for(int j = 0; j < cols; j++) {

   			int sum = 0;

			
   			// iterate over mask
   			for (int k = -mask_y/2; k < mask_y/2; k++) //(int k = -mask_y/2; k < mask_y/2; k++)
   			{
   				for (int l = -mask_x/2; l < mask_x/2; l++) //(int l = -mask_x/2; l < mask_x/2; l++)
   				{

   					//check bounds
   					if(i + k < 0 || i + k >= rows || j + l < 0 || j + l >= cols){
   						sum += 0;
					}
   					else{
						int num1 = i+k;
						int num2 = j+l;
						image.getPixelVal(num1,num2, temp1); //temp1 = image[i+k][j+l];
						int num3 = k+mask_y/2;
						int num4 = l+mask_x/2;
						mask.getPixelVal(num3, num4, temp2);
						sum = sum + (temp1 * temp2);
						//sum += originalImage[i + k][j + l] * mask[k + size_y/2][l + size_x/2];
					}
				
   				}

   			}

   			corr_values.push_back(sum);
			resizedImage.setPixelVal(i, j, sum);
   		}
   	}

	rescale(rows, cols, resizedImage, writePath);
}


// a function to perform averaging on a given image
void averaging(ImageType& imageToAverage, int maskSize, char * writePath){

	int orig_rows, orig_cols, orig_qlvl;
	imageToAverage.getImageInfo(orig_rows, orig_cols, orig_qlvl);

	// Iterate through image pixels
  	for(int i = 0; i < orig_cols; i++){
   		for(int j = 0; j < orig_rows; j++) {

   			int average = 0;

   			//iterate through mask
   			for (int k = -maskSize/2; k < maskSize/2; k++){
   				for (int l = -maskSize/2; l < maskSize/2; l++){
   					// calcualte average
   					if(i + k < 0 || i + k >= orig_cols || j + l < 0 || j + l >= orig_rows){
   						average += 0;
					   
					}
   					else{
   						int temp;
						imageToAverage.getPixelVal((i+k), (j+l), temp);
						average = average + temp;
					}
   				}
   			}
		// Compute Average
		int temp2=average / (maskSize * maskSize);
		imageToAverage.setPixelVal(i, j, temp2);
   		}
   	}
	rescale(orig_rows, orig_cols, imageToAverage, writePath);
}


// a function to perform gaussian filtering on a given image
void gaussian(ImageType& imageToGaussian, int maskSize, char * writePath){
	// 7x7 Gaussian mask
	const int mask_7x7[7][7] = {

		{1, 1, 2, 2, 2, 1, 1},
		{1, 2, 2, 4, 2, 2, 1},
		{2, 2, 4, 8, 4, 2, 2},
		{2, 4, 8, 16, 8, 4, 2},
		{2, 2, 4, 8, 4, 2, 2},
		{1, 2, 2, 4, 2, 2, 1},
		{1, 1, 2, 2, 2, 1, 1}
	};

	//15x15 Gaussian mask
	const int mask_15x15[15][15] = {
		{2, 2,  3,  4,  5,  5,  6,  6,  6,  5,  5,  4,  3, 2, 2},
		{2, 3,  4,  5,  7,  7,  8,  8,  8,  7,  7,  5,  4, 3, 2},
		{3, 4,  6,  7,  9, 10, 10, 11, 10, 10,  9,  7,  6, 4, 3},
		{4, 5,  7,  9, 10, 12, 13, 13, 13, 12, 10,  9,  7, 5, 4},
		{5, 7,  9, 11, 13, 14, 15, 16, 15, 14, 13, 11,  9, 7, 5},
		{5, 7, 10, 12, 14, 16, 17, 18, 17, 16, 14, 12, 10, 7, 5},
		{6, 8, 10, 13, 15, 17, 19, 19, 19, 17, 15, 13, 10, 8, 6},
		{6, 8, 11, 13, 16, 18, 19, 20, 19, 18, 16, 13, 11, 8, 6},
		{6, 8, 10, 13, 15, 17, 19, 19, 19, 17, 15, 13, 10, 8, 6},
		{5, 7, 10, 12, 14, 16, 17, 18, 17, 16, 14, 12, 10, 7, 5},
		{5, 7,  9, 11, 13, 14, 15, 16, 15, 14, 13, 11,  9, 7, 5},
		{4, 5,  7,  9, 10, 12, 13, 13, 13, 12, 10,  9,  7, 5, 4},
		{3, 4,  6,  7,  9, 10, 10, 11, 10, 10,  9,  7,  6, 4, 3},
		{2, 3,  4,  5,  7,  7,  8,  8,  8,  7,  7,  5,  4, 3, 2},
		{2, 2,  3,  4,  5,  5,  6,  6,  6,  5,  5,  4,  3, 2, 2},
	};


	int orig_rows, orig_cols, orig_qlvl;
	imageToGaussian.getImageInfo(orig_rows, orig_cols, orig_qlvl);

	// calculate normalizion factor
	int factor = 0; 
	for(int i = 0; i < maskSize; i++){
		for (int j = 0; j < maskSize; j++){

			if(maskSize == 7){
				factor += mask_7x7[i][j];
			}
			else{
				factor += mask_15x15[i][j];
			}
		}
	}

	// iterate through image pixels
    for(int i = 0; i < orig_cols; i++){
   		for(int j = 0; j < orig_rows; j++) {

   			int pixelVal = 0;

   			// iterate through mask
   			for (int k = -maskSize/2; k < maskSize/2; k++){
   				for (int l = -maskSize/2; l < maskSize/2; l++){
   					
					// check bounds
   					if(i + k < 0 || i + k >= orig_cols || j + l < 0 || j + l >= orig_rows){
   						pixelVal += 0;
   					}
					
   					// calculate output for desired mask size
   					else if(maskSize == 7){
						int temp;
						imageToGaussian.getPixelVal((i+k), (j+l), temp);
						pixelVal = pixelVal + (temp * mask_7x7[k + maskSize/2][l + maskSize/2]);
					}

   					else if(maskSize == 15){
						int temp;
						imageToGaussian.getPixelVal((i+k), (j+l), temp);
						pixelVal = pixelVal + (temp * mask_15x15[k + maskSize/2][l + maskSize/2]);

					}
   				}
   			}

   			// update image pixel
			int temp2 = pixelVal/factor;
			imageToGaussian.setPixelVal(i, j, temp2);
   		}
   	}
	
	rescale(orig_rows, orig_cols, imageToGaussian, writePath);

}


// a function to perform median filtering on a given image
void medianFilter(ImageType& imageToFilter, int maskSize, char * writePath){

	int orig_rows, orig_cols, orig_qlvl;
	imageToFilter.getImageInfo(orig_rows, orig_cols, orig_qlvl);
	// resized image after calculating median
	ImageType resizedImage(orig_rows, orig_cols, orig_qlvl);

    // Iterate through pixel values
    for(int i = 0; i < orig_cols; i++){
   		for(int j = 0; j < orig_rows; j++) {

   		vector<int> pixelVal;

        // Iterate though mask
   			for (int k = i; k < i+maskSize; k++)
   			{
   				for (int l = j; l < j+maskSize; l++)
   				{
   					//check bounds
   					if(i + k < 0 || i + k >= orig_cols || j + l < 0 || j + l >= orig_rows){
   						pixelVal.push_back(0);
					}
   					else{
						int temp;
						imageToFilter.getPixelVal((i+k), (j+l), temp);
						pixelVal.push_back(temp);
					}	   
   					
   				}
   			}


        	// Sort pixelVal
   			sort(pixelVal.begin(), pixelVal.end());

			// Assign median value to the desired location i,j
			int temp2 = pixelVal[(maskSize * maskSize) / 2];
			imageToFilter.setPixelVal(i, j, temp2);
   		}
	
	// reused from proj 1 to resize image
	int factor = 2; // factor is 2 becuase the image was 1/4 of the expected size
	for (int x = 1; x < orig_cols; x++)
	{
		for (int y = 1; y < orig_rows; y++)
		{
			int temp;
			imageToFilter.getPixelVal(x / factor, y / factor, temp);
			//ImageType
			resizedImage.setPixelVal(x, y, temp);
		}
	}
}

	//rescale(orig_rows, orig_cols, imageToFilter, writePath);
	writeImage(writePath, resizedImage);


}


// a function to artifically add salt and pepper noise to an image
void saltAndPepper(ImageType& imageToCorrupt, int noise_pct, char * writePath){

	int rows, cols, qlvl;
	imageToCorrupt.getImageInfo(rows, cols, qlvl);
	vector<int> idx; // Location of corrupt pixels
	int corruptPixels = cols * rows * (noise_pct / 100.0); // # of corrupt pixels

    // save each image location
    for (int i = 0; i < cols * rows; i++)
    {
      idx.push_back(i);
    }

    // distribute salt and pepper pixels
    random_shuffle(idx.begin(), idx.end());

    for (int k = 0; k < corruptPixels; k++){
    	int i = idx[k] / cols;
    	int j = idx[k] % cols;

    	// generate noise value
    	int noisePixelVal = 0; //initialize to all black
		int temp;
		// make half of the noise pixels white
    	if(rand() % 100 < 50){
        	noisePixelVal = 255;
		}

	  imageToCorrupt.getPixelVal(i, j, temp);
	  imageToCorrupt.setPixelVal(i, j, noisePixelVal);
    }

	rescale(cols, rows, imageToCorrupt, writePath); // rescale back to 0-255
}


// this function finds the gradient of an image given a specific mask
void gradient(ImageType& originalImage, int option, char * writePath){
	
	// Prewitt mask in x
	const int prewitt_x[3][3] = {
    	{-1, 0, 1},
    	{-1, 0, 1},
    	{-1, 0, 1}
	};

	// Prewitt mask in y
	const int prewitt_y[3][3] = {
    	{-1, -1, -1},
    	{0, 0, 0},
    	{1, 1, 1}
	};

	// Sobel mask in x
	const int sobel_x[3][3] = {
    	{-1, 0, 1},
    	{-2, 0, 2},
    	{-1, 0, 1}
	};

	// Sobel mask in y
	const int sobel_y[3][3] = {
    	{-1, -2, -1},
    	{0, 0, 0},
    	{1, 2, 1}
	};
	
	int orig_rows, orig_cols, orig_qlvl;
	int maskSize = 3;
   	vector<int> pixelVal_x;
	vector<int> pixelVal_y;
	ImageType prewittX(256, 256, 255);
	ImageType prewittY(256, 256, 255);
	char *writePathPrewittX = "./Images/lenna_prewittX.pgm";     // output image
	char *writePathPrewittY = "./Images/lenna_prewittY.pgm";     // output image

	ImageType SobelX(256, 256, 255);
	ImageType SobelY(256, 256, 255);
	char *writePathSobelX = "./Images/lenna_sobelX.pgm";     // output image
	char *writePathSobelY = "./Images/lenna_SobelY.pgm";     // output image

	/*ImageType prewittImage(256, 256, 255);
	ImageType sobelImage(256, 256, 255);
	char *writePathPrewitt = "./Images/lenna_prewitt.pgm";     // output image
	char *writePathSobel = "./Images/lenna_sobel.pgm";     // output image*/


	originalImage.getImageInfo(orig_rows, orig_cols, orig_qlvl);

	if(option==1){ // prewitt 
		// prewitt implementation attempt
		int temp, temp2;
		// prewitt in x direction
		for(int i=0; i<orig_rows; i++){
			for(int j=0; j<orig_cols; j++){
				// iterate through mask
				for (int k = -maskSize/2; k < maskSize/2; k++){
					for (int l = -maskSize/2; l < maskSize/2; l++){
						// check bounds
						if(i + k < 0 || i + k >= orig_cols){
							pixelVal_x.push_back(0);
						}
						else{
							originalImage.getPixelVal((i+k), (j+l), temp);
							temp2 += temp * prewitt_x[k + maskSize/2][l + maskSize/2];
							pixelVal_x.push_back(temp2);
							prewittX.setPixelVal(i, j, temp2);
						}

					}
				}

			}
		
		}
		rescale(256, 256, prewittX, writePathPrewittX);

		// prewitt in y
		for(int i=0; i<orig_rows; i++){
			for(int j=0; j<orig_cols; j++){
				// iterate through mask
				for (int k = -maskSize/2; k < maskSize/2; k++){
					for (int l = -maskSize/2; l < maskSize/2; l++){
						// check bounds
						if(i + k < 0 || i + k >= orig_cols){
							pixelVal_y.push_back(0);
						}
						else{
							originalImage.getPixelVal((i+k), (j+l), temp);
							temp2 += temp * prewitt_y[k + maskSize/2][l + maskSize/2];
							pixelVal_y.push_back(temp2);
							prewittY.setPixelVal(i, j, temp2);


						}

					}
				}	

			}
		}
		rescale(256, 256, prewittY, writePathPrewittY);

	}

	else if(option==2){ //sobel
		// sobel implementation attempt
		int temp, temp2;
		// sobel in x direction
		for(int i=0; i<orig_rows; i++){
			for(int j=0; j<orig_cols; j++){
				// iterate through mask
				for (int k = -maskSize/2; k < maskSize/2; k++){
					for (int l = -maskSize/2; l < maskSize/2; l++){
						// check bounds
						if(i + k < 0 || i + k >= orig_cols){
							pixelVal_x.push_back(0);
						}
						else{
							originalImage.getPixelVal((i+k), (j+l), temp);
							temp2 += temp * sobel_x[k + maskSize/2][l + maskSize/2];
							pixelVal_x.push_back(temp2);
							SobelX.setPixelVal(i, j, temp2);

							}

					}
				}
			}
		}
		rescale(256, 256, SobelX, writePathSobelX);


		// sobel in y direction
		for(int i=0; i<orig_rows; i++){
			for(int j=0; j<orig_cols; j++){
				// iterate through mask
				for (int k = -maskSize/2; k < maskSize/2; k++){
					for (int l = -maskSize/2; l < maskSize/2; l++){
						// check bounds
						if(i + k < 0 || i + k >= orig_cols){
							pixelVal_y.push_back(0);
						}
						else{
							originalImage.getPixelVal((i+k), (j+l), temp);
							temp2 += temp * sobel_y[k + maskSize/2][l + maskSize/2];
							pixelVal_y.push_back(temp2);
							SobelY.setPixelVal(i, j, temp2);

						}

					}
				}	

			}
		}	
			rescale(256, 256, SobelY, writePathSobelY);

	}



}


// this function finds the laplacian of an image given the laplacian mask
void laplacian(ImageType& originalImage, char* writePath){
	int maskSize = 3; // size of the laplacian mask

	// the laplacian mask
	const int laplacianMask[3][3] = {
    	{0, 1, 0},
    	{1, -4, 1},
    	{0, 1, 0}
	};

	int factor = 0;
	int orig_rows, orig_cols, orig_qlvl;
	originalImage.getImageInfo(orig_rows, orig_cols, orig_qlvl);

	// calculate normalizion factor
	for(int i = 0; i < maskSize; i++){
		for (int j = 0; j < maskSize; j++){
			factor += laplacianMask[i][j];
		}
	}

	// iterate through image pixels
    for(int i = 0; i < orig_cols; i++){
   		for(int j = 0; j < orig_rows; j++) {

   			int pixelVal = 0;

   			// iterate through mask
   			for (int k = -maskSize/2; k < maskSize/2; k++){
   				for (int l = -maskSize/2; l < maskSize/2; l++){
   					
					// check bounds
   					if(i + k < 0 || i + k >= orig_cols || j + l < 0 || j + l >= orig_rows){
   						pixelVal += 0;
   					}
					
   					// calculate output value depending on mask size
   					else{
						int temp;
						originalImage.getPixelVal(i, j, temp);
						//cout << temp << endl;
						pixelVal = pixelVal + (temp * laplacianMask[k + maskSize/2][l + maskSize/2]);
					}
   				}
   			}

   			// update image pixel
			int temp2 = pixelVal;
			originalImage.setPixelVal(i, j, temp2);
   		}
   	}
	
	rescale(orig_rows, orig_cols, originalImage, writePath);
}




// QUESTION 1: CORRELATION
void Q1()
{
	// the size (m x n) and the quantization level for the original image
	int m, n, q;
	// the size (m x n) and the quantization level for the mask image
	int m_mask, n_mask, q_mask;

	int mask_x, mask_y;

	bool type;
	bool maskType;

	char *readPath = "./Images/Image.pgm"; // reading the original image
	char *readMaskPath = "./Images/Pattern.pgm"; // reading the masked image
	char *writePath = "./Images/Image_correlated.pgm"; // the path where the correlated image will be writen

	readImageHeader(readPath, n, m, q, type);
	readImageHeader(readMaskPath, n_mask, m_mask, q_mask, maskType);

	ImageType OriginalImage(n, m, q); // original image
	ImageType MaskImage(n_mask, m_mask, q_mask); // mask image

	OriginalImage.getImageInfo(n, m, q);
	MaskImage.getImageInfo(n_mask, m_mask, q_mask);

	int x = m;
	int y = q;
	int z = n;

	// final (correalted) image
	ImageType CorrelatedImage(z, x, y);

	readImage(readPath, OriginalImage);

	readImage(readMaskPath, MaskImage);

	correlation(OriginalImage, MaskImage, m_mask, n_mask, writePath);

	cout<< "File written to: " << writePath << endl;
}


// QUESTION 2: AVERAGING AND GAUSSIAN SMOOTHING
void Q2()
{
	int q2option;
	cout<<"\n\nSelect from the following options:\n1. Averaging\n2. Gaussian Smoothing\n";
	cout << "Please enter an option: ";
	cin>>q2option;

	if(q2option==1){
		int avgOption;
		int n, m, q;
		bool type;
		cout<<"\n\nSelect from the following options:\n1. 7x7 averaging filter \n2. 15x15 averaging filter\n";
		cout << "Please enter an option: ";
		cin >> avgOption;

		char *readPath = "./Images/lenna.pgm";			// input image
		char *writePath = "./Images/lenna_averaged.pgm"; // output image

		readImageHeader(readPath, n, m, q, type);
		ImageType imageToAverage(n, m, q);
		readImage(readPath, imageToAverage);

		switch (avgOption)
		{
		case 1:
			averaging(imageToAverage, 7, writePath);
			break;
		case 2:
			averaging(imageToAverage, 15, writePath);
			break;
		}
		cout<< "File written to: " << writePath << endl;


	}

	else if(q2option != 1){
		int n, m, q;
		bool type;
		int gaussianOption;
		cout<<"\n\nSelect from the following options:\n1. 7x7 gaussian filter \n2. 15x15 gaussian filter\n";
		cout << "Please enter an option: ";
		cin>>gaussianOption;
		
		char *readPath = "./Images/lenna.pgm";			// input image
		char *writePath = "./Images/lenna_gaussian.pgm"; // output image

		readImageHeader(readPath, n, m, q, type);
		ImageType imageToGaussian(n, m, q);
		readImage(readPath, imageToGaussian);

		switch (gaussianOption)
		{
		case 1:
			gaussian(imageToGaussian, 7, writePath);
			break;
		case 2:
			gaussian(imageToGaussian, 15, writePath);
			break;
		}


		cout<< "File written to: " << writePath << endl;

	}

}


// QUESTION 3: MEDIAN FILTERING
void Q3()
{
	int n, m, q;
	bool type;
	// Variables for mask size and noise percentage
	int maskSize;
	int noise_pct;

	// Variables to get the option selected
	int ms_option;
	cout<<"\n\nSelect from the following Mask Size options:\n1. 7x7 \n2. 15x15\n";
	cout << "Please enter an option: ";
	cin >> ms_option;

	cout << "\nPlease Enter a Noise Percentage (0-100): ";
	cin >> noise_pct;

	// Read and write paths for images
	char *readPath = "./Images/boat.pgm";			// input image
	char *writePath = "./Images/boat_median.pgm"; 	// output image
	char *writePathSP = "./Images/boat_sp.pgm";
	char *readPathSP = "./Images/boat_sp.pgm";
	char *writePathAverage = "./Images/boat_med_avg.pgm";

	// Read original image
	readImageHeader(readPath, n, m, q, type);
	ImageType imageToCorrupt(n, m, q);
	readImage(readPath, imageToCorrupt);


	saltAndPepper(imageToCorrupt, noise_pct, writePathSP);

	// Create a copy of corrupted Image
	int x = m;
	int y = q;
	int z = n;


	// Read Corrupted Image
	readImageHeader(readPathSP, z, x, y, type);
	ImageType corruptedImage(z, x, y);
	readImage(readPathSP, corruptedImage);



	// Uncomment and Comment for whichever filtering you want 
	switch (ms_option)
		{
		case 1:
			medianFilter(corruptedImage, 7, writePath);
			//averaging(corruptedImage, 7, writePathAverage);
			break;
		case 2:
			medianFilter(corruptedImage, 15, writePath);
			//averaging(corruptedImage, 15, writePathAverage);
			break;
		}

	// if doing averaging, change to writePathAverage
	cout<< "File written to: " << writePath << endl;

}


// QUESTION 4: UNSHARP MASKING AND HIGH BOOST FILTERING
void Q4()
{
	int n, m, q;
    bool type;

	int k; // Modify K depending on unsharp or high boost
	cout<<"\nPlease enter desired constant k"<<endl;
	cout<<"k=1 --> Unsharp Masking\nk>1 --> High Boost Filtering"<<endl;
	cout<<"Enter k: ";
	cin>>k;

    char *readPath = "./Images/lenna.pgm";            // input image
    char *writePath = "./Images/lenna_k=1.pgm";     // output image
    char *tmpWritePath = "./Images/lenna_lpSmooth.pgm";
	char *tmpRescale = "./Images/lenna_tmp.pgm";

    readImageHeader(readPath, n, m, q, type);
    ImageType originalImage(n, m, q);
    readImage(readPath, originalImage);

    // Create a copy of original Image
    int x = m;
    int y = q;
    int z = n;

    ImageType lpImage(z, x, y);
    readImage(readPath, lpImage);

    gaussian(lpImage, 7, tmpWritePath);

	// Calculating the gmask
	ImageType gMask(256, 256, 255);
    vector<int> gMaskVal;
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            int temp1, temp2, temp3;
            originalImage.getPixelVal(i, j, temp1);
            lpImage.getPixelVal(i, j, temp2);
            temp3 = temp1 - temp2;
            gMaskVal.push_back(temp3);
            gMask.setPixelVal(i, j, temp3);
        }
    }
	rescale(n, m, gMask, tmpRescale);

	//Calculate final image
	ImageType sharpImage(256, 256, 255);
	vector<int> sharpImageVal;
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            int temp1, temp2, temp3;
            originalImage.getPixelVal(i, j, temp1);
            gMask.getPixelVal(i, j, temp2);
			temp3 = temp1 + (k * temp2);
            sharpImageVal.push_back(temp3);
            sharpImage.setPixelVal(i, j, temp3);
			

        }
    }

	rescale(n, m, sharpImage, writePath);

	cout<< "File written to: " << writePath << endl;


}


// QUESTION 5: GRADIENT AND LAPLACIAN
void Q5()
{
	int n, m, q;
    bool type;

	char *readPath = "./Images/lenna.pgm";            // input image
    char *writePath = "./Images/lenna_sharpened.pgm";     // output image

    readImageHeader(readPath, n, m, q, type);
    ImageType originalImage(n, m, q);
    readImage(readPath, originalImage);

	int gradOrLap; // variable to choose gradient or laplacian
	cout<<"\nSelect from the following options:\n1. Gradient\n2. Laplacian\n";
	cout << "Please enter an option: ";
	cin>> gradOrLap;

	if (gradOrLap==1){
		int maskOption;
		cout<<"\nSelect from the following mask:\n1. Prewitt\n2. Sobel\n";
		cout << "Please enter an option: ";
		cin>> maskOption;

		switch (maskOption){
			case 1:
				gradient(originalImage, 1, writePath); // prewitt
				break;
			case 2:
				gradient(originalImage, 2, writePath); // sobel
				break;
		}

	}

	else if(gradOrLap!=1){
		laplacian(originalImage, writePath);
	}

	cout<< "File written to: " << writePath << endl;

}


int main()
{
	cout << "CS 474 Programming Assignment 1: Team 8\n\n";

	int option;
	cout << "Please select one of the following options:\n";
	cout << "1. Correlation\n2. Averaging and Gaussian Smoothing\n3. Median Filtering\n4. Unsharp Masking and High Boost Filtering\n5. Gradient and Laplacian\n6. Exit\n\n";
	cout << "Please enter an option: ";
	cin >> option;

	switch (option)
	{
	case 1:
		Q1();
		break;
	case 2:
		Q2();
		break;
	case 3:
		Q3();
		break;
	case 4:
		Q4();
		break;
	case 5:
		Q5();
		break;
	}


	return 0;
}