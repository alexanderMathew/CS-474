//CS 474 Project 4
//Team #8: Alexander Mathew & Tal Shahar Zemach

#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <math.h>
#include "image.h"

using namespace std;

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

// extern float rand();         /* ranf() is uniform in 0..1 */ WAS IN BOX-MULLER.c


// Function declarations
int readImage(char fname[], ImageType &image);
int writeImage(char fname[], ImageType &image);
int readImageHeader(char fname[], int &N, int &M, int &Q, bool &type);
void fft(float data[], unsigned long nn, int isign);

//this function rescales values to [0-255]
void rescale(int n, int m, ImageType OriginalImage, char* writePath){

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

// this function rescales values to [0-255] for a float array
void rescaleFloatArr(int n, int m, float array[])
{

	float max = -10000.0;
	float min = 100000.0;

	for (int i = 0; i < n*m; i++){
		if (array[i] > max){
			max=array[i];
		}
		if (array[i] < min){
			min=array[i];
		}
	}

	//cout << min << endl << "\t"<< max << endl;
	// map to [0,255] equation
	for (int i = 0; i < n*m; i++){
		array[i] = 255.0 * ((array[i]) - min) / (max - min);
	}
}

// This function computes the Fast Fourier Transform
// Parameters: elements of the array, the size of the array, and a bool to indicate forward or inverse FFT
// nn must be a power of 2
void fft(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP

// This function compute the 2D DFT using the FFT implementation above
void fft2D(int N, int M, float real_Fuv[], float imag_Fuv[], int isign){
	
	float data[N][M][2];

	// put the real and imaginary parts into the array
	//int index = 0;
	for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
			data[i][j][0] = real_Fuv[i*M+j];
			data[i][j][1] = imag_Fuv[i*M+j];
			//index++;   
        }
    }

	// FFT along the rows 
	for (int i = 0; i < N; i++) { 
		fft(data[i][0]-1, M, isign); 
	}

	// Flip the rows and columns 
    for (int i = 0; i < N; i++){
        for (int j = i; j < M; j++){
			swap(data[i][j][0], data[j][i][0]);
			swap(data[i][j][1], data[j][i][1]);     
        }
    }

	// FFT along the columns
	for (int i = 0; i < M; i++) { 
		fft(data[i][0]-1, N, isign); 
	}

	// Multiply by correlation factor
	float corrFactor = 1 / sqrt(N * M);
	for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
  			data[i][j][0] *= corrFactor; 
			data[i][j][1] *= corrFactor; 
        }
    }

	// Flip the columns and rows 
    for (int i = 0; i < N; i++){
        for (int j = i; j < M; j++){
			swap(data[i][j][0], data[j][i][0]);
			swap(data[i][j][1], data[j][i][1]);     
        }
    }

	// put the DFT values back into the real and imaginary arrays
	for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
			real_Fuv[i*M+j] = data[i][j][0];
			imag_Fuv[i*M+j] = data[i][j][1];
        }
    }

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


// Function to perform convolution with a sobel mask in the spatial domain
void sobelConvolution(ImageType& image, char * writePath){

	int maskSize = 3;

	// The sobel mask (3x3)
	const int sobel[3][3] = {
		{ -1, 0, 1},
		{ -2, 0, 2},
		{ -1, 0, 1}
	};

	//Image originalImage = Image(image);
	int rows, cols, qlvl;

	image.getImageInfo(rows, cols, qlvl);

	// Vector to store values
	vector<int> values;
	int temp1;
	int temp2;

	ImageType finalImage(rows, cols, qlvl);

	// iterate through image pixels
    for(int i = 0; i < rows; i++){
   		for(int j = 0; j < cols; j++) {

   			int sum = 0;

   			// iterate over mask
   			for (int k = -maskSize/2; k <= maskSize/2; k++){ //(int k = -maskSize/2; k < maskSize/2; k++)
   				for (int l = -maskSize/2; l <= maskSize/2; l++){ //(int l = -maskSize/2; l < maskSize/2; l++)
   				
   					//check bounds
   					if(i + k < 0 || i + k >= rows || j + l < 0 || j + l >= cols){
   						sum += 0;
					}

   					else{
						int num1 = i+k;
						int num2 = j+l;
						image.getPixelVal(num1,num2, temp1);
						temp2 = sobel[maskSize/2 - k][maskSize/2 - l];
						//IT WAS THIS: temp2 = sobel[k+maskSize/2][l+maskSize/2];
						sum = sum + (temp1 * temp2);
					}
   				}
   			}

   			values.push_back(sum);
			finalImage.setPixelVal(i, j, sum);
   		}
   	}

	rescale(rows, cols, finalImage, writePath);
}

// a function to find the spectrum of the image starting at the spatial domain
// this function will fetch pixel values from original image, center, take the FFT,
void findSpectrum(int n, int m, float realArray[], float imagArray[], ImageType& originalImage, ImageType&spectrumImage, char * spectrumImageWritePath){
	// get the pixel value from the original image
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			int temp;
			originalImage.getPixelVal(i, j, temp); 
			realArray[i*m+j] = temp;
		}
	}
	
	// center using the pow eq
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			realArray[i*m+j] = realArray[i*m+j] * pow(-1, i+j); // centering
		}
	}

	// FFT
	fft2D(n, m, realArray, imagArray, -1);

	float magnitudeSpectrum[n*m]; // temporary magnitude array used to find the magnitude of the spectrum of the original image

	// calculate magnitude
	for(int i=0; i<n*m; i++){
		magnitudeSpectrum[i] = sqrt (realArray[i] * realArray[i] + imagArray[i] * imagArray[i]); // magnitude = sqrt (real^2 + imaginary^2)
	}

	// set the pixel values to the magnitude
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){	
			int temp;
			spectrumImage.setPixelVal(i, j, magnitudeSpectrum[i*m+j]);
		}
	}

	writeImage(spectrumImageWritePath, spectrumImage);
}		



// The band reject filter used for part 1a
int bandRejectFilter(int u, int v){
	float d0, w, D;
	d0 = 35.777;
	w = 2.0;
	D = sqrt((u*u)+(v*v));

	// uncomment for band or notch 
	if((d0 - w/2) <= D && D <= d0 + w/2){
		return 0;
	}
	else{
		return 1;
	}
	 
}

// The band reject filter used for part 1b
int bandRejectFilterNoise(int u, int v){
	float d0, w, D;
	d0 = 35.777;
	w = 2.0;
	D = sqrt((u*u)+(v*v));

	// uncomment for band or notch 
	if((d0 - w/2) <= D && D <= d0 + w/2){
		return 1;
	}
	else{
		return 0;
	}
	 
}

// The notch filter used for part 1a
int notchFilter(int u, int v){

	if(abs(u) == 16 && abs(v) == 32){
		return 0;
	}
	else{
		return 1;
	}
	 
}

// The notch filter used for part 1b
int notchFilterNoise(int u, int v){

	if(abs(u) == 16 && abs(v) == 32){
		return 1;
	}
	else{
		return 0;
	}
	 
}

// the mtion blur degredation filter used for part 3
complex<float> motionBlur(int u, int v, float T, float a, float b){
	float temp;
	if(u+v==0){
		temp = T;
	}
	else{	
	 temp = (T/(M_PI*((u*a)+(v*b))))*sin(M_PI*((u*a)+(v*b)));
	}
	/*
	using Euler's Formula: e^-j(theta)=cos(theta)-jsin(theta)
	*/
	//float temp2 = M_PI*((u*a)+(v*b));
	//temp = temp * temp2;

	return temp * exp(complex<float>(0, -M_PI*((u*a)+(v*b))));
}

// the homomorphic filter used for part 4
float homomorphicFilter(int u, int v, float gamma_H, float gamma_L, float c, float D_0){
	float temp = gamma_H - gamma_L;
	float temp2 = -c*((u*u+v*v)/(D_0*D_0));
	float temp3 = 1 - exp(temp2);
	float result = temp * temp3 + gamma_L;
	return result;
}




// Experiment 1: attempting to remove cosine (periodic) noise and isolating the noise pattern
void exp1(int experimentPart){

	// rows and cols of input image
	int n, m, q;
    bool type;

	char *readPath = "./Images/boy_noisy.pgm";          // input image
    char *writePath = "./Images/part1/boy_noiseFree.pgm";     // output image
	char *writePath2 = "./Images/part1/boy_Gaussian.pgm";	// output after conducting gaussian filter
	char *writePathSpectrumFilter = "./Images/part1/boy_SpectrumFilter.pgm";	// output the filter of boy_noisy


	// select whether band-reject for notch filter
	int filterOption;
	cout << "Please select from one of the following options: \n";
	cout << "1. Band-reject filter\n";
	cout << "2. Notch Filter\n";
	cout << "Please enter an option: ";
	cin >> filterOption;

	// select whether 7x7 or 15x15 gaussian noise
	int gaussianOption;
	cout << "Please select from one of the following options: \n";
	cout << "1. Attempt with 7x7 gaussian noise\n";
	cout << "2. Attempt with 15x15 gaussian noise\n";
	cout << "Please enter an option: ";
	cin >> gaussianOption;

	// read input image 
    readImageHeader(readPath, n, m, q, type);
    ImageType originalImage(n, m, q);
	ImageType originalImageSpectrumFilter(n, m, q);
    readImage(readPath, originalImage);

	float realArray[n*m];
	float imagArray[n*m];

	// initialize the imaginary array to 0
	for(int i=0; i<n*m; i++){
		imagArray[i] = 0;
	}

	// initialize the real array according to pixel values from originalImage
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			int temp;
			originalImage.getPixelVal(i, j, temp); 
			realArray[i*m+j] = temp;
			realArray[i*m+j] = realArray[i*m+j]*pow(-1, i+j); // shift to the center
		}
	}

	// find the spectrum of boy
	ImageType originalImageSpectrum(n, m, q);
    char *writePathBoySpectrum = "./Images/part1/boy_spectrum.pgm";     // output image
	float realSpectrumArray[n*m]; // a real array just for calculating the spectrum of boy
	float imagSpectrumArray[n*m]; // an imaginary array just for calculating the spectrum of boy
	findSpectrum(n, m, realSpectrumArray, imagSpectrumArray, originalImage, originalImageSpectrum, writePathBoySpectrum);

	// Conduct 2D FFT
	fft2D(n, m, realArray, imagArray, -1);


	
	// Part 1a: Noise-free
	if(experimentPart==1){
		// band-reject filter
		if(filterOption==1){
			// Create the band reject filter in the frequency domain
			for(int i=0; i<n; i++){
				for(int j=0; j<m; j++){
					int h = bandRejectFilter(i-n/2, j-m/2);
					realArray[i*m+j] = realArray[i*m+j] * h;
					imagArray[i*m+j] = imagArray[i*m+j] * h;
				}
			}
		}

		// notch filter
		if(filterOption==2){
			// Create the notch filter in the frequency domain
			for(int i=0; i<n; i++){
				for(int j=0; j<m; j++){
					int notch = notchFilter(i-n/2, j-m/2);
					realArray[i*m+j] = realArray[i*m+j] * notch;
					imagArray[i*m+j] = imagArray[i*m+j] * notch;
				}
			}
		}
	}

	// Part 1b: Isolate noise
	if(experimentPart==2){
		// band-reject filter
		if(filterOption==1){
			// Create the band reject filter in the frequency domain
			for(int i=0; i<n; i++){
				for(int j=0; j<m; j++){
					int h = bandRejectFilterNoise(i-n/2, j-m/2);
					realArray[i*m+j] = realArray[i*m+j] * h;
					imagArray[i*m+j] = imagArray[i*m+j] * h;
				}
			}
		}

		// notch filter
		if(filterOption==2){
			// Create the notch filter in the frequency domain
			for(int i=0; i<n; i++){
				for(int j=0; j<m; j++){
					int notch = notchFilterNoise(i-n/2, j-m/2);
					realArray[i*m+j] = realArray[i*m+j] * notch;
					imagArray[i*m+j] = imagArray[i*m+j] * notch;
				}
			}
		}
	}


	float magnitude[n*m]; // magnitude array for visualization
	// calculate magnitude
	for(int i=0; i<n*m; i++){
		magnitude[i] = sqrt (realArray[i] * realArray[i] + imagArray[i] * imagArray[i]); // magnitude = sqrt (real^2 + imaginary^2)
	}

	// set the pixel values to the magnitude
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){	
			int temp;
			originalImageSpectrumFilter.setPixelVal(i, j, magnitude[i*m+j]);
		}
	}

	// write spectrum of the original image with the filter
	writeImage(writePathSpectrumFilter, originalImageSpectrumFilter);


	// Inverse DFT
	fft2D(n, m, realArray, imagArray, 1);

	// initialize the real array according to pixel values from originalImage
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			realArray[i*m+j] = realArray[i*m+j]*pow(-1, i+j);
			originalImage.setPixelVal(i, j, realArray[i*m+j]); // conducted centering before-now need to 'undo' it
	
		}
	}
	
	writeImage(writePath, originalImage); // writing boy without noise
	
	
	// utilize gaussian filtering in the spatial domain for comparison purposes
	if(gaussianOption==1){
		gaussian(originalImage, 7, writePath2); // 7x7 gaussian filter
	}
	else if(gaussianOption == 2){
		gaussian(originalImage, 15, writePath2); // 15x15 gaussian filter
	}


}



// Experiment 2: Obtain a frequency domain transfrer function from a spatial filter
void exp2(){

	int n, m, q;
    bool type;

	char *readPath = "./Images/lenna.pgm";            // input image
    char *writePath = "./Images/part2/lenna_convoluted.pgm";     // output image

    readImageHeader(readPath, n, m, q, type);
    ImageType originalImage(n, m, q);
    readImage(readPath, originalImage);

	int spatialOrFreq;
	cout<<"Do you want to perform the Convolution in the spatial domain or the frequency domain?\n";
	cout<<"1. Spatial Domain\n2. Frequency Domain\n\nPlease enter an option: ";
	cin>>spatialOrFreq;

	if(spatialOrFreq==1){
		sobelConvolution(originalImage, writePath);
	}

	if(spatialOrFreq==2){
		int p_size = n*2;
		int q_size = m*2;
		ImageType paddedImage(p_size, q_size, q);


		float realArray[n*m];
		float imagArray[n*m];
		float paddedrealArray[p_size*q_size]; // this is a padded realArray
		float paddedimagArray[p_size*q_size]; // this is a padded imagArray

		// conduct padding to the image
    	for (int i = 0; i < n; i++){
        	for (int j = 0; j < m; j++){
           		int temp;
				originalImage.getPixelVal(i, j, temp);
            	paddedImage.setPixelVal(i, j, temp);
        	}
    	}

    	for (int i = n; i < p_size; i++){
        	for (int j = m; j < q_size; j++){
        	    paddedImage.setPixelVal(i, j, 0);
        	}
    	}

		
		// initialize the imaginary array to 0
		for(int i=0; i<n*m; i++){
			paddedrealArray[i] = 0;
		}

		// initialize the real array according to pixel values from paddedImage
		for(int i=0; i < p_size; i++){
			for(int j=0; j < q_size; j++){
				int temp;
				paddedImage.getPixelVal(i, j, temp); 
				paddedrealArray[i*q_size+j] = temp;
			}
		}


		// find the spectrum of lenna
		ImageType paddedLennaSpectrum(n, m, q);
    	char *writePathPaddedLenna = "./Images/part2/paddedLenna_spectrum.pgm";     // output image
		findSpectrum(n, m, realArray, imagArray, originalImage, paddedLennaSpectrum, writePathPaddedLenna);
		
		
		
		// center using the pow eq
		for(int i=0; i < p_size; i++){
			for(int j=0; j < q_size; j++){
				paddedrealArray[i*q_size+j] = paddedrealArray[i*q_size+j] * pow(-1, i+j); // centering
			}
		}


		
		// specifiy the sobel filter in the spatial domain
		// the sobel filter is padded with zeros to preserve odd symmetry
		const int filter_size = 4;

		// sobel filter while preserving the odd symmetry
		const int sobel[filter_size][filter_size] = {
			{0, 0, 0, 0},
			{0, -1, 0, 1},
			{0, -2, 0, 2},
			{0, -1, 0, 1}
		};

		float sobelrealArray[filter_size*filter_size];
		float sobelimagArray[filter_size*filter_size];
		float sobelpaddedrealArray[p_size*q_size]; // this is a padded realArray
		float sobelpaddedimagArray[p_size*q_size]; // this is a padded imagArray

		ImageType paddedSobel(p_size, q_size, q);


		// conduct padding to the sobel filter
    	for (int i = 0; i < filter_size; i++){
        	for (int j = 0; j < filter_size; j++){
           		int temp = sobel[i][j];
            	paddedSobel.setPixelVal(i, j, temp);
        	}
    	}

    	for (int i = filter_size; i < p_size; i++){
        	for (int j = filter_size; j < q_size; j++){
        	    paddedSobel.setPixelVal(i, j, 0);
        	}
    	}

		// initialize the imaginary array to 0
		for(int i=0; i<p_size*q; i++){
			sobelpaddedrealArray[i] = 0;
		}

		// initialize real array to 0
		for(int i=0; i<p_size*q_size; i++){
			sobelpaddedrealArray[i] = 0;
		}
		
		// get the pixel values from paddedSobel
		for (int i=0; i<filter_size; i++) {
			for (int j=0; j<filter_size; j++) {
				int temp;
				paddedSobel.getPixelVal(i, j, temp); 
				sobelpaddedrealArray[i*q_size+j] = temp;
			}
		}
		
		
		// Center the padded sobel mask using the centering equation
		for(int i=0; i < p_size; i++){
			for(int j=0; j < q_size; j++){
				sobelpaddedrealArray[i*q_size+j] = sobelpaddedrealArray[i*q_size+j]*pow(-1, i+j); // centering
			}
		}

		// take the forward FFT of Lenna
		fft2D(p_size, q_size, paddedrealArray, paddedimagArray, -1);

		// take the forward FFT of Sobel Filter
		fft2D(p_size, q_size, sobelpaddedrealArray, sobelpaddedimagArray, -1);

		/*
		////////////////////////////////////////////
		// Finding the spectrum of padded sobel and outputting it
		ImageType sobelSpectrum(p_size, q_size, q);
		char *writePathSobelSpectrum = "./Images/paddedsobel_Spectrum.pgm";

		float magnitudeSpectrumSobel[p_size*q_size]; // used to find the magnitude of the spectrum of the original image

		// calculate magnitude
		for(int i=0; i<p_size*q_size; i++){
			magnitudeSpectrumSobel[i] = sqrt(sobelpaddedrealArray[i] * sobelpaddedrealArray[i] + sobelpaddedimagArray[i] * sobelpaddedimagArray[i]); // magnitude = sqrt (real^2 + imaginary^2)
		}

		// rescale and writing the spectrum of the original image
		rescaleFloatArr(p_size, q_size, magnitudeSpectrumSobel);

		// set the pixel values to the magnitude
		for(int i=0; i < p_size; i++){
			for(int j=0; j < q_size; j++){	

				sobelSpectrum.setPixelVal(i, j, magnitudeSpectrumSobel[i*q_size+j]);
			}
		}

		writeImage(writePathSobelSpectrum, sobelSpectrum);
		////////////////////////////////////////////

		
		////////////////////////////////////////////
		// Finding the spectrum of lenna and outputting it
		ImageType originalImageSpectrum(n, m, q);
		char *writePathLennaSpectrum = "./Images/Lenna_Spectrum.pgm";
		
		float magnitudeSpectrum[n*m]; // used to find the magnitude of the spectrum of the original image

		// calculate magnitude
		for(int i=0; i<n*m; i++){
			magnitudeSpectrum[i] = sqrt(realArray[i] * realArray[i] + imagArray[i] * imagArray[i]); // magnitude = sqrt (real^2 + imaginary^2)
		}

		// set the pixel values to the magnitude
		for(int i=0; i < n; i++){
			for(int j=0; j < m; j++){	
				int temp;
				originalImageSpectrum.setPixelVal(i, j, magnitudeSpectrum[i*m+j]);
			}
		}

		// writing the spectrum of the original image
		writeImage(writePathLennaSpectrum, originalImageSpectrum);	
		////////////////////////////////////////////
		*/
		
		
		// complex multiplication with the sobelFilter
		for(int i=0; i<p_size; i++){ //rows of the image
			for(int j=0; j<q_size; j++){ //columns of the image
				/*
				This is the complex multipication equation I need to utilize
				(a+bi) * (c+di) = (ac-bd) + (ad+bc)i
				*/
				paddedrealArray[i*q_size+j] = paddedrealArray[i*q_size+j] * sobelpaddedrealArray[i*q_size+j] -  paddedimagArray[i*q_size+j] * sobelpaddedimagArray[i*q_size+j];
            	paddedimagArray[i*q_size+j] = paddedrealArray[i*q_size+j] * sobelpaddedimagArray[i*q_size+j] + paddedimagArray[i*q_size+j] * sobelpaddedrealArray[i*q_size+j];
			}
		}

		///////////////////////////////////////
		float magnitude[p_size*q_size]; // magnitude array for visualization
		ImageType originalImageSpectrumConvoluted(p_size, q_size, q);

		// visualizing the pixel-by-pixel multipiclation
		// calculate magnitude
		for(int i=0; i<q_size*q_size; i++){
			// commented line below can be used to visualize spectrum of sobel
			// magnitude[i] = sqrt (sobelpaddedrealArray[i] * sobelpaddedrealArray[i] + sobelpaddedimagArray[i] * sobelpaddedimagArray[i]); // magnitude = sqrt (real^2 + imaginary^2)
			magnitude[i] = sqrt (paddedrealArray[i] * paddedrealArray[i] + paddedimagArray[i] * paddedimagArray[i]); // magnitude = sqrt (real^2 + imaginary^2)
		}

		rescaleFloatArr(p_size, q_size, magnitude);

		// set the pixel values to the magnitude
		for(int i=0; i < p_size; i++){
			for(int j=0; j < q_size; j++){	
				int temp;
				originalImageSpectrumConvoluted.setPixelVal(i, j, magnitude[i*q_size+j]);
			}
		}
		char *writePathLennaSpectrumConvoluted = "./Images/part2/Lenna_convolutedFreq.pgm";
		writeImage(writePathLennaSpectrumConvoluted, originalImageSpectrumConvoluted);
		///////////////////////////////////////
		

		// take inverse DFT
		fft2D(p_size, q_size, paddedrealArray, paddedimagArray, 1);
		
		// remove imaginary part
		for(int i=0; i<p_size*q_size; i++){
			paddedimagArray[i]=0;
		}
		
		// undo centering
		for(int i=0; i < p_size; i++){
			for(int j=0; j < q_size; j++){
				paddedrealArray[i*p_size+j] = paddedrealArray[i*p_size+j]*pow(-1, i+j);
			}
		}

		rescaleFloatArr(p_size, q_size, paddedrealArray);

		ImageType finalImage(n, m, q); // the final image will only be 256*256

		// set the real part to the new array
		for(int i=0; i < p_size/2; i++){ // p_size/2 because the final image will only be 256*256
			for(int j=0; j < q_size/2; j++){
				finalImage.setPixelVal(i, j, paddedrealArray[i*p_size+j]);
			}
		}
		

		char *writePathFinal = "./Images/part2/Lenna_final_convolutedFreqFinal.pgm";		


		// write spectrum of the original image with the filter
		writeImage(writePathFinal, finalImage);
	}
	
}



// Experiment 3: apply motion blur and then restore Lenna
void exp3(int experimentPart){

	int mu = 0; // mu used for box_muller
	int sigma = 1; // sigma used for box_muller

	int n, m, q;
    bool type;

	char *readPath = "./Images/lenna.pgm";            // input image
    char *writePath = "./Images/part3/lenna_part3.pgm";     // output image

    readImageHeader(readPath, n, m, q, type);
    ImageType originalImage(n, m, q);
    readImage(readPath, originalImage);

	
	////////////////////////////////////////////////////////////////
	// Generate the noise in the spatial domain, take its FT
	float noiserealArray[n*m];
	float noiseimagArray[n*m];

	// setting the imaginary part equal to 0
	for(int i=0; i<n*m; i++){
		noiseimagArray[i] = 0;
	}

	// setting the real part according to the pixel values from the array
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution(mu, sigma);
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			noiserealArray[i*m+j] = distribution(generator) * pow(-1, i+j); // used normal distribution and centered the array
		}
	}


	// take the FT of the noise
	fft2D(n, m, noiserealArray, noiseimagArray, -1);
	////////////////////////////////////////////////////////////////
	


	/////////////////////////////////////////////
	// setup the original image and take the FT
	float realArray[n*m];
	float imagArray[n*m];

	// setting the imaginary part equal to 0
	for(int i=0; i<n*m; i++){
		imagArray[i] = 0;
	}

	// setting the real part according to the pixel values from the array
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			int temp;
			originalImage.getPixelVal(i, j, temp); 
			realArray[i*m+j] = temp;
		}
	}

	// conduct centering
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			realArray[i*m+j] = realArray[i*m+j]*pow(-1, i+j);
		}
	}
	
	// take the forward FFT
	fft2D(n, m, realArray, imagArray, -1);
	/////////////////////////////////////////////

	// specify array according to the original image
	complex<float> gArray[n*m];
	for(int i=0; i<n*m; i++){
		gArray[i] = complex<float>(realArray[i], imagArray[i]);
	}

	// specify nArray according to the noise
	complex<float> nArray[n*m];
	for(int i=0; i<n*m; i++){
		nArray[i] = complex<float>(noiserealArray[i], noiseimagArray[i]);
	}

	// complex multipication between image and motion blur G(u,v) = H(u,v)F(u,v)
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			gArray[i*m+j] *= motionBlur(i-n/2, j-m/2, 1.0, 0.1, 0.1);
		}		
	}

	///////////////////////////////////////////////////////////////
	/*// this code is for visualize the motion blur WITHOUT NOISE
	// If you do not want to visualize this, then comment this whole section.
	ImageType motionBlurImage(n, m, q);
	char *writePathMB = "./Images/part3/lenna_MB.pgm"; // noise generated in spatial domain for part 3
	float realArrayMB[n*m]; 
	float imagArrayMB[n*m];

	for(int i=0; i<n*m; i++){
		realArrayMB[i] = gArray[i].real();
		imagArrayMB[i] = gArray[i].imag();
	}

	fft2D(n, m, realArrayMB, imagArrayMB, 1);

	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			realArrayMB[i*m+j] = realArrayMB[i*m+j]*pow(-1, i+j); //undo centering
		}
	}

	rescaleFloatArr(n, m, realArrayMB);

	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			motionBlurImage.setPixelVal(i, j, realArrayMB[i*m+j]);
		}
	}

	writeImage(writePathMB, motionBlurImage);*/
	///////////////////////////////////////////////////////////////

	// adding noise to G: G(u,v) = H(u,v)F(u,v) + N(u,v)
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			gArray[i*m+j] += nArray[i*m+j];
		}		
	}

	// need to set fHatArray from complex back to what it was before
	for(int i=0; i<n*m; i++){
		realArray[i] = gArray[i].real();
		imagArray[i] = gArray[i].imag();
	}

	///////////////////////////////////////////////////////////////
	//this code is for visualize the motion blur WITH noise
	//If you do not want to visualize this, then comment this whole section.
	ImageType noisyImage(n, m, q);
	char *writePathNoise = "./Images/part3/Lenna_MB_Noise.pgm"; // noise generated in spatial domain for part 3
	float realArrayMBNoise[n*m]; 
	float imagArrayMBNoise[n*m];

	for(int i=0; i<n*m; i++){
		realArrayMBNoise[i] = realArray[i];
		imagArrayMBNoise[i] = imagArray[i];
	}

	fft2D(n, m, realArrayMBNoise, imagArrayMBNoise, 1);

	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			realArrayMBNoise[i*m+j] = realArrayMBNoise[i*m+j]*pow(-1, i+j); //undo centering
		}
	}

	rescaleFloatArr(n, m, realArrayMBNoise);

	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			noisyImage.setPixelVal(i, j, realArrayMBNoise[i*m+j]);
		}
	}

	writeImage(writePathNoise, noisyImage);
	///////////////////////////////////////////////////////////////


	// implement inverse filtering F^(u,v)=G(u,v)/H(u,v)
	complex<float> fHatArray[n*m];
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			fHatArray[i*m+j] = gArray[i*m+j] / motionBlur(i-n/2, j-m/2, 1.0, 0.1, 0.1); //F^ = G / H
		}		
	}

	// need to set fHatArray from complex back to what it was before
	for(int i=0; i<n*m; i++){
		realArray[i] = fHatArray[i].real();
		imagArray[i] = fHatArray[i].imag();
	}


	// Inverse DFT
	fft2D(n, m, realArray, imagArray, 1);



	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			realArray[i*m+j] = realArray[i*m+j]*pow(-1, i+j); // conducted centering before, now need to 'undo' it
	
		}
	}

	rescaleFloatArr(n, m, realArray);

	// initialize the real array according to pixel values from originalImage
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			originalImage.setPixelVal(i, j, realArray[i*m+j]);
	
		}
	}

	writeImage(writePath, originalImage);
}


// Experiment 4: conduct homomorphic filtering on a given image
void exp4(){

	int n, m, q;
	bool type;

	char *readPath = "./Images/girl.pgm";            // input image
	char *writePath = "./Images/part4/girl_filtered.pgm";     // output image

	readImageHeader(readPath, n, m, q, type);
	ImageType originalImage(n, m, q);
	readImage(readPath, originalImage);

	float realArray[n*m];
	float imagArray[n*m];

	// initialize the imaginary array to 0
	for(int i=0; i<n*m; i++){
		imagArray[i] = 0;
	}

	// initialize the real array according to pixel values from the padded image
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			int temp;
			originalImage.getPixelVal(i, j, temp); 
			realArray[i*m+j] = temp;
		}
	}

	// applying log transformation for real part to ensure illuminance and reflection are manipulated separately
	for(int i=0; i<n*m; i++){
		realArray[i] = log(1+realArray[i]);                   
	}

	// shift to the center
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			realArray[i*m+j] = realArray[i*m+j]*pow(-1, i+j);
		}
	}

	// Conduct 2D FFT
	fft2D(n, m, realArray, imagArray, -1);

	// Find the spectrum of the original image
	ImageType originalImageSpectrum(n, m, q);
	char *writePathSpectrum = "./Images/part4/girl_Spectrum.pgm";	// the spectrum of the original image
	float realArraySpectrum[n*m];
	float imagArraySpectrum[n*m];
	findSpectrum(n, m, realArraySpectrum, imagArraySpectrum, originalImage, originalImageSpectrum, writePathSpectrum);


	// create the homomorphic filter and pixel-by-pixel multiply
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			float f = homomorphicFilter(i-n/2, j-m/2, 1.5, 0.5, 1.0, 1.8);
			//cout<<f<<endl;
			realArray[i*m+j] = realArray[i*m+j] * f;
			imagArray[i*m+j] = imagArray[i*m+j] * f;
		}
	}

	//////////////////////////////////////////////
	// pixel-by-pixel visualization
	ImageType finalImage(n, m, q);
	char* writePathFinal = "./Images/part4/girl_freq.pgm";
	float magnitude[n*m]; // magnitude array for visualization
	// calculate magnitude
	for(int i=0; i<n*m; i++){
		magnitude[i] = sqrt (realArray[i] * realArray[i] + imagArray[i] * imagArray[i]); // magnitude = sqrt (real^2 + imaginary^2)
	}

	// set the pixel values to the magnitude
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){	
			finalImage.setPixelVal(i, j, magnitude[i*m+j]);
		}
	}

	// write spectrum of the original image with the filter
	writeImage(writePathFinal, finalImage);
	//////////////////////////////////////////////

	// Inverse DFT
	fft2D(n, m, realArray, imagArray, 1);

	// initialize the real array according to pixel values from originalImage
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			realArray[i*m+j] = realArray[i*m+j] * pow(-1, i + j); // undo centering
			realArray[i*m+j] = exp(realArray[i*m+j]) - 1; // conducting an exponential transformation
		}
	}

	rescaleFloatArr(n, m, realArray);

	// initialize the real array according to pixel values from originalImage
	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			originalImage.setPixelVal(i, j, realArray[i*m+j]);
		}
	}

	writeImage(writePath, originalImage); // writing girl without noise
}


int main()
{
	cout << "CS 474 Programming Assignment 4: Team 8\n\n";

	int option;
	cout << "Please select one of the following options:\n";
	cout << "1. Experiment 1a: Remove noise of a given image\n";
	cout << "2. Experiment 1b: Extract the noise pattern of a given image\n\n";
	cout << "3. Experiment 2: Conduct convolution in the frequency domain for Lenna.pgm (filter specified in spatial domain)\n\n";
	cout << "4. Experiment 3a: Generate motion blur and gaussian noise to Lenna.pgm and apply Inverse Filtering\n";
	cout << "5. Experiment 3b: Generate motion blur and gaussian noise to Lenna.pgm and apply Wiener Filtering\n\n";
	cout << "6. Experiment 4: Conduct homomorphic filtering on a given image\n\n";
	cout << "7. Exit\n\n\n";
	cout << "Please enter an option: ";
	cin >> option;

	switch (option)
	{
	case 1:
		exp1(1);
		break;
	case 2:
		exp1(2);
		break;
	case 3:
		exp2();
		break;
	case 4:
		exp3(1);
		break;
	case 5:
		exp3(2);
		break;
	case 6:
		exp4();
		break;
	}


	return 0;
}