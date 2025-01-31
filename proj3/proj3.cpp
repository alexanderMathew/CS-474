//CS 474 Project 2
//Team #8: Alexander Mathew & Tal Shahar Zemach

#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "image.h"

using namespace std;

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

// Function declarations
int readImage(char fname[], ImageType &image);
int writeImage(char fname[], ImageType &image);
int readImageHeader(char fname[], int &N, int &M, int &Q, bool &type);
void fft(float data[], unsigned long nn, int isign);


// this function rescales values to [0-255] for a float array
void rescale(int n, int m, float array[])
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



// Example 1a: compute FFT of {2,3,4,4}
void exp1a(){
	// f[0] isn't being used
	// odd is the real part
	// even is the imaginary part
	float f[] = {0, 2, 0, 3, 0, 4, 0, 4, 0, 0}; // given signal
	float n = 4.0; // length of the input

	// Perform the FFT
	fft(f, 4, -1);

	for (int i = 0; i < 9; i++) {
		f[i] = f[i] / n; // divide by the normalization factor
	}

	int printOption;
	cout<<"\nSelect which values you want printed:";
	cout<< "\n1. DFT of [2, 3, 4, 4]\n2. Real part of [2, 3, 4, 4]\n3. Imaginary Part of [2, 3, 4, 4]\n";
	cout << "4. Magnitude of the [2, 3, 4, 4]\n5. Use Inverse DFT to Verify Results\n\nPlease Enter a value: ";
	cin>>printOption;

	if(printOption==1){
		// Output results from DFT
		cout<<"\nDFT of {2, 3, 4, 4}:\n";
		for (int i = 0; i < 9; i++) {
			cout << "F("<< i <<") = " << f[i] << endl;
		}
	}

	if(printOption==2){
		//OUTPUT REAL
		cout<<"\n\nREAL PART\n";
		for(int i = 0; i < 4; i++){
			cout<<f[2*i+1]<<endl;
		}
	}

	if(printOption==3){
		//OUTPUT IMAGINARY
		cout<<"\n\nIMAGINARY PART\n";
		for(int i = 0; i < 4; i++){
			cout<<f[2*(i+1)]<<endl;
		}
	}

	if(printOption==4){
		//OUTPUT MAGNITUDE
		cout<<"\n\nMAGNITUDE PART\n";
		for(int i = 0; i < 4; i++){
			float real = f[2*i+1];
			float imaginary = f[2*(i+1)];
			real = real * real; // real^2
			imaginary = imaginary * imaginary; // imaginary^2
			float total = real + imaginary; // real^2 + imaginary^2
			total = sqrt (total); // magnitude = sqrt (real^2 + imaginary^2)
			cout << total <<endl;
		}
	}

	if(printOption==5){
	//  Inverse DFT to verify results
		fft(f, 4, 1);
		cout<<"\n\n\nUsing Inverse DFT to verify forward DFT:\n";
		for (int i = 0; i < 9; i++) {
			cout << "f(F("<<i<<")) =  " << f[i] << endl;
		}	
	}
}


// Example 1b: compute FFT of cosine fucntion with 128 samples
void exp1b(){
	float step = 1.0 / 128.0;
	float cosine[128];
	for (int i = 0; i <= 128; i++) {
		cosine[i] = 0;
	}
	float n = 128.0;

	// Calculate cosine function using cos(2*pi*u0*x)
	for (int i = 0; i < 128; i++) {
		cosine[i] = cos( 2.0 * M_PI * 8.0 * (float)(i * step) );

		cosine[i] = cosine[i] * pow(-1, i); //  shift magnitude to center of frequency domain
	}


	// intialize complexCosine, which is the 1D array with real and imaginary parts
	float complexCosine[257];
	for (int i = 0; i < 257; i++) {
		complexCosine[i] = 0; // initialization
	}


	// making the cos into an array with real and imaginary parts
	for(int i=0; i < 128; i++){
		complexCosine[2*i+1]=cosine[i]; // putting all of the real part into the complexCosine
	}
	

	// Perform the DFT on the Cosine samples
	fft(complexCosine, 128, -1);

	for (int i = 0; i < 257; i++){
		complexCosine[i] = complexCosine[i] / n; // divide by the normalization factor
	}

	int printOption;
	cout<<"\nSelect which values you want printed:";
	cout<< "\n1. DFT of the Cosine function\n2. Real part of the cosine function\n3. Imaginary Part of the cosine function\n";
	cout<< "4. Magnitude of the cosine funciton\n5. Phase of the cosine function\n6. Use inverse DFT to verify results\n\nPlease Enter a value: ";
	cin>>printOption;

	if(printOption==1){
		// Ouput DFT cosine results
		cout << "\nDFT of the cosine function" << endl;
		cout<<"\nDFT of cosine:\n";
		for (int i = 0; i < 257; i++) {
			cout<<complexCosine[i]<<endl;
		}
	}

	if(printOption==2){
		//OUTPUT REAL
		cout<<"\n\nREAL PART\n";
		for(int i = 0; i < 128; i++){
			cout<<complexCosine[2*i+1]<<endl;
		}
	}

	if(printOption==3){
		//OUTPUT IMAGINARY
		cout<<"\n\nIMAGINARY PART\n";
		for(int i = 0; i < 128; i++){
			cout<<complexCosine[2*(i+1)]<<endl;
		}
	}
	
	
	if(printOption==4){
		//OUTPUT MAGNITUDE
		cout<<"\n\nMAGNITUDE PART\n";
		for(int i = 0; i < 128; i++){
			float real = complexCosine[2*i+1];
			float imaginary = complexCosine[2*(i+1)];
			real = real * real; // real^2
			imaginary = imaginary * imaginary; // imaginary^2
			float total = real + imaginary; // real^2 + imaginary^2
			total = sqrt (total); // magnitude = sqrt (real^2 + imaginary^2)
			cout << total <<endl;
		}
	}
	
	if(printOption==5){
		//OUTPUT PHASE
		cout<<"\n\nPHASE PART\n";
		for(int i = 0; i < 128; i++){
			float real = complexCosine[2*i+1];
			float imaginary = complexCosine[2*(i+1)];
			float phase = atan2(imaginary, real);
			cout << phase <<endl;
		}
	}

	if(printOption==6){

		fft(complexCosine, 128, 1);

		cout<<"\n\n\nUsing Inverse DFT to verify forward DFT:\n";
		float originalCosine[128];

		for(int i=0; i < 128; i++){
			originalCosine[i]=0; // initialization
			originalCosine[i]=complexCosine[2*i+1]; // take the real values from complexCosine
			cout << originalCosine[i] << endl; // outputting the cosine funtion after taking the DFT and inverse DFT
		}

	}

}


// Example 1c: compute FFT of rectangular fucntion
void exp1c(){
	ifstream inRectFile("Rect_128.dat");
	float Rect_data[128];
	float n = 128.0;

	
	// Read Rect data from file
	int idx = 0;
	string line;
	while (getline(inRectFile, line)) {
		Rect_data[idx] = stof(line)*pow(-1, idx);
		idx++;
	}

	float complexRect_data[257];
	for (int i = 0; i < 257; i++) {
		complexRect_data[i] = 0;
	}

	//making the array with real and imaginary parts
	for(int i=0; i< 128; i++){
		complexRect_data[2*i+1]=Rect_data[i];
	}


	// Perform the DFT on Rect data
	fft(complexRect_data, 128, -1);

	for (int i = 0; i < 257; i++) {
		complexRect_data[i] = complexRect_data[i] / n; // divide by the normalization factor
	}

	int printOption;
	cout<<"\nSelect which values you want printed:";
	cout<< "\n1. DFT of the Cosine function\n2. Real part of the cosine function\n3. Imaginary Part of the cosine function\n4. Magnitude of the cosine funciton\n5. Phase of the cosine function\n\nPlease Enter a value: ";
	cin>>printOption;

	if(printOption==1){
		// Output results of DFT on Rect data
		cout << "\nDFT of the rectangular function (Rect_128.dat)" << endl;
		for (int i = 0; i < 257; i++) {
			cout << complexRect_data[i] << endl;

		}
	}
	
	if(printOption==2){
		//OUTPUT REAL
		cout<<"\n\nREAL PART\n";
		for(int i = 0; i < 128; i++){
			cout<<complexRect_data[2*i+1]<<endl;
		}
	}

	if(printOption==3){
		//OUTPUT IMAGINARY
		cout<<"\n\nIMAGINARY PART\n";
		for(int i = 0; i < 128; i++){
			cout<<complexRect_data[2*(i+1)]<<endl;
		}
	}
	
	if(printOption==4){
		//OUTPUT MAGNITUDE
		cout<<"\n\nMAGNITUDE PART\n";
		for(int i = 0; i < 128; i++){
			float real = complexRect_data[2*i+1];
			float imaginary = complexRect_data[2*(i+1)];
			real = real * real; // real^2
			imaginary = imaginary * imaginary; // imaginary^2
			float total = real + imaginary; // real^2 + imaginary^2
			total = sqrt (total); // magnitude = sqrt (real^2 + imaginary^2)
			cout << total <<endl;
		}
	}	
	
	if(printOption==5){
		//OUTPUT PHASE
		cout<<"\n\nPHASE PART\n";
		for(int i = 0; i < 128; i++){
			float real = complexRect_data[2*i+1];
			float imaginary = complexRect_data[2*(i+1)];
			float phase = atan2(imaginary, real);
			cout << phase <<endl;
		}
	}
	
}



// Compute 2D FFT given a black image with a square in the middle (square size varies)
void exp2(int squareSize){

	//generate 512x512 image
	int n = 512;
	int m = 512;
	int q = 255;

	char *originalImageWritePath = "./Images/test_original.pgm";
	char *writePath = "./Images/test.pgm";
	ImageType generatedImage(n, m, q);

	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			generatedImage.setPixelVal(i, j, 0); // setting the pixel value to black
		}
	}

	// calculate the location of the square so that it's in the middle of the image
	int lowerBound = (n - squareSize) / 2;
	int higherBound = ((n - squareSize) / 2 ) + squareSize;
	for(int i=lowerBound; i <= higherBound; i++){
		for(int j=lowerBound; j <= higherBound; j++){
			int temp;
			generatedImage.getPixelVal(i,j, temp);
			generatedImage.setPixelVal(i,j,255); // setting the pixel vlaue to white
		}
	}

	writeImage(originalImageWritePath, generatedImage); // write the original image with the white square

	// real and imaginary arrays
	float realArray[n*m];
	float imagArray[n*m];

	for(int i=0; i<n*m; i++){
		imagArray[i] = 0; // setting the imaginary array equal to 0
	}

	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			int temp;
			generatedImage.getPixelVal(i, j, temp); 
			realArray[i*m+j] = temp; // setting the values of the realArray from the original image
		}
	}
		
	writeImage(writePath, generatedImage);	

	
	fft2D(n, m, realArray, imagArray, -1);

	float total[n*m];
	for(int i=0; i<n*m; i++){
		total[i] = sqrt (realArray[i] * realArray[i] + imagArray[i] * imagArray[i]); // magnitude = sqrt (real^2 + imaginary^2)
	}

	// allowing user to shift magnitude to the center of the frequency domain
	int shiftMagnitude;
	cout << "\nWould you like to shift the magnitude to the center of the frequency domain?\n1. No\n2. Yes\n\nPlease enter your option: ";
	cin>> shiftMagnitude;
	if(shiftMagnitude==2){
		for (int i = 0; i < n*m; i++) {
			total[i] = total[i] * pow(-1, i); //  shift magnitude to center of frequency domain
		}
	}

	// allowing user to visualize the image with a log tranformation
	int logTransformation;
	cout << "\nWould you like to visualize the image with a log transformation?\n1. No\n2. Yes\n\nPlease enter your option: ";
	cin>> logTransformation;
	if(logTransformation==2){
		// Applying log transformation for visualization
		for(int i=0; i<n*m; i++){
			total[i] = log(1+total[i]);                       
		}
	}

	rescale(n, m, total);


	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){	
			int temp;
			generatedImage.getPixelVal(i, j, temp); 
			generatedImage.setPixelVal(i, j, total[i*m+j]);
		}
	}
	

	writeImage(writePath, generatedImage);	


}



// Compute the 2D FFT of Lenna
void exp3(int experimentPart){

	int n, m, q;
    bool type;

	char *readPath = "./Images/lenna.pgm";            // input image
    char *writePath = "./Images/lenna_dft.pgm";     // output image

    readImageHeader(readPath, n, m, q, type);
    ImageType originalImage(n, m, q);
    readImage(readPath, originalImage);


	float realArray[n*m];
	float imagArray[n*m];

	for(int i=0; i<n*m; i++){
		imagArray[i] = 0;
	}

	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			int temp;
			originalImage.getPixelVal(i, j, temp); 
			realArray[i*m+j] = temp;
		}
	}
	
	fft2D(n, m, realArray, imagArray, -1);

	float magnitude[n*m];
	float phase[n*m];

	// Part 3a 
	if(experimentPart==1){
		// calculate magnitude
		for(int i=0; i<n*m; i++){
			magnitude[i] = sqrt (realArray[i] * realArray[i] + imagArray[i] * imagArray[i]); // magnitude = sqrt (real^2 + imaginary^2)
		}
		// Set phase to 0 by setting real to magnitude and imaginary to 0
		for(int i=0; i<n*m; i++){
			realArray[i]=magnitude[i]; // set real to magnitude                           
		}
		for(int i=0; i<n*m; i++){
			imagArray[i]=0; // set imaginary to 0
		}
	}


	// Part 3b
	if(experimentPart==2){
		//set magnitude to 1
		for(int i=0; i<n*m; i++){
			magnitude[i] = 1; // magnitude = 1
		}
		//setting phase to original by setting the realArray equal to cos(theta) and imagArray equal to sin(theta)
		for(int i=0; i<m*n; i++){
			phase[i]=atan2(imagArray[i], realArray[i]); // calculate the phase
			realArray[i]=cos(phase[i]); // setting the realArray equal to cos(theta)
			imagArray[i]=sin(phase[i]); // setting the imaginary arrary equal to sin(theta)
		}
	}

	// Inverse DFT
	fft2D(n, m, realArray, imagArray, 1);

	// allowing user to decide whether to include log tranformation or not
	int logTranformationOption;
	cout << "\nWould you like to visualize the image with a log transformation?\n1. No\n2. Yes\n\nPlease enter your option: ";
	cin >> logTranformationOption;
	if(logTranformationOption==2){
		// Applying log transformation for visualization
		for(int i=0; i<n*m; i++){
			realArray[i] = log(1+realArray[i]);                       
		}
	}

	rescale(n, m, realArray);

	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){	
			int temp;
			originalImage.getPixelVal(i, j, temp); 
			originalImage.setPixelVal(i, j, realArray[i*m+j]);
		}
	}
	
	writeImage(writePath, originalImage);

}



int main()
{
	cout << "CS 474 Programming Assignment 1: Team 8\n\n";

	int option;
	cout << "Please select one of the following options:\n";
	cout << "1. Experiment 1a: DFT of the signal [2, 3, 4, 4]\n";
	cout << "2. Experiment 1b: DFT of cosine function with 128 samples\n";
	cout << "3. Experiment 1c: DFT of the rectangular function (Rect128.dat)\n\n";
	cout << "4. Experiment 2a: 2D DFT of image with 32x32 square inside\n";
	cout << "5. Experiment 2b: 2D DFT of image with 64x64 square inside\n";
	cout << "6. Experiment 2c: 2D DFT of image with 128x128 square inside\n\n";
	cout << "7. Experiment 3a: Visualizing magnitude of Lenna.pgm\n";
	cout << "8. Experiment 3b: Visualizing phase of Lenna.pgm\n\n";
	cout<< "9. Exit\n\n\n";
	cout << "Please enter an option: ";
	cin >> option;

	switch (option)
	{
	case 1:
		exp1a();
		break;
	case 2:
		exp1b();
		break;
	case 3:
		exp1c();
		break;
	case 4:
		exp2(32);
		break;
	case 5:
		exp2(64);
		break;
	case 6:
		exp2(128);
		break;
	case 7:
		exp3(1);
		break;
	case 8:
		exp3(2);
		break;

	}


	return 0;
}