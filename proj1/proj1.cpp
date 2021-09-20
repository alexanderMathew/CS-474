#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
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

			// RESIZE ATTEMPT
			/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for(int x = 0; x < orig_cols; x++){
				for(int y = 0; y < orig_rows; y++){
					int temp;
					if((x % factor) != 0 || (y % factor) != 0){
						newImage.resizePixel(x, y, pixelVal); //newImage[x][y] = pixelVal;
						//pixelVal[x][y] = temp;
						
					}
					else{
						newImage.assignPixelVal(x, y, pixelVal); //assign pixel values
						//temp = pixelVal[x][y];

					}
						newImage.printPixelVal(x, y); //print pixel values
				}				
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		sampled_rows++;

		}
	sampled_cols++;
	}

}


void quantizeImage(int qLevel, int n, int m, ImageType originalImage, ImageType newImage){
    int val;
    //originalImage.getImageInfo(n, m, qLevel); //Update the image so the quantum level matches Q ********
    qLevel = 256 / qLevel;

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            originalImage.getPixelVal(i, j, val);
			//originalImage.printPixelVal(i, j);
            val /= qLevel;
            newImage.setPixelVal(i, j, val);
			//newImage.printPixelVal(i, j);
        }
    }

}

/*void rescale(int n, int m, int q, ImageType OriginalImage){
	OriginalImage.getImageInfo(n, m, q);

	int max = 0;
	//OriginalImage.setPixelVal(0, 0, max); //initialization - point 0,0 is our original max
	int min = 0;
	//OriginalImage.setPixelVal(0, 0, min); //initialization - point 0,0 is our original min
	
	/////////////PRINT IMAGE////////////////////
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			//OriginalImage.printPixelVal(i,j);
		}
	}	
	///////////////////////////////////////////

	int temp;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			OriginalImage.getPixelVal(i,j,temp);
			if (temp > max){
				//OriginalImage.setPixelVal(i, j, max); //max=OriginalImage[i][j];
			}
			if(temp < min){
				//OriginalImage.setPixelVal(i, j, min); //min=OriginalImage[i][j];
			}
			cout<<"MAX: "<<max<<"\tMIN: "<<min<<endl;
		}
	}
}*/

// QUESTION 1: IMAGE SAMPLING
void Q1() {
	int m, n, q;
	bool type;
	int ss_factor; // subsampling factor

	cout << "Enter subsampling factor: ";
	cin >> ss_factor;

	char* readPath = "./Images/lenna.pgm"; // input image
	char *writePath = "./Images/lenna_sampled.pgm"; // output image

	readImageHeader(readPath, n, m, q, type);
	ImageType ImageToSample(n, m, q); // input image created as an ImageType object
	ImageToSample.getImageInfo(n, m, q); // gets image info of ImageToSample

	int x = m;
	int y = q;
	int z = n;

	ImageType ssImage(z, x, y); // the subsampled image

	readImage(readPath, ImageToSample); //read in the ImageToSample

	// sample image
	subSampleImage(ImageToSample, ssImage, ss_factor);

	// write image
	writeImage(writePath, ssImage);

	cout << "File written to: " << writePath <<endl;

}

// QUESTION 2: IMAGE QUANTIZATION
void Q2(){
    int q_level; //our desired q-level
    int m, n, q; //n, m, and q for the input image
	bool type;
    char* readPath = "./Images/lenna.pgm";
	char *writePath = "./Images/lenna_quantized.pgm"; 
	

	readImageHeader(readPath, n, m, q, type);
	ImageType ImageToQuantize(n, m, q);
	ImageToQuantize.getImageInfo(n, m, q);

	int x = n;
	int y = m;
	int z = q;

	readImage(readPath, ImageToQuantize);


    cout << "Enter quantization level: ";
    cin >> q_level;

	ImageType resizedImage(x, y, q_level);
	
	/*for(int i = 0; i < z; i++){
		for(int j = 0; j < x; j++){
			resizedImage.printPixelVal(i, j);
		}
	}*/

	//PROBLEM WITHH FUNCTION
    quantizeImage(q_level, x, y, ImageToQuantize, resizedImage);
	for(int i = 0; i < y; i++){
		for(int j = 0; j < x; j++){
			resizedImage.printPixelVal(i, j);
		}
	}
    writeImage(writePath, resizedImage);

	cout << "File written to: " << writePath <<endl;

}

// QUESTION 3: HISTOGRAM EQUALIZATION
//////////////////////////////////////////////////////////////////////////////////////////CHANGE COMMENTS!!!!!!!!!!!!
void Q3()
{
	int m, n, q; //n, m, and q for the input image
	bool type;
	char* readPath = "./Images/lenna.pgm";
	char *writePath = "./Images/lenna_hist_eq.pgm"; 

	readImageHeader(readPath, n, m, q, type);

    ImageType originalImage(n, m, q);
	originalImage.getImageInfo(n, m, q);
    readImage(readPath, originalImage);

    //Declaration for storing the two histograms
    vector<int> originalHist(256, 0);
    vector<int> equalizedHist(256, 0);

    //Calculate Histogram
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            int temp;
            originalImage.getPixelVal(i, j, temp);
            originalHist[temp]++;
        }
    }

    //Calculate total amount of pixels
    int total = m * n;
    int frequency = 0;

    // calculating cumulative frequency and new gray levels
    for (int i = 0; i < 256; i++) {
        // calculating cumulative frequency
        frequency = frequency + originalHist[i];
  
		// calculating equalized histogram's gray level
        equalizedHist[i] = round((((float)frequency) * 255) / total);
    }

    ImageType equalizedImage(n, m, q);

    // performing histogram equalisation by mapping new gray levels
    for (int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++){
            // mapping to new gray level values
            int temp, temp2;
            originalImage.getPixelVal(i, j, temp);
            equalizedImage.setPixelVal(i, j, equalizedHist[temp]);
        }
    }

    //Output equalized image
    writeImage(writePath, equalizedImage);

	cout << "File written to: " << writePath <<endl;


}

int main(){
	cout<<"CS 474 Programming Assignment 1: Team 8\n\n";

	int option;
	cout<<"Please select one of the following options:\n";
	cout<<"1. Image Sampling\n2. Image Quantization\n3. Histogram Equalization\n4. Histogram Specification\n5. Exit\n\n";
	cout<<"Please enter an option: ";
	cin>>option;

	//while (option !=5){
		switch (option){
			case 1:
				Q1();
				break;
			case 2:
				Q2();
				break;
			case 3:
				Q3();
				break;
		}
	//}

	return 0;
}

