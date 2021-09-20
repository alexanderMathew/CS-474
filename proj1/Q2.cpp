// CS 474 Project 1, Question 2
// Team #8: Alexander Mathew and Tal Shahar Zemach

#include <iostream>
#include <fstream>
#include <sstream>
#include "image.h"
using namespace std;


// Function declarations
int readImage(char fname[], ImageType& image);
int writeImage(char fname[], ImageType& image);
int readImageHeader(char fname[], int& N, int& M, int& Q, bool& type);

void quantizeImage(int qLevel, int n, int m, ImageType originalImage/*, ImageType& newImage*/){
    int val;
    originalImage.setImageInfo(n, m, qLevel); //Update the image so the quantum level matches Q ********
    qLevel = 256 / qLevel;

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            originalImage.getPixelVal(i, j, val);
            val = val / qLevel;
            originalImage.setPixelVal(i, j, val);
        }
    }

}


int main(){
    int q_level;
    int m, n, q;
	bool type;
	int ss_factor = 2;
    char* readPath = "./Images/lenna.pgm";
	char *writePath = "./Images/lenna_quantized.pgm"; 
	

	readImageHeader(readPath, n, m, q, type);
	ImageType ImageToQuantize(n, m, q);
	ImageToQuantize.getImageInfo(n, m, q);

	int x = m;
	int y = q;
	int z = n;

	ImageType resizedImage(z, x, y);

	readImage(readPath, ImageToQuantize);
    cout << "Enter quantization level: ";
    cin >> q_level;


    quantizeImage(q_level, z, x, ImageToQuantize);
    writeImage(writePath, resizedImage);


    return 0;
}