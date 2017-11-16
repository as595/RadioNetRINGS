#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

// my casacore interface.
#include "./cppCasaImageInterface.cpp"

//
//	main()
//
//	CJS: 15/11/2017
//
//	Demonstrates how to use the CASA Image Interface. This code just loads an image of the Lovell telescope from a Bitmap file, and saves it as a CASA image.
//

int main( int pArgc, char ** pArgv )
{

	const double RA = 140.0;
	const double DEC = 60.0;
	const double PIXEL_SIZE = 0.000002778; // in degrees
	const double FREQUENCY = 1400000000.0; // in Hz

	// store the image width and height.
	int imageWidth = 0, imageHeight = 0;

	// create a memory area to store the image.
	complex<double> * myImage = NULL;
	
	// create a casa image interface.
	CasaImageInterface<double> casaImageInterface;

	// filename of a bitmap image to load. the bitmap needs to be 256-colour greyscale.
	const char inputFilename[] = "./lovell.bmp";

	// load an image from a bitmap file. the bitmap needs to be 256-colour greyscale.
	casaImageInterface.LoadBitmap( inputFilename, &myImage, &imageWidth, &imageHeight );

	// filename of the CASA image.
	const char outputFilename[] = "./lovell.image";

	// now save this image as a CASA image.
	casaImageInterface.WriteCasaImage( outputFilename, imageWidth, imageHeight, RA, DEC, PIXEL_SIZE, myImage, FREQUENCY );

} // main
