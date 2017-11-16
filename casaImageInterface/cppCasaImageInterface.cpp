#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

// include casacore libraries.
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/ImageUtilities.h>
#include <casacore/lattices/Lattices/TiledShape.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/casa/Arrays/ArrayIO.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/coordinates/Coordinates/Coordinate.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/Projection.h>
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/tables/TaQL/TableParse.h>
#include <casacore/tables/Tables/ScalarColumn.h>

using namespace std;
using namespace casacore;
using namespace casa;

//
//	CasaImageInterface
//
//	CJS: 08/03/2017
//
//	Provides an interface for writing images to CASA or bitmap files. We can also load images from bitmap, provided
//	the input files are 256-colour grayscale.
//

template<class T>
class CasaImageInterface
{
	
	public:
		
		// default constructor
		CasaImageInterface();
   
		// destructor
		~CasaImageInterface();
	
		// load bitmap
		bool LoadBitmap( const char * pFilename, complex<T> ** pImageData );
		bool LoadBitmap( const char * pFilename, complex<T> ** pImageData, int * pWidth, int * pHeight );
		bool LoadBitmap( const char * pFilename, complex<T> ** pImageData, double * pImageScale );
		bool LoadBitmap( const char * pFilename, complex<T> ** pImageData, int * pWidth, int * pHeight,
					double * pImageScale );
						
		// save bitmap
		bool SaveBitmap( const char * pFilename, complex<T> * pImageData, int pWidth, int pHeight );
		bool SaveBitmap( const char * pFilename, complex<T> * pImageData, int pWidth, int pHeight,
					double * pImageScale );
		
		// write a casa image
		void WriteCasaImage( const char * pFilename, int pWidth, int pHeight, double pRA, double pDec,
					double pPixelSize, complex<T> * pImage, int pFrequency );
	
	private:

		static const double PI = 3.141592654;
		
		ImageUtilities * _imageUtilities = NULL;

		// bitmap file header positions.
		static const int BIT_CONST = 0x00;
		static const int MAP_CONST = 0x01;
		static const int IMAGE_SIZE = 0x02;
		static const int RESERVED = 0x06;
		static const int FILE_HEADER_SIZE = 0x0A;
		static const int BITMAP_INFO_HEADER = 0x0E;
		static const int IMAGE_WIDTH = 0x12;
		static const int IMAGE_HEIGHT = 0x16;
		static const int COLOUR_PLANES = 0x1A;
		static const int BIT_COUNT = 0x1C;
		static const int COMPRESSION_TYPE = 0x1E;
		static const int COLOURS_USED = 0x2E;
		static const int SIGNIFICANT_COLOURS = 0x32;
		
		// load the bitmap.
		bool loadBitmap( const char * pFilename, complex<T> ** pImageData, int * pWidth, int * pHeight,
					double * pImageScale );
						
		// save the bitmap.
		bool saveBitmap( const char * pFilename, complex<T> * pImageData, int pWidth, int pHeight,
					double * pImageScale );
		
		// write a casa image.
		void writeCasaImage( const char * pFilename, int pWidth, int pHeight, double pRA, double pDec,
					double pPixelSize, complex<T> * pImage, int pFrequency );
    
}; // CasaImageInterface

//
//	P U B L I C   C L A S S   M E M B E R S
//

//
//	CasaImageInterface<T>::CasaImageInterface()
//
//	CJS: 08/03/2017
//
//	The constructor.
//

template<class T>
CasaImageInterface<T>::CasaImageInterface()
{
	
	// create a new image utilities object.
	_imageUtilities = new ImageUtilities();
	
} // CasaImageInterface<T>::CasaImageInterface

//
//	CasaImageInterface<T>::~CasaImageInterface()
//
//	CJS: 08/03/2017
//
//	The destructor.
//

template<class T>
CasaImageInterface<T>::~CasaImageInterface()
{
	
	// delete the image utilities object.
	delete _imageUtilities;
	
} // CasaImageInterface<T>::~CasaImageInterface

//
//	CasaImageInterface<T>::LoadBitmap()
//
//	CJS: 08/03/2017
//
//	public interface to load bitmap files, with various overloads.
//

template<class T>
bool CasaImageInterface<T>::LoadBitmap( const char * pFilename, complex<T> ** pImageData )
{
	
	int width = 0, height = 0;
	double imageScale = 0;
	return loadBitmap( pFilename, pImageData, &width, &height, &imageScale );
	
} // CasaImageInterface<T>::LoadBitmap

template<class T>
bool CasaImageInterface<T>::LoadBitmap( const char * pFilename, complex<T> ** pImageData, int * pWidth,
					int * pHeight )
{
	
	double imageScale = 0;
	return loadBitmap( pFilename, pImageData, pWidth, pHeight, &imageScale );
	
} // CasaImageInterface<T>::LoadBitmap

template<class T>
bool CasaImageInterface<T>::LoadBitmap( const char * pFilename, complex<T> ** pImageData, double * pImageScale )
{
	
	int width = 0, height = 0;
	return loadBitmap( pFilename, pImageData, &width, &height, pImageScale );
	
} // CasaImageInterface<T>::LoadBitmap

template<class T>
bool CasaImageInterface<T>::LoadBitmap( const char * pFilename, complex<T> ** pImageData, int * pWidth, int * pHeight,
					double * pImageScale )
{
	
	return loadBitmap( pFilename, pImageData, pWidth, pHeight, pImageScale );
	
} // CasaImageInterface<T>::LoadBitmap

//
//	CasaImageInterface::SaveBitmap()
//
//	CJS: 13/03/2017
//
//	public interface to save bitmap files, with various overloads.
//

template<class T>
bool CasaImageInterface<T>::SaveBitmap( const char * pFilename, complex<T> * pImageData, int pWidth, int pHeight )
{
	
	double imageScale = 0;
	return saveBitmap( pFilename, pImageData, pWidth, pHeight, &imageScale );
	
} // CasaImageInterface<T>::SaveBitmap

template<class T>
bool CasaImageInterface<T>::SaveBitmap( const char * pFilename, complex<T> * pImageData, int pWidth, int pHeight,
					double * pImageScale )
{
	
	return saveBitmap( pFilename, pImageData, pWidth, pHeight, pImageScale );
	
} // CasaImageInterface<T>::SaveBitmap

//
//	CasaImageInterface<T>::WriteCasaImage()
//
//	CJS: 08/03/2017
//
//	public interface to writing casa images.
//

template<class T>
void CasaImageInterface<T>::WriteCasaImage( const char * pFilename, int pWidth, int pHeight, double pRA, double pDec,
						double pPixelSize, complex<T> * pImage, int pFrequency )
{
	
	writeCasaImage( pFilename, pWidth, pHeight, pRA, pDec, pPixelSize, pImage, pFrequency );
				
} // CasaImageInterface<T>::WriteCasaImage

//
//	P R I V A T E   C L A S S   M E M B E R S
//

//
//	CasaImageInterface<T>::loadBitmap()
//
//	CJS: 07/07/2015
//
//	load a bitmap file, and return the image size and a boolean indicating success.
//	the image must be 8-bit greyscale.
//

template<class T>
bool CasaImageInterface<T>::loadBitmap( const char * pFilename, complex<T> ** pImageData, int * pWidth, int * pHeight,
					double * pImageScale )
{
	
	bool ok = true;
	unsigned char * fileInfo = NULL, * fileHeader = NULL, * tmpData = NULL;
	
	// open the bitmap file.
	FILE * inputFile = fopen( pFilename, "r" );
	if (inputFile == NULL)
	{
		printf("Could not open file \"%s\".\n", pFilename);
		ok = false;
	}
	else
	{
		
		// reserve memory for the start of the file header, and read it from the file. we only
		// read the first 18 bytes, because these contain information about how large the header is. once we
		// know this we can read the rest of the header.
		fileInfo = (unsigned char *) malloc( 18 );
		size_t num_read = fread( fileInfo, sizeof( unsigned char ), 18, inputFile );
				
		// ensure we've read the correct number of bytes.
		if (num_read != 18)
		{
			printf( "Error: read only %lu values from the file header.\n", num_read );
			ok = false;
		}

		// make sure this is a bitmap file by checking that the first two bytes are ASCII codes 'B' and 'M'.
		if (ok == true)
			if ((fileInfo[BIT_CONST] != 'B') || (fileInfo[MAP_CONST] != 'M'))
			{
				printf( "Error: this is not a bitmap file.\n" );
				ok = false;
			}
			
		// get the size of the file header (i.e. a pointer to the start of the actual image).
		int fileHeaderSize = 0;
		if (ok == true)
			memcpy( &fileHeaderSize, &fileInfo[FILE_HEADER_SIZE], 4 );
			
		// get the size of the bitmap info header (the bitmap info header is followed by the colour table,
		// so we need to know the offset in order to read the colours).
		int bitmapInfoHeaderSize = 0;
		if (ok == true)
			memcpy( &bitmapInfoHeaderSize, &fileInfo[BITMAP_INFO_HEADER], 4 );
		
		// need to add 14 because the bitmap info header size does not include the first 14 bytes of the file (which
		// technically are part of the file header but not the bitmap header; we lump everything in together so that
		// all of our offsets are from the start of the file - less confusing this way).
		bitmapInfoHeaderSize = bitmapInfoHeaderSize + 14;
			
		// get the rest of the file header now we know how big it is. we already have the first 18 bytes,
		// which should be copied to the start of the new memory area.
		if (ok == true)
		{
			fileHeader = (unsigned char *) malloc( fileHeaderSize );
			memcpy( fileHeader, fileInfo, 18 );
			num_read = fread( &fileHeader[18], sizeof( unsigned char ), fileHeaderSize - 18, inputFile );
			if (num_read != (fileHeaderSize - 18))
			{
				printf( "Error: read only %lu values from the file header.\n", num_read + 18 );
				ok = false;
			}
		}
		
		// get the input image flux scale. this value may be stored in the reserved part of the bitmap file header
		// (0x06 -> 0x09), and will not be supplied if the input image has been saved using something like GIMP or
		// Photoshop. if it is zero, then we assume a scale of 1 Jy/PIXEL. this value gets re-scaled along with our
		// image, and is then written back to the output file.
		if (ok == true)
			memcpy( pImageScale, &fileHeader[RESERVED], 4 );
		
		if (*pImageScale == 0)
			*pImageScale = 1;
			
		// ensure we have an 8-bit image.
		if (ok == true)
		{
			short bitCount;
			memcpy( &bitCount, &fileHeader[BIT_COUNT], 2 );
			if (bitCount != 8)
			{
				printf( "Error: expecting an 8-bit greyscale image. This one is %hi bit.\n", bitCount );
				ok = false;
			}
		}
			
		// ensure the image in not compressed.
		if (ok == true)
		{
			int compressionMethod;
			memcpy( &compressionMethod, &fileHeader[COMPRESSION_TYPE], 4 );
			if (compressionMethod != 0)
			{
				printf( "Error: can only handle uncompressed bitmaps." );
				ok = false;
			}
		}
			
		if (ok == true)
		{
			
			// get the width and height of the image.
			memcpy( pWidth, &fileHeader[IMAGE_WIDTH], 4 );
			memcpy( pHeight, &fileHeader[IMAGE_HEIGHT], 4 );
		
			// ensure width and height are greater than zero.
			if (*pWidth <= 0 || *pHeight <= 0)
			{
				printf( "Error: invalid image size (%i x %i).\n", *pWidth, *pHeight );
				ok = false;
			}
			
		}
		
		if (ok == true)
		{
			
			// ensure the number of colours used is 256.
			int coloursUsed = 0;
			memcpy( &coloursUsed, &fileHeader[COLOURS_USED], 4 );
			if (coloursUsed != 256)
			{
				printf( "ERROR: Can only handle 256 colours in pallette.\n" );
				ok = false;
			}
			
		}
		
		// get the number of significant colours used. this value can (theoretically) be less than COLOURS_USED
		// if an image is only using (e.g.) 37 shades rather than all 256. in practice, this is never implemented, and
		// this value will either be 0 (= all colours) or will match COLOURS_USED. however, only SIGNIFICANT_COLOURS are
		// written to the pallette, so we have to handle this parameter just in case.
		int significantColours = 0;
		if (ok == true)
			memcpy( &significantColours, &fileHeader[SIGNIFICANT_COLOURS], 4 );
		
		// if significant colours = 0, then they are ALL significant so set to 256.
		if (significantColours == 0)
			significantColours = 256;
			
		unsigned int colour[256];
		if (ok == true)
		{
				
			// load colour table from bmp.
			for ( unsigned int i = 0; i < significantColours; ++i )
			{
				
				memcpy( &colour[i], &fileHeader[bitmapInfoHeaderSize + (i * 4)], 4 );
				
				// convert pallette colour to greyscale, using 0.2990, 0.5870, 0.1140 RGB weighting. add 0.5
				// to round to nearest integer (since C only rounds down).
				unsigned char red = colour[i] >> 16;
				unsigned char green = (colour[i] >> 8) - (red << 8);
				unsigned char blue = colour[i] - (red << 16) - (green << 8);
				colour[i] = (unsigned int) ((((double)red * 0.2990) + ((double)green * 0.5870) +
								((double)blue * 0.1140)) + 0.5);
				
			}
				
			// reserve some memory for the image, and read it from the file.
			tmpData = (unsigned char *) malloc( *pWidth * *pHeight );
			*pImageData = (complex<T> *) malloc( *pWidth * *pHeight * sizeof( complex<T> ) );
			num_read = fread( tmpData, sizeof( unsigned char ), *pWidth * *pHeight, inputFile );
				
			// ensure we've read the correct number of bytes.
			if (num_read != *pWidth * *pHeight)
			{
				printf( "Error: read only %lu values from the image.\n", num_read );
				ok = false;
			}
				
		}
			
		if (ok == true)
		{
				
			// update image values using the values from the colour table.
			complex<T> * complexData = *pImageData;
			for ( int i = 0; i < *pWidth * *pHeight; i++ )
				complexData[i] = complex<T>( (double)colour[tmpData[i]] * *pImageScale, 0 );
				
		}
		
		// close file.
		fclose( inputFile );
	
	}
	
	// tidy up memory.
	if (fileInfo != NULL)
		free( (void *) fileInfo );
	if (fileHeader != NULL)
		free( (void *) fileHeader );
	if (tmpData != NULL)
		free( (void *) tmpData );
	
	// return success flag.
	return ok;
	
} // CasaImageInterface<T>::loadBitmap

//
//	CasaImageInterface<T>::saveBitmap()
//
//	CJS: 10/08/2015
//
//	Save a bitmap file.
//

template<class T>
bool CasaImageInterface<T>::saveBitmap( const char * pFilename, complex<T> * pImageData, int pWidth, int pHeight,
					double * pImageScale )
{
	
	unsigned char * image = NULL;
	
	const int HEADER_SIZE = 1078;
	
	// allocate and build the header.
	unsigned char * fileHeader = (unsigned char *) malloc( HEADER_SIZE );
	memset( fileHeader, 0, HEADER_SIZE );

	// file header.
	fileHeader[BIT_CONST] = 'B'; fileHeader[MAP_CONST] = 'M';					// bfType
	int size = (pWidth * pHeight) + HEADER_SIZE; memcpy( &fileHeader[IMAGE_SIZE], &size, 4 );	// bfSize
	int offBits = HEADER_SIZE; memcpy( &fileHeader[FILE_HEADER_SIZE], &offBits, 4 );		// bfOffBits

	// image header.
	size = 40; memcpy( &fileHeader[BITMAP_INFO_HEADER], &size, 4 );					// biSize
	memcpy( &fileHeader[IMAGE_WIDTH], &pWidth, 4 );							// biWidth
	memcpy( &fileHeader[IMAGE_HEIGHT], &pHeight, 4 );						// biHeight
	short planes = 1; memcpy( &fileHeader[COLOUR_PLANES], &planes, 2 );				// biPlanes
	short bitCount = 8; memcpy( &fileHeader[BIT_COUNT], &bitCount, 2 );				// biBitCount
	int coloursUsed = 256; memcpy( &fileHeader[COLOURS_USED], &coloursUsed, 4 );			// biClrUsed

	// colour table.
	for (unsigned int i = 0; i < 256; ++i)
	{
		unsigned int colour = (i << 16) + (i << 8) + i;
		memcpy( &fileHeader[54 + (i * 4)], &colour, 4 );
	}
	
	bool ok = true;

	// open file.
	FILE * outputFile = fopen( pFilename, "w" );
	if (outputFile == NULL)
	{
		printf( "Could not open file \"%s\".\n", pFilename );
		ok = false;
	}
	else
	{

		// write the file header.
		size_t num_written = fwrite( fileHeader, 1, 1078, outputFile );
		if (num_written != 1078)
		{
			printf( "Error: cannot write to file.\n" );
			ok = false;
		}
		
		// find the maximum and minimum pixel values.
		T min = real( pImageData[ 0 ] );
		T max = real( pImageData[ 0 ] );
		for ( int i = 1; i < pWidth * pHeight; i++ )
		{
			if (real( pImageData[ i ] ) < min)
				min = real( pImageData[ i ] );
			if (real( pImageData[ i ] ) > max)
				max = real( pImageData[ i ] );
		}
		
		printf("min: %f, max: %f\n", min, max );

		// add 1% allowance to max - we don't want saturation.
		max = ((max - min) * 1.1) + min;
		
		// construct the image.
		image = (unsigned char *) malloc( pWidth * pHeight * sizeof( unsigned char ) );
		for ( int i = 0; i < pWidth * pHeight; i++ )
			image[i] = (unsigned char)( (real( pImageData[ i ] ) - min) * ((double)256 / (max - min)) );
		
		// write the data.
		if (ok == true)
		{
			
			size_t num_written = fwrite( image, 1, pWidth * pHeight, outputFile );
			if (num_written != (pWidth * pHeight))
			{
				printf( "Error: cannot write to file.\n" );
				ok = false;
			}
			
		}

		// close file.
		fclose( outputFile );
		
	}

	// cleanup memory.
	free( (void *) fileHeader );
	if (image != NULL)
		free( image );
	
	// return success flag.
	return ok;
	
} // CasaImageInterface<T>::saveBitmap

//
//	CasaImageInterface<T>::writeCasaImage()
//
//	CJS: 08/03/2017
//
//	write an image to a casa image.
//
		
template<class T>
void CasaImageInterface<T>::writeCasaImage( const char * pFilename, int pWidth, int pHeight, double pRA,
						double pDec, double pPixelSize, complex<T> * pImage,
						int pFrequency )
{
	
	// image size.
	IPosition imageSize = IPosition( 3, pWidth, pHeight, 1 );
	
	// create a tiled shape.
	TiledShape outShape( imageSize );
	
	// the transformation matrix is diagonal so that the image is aligned with ra and dec (casa doesn't
	// like rotated coordinate systems).
	Matrix<Double> dirTransform( 2, 2, 2 );
	dirTransform( 0, 0 ) = 1;
	dirTransform( 0, 1 ) = 0;
	dirTransform( 1, 0 ) = 0;
	dirTransform( 1, 1 ) = 1;
	
	// create a new direction coordinate.
	DirectionCoordinate outDirectionCoordinate(	MDirection::J2000, 
							Projection( Projection::SIN ),
							pRA * PI / (double)180.0, pDec * PI / (double)180.0,
							-pPixelSize * PI / ((double)180.0 * (double)3600.0),
							pPixelSize * PI / ((double)180.0 * (double)3600.0),
							dirTransform,
							pWidth / 2, pHeight / 2 );
	
	// to change the axis units to degrees:
	Vector<String> outUnits( 2 ); outUnits = "deg";
	outDirectionCoordinate.setWorldAxisUnits( outUnits );
							
	// create a new spectral coordinate.
	SpectralCoordinate outSpectralCoordinate( MFrequency::REST, pFrequency, 1, 0, pFrequency );
	
	// create a new coordinate system.
	CoordinateSystem outCoordinateSystem;
	
	// add the direction and spectral coordinates.
	outCoordinateSystem.addCoordinate( outDirectionCoordinate );
	outCoordinateSystem.addCoordinate( outSpectralCoordinate );
	
	printf( "Num Pixel Axes: %i\n", outCoordinateSystem.nPixelAxes() );
	int coordinate = -1; int& coordinateRef = coordinate;
	int axis = -1; int& axisRef = axis;
	outCoordinateSystem.findPixelAxis( coordinateRef, axisRef, 0 );
	printf( "ra: %i, %i, 0\n", coordinate, axis );
	outCoordinateSystem.findPixelAxis( coordinateRef, axisRef, 1 );
	printf( "ra: %i, %i, 1\n", coordinate, axis );
	outCoordinateSystem.findPixelAxis( coordinateRef, axisRef, 2 );
	printf( "ra: %i, %i, 2\n", coordinate, axis );
	
	// create memory for the image pixels and the mask.
	Array<Float> outPixels( imageSize );
	Array<Bool> outMask( imageSize );
	
	// copy the image into the arrays.
	for ( int i = 0; i < pWidth; i++ )
		for ( int j = 0; j < pHeight; j++ )
		{
			
			IPosition arrayPos( 3, i, j, 0 );
			outPixels( arrayPos ) = (float)real( pImage[ (j * pWidth) + i ] );
			outMask( arrayPos ) = true;
			
		}
	
	// create an IO object?
	LogIO outIO;

	// write the image to a casa image directory (AIPSPP).
	_imageUtilities->writeImage(	outShape,
					outCoordinateSystem,
					pFilename,
					outPixels,
					outIO,
					outMask );
						
} // CasaImageInterface<T>::writeCasaImage
