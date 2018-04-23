// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
//
// This header contains constant definitions and auxiliary functions used in many places in the lightfield compressor.

#ifndef GEN_TYPES_HH
#define GEN_TYPES_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sys/stat.h>


#define NBIT_GR 32

const bool verbose = false;

// number of pixels in the original lenslet image.
const int num_pixels_in_lenslet = 41483904;

int getMedian(std::vector<int> scores)
{
	int median;
	size_t size = scores.size();

	std::sort(scores.begin(), scores.end());

	if (size % 2 == 0)
	{
		median = floor( ( (float)scores[size / 2 - 1] + (float)scores[size / 2]) / 2 );
	}
	else
	{
		median = scores[size / 2];
	}

	return median;
}

unsigned short getMedian(std::vector<unsigned short> scores)
{
	unsigned short median;
	size_t size = scores.size();

	std::sort(scores.begin(), scores.end());

	if (size % 2 == 0)
	{
		median = floor(((float)scores[size / 2 - 1] + (float)scores[size / 2]) / 2);
	}
	else
	{
		median = scores[size / 2];
	}

	return median;
}


float getMedian(std::vector<float> scores)
{
	float median;
	size_t size = scores.size();

	std::sort(scores.begin(), scores.end());

	if (size % 2 == 0)
	{
		median = (((float)scores[size / 2 - 1] + (float)scores[size / 2]) / 2);
	}
	else
	{
		median = scores[size / 2];
	}

	return median;
}

void medfilt2D(int* input, int* output, int SZ, int nr, int nc)
{
	int dsz = floor(SZ / 2);
	std::vector<int> scores;
	for (int y = 0; y < nr; y++){
		for (int x = 0; x < nc; x++){
			scores.clear();
			for (int dy = -dsz; dy < dsz; dy++){
				for (int dx = -dsz; dx < dsz; dx++){
					if ((y + dy) >= 0 && (y + dy) < nr
						&& (x + dx) >= 0 && (x + dx) < nc)
						scores.push_back(input[y + dy + (x + dx)*nr]);
				}
			}
			output[y + x*nr] = getMedian(scores);
		}
	}
}

void medfilt2D(unsigned short* input, unsigned short* output, int SZ, int nr, int nc)
{
	int dsz = floor(SZ / 2);
	std::vector<unsigned short> scores;
	for (int y = 0; y < nr; y++){
		for (int x = 0; x < nc; x++){
			scores.clear();
			for (int dy = -dsz; dy < dsz; dy++){
				for (int dx = -dsz; dx < dsz; dx++){
					if ((y + dy) >= 0 && (y + dy) < nr
						&& (x + dx) >= 0 && (x + dx) < nc)
						scores.push_back(input[y + dy + (x + dx)*nr]);
				}
			}
			output[y + x*nr] = getMedian(scores);
		}
	}
}


void medfilt2D(float* input, float* output, int SZ, int nr, int nc)
{
	int dsz = floor(SZ / 2);
	std::vector<float> scores;
	for (int y = 0; y < nr; y++){
		for (int x = 0; x < nc; x++){
			scores.clear();
			for (int dy = -dsz; dy < dsz; dy++){
				for (int dx = -dsz; dx < dsz; dx++){
					if ((y + dy) >= 0 && (y + dy) < nr
						&& (x + dx) >= 0 && (x + dx) < nc)
						scores.push_back(input[y + dy + (x + dx)*nr]);
				}
			}
			output[y + x*nr] = getMedian(scores);
		}
	}
}

// write a metadata file about sizes of images
inline void aux_write_header_file(int nr, int nc, int nvr, int nvc, const char* filename) {
	FILE* f_main_header = fopen(filename, "wb");
	fwrite(&nr, sizeof(int), 1, f_main_header);
	fwrite(&nc, sizeof(int), 1, f_main_header);
	fwrite(&nvr, sizeof(int), 1, f_main_header);
	fwrite(&nvc, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

inline void aux_read_header_file(int* nr, int* nc, int* nvr, int* nvc, const char* filename) {
	FILE* f_main_header = fopen(filename, "rb");
	fread(nr, sizeof(int), 1, f_main_header);
	fread(nc, sizeof(int), 1, f_main_header);
	fread(nvr, sizeof(int), 1, f_main_header);
	fread(nvc, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

inline void aux_read_file_float(const int nr, const int nc, const int ncomponents, const char* filename, float *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(float), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_read_file_uint32(const int nr, const int nc, const int ncomponents, const char* filename, unsigned int *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(unsigned int), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_read_file_int32(const int nr, const int nc, const int ncomponents, const char* filename, int *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(int), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_read_file_uint8(const int nr, const int nc, const int ncomponents, const char* filename, unsigned char *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(unsigned char), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_read_file_uint16(const int nr, const int nc, const int ncomponents, const char* filename, unsigned short *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(unsigned short), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_readDSP_uint16(const int nr, const int nc, const int ncomponents, const char* filename, unsigned short *data) {

	FILE* f_file = fopen(filename, "rb");

	//fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(unsigned short), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

// write a metadata file about sizes of images
inline void aux_write_header(int nr, int nc, int nvr, int nvc) {
	FILE* f_main_header = fopen("HDR", "wb");
	fwrite(&nr, sizeof(int), 1, f_main_header);
	fwrite(&nc, sizeof(int), 1, f_main_header);
	fwrite(&nvr, sizeof(int), 1, f_main_header);
	fwrite(&nvc, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

// read the metadata file
inline void aux_read_header(int* nr, int* nc, int* nvr, int* nvc) {
	FILE* f_main_header = fopen("HDR", "rb");
	fread(nr, sizeof(int), 1, f_main_header);
	fread(nc, sizeof(int), 1, f_main_header);
	fread(nvr, sizeof(int), 1, f_main_header);
	fread(nvc, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

// write a header about prediction order and maximum region number
inline void aux_write_pred_header(int Ms, int maxiS) {
	FILE* f_main_header = fopen("HDRP", "wb");
	fwrite(&Ms, sizeof(int), 1, f_main_header);
	fwrite(&maxiS, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

// read the header about prediction order and maximum region number
inline void aux_read_pred_header(int* Ms, int* maxiS) {
	FILE* f_main_header = fopen("HDRP", "rb");
	fread(Ms, sizeof(int), 1, f_main_header);
	fread(maxiS, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}


// read an 8 bit ppm file (here only one channel)
inline void aux_read8ppm(FILE *filept, int width, int height, int *img){
	int  max, x, y;
	int red;
	char dummy[100];

	unsigned char *Image8bit = NULL;

	/*--< Read header information of 16bit ppm image from filept >--*/
	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);
	printf("%s %d %d %d\n", dummy, max, width, height);


	Image8bit = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image8bit, sizeof(unsigned char), width*height * 3, filept);

	int i = 0;
	/*--< Find maximum value in the image >--*/
	for (x = 0; x < width; x++){
		for (y = 0; y < height; y++){
			red = (int)Image8bit[(x + y*width) * 3];
			img[i] = red;
			i++;
		}
	}

	free(Image8bit);

}

// read a 16 bit color pgm image
void aux_read16pgm(FILE *filept, int *img)
{
	int  max, x, y;
	int red, green, blue, pixelmax;
	char dummy[100];

	bool divd = 0;

	unsigned short int *Image16bit = NULL;

	int width, height;

	//FILE* filept = fopen("007_007.ppm", "r");
	/*--< Read header information of 16bit ppm image from filept >--*/
	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);
	//printf("%s\n", dummy);

	/*--< Really 16bit ppm? Check it >--*/
	if (strncmp(dummy, "P5", 2) != 0) {
		fprintf(stderr, "Error: The input data is not PGM.\n");
		exit(1);
	}

	Image16bit = (unsigned short int *)malloc(width*height* sizeof(unsigned short int));

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image16bit, sizeof(unsigned short int), width*height, filept);

	int i = 0;
	/*--< Find maximum value in the image >--*/
	pixelmax = 0;
	/* UNSW inverse depth already in 4K resolution, just crop to 1080p */
	for (x = 0; x < width; x++){
		for (y = 0; y < height; y++){
			red = Image16bit[(x + y*width)];

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);

			if (pixelmax < red) pixelmax = red;


			img[i] = red;


			if (0){ //debug stuff
				fprintf(stderr, "sample %i\n", img[i]);
			}

			i++;


		}
	}
	free(Image16bit);
}


// read a 16 bit color pgm image
inline void aux_read16pgm_1080p(FILE *filept, int *img)
{
	int  max, x, y;
	int red, green, blue, pixelmax;
	char dummy[100];

	bool divd = 0;

	unsigned short int *Image16bit = NULL;

	int width, height;

	//FILE* filept = fopen("007_007.ppm", "r");
	/*--< Read header information of 16bit ppm image from filept >--*/
	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);
	//printf("%s\n", dummy);

	/*--< Really 16bit ppm? Check it >--*/
	if (strncmp(dummy, "P5", 2) != 0) {
		fprintf(stderr, "Error: The input data is not PGM.\n");
		exit(1);
	}
	//if (max == 65535){
	//	divd = 1;
	//}

	Image16bit = (unsigned short int *)malloc(width*height* sizeof(unsigned short int));

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image16bit, sizeof(unsigned short int), width*height, filept);

	int i = 0;
	/*--< Find maximum value in the image >--*/
	pixelmax = 0;
	/* UNSW inverse depth already in 4K resolution, just crop to 1080p */
	for (x = 960; x < 960 + 1920; x++){
		for (y = 540; y < 540 + 1080; y++){
			//if (y>209 + 540 && x > 85 + 960 && y < 209 + 540 + 1080 && x < 85 + 960 + 1920)
				{
					red = Image16bit[(x + y*width)];
					//green = Image16bit[(x + y*width) * 3 + 1];
					//blue = Image16bit[(x + y*width) * 3 + 2];

					// Exhange upper 8bit and lower 8bit for Intel x86
					red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
					//green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
					//blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);

					if (pixelmax < red) pixelmax = red;
					//if (pixelmax < green) pixelmax = green;
					//if (pixelmax < blue) pixelmax = blue;


					//if (divd) // fix for 16bit to 10bit
					//{
					//	red = red >> 6;
					//	//green = green >> 6;
					//	//blue = blue >> 6;
					//}

					img[i] = red;
					//img[i + height*width] = green;
					//img[i + 2 * height*width] = blue;

					if (0){ //debug stuff
						fprintf(stderr, "sample %i\n", img[i]);
						/*fprintf(stderr, "sample %i\t%i\t%i\n", img[i] >> 6, img[i + height*width] >> 6, img[i + 2 * height*width] >> 6);*/
					}

					i++;

				}
		}
	}
	free(Image16bit);
}

inline void aux_read16ppm(FILE *filept, int width, int height, unsigned short *img)
{
	int  max, x, y;
	int red, green, blue, pixelmax;
	char dummy[100];

	unsigned short int *Image16bit = NULL;

	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);

	if (strncmp(dummy, "P6", 2) != 0) {
		fprintf(stderr, "Error: The input data is not binary PPM.\n");
		exit(1);
	}

	Image16bit = (unsigned short int *)malloc(width*height * 3 * sizeof(unsigned short int));

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image16bit, sizeof(unsigned short int), width*height * 3, filept);

	int i = 0;
	/*--< Find maximum value in the image >--*/
	pixelmax = 0;
	for (x = 0; x < width; x++){
		for (y = 0; y < height; y++){
			red = Image16bit[(x + y*width) * 3];
			green = Image16bit[(x + y*width) * 3 + 1];
			blue = Image16bit[(x + y*width) * 3 + 2];

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
			green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
			blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);

			if (pixelmax < red) pixelmax = red;
			if (pixelmax < green) pixelmax = green;
			if (pixelmax < blue) pixelmax = blue;


			img[i] = red;
			img[i + height*width] = green;
			img[i + 2 * height*width] = blue;

			if (0){ //debug stuff
				fprintf(stderr, "sample %i\t%i\t%i\n", img[i], img[i + height*width], img[i + 2 * height*width]);
				fprintf(stderr, "sample %i\t%i\t%i\n", img[i] >> 6, img[i + height*width] >> 6, img[i + 2 * height*width] >> 6);
			}

			i++;


		}
	}

	free(Image16bit);

}

void aux_read16ppm_1(FILE *filept, int &width, int &height, unsigned short *&img)
{
	int  max, x, y;
	int red, green, blue;
	char dummy[100];

	unsigned short *Image16bit = NULL;

	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);

	if (strncmp(dummy, "P6", 2) != 0) {
		fprintf(stderr, "Error: The input data is not binary PPM.\n");
		exit(1);
	}

	std::cout << width << "\t" << height << "\n";

	img = new unsigned short[width*height * 3]();

	Image16bit = new unsigned short[width*height * 3]();

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image16bit, sizeof(unsigned short int), width*height * 3, filept);

	int i = 0;

	for (x = 0; x < width; x++){
		for (y = 0; y < height; y++){

			red = Image16bit[(x + y*width) * 3];
			green = Image16bit[(x + y*width) * 3 + 1];
			blue = Image16bit[(x + y*width) * 3 + 2];

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
			green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
			blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);

			img[i] = red;
			img[i + height*width] = green;
			img[i + 2 * height*width] = blue;

			i++;

		}
	}

	delete[](Image16bit);

}

// read a 16 bit color ppm image
inline void aux_read16ppm(FILE *filept, int width, int height, int *img)
{
	int  max, x, y;
	int red, green, blue, pixelmax;
	char dummy[100];

	bool divd = 0;

	unsigned short int *Image16bit = NULL;

	//FILE* filept = fopen("007_007.ppm", "r");
	/*--< Read header information of 16bit ppm image from filept >--*/
	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);
	//printf("%s\n", dummy);

	/*--< Really 16bit ppm? Check it >--*/
	if (strncmp(dummy, "P6", 2) != 0) {
		fprintf(stderr, "Error: The input data is not binary PPM.\n");
		exit(1);
	}
	if (max == 65535){
		//fprintf(stderr, "Warning: The input data is 16bit PPM.\n");
		divd = 1;
		//exit(1);
	}
	//if (max != 1023){
	//    fprintf(stderr, "Error: The input data is not 10bit PPM.\n");
	//    exit(1);
	//}

	Image16bit = (unsigned short int *)malloc(width*height * 3 * sizeof(unsigned short int));

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image16bit, sizeof(unsigned short int), width*height * 3, filept);

	int i = 0;
	/*--< Find maximum value in the image >--*/
	pixelmax = 0;
	for (x = 0; x < width; x++){
		for (y = 0; y < height; y++){
			red = Image16bit[(x + y*width) * 3];
			green = Image16bit[(x + y*width) * 3 + 1];
			blue = Image16bit[(x + y*width) * 3 + 2];

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
			green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
			blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);

			if (pixelmax < red) pixelmax = red;
			if (pixelmax < green) pixelmax = green;
			if (pixelmax < blue) pixelmax = blue;


			if (divd) // fix for 16bit to 10bit
			{
				red = red >> 6;
				green = green >> 6;
				blue = blue >> 6;
			}

			img[i] = red;
			img[i + height*width] = green;
			img[i + 2 * height*width] = blue;

			if (0){ //debug stuff
				fprintf(stderr, "sample %i\t%i\t%i\n", img[i], img[i + height*width], img[i + 2 * height*width]);
				fprintf(stderr, "sample %i\t%i\t%i\n", img[i] >> 6, img[i + height*width] >> 6, img[i + 2 * height*width] >> 6);
			}

			i++;


		}
	}

	free(Image16bit);

}

inline void aux_write16pgm(FILE *fp, int width, int height, int *img)
{
	unsigned char *p;
	int i, tmp, j;

	unsigned short* img16bit = (unsigned short*)malloc(height*width*sizeof(unsigned short));
	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = (unsigned short)img[j + i*height];
			lin_ind = lin_ind + 1;
		}
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++){
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		//tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		//tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
	}

	fprintf(fp, "P5\n%d %d\n65535\n", width, height);
	fwrite(img16bit, sizeof(unsigned short int), width*height, fp);
	free(img16bit);
}

inline void aux_write16ppm_16(FILE *fp, int width, int height, unsigned short int *img){
	unsigned char *p;
	int i, tmp, j;

	//unsigned short int maxi = 0;

	unsigned short* img16bit = (unsigned short*)malloc(height*width * 3 * sizeof(unsigned short));
	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = img[j + i*height];
			img16bit[lin_ind + 1] = img[j + i*height + width*height];
			img16bit[lin_ind + 2] = img[j + i*height + 2 * width*height];
			lin_ind = lin_ind + 3;
		}
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++){
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
	}
	fprintf(fp, "P6\n%d %d\n65535\n", width, height);
	fwrite(img16bit, sizeof(unsigned short int), width*height * 3, fp);
	free(img16bit);
}

// write a 16 bit color ppm image
inline void aux_write16ppm(FILE *fp, int width, int height, unsigned short int *img){
	unsigned char *p;
	int i, tmp, j;

	unsigned short* img16bit = (unsigned short*)malloc(height*width * 3 * sizeof(unsigned short));
	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = img[j + i*height];
			img16bit[lin_ind + 1] = img[j + i*height + width*height];
			img16bit[lin_ind + 2] = img[j + i*height + 2 * width*height];
			lin_ind = lin_ind + 3;
		}
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++){
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
	}
	fprintf(fp, "P6\n%d %d\n1023\n", width, height);
	fwrite(img16bit, sizeof(unsigned short int), width*height * 3, fp);
	free(img16bit);
}

// allocates an int vector of size m
inline int* alocaVector(int m)
{
	int *vector, i;
	vector = (int*)malloc(m*sizeof(int));
	for (i = 0; i < m; i++)
	{
		vector[i] = 0;
	}
	return vector;
}

// allocates a double vector of size m
inline double* alocaDoubleVector(int m)
{
	double *vector;
	int i;
	vector = (double*)malloc(m*sizeof(double));
	for (i = 0; i < m; i++)
	{
		vector[i] = 0;
	}
	return vector;
}

inline long aux_GetFileSize(char* filename)
{
	struct stat stat_buf;
	int rc = stat(filename, &stat_buf);
	return rc == 0 ? stat_buf.st_size : -1;
}

int FastOLS_new(const double *AA, const double *Yd, int *PredRegr0, double *PredTheta0, const int Ms, const int MT, const int MPHI, const int N)
{
	int mTheta, M, iM, iM1;
	double *B, *C, sigerr, *Ag, *g;
	double C1, valm1, temp, crit, sabsval;
	int p, j_p, i, j, k, itemp;


	double *PHI = new double[MT*MT]();
	double *PSI = new double[MT]();

	/* make ATA */
	for (int i1 = 0; i1 < MT; i1++){
		for (int j1 = 0; j1 < MT; j1++){
			for (int ii = 0; ii < N; ii++){
				*(PHI + i1 + j1*MT) += (*(AA + ii + i1*N))*(*(AA + ii + j1*N));
			}
		}
	}

	//for (int i1 = 0; i1 < MT; i1++){
	//	for (int j1 = 0; j1 < MT; j1++){
	//		std::cout << *(ATA + i1 + j1*MT) << "\n";
	//	}
	//}

	//std::cout << "------------------------------------------------\n";

	/* make ATYd */
	for (int i1 = 0; i1 < MT; i1++){
		for (int ii = 0; ii < N; ii++){
			*(PSI + i1) += (*(AA + ii + i1*N))*(*(Yd + ii));
		}
	}

	/* YdTYd */
	double yd2 = 0;
	for (int ii = 0; ii < N; ii++)
		yd2 += (*(Yd + ii))*(*(Yd + ii));

	// Usage example: Ms= 3 says the sparsity (length of final predictor) and MT =42 tells how many regressors are available
	// Finally, MPHI = 63 tells the dimensions of the matrices, for getting linear indices in PHI
	M = MT + 1;
	B = alocaDoubleVector(M*M);// ((M+2)*(M+2));
	C = alocaDoubleVector(M*M);// ((M+2)*(M+2));
	Ag = alocaDoubleVector(M*M);// ((M+2)*(M+2));
	g = alocaDoubleVector(M);


	// Inputs: PHI is MTxMT, PSI is MTx1;
	// Outputs: PredRegr and PredTheta are also MTx1 but only first Ms entries are needed
	// Internal variables: B and C are (MT+1)x(MT+1) i.e. MxM

	B[MT + MT*M] = yd2; //B[MT,MT] = yd2; // we start from B[0,0]
	for (iM = 0; iM<MT; iM++)
	{
		PredRegr0[iM] = iM;
		B[iM + MT*M] = PSI[iM];// B[iM,MT] = PSI[iM];
		B[MT + iM*M] = PSI[iM];//B[MT,iM] = PSI[iM];
		for (iM1 = 0; iM1<MT; iM1++)
		{
			B[iM + iM1*M] = PHI[iM + iM1*MPHI];//B[iM,iM1]=PHI[iM,iM1];
		}
	}
	for (iM = 0; iM<M; iM++)
		for (iM1 = 0; iM1<M; iM1++)
			C[iM + iM1*M] = 0;//C[iM,iM1] = 0
	for (iM = 0; iM<MT; iM++)
		C[iM + iM*M] = 1;// C[iM,iM] = 1;
	crit = B[MT + MT*M];//crit = B[MT,MT];
	if (crit < 0.0000001)
	{
		printf("crir, yd2 [%f] [%f] ", crit, yd2);
		i = 0;
		return i;
	}


	for (p = 0; p<Ms; p++)
	{
		valm1 = 0; j_p = 0; // pick the max value in next loop
		for (j = p; j<MT; j++)
		{
			//if(B[j+j*M] > 0.00000000000000001)
			sigerr = B[j + MT*M] * B[j + MT*M] / B[j + j*M];//sigerr  = B[j,MT]*B[j,MT]/B[j,j];
			//else
			//	sigerr = 0;
			if (sigerr > valm1)
			{
				valm1 = sigerr;
				j_p = j;
			}
		} // j_p is the index of maximum
		crit = crit - valm1;
		itemp = PredRegr0[j_p]; PredRegr0[j_p] = PredRegr0[p]; PredRegr0[p] = itemp;
		for (j = p; j<M; j++)
		{
			//% interchange B(p:end,j_p) with B(p:end,p)
			temp = B[j + j_p*M]; //temp = B[j,j_p];
			B[j + j_p*M] = B[j + p*M]; //B[j,j_p] = B[j,p];
			B[j + p*M] = temp;// B[j,p] = temp;
		}
		for (j = p; j<M; j++)
		{
			//% interchange B(j_p,p:end) with B(p,p:end)
			temp = B[j_p + j*M]; // temp = B[j_p,j];
			B[j_p + j*M] = B[p + j*M]; //B[j_p,j] = B[p,j];
			B[p + j*M] = temp;//B[p,j] = temp;
		}

		//% Fast
		for (j = 0; j <= p - 1; j++)
		{
			// % interchange C(1:p-1,j_p) with C(1:p-1,p)
			temp = C[j + j_p*M]; // temp = C[j,j_p];
			C[j + j_p*M] = C[j + p*M]; //C[j,j_p] = C[j,p];
			C[j + p*M] = temp; //C[j,p] = temp;
		}
		for (j = (p + 1); j<M; j++)
		{
			//if(B[p+p*M] > 0.00000000000000000000001)
			C[p + j*M] = B[p + j*M] / B[p + p*M];//C[p,j] = B[p,j]/B[p,p];
			//else
			//C[p+j*M] = 0;
		}

		for (j = (p + 1); j<MT; j++)
			for (k = j; k <= MT; k++)
			{
			B[j + k*M] = B[j + k*M] - C[p + j*M] * C[p + k*M] * B[p + p*M];//B[j,k] = B[j,k]-C[p,j]*C[p,k]*B[p,p];
			}
		for (j = (p + 1); j<MT; j++)
			for (k = j; k <= MT; k++)
			{
			B[k + j*M] = B[j + k*M];//B[k,j] = B[j,k];
			}
		//for j = (p+1):M
		//    for k = j:(M+1)
		//        B(j,k) = B(j,k)-C(p,j)*C(p,k)*B(p,p);
		//        %B(j,k) = B(j,k)-C(p,j)*B(p,k);
		//        B(k,j) = B(j,k);
		//    end
		//end
		//        for( iM=0; iM<Ms; iM++)
		//        {
		//        for( iM1=0; iM1<Ms; iM1++)
		//        	printf("C[%f] ",C[iM+iM1*M]);
		//        printf(" \n" );
		//        }
		// scanf("%d",&i);
	} //% for( p=0; p<M; p++ )


	// final triangular backsolving
	for (i = 0; i<Ms; i++)
	{
		g[i] = C[i + MT*M];//g[i] = C[i,MT];
		for (j = 0; j<Ms; j++)
			Ag[i + j*M] = C[i + j*M];//Ag[i,j] = C[i,j];
	}
	PredTheta0[Ms - 1] = g[Ms - 1];
	for (i = Ms - 2; i >= 0; i--)
	{
		PredTheta0[i] = g[i];
		for (j = i + 1; j<Ms; j++)
			PredTheta0[i] = PredTheta0[i] - Ag[i + j*M] * PredTheta0[j];//PredTheta[i] = PredTheta[i]- Ag[i,j]*PredTheta[j];
	}
	//printf("pred FASTOLS [%f][%f][%f][%f][%f][%f]\n",PredTheta0[0], PredTheta0[1], PredTheta0[2], PredTheta0[3], PredTheta0[4], PredTheta0[5]);
	// printf("pred [%d][%d][%d][%d][%d]\n",PredRegr0[0], PredRegr0[1], PredRegr0[2], PredRegr0[3], PredRegr0[4] );

	if (PredTheta0[0] != PredTheta0[0])
	{// if is nan
		printf("PredTheta0[0]  is NaN\n");
		PredTheta0[0] = 1.0;
		for (i = 1; i < Ms; i++)
		{
			PredTheta0[i] = 0.0;
		}
	}

	sabsval = 0;
	for (i = 0; i<Ms; i++)
	{

		if (PredTheta0[i] != PredTheta0[i])
		{// if is nan
			PredTheta0[i] = 0.0;
			printf("PredTheta0[%i]  is NaN\n", i);
		}

		if (PredTheta0[i] > 0)
			sabsval = sabsval + PredTheta0[i];
		else
			sabsval = sabsval - PredTheta0[i];
	}
	//printf("%f\n", sabsval);
	if (sabsval > 2 * Ms) // if average coefficients are too high forget about intrpolation
	{
		PredTheta0[0] = C[0 + MT*M];//g[0] = C[0,MT];

		//PredTheta0[0] = 1;

		if (PredTheta0[0] != PredTheta0[0])
		{// if is nan
			printf("C[0+MT*M]  is NaN\n", i);
			PredTheta0[0] = 1.0;
		}

		// fix ???
		//if (abs(PredTheta0[0]) > 100)
		//PredTheta0[0] = PredTheta0[0] / abs(PredTheta0[0]);

		for (i = 1; i<Ms; i++) {
			PredTheta0[i] = 0.0;
		}

		i = 1;
	}
	else {
		i = Ms;
	}

	free(B);
	free(C);
	free(Ag);
	free(g);

	delete[](PSI);
	delete[](PHI);

	//printf("pred Ms [%d]",Ms);
	return i;

}

int FastOLS(const double *PHI, const double *PSI, const double yd2, int *PredRegr0, double *PredTheta0, const int Ms, const int MT, const int MPHI)
{
	int mTheta, M, iM, iM1;
	double *B, *C, sigerr, *Ag, *g;
	double C1, valm1, temp, crit, sabsval;
	int p, j_p, i, j, k, itemp;

	// Usage example: Ms= 3 says the sparsity (length of final predictor) and MT =42 tells how many regressors are available
	// Finally, MPHI = 63 tells the dimensions of the matrices, for getting linear indices in PHI
	M = MT + 1;
	B = alocaDoubleVector(M*M);// ((M+2)*(M+2));
	C = alocaDoubleVector(M*M);// ((M+2)*(M+2));
	Ag = alocaDoubleVector(M*M);// ((M+2)*(M+2));
	g = alocaDoubleVector(M);


	// Inputs: PHI is MTxMT, PSI is MTx1;
	// Outputs: PredRegr and PredTheta are also MTx1 but only first Ms entries are needed
	// Internal variables: B and C are (MT+1)x(MT+1) i.e. MxM

	B[MT + MT*M] = yd2; //B[MT,MT] = yd2; // we start from B[0,0]
	for (iM = 0; iM<MT; iM++)
	{
		PredRegr0[iM] = iM;
		B[iM + MT*M] = PSI[iM];// B[iM,MT] = PSI[iM];
		B[MT + iM*M] = PSI[iM];//B[MT,iM] = PSI[iM];
		for (iM1 = 0; iM1<MT; iM1++)
		{
			B[iM + iM1*M] = PHI[iM + iM1*MPHI];//B[iM,iM1]=PHI[iM,iM1];
		}
	}
	for (iM = 0; iM<M; iM++)
		for (iM1 = 0; iM1<M; iM1++)
			C[iM + iM1*M] = 0;//C[iM,iM1] = 0
	for (iM = 0; iM<MT; iM++)
		C[iM + iM*M] = 1;// C[iM,iM] = 1;
	crit = B[MT + MT*M];//crit = B[MT,MT];
	if (crit < 0.0000001)
	{
		printf("crir, yd2 [%f] [%f] ", crit, yd2);
		i = 0;
		return i;
	}


	for (p = 0; p<Ms; p++)
	{
		valm1 = 0; j_p = 0; // pick the max value in next loop
		for (j = p; j<MT; j++)
		{
			//if(B[j+j*M] > 0.00000000000000001)
			sigerr = B[j + MT*M] * B[j + MT*M] / B[j + j*M];//sigerr  = B[j,MT]*B[j,MT]/B[j,j];
			//else
			//	sigerr = 0;
			if (sigerr > valm1)
			{
				valm1 = sigerr;
				j_p = j;
			}
		} // j_p is the index of maximum
		crit = crit - valm1;
		itemp = PredRegr0[j_p]; PredRegr0[j_p] = PredRegr0[p]; PredRegr0[p] = itemp;
		for (j = p; j<M; j++)
		{
			//% interchange B(p:end,j_p) with B(p:end,p)
			temp = B[j + j_p*M]; //temp = B[j,j_p];
			B[j + j_p*M] = B[j + p*M]; //B[j,j_p] = B[j,p];
			B[j + p*M] = temp;// B[j,p] = temp;
		}
		for (j = p; j<M; j++)
		{
			//% interchange B(j_p,p:end) with B(p,p:end)
			temp = B[j_p + j*M]; // temp = B[j_p,j];
			B[j_p + j*M] = B[p + j*M]; //B[j_p,j] = B[p,j];
			B[p + j*M] = temp;//B[p,j] = temp;
		}

		//% Fast
		for (j = 0; j <= p - 1; j++)
		{
			// % interchange C(1:p-1,j_p) with C(1:p-1,p)
			temp = C[j + j_p*M]; // temp = C[j,j_p];
			C[j + j_p*M] = C[j + p*M]; //C[j,j_p] = C[j,p];
			C[j + p*M] = temp; //C[j,p] = temp;
		}
		for (j = (p + 1); j<M; j++)
		{
			//if(B[p+p*M] > 0.00000000000000000000001)
			C[p + j*M] = B[p + j*M] / B[p + p*M];//C[p,j] = B[p,j]/B[p,p];
			//else
			//C[p+j*M] = 0;
		}

		for (j = (p + 1); j<MT; j++)
			for (k = j; k <= MT; k++)
			{
			B[j + k*M] = B[j + k*M] - C[p + j*M] * C[p + k*M] * B[p + p*M];//B[j,k] = B[j,k]-C[p,j]*C[p,k]*B[p,p];
			}
		for (j = (p + 1); j<MT; j++)
			for (k = j; k <= MT; k++)
			{
			B[k + j*M] = B[j + k*M];//B[k,j] = B[j,k];
			}
		//for j = (p+1):M
		//    for k = j:(M+1)
		//        B(j,k) = B(j,k)-C(p,j)*C(p,k)*B(p,p);
		//        %B(j,k) = B(j,k)-C(p,j)*B(p,k);
		//        B(k,j) = B(j,k);
		//    end
		//end
		//        for( iM=0; iM<Ms; iM++)
		//        {
		//        for( iM1=0; iM1<Ms; iM1++)
		//        	printf("C[%f] ",C[iM+iM1*M]);
		//        printf(" \n" );
		//        }
		// scanf("%d",&i);
	} //% for( p=0; p<M; p++ )


	// final triangular backsolving
	for (i = 0; i<Ms; i++)
	{
		g[i] = C[i + MT*M];//g[i] = C[i,MT];
		for (j = 0; j<Ms; j++)
			Ag[i + j*M] = C[i + j*M];//Ag[i,j] = C[i,j];
	}
	PredTheta0[Ms - 1] = g[Ms - 1];
	for (i = Ms - 2; i >= 0; i--)
	{
		PredTheta0[i] = g[i];
		for (j = i + 1; j<Ms; j++)
			PredTheta0[i] = PredTheta0[i] - Ag[i + j*M] * PredTheta0[j];//PredTheta[i] = PredTheta[i]- Ag[i,j]*PredTheta[j];
	}
	//printf("pred FASTOLS [%f][%f][%f][%f][%f][%f]\n",PredTheta0[0], PredTheta0[1], PredTheta0[2], PredTheta0[3], PredTheta0[4], PredTheta0[5]);
	// printf("pred [%d][%d][%d][%d][%d]\n",PredRegr0[0], PredRegr0[1], PredRegr0[2], PredRegr0[3], PredRegr0[4] );

	if (PredTheta0[0] != PredTheta0[0])
	{// if is nan
		printf("PredTheta0[0]  is NaN\n");
		PredTheta0[0] = 1.0;
		for (i = 1; i < Ms; i++)
		{
			PredTheta0[i] = 0.0;
		}
	}

	sabsval = 0;
	for (i = 0; i<Ms; i++)
	{

		if (PredTheta0[i] != PredTheta0[i])
		{// if is nan
			PredTheta0[i] = 0.0;
			printf("PredTheta0[%i]  is NaN\n", i);
		}

		if (PredTheta0[i] > 0)
			sabsval = sabsval + PredTheta0[i];
		else
			sabsval = sabsval - PredTheta0[i];
	}
	//printf("%f\n", sabsval);
	if (sabsval > 2 * Ms) // if average coefficients are too high forget about intrpolation
	{
		PredTheta0[0] = C[0 + MT*M];//g[0] = C[0,MT];

		//PredTheta0[0] = 1;

		if (PredTheta0[0] != PredTheta0[0])
		{// if is nan
			printf("C[0+MT*M]  is NaN\n", i);
			PredTheta0[0] = 1.0;
		}

		// fix ???
		//if (abs(PredTheta0[0]) > 100)
		//PredTheta0[0] = PredTheta0[0] / abs(PredTheta0[0]);

		for (i = 1; i<Ms; i++) {
			PredTheta0[i] = 0.0;
		}

		i = 1;
	}
	else {
		i = Ms;
	}

	free(B);
	free(C);
	free(Ag);
	free(g);

	//printf("pred Ms [%d]",Ms);
	return i;

}

#endif
