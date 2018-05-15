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

#define IO_V true
#define NBIT_GR 32
#define SYSTEM_VERBOSE 0

#define BIT_DEPTH 10
#define D_DEPTH 14
#define BIT_DEPTH_RESIDUAL 16


const bool verbose = false;

// number of pixels in the original lenslet image.
const int num_pixels_in_lenslet = 41483904;

int system_1(char *str){

	std::string sys_call_str(str);

#ifdef SYSTEM_VERBOSE

#ifdef _WIN32 || _WIN64 
	sys_call_str.append("> nul");
#endif
#ifdef __unix__
	sys_call_str.append("> /dev/null");
#endif

#endif

	return system(sys_call_str.c_str());

}


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
	int dsz = (SZ / 2);
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
	int dsz = (SZ / 2);
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
	int dsz = (SZ / 2);
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

bool aux_read16PGMPPM(const char* filename, int &width, int &height, int &ncomp, unsigned short *&img)
{
	if (IO_V)
		printf("Reading %s\n", filename);

	int  max, x, y;
	int red, green, blue;
	char dummy[100];

	unsigned short *Image16bit = NULL;

	FILE *filept;

	filept = fopen(filename, "rb");

	if (filept == NULL){
		printf("%s does not exist\n", filename);
		return false;
	}


	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);

	if (!strncmp(dummy, "P6", 2)) {
		ncomp = 3;
	}
	else if (!strncmp(dummy, "P5", 2)){ 
		ncomp = 1;  
	}
	else{ printf("ERROR NOT PGM OR PPM\n"); return false; }


	//std::cout << width << "\t" << height << "\n";

	img = new unsigned short[width*height * ncomp]();

	Image16bit = new unsigned short[width*height * ncomp]();

	/*--< Read 16bit ppm image from filept >--*/
	int nread = fread(Image16bit, sizeof(unsigned short), width*height * ncomp, filept);

	if (nread != width*height * ncomp)
	{
		fprintf(stderr, "READ ERROR aux_read16ppm() %s\n", filename);
		return false;
	}

	fclose(filept);

	int i = 0;

	for (x = 0; x < width; x++){
		for (y = 0; y < height; y++){

			red = Image16bit[(x + y*width) * ncomp];
			if (ncomp == 3){
				green = Image16bit[(x + y*width) * ncomp + 1];
				blue = Image16bit[(x + y*width) * ncomp + 2];
			}

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
			if (ncomp == 3){
				green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
				blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);
			}

			img[i] = red;
			if (ncomp == 3){
				img[i + height*width] = green;
				img[i + 2 * height*width] = blue;
			}

			i++;

		}
	}

	delete[](Image16bit);

	return true;

}

bool aux_read16pgm(const char* filename, int &width, int &height, unsigned short *&img)
{
	if (IO_V)
		printf("Reading %s\n", filename);

	int  max, x, y;
	int red, green, blue, pixelmax;
	char dummy[100];

	//int width, height;

	FILE *filept;

	filept = fopen(filename, "rb");

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
		return false;
	}

	unsigned short *Image16bit = new unsigned short[width*height]();
	img = new unsigned short[width*height]();

	/*--< Read 16bit ppm image from filept >--*/
	int nread = fread(Image16bit, sizeof(unsigned short), width*height, filept);

	if (nread != width*height){
		fprintf(stderr, "READ ERROR aux_read16pgm() %s\n", filename);
		return false;
	}


	fclose(filept);

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
	
	delete[](Image16bit);

	return true;
}

bool aux_read16ppm(const char* filename, int &width, int &height, unsigned short *&img)
{
	if (IO_V)
		printf("Reading %s\n", filename);

	int  max, x, y;
	int red, green, blue;
	char dummy[100];

	unsigned short *Image16bit = NULL;

	FILE *filept;

	filept = fopen(filename, "rb");

	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);

	if (strncmp(dummy, "P6", 2) != 0) {
		fprintf(stderr, "Error: The input data is not binary PPM.\n");
		return false;
	}

	//std::cout << width << "\t" << height << "\n";

	img = new unsigned short[width*height * 3]();

	Image16bit = new unsigned short[width*height * 3]();

	/*--< Read 16bit ppm image from filept >--*/
	int nread = fread(Image16bit, sizeof(unsigned short), width*height * 3, filept);

	if (nread != width*height * 3)
	{
		fprintf(stderr,"READ ERROR aux_read16ppm() %s\n",filename);
		return false;
	}

	fclose(filept);

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

	return true;

}

bool aux_write16PGMPPM(const char* filename, const int width, const int height, const int ncomp, unsigned short *img)
{

	if (IO_V)
		printf("Writing %s\n", filename);

	unsigned char *p;
	int i, tmp, j;

	unsigned short maxi = 0;

	FILE *filept;

	filept = fopen(filename, "wb");

	if (filept == NULL){
		printf("Cannot open %s\n", filename);
		return false;
	}

	unsigned short* img16bit = new unsigned short[height*width * ncomp]();

	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = img[j + i*height];
			if (ncomp == 3){
				img16bit[lin_ind + 1] = img[j + i*height + width*height];
				img16bit[lin_ind + 2] = img[j + i*height + 2 * width*height];
				lin_ind = lin_ind + 3;
			}
			else{
				lin_ind = lin_ind + 1;
			}
			
		}
	}

	for (i = 0; i < width*height * ncomp; i++)
	{
		if (*(img16bit + i) > maxi)
			maxi = *(img16bit + i);
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++){
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		if (ncomp == 3){
			tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
			tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		}
	}


	if (ncomp == 3){
		fprintf(filept, "P6\n%d %d\n%d\n", width, height, maxi);
	}
	else{
		fprintf(filept, "P5\n%d %d\n%d\n", width, height, maxi);
	}

	fwrite(img16bit, sizeof(unsigned short), width*height * ncomp, filept);

	fclose(filept);

	delete[](img16bit);

	return true;
}

void aux_write16pgm(const char* filename, int width, int height, unsigned short *img)
{
	if (IO_V)
		printf("Writing %s\n", filename);

	unsigned char *p;
	int i, tmp, j;

	unsigned int maxi = 0;


	FILE *filept;

	filept = fopen(filename, "wb");

	unsigned short* img16bit = new unsigned short[height*width]();

	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = (unsigned short)img[j + i*height];
			lin_ind = lin_ind + 1;
		}
	}

	for (i = 0; i < width*height; i++)
	{
		if (*(img16bit + i) > maxi)
			maxi = *(img16bit + i);
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++){
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
	}

	fprintf(filept, "P5\n%d %d\n%d\n", width, height, maxi);
	fwrite(img16bit, sizeof(unsigned short int), width*height, filept);
	
	fclose(filept);

	delete[](img16bit);
}

void aux_write16ppm(const char* filename, int width, int height, unsigned short *img)
{
	if (IO_V)
		printf("Writing %s\n", filename);

	unsigned char *p;
	int i, tmp, j;

	unsigned short int maxi = 0;

	FILE *filept;

	filept = fopen(filename, "wb");

	unsigned short* img16bit = new unsigned short[height*width*3]();

	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = img[j + i*height];
			img16bit[lin_ind + 1] = img[j + i*height + width*height];
			img16bit[lin_ind + 2] = img[j + i*height + 2 * width*height];
			lin_ind = lin_ind + 3;
		}
	}

	for (i = 0; i < width*height * 3; i++)
	{
		if (*(img16bit + i) > maxi)
			maxi = *(img16bit + i);
	}

	p = (unsigned char *) img16bit;
	for (i = 0; i < width*height; i++){
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
	}



	fprintf(filept, "P6\n%d %d\n%d\n", width, height, maxi);
	fwrite(img16bit, sizeof(unsigned short int), width*height * 3, filept);

	fclose(filept);

	delete[](img16bit);
}

void decodeResidualJP2(unsigned short *ps, const char *kdu_expand_path, const char *jp2_residual_path_jp2, const char *ppm_residual_path, int ncomp, const int offset, const int maxvali)
{
	/* decode residual with kakadu */
	char kdu_expand_s[256];
	sprintf(kdu_expand_s, "\"%s\"%s%s%s%s", kdu_expand_path, " -i ", jp2_residual_path_jp2, " -o ", ppm_residual_path);

	//std::cout << kdu_expand_s << "\n";

	int status = system_1(kdu_expand_s);

	/* apply residual */

	unsigned short* jp2_residual;

	int nc1, nr1;

	if (aux_read16PGMPPM(ppm_residual_path, nc1, nr1, ncomp, jp2_residual))
	{

		for (int iir = 0; iir < nc1*nr1 * ncomp; iir++)
		{
			signed int val = ((signed int)*(ps + iir)) + ((signed int)jp2_residual[iir]) - offset;
			if (val < 0)
				val = 0;
			if (val > maxvali)
				val = maxvali;
			*(ps + iir) = (unsigned short)(val);
		}

		delete[](jp2_residual);
	}
}

void encodeResidualJP2(const int nr, const int nc, unsigned short *original_intermediate_view, unsigned short *ps, const char *ppm_residual_path,
	const char *kdu_compress_path, const char *jp2_residual_path_jp2, const float residual_rate, const int ncomp, const int offset)
{
	/*establish residual*/
	unsigned short *residual_image = new unsigned short[nr*nc * ncomp]();

	for (int iir = 0; iir < nr*nc*ncomp; iir++){
		signed int res_val = (((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + offset);
		if (res_val>pow(2, 16) - 1)
			res_val = pow(2, 16) - 1;
		if (res_val < 0)
			res_val = 0;
		*(residual_image + iir) = (unsigned short)res_val;
	}

	aux_write16PGMPPM(ppm_residual_path, nc, nr, ncomp, residual_image);

	delete[](residual_image);

	/* here encode residual with kakadu */

	char kdu_compress_s[256];
	sprintf(kdu_compress_s, "\"%s\"%s%s%s%s%s%f", kdu_compress_path, " -i ", ppm_residual_path, " -o ", jp2_residual_path_jp2, " -no_weights -precise -full -rate ", residual_rate);

	//std::cout << kdu_compress_s << "\n";

	int status = system_1(kdu_compress_s);
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

#endif
