#include <iostream>
#include <stdio.h>

#include "gen_types.hh"
#include "msImageProcessor.h"

#ifdef _MSC_VER
typedef __int8 int8_t;
typedef __int32 int32_t;
typedef unsigned __int8 uint8_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

int main(const int argc, const char** argv) {


	// Read in the input variables
	int sigmaS = atoi( argv[2] );
	float sigmaR = atof(argv[3]);
	int minRegion = atoi(argv[4]);

	msImageProcessor im_proc;

	int nr, nc, ncomponents, nvc;

	aux_read_header_file(&nr, &nc, &ncomponents, &nvc, argv[1]);

	printf("nr %d \t nc %d \t ncomponents %d \n", nr, nc, ncomponents);

	unsigned char *rgb_data_uint8 = new unsigned char[nr*nc*ncomponents];

	aux_read_file_uint8(nr, nc, ncomponents, argv[1], rgb_data_uint8);

	//printf("debug at 5678 \t %d", rgb_data_uint8[5678 - 1]);

	unsigned char *rgb_image = new unsigned char[nr*nc*ncomponents];

	unsigned char *B = rgb_image;

	for (int row = 0; row < nr; row++){
		for (int col = 0; col < nr*nc; col += nr){
			*B++ = rgb_data_uint8[row + col];
			*B++ = rgb_data_uint8[row + col + nr*nc];
			*B++ = rgb_data_uint8[row + col + 2 * nr*nc];
		}
	}

	im_proc.DefineImage(rgb_image, COLOR, nr, nc);

	im_proc.Segment(sigmaS, sigmaR, minRegion, HIGH_SPEEDUP);

	int *labels = im_proc.GetLabels();

	uint32_t *labelIm = new uint32_t[nr*nc];

	for (int row = 0; row < nr; row++){
		for (int col = 0; col < nr*nc; col += nr){
			labelIm[row + col] = (uint32_t)(*labels++) + 1;
		}
	}

	aux_write_header_file(nr, nc, 1, 1, argv[5]);

	FILE* f_file = fopen(argv[5], "a+b");

	fwrite(labelIm, sizeof(uint32_t), nr*nc, f_file);

	return 0;

}