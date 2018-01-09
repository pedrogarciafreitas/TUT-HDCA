/* Reads Ostendo disparity (.dsp) and RGB image (.ppm), and creates quantized disparity of 512 levels, based on aggregation of disparity+color into consistent regions. */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "gen_types.hh"

#include "../CERV/cerv.h"

int main(const int argc, const char** argv) {

	bool encode = 0;

	if (encode){

		int nr, nc, ncomponents, nvc;

		aux_read_header_file(&nr, &nc, &ncomponents, &nvc, argv[1]);

		printf("nr %d \t nc %d \t ncomponents %d \n", nr, nc, ncomponents);

		int *segmentation_data = new int[nr*nc*ncomponents];

		aux_read_file_int32(nr, nc, ncomponents, argv[1], segmentation_data);

		/* encode quantized labels with cerv */
		int** SEGM2D = (int**)malloc(nr*sizeof(int*));
		for (int i = 0; i < nr; ++i) { SEGM2D[i] = (int*)malloc(nc*sizeof(int)); }
		for (int i = 0; i < nr; ++i) {
			for (int j = 0; j < nc; ++j) {
				SEGM2D[i][j] = segmentation_data[i + j*nr];
			}
		}

		cerv_encode(SEGM2D, nr, nc, argv[2]);

	}
	
	if (!encode){

		/* decode quantized labels with cerv */

		int nr = atoi(argv[2]);
		int nc = atoi(argv[3]);

		int* SEGMFINAL = alocaVector(nr*nc);

		int** SEGM2D = (int**)malloc(nr*sizeof(int*));
		for (int i = 0; i < nr; ++i) { SEGM2D[i] = (int*)malloc(nc*sizeof(int)); }
		for (int i = 0; i < nr; ++i) {
			for (int j = 0; j < nc; ++j) {
				SEGM2D[i][j] = 0;
			}
		}

		cerv_decode(SEGM2D, nr, nc, argv[1]);

		int number_of_regions = 0;

		for (int i = 0; i < nr; ++i) {
			for (int j = 0; j < nc; ++j) {
				SEGMFINAL[i + j*nr] = SEGM2D[i][j];
				if (SEGMFINAL[i + j*nr] > number_of_regions)
					number_of_regions++;
			}
		}

		/* WRITE QUANTIZED DISPARITY AND LABELS TO DISK */
		aux_write_header_file(nr, nc, 1, 1, argv[4]);

		FILE *f_file = fopen(argv[4], "a+b");

		fwrite(SEGMFINAL, sizeof(int), nr*nc, f_file);

		fclose(f_file);

	}

	return 0;

}