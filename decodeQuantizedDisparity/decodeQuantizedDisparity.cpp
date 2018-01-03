/* Reads Ostendo disparity (.dsp) and RGB image (.ppm), and creates quantized disparity of 512 levels, based on aggregation of disparity+color into consistent regions. */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "gen_types.hh"
#include "msImageProcessor.h"

#include "../CERV/cerv.h"

int main(const int argc, const char** argv) {

	int nr, nc, nregs, tmp;

	aux_read_header_file(&nr, &nc, &nregs, &tmp, argv[2]);

	printf("%d\t%d\t%d\n", nr, nc, nregs);

	/* encode quantized labels with cerv */
	int** SEGM2D = (int**)malloc(nr*sizeof(int*));
	for (int i = 0; i < nr; ++i) {
		SEGM2D[i] = (int*)malloc(nc*sizeof(int));
	}
	/* decode quantized labels with cerv */

	int* SEGMFINAL = alocaVector((int)nr*nc);

	for (int i = 0; i < nr; ++i) {
		for (int j = 0; j < nc; ++j) {
			SEGM2D[i][j] = 0;
		}
	}

	cerv_decode(SEGM2D, nr, nc, argv[1]);
	for (int i = 0; i < nr; ++i) {
		for (int j = 0; j < nc; ++j) {
			SEGMFINAL[i + j*nr] = SEGM2D[i][j];
		}
	}

	/* these disparity values corresponding to different labels need to be transmitted also */
	float *quantDisparities_labels = new float[nregs];

	FILE *f_file;
	f_file = fopen(argv[2], "rb");
	fseek(f_file, 4 * sizeof(int), SEEK_SET);
	fread(quantDisparities_labels, sizeof(float), nregs, f_file);

	float *quantDM = new float[nr*nc];

	/* reassign labels based on cerv labeling */
	for (int ij = 0; ij < nr*nc; ij++){
		quantDM[ij] = quantDisparities_labels[ SEGMFINAL[ij] ];
	}


	/* WRITE QUANTIZED DISPARITY TO DISK */
	aux_write_header_file(nr, nc, 1, 1, argv[3]);

	f_file = fopen(argv[3], "a+b");

	fwrite(quantDM, sizeof(float), nr*nc, f_file);

	fclose(f_file);

	return 0;

}