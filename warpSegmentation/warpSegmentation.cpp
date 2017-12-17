/* View warping. */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "gen_types.hh"

int main(const int argc, const char** argv) {

	int nr, nc, ncomponents, nvc;

	int ii0 = atoi(argv[3]), jj0 = atoi(argv[4]), ii1 = atoi(argv[5]), jj1 = atoi(argv[6]);

	aux_read_header_file(&nr, &nc, &ncomponents, &nvc, argv[1]);


	/* reference */
	unsigned int *St = new unsigned int[nr*nc*1];
	aux_read_file_uint32(nr, nc, 1, argv[1], St);

	/* quantized disparity */
	float *quantDM = new float[nr*nc];
	aux_read_file_float(nr, nc, 1, argv[2], quantDM);


	int *ROWS = new int[nr*nc];
	int *COLS = new int[nr*nc];

	aux_read_file_int32(nr, nc, 1, argv[7], ROWS);
	aux_read_file_int32(nr, nc, 1, argv[8], COLS);


	float *DispTarg = new float[nr*nc];

	unsigned int *WarpedSt = new unsigned int[nr*nc];

	for (int ij = 0; ij < nr*nc; ij++)
		DispTarg[ij] = -1;


	for (int ij = 0; ij < nr*nc; ij++)
	{
		float disp0 = quantDM[ij];

		int ix = ij % nr; //row
		int iy = (ij - ix) / nr; //col

		int iynew = iy + (int)round((float)(ii1 - ii0)*quantDM[ij] + (float)COLS[ij]);
		int ixnew = ix + (int)round((float)(jj1 - jj0)*quantDM[ij] * (60.0 / 40.0) + (float)ROWS[ij]);

		if (iynew >= 0 && iynew<nc && ixnew >= 0 && ixnew < nr){
			int indnew = ixnew + iynew*nr;
			if (DispTarg[indnew] < disp0){
				DispTarg[indnew] = disp0;
				WarpedSt[indnew] = St[ij];
			}
		}

	}

	aux_write_header_file(nr, nc, 1, 1, argv[9]);

	FILE *f_file = fopen(argv[9], "a+b");

	fwrite(WarpedSt, sizeof(unsigned int), nr*nc, f_file);

	fclose(f_file);

	/*
	aux_write_header_file(nr, nc, 1, 1, argv[10]);

	f_file = fopen(argv[10], "a+b");

	fwrite(DispTarg, sizeof(float), nr*nc, f_file);

	fclose(f_file);
	*/

	return 0;

}