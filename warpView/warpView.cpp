/* View warping. */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "gen_types.hh"

int main(const int argc, const char** argv) {

	int nr, nc, ncomponents, nvc;

	int ii0 = atof(argv[3]), jj0 = atof(argv[4]), ii1 = atoi(argv[5]), jj1 = atoi(argv[6]);

	aux_read_header_file(&nr, &nc, &ncomponents, &nvc, argv[1]);


	/* reference */
	unsigned short *AA1 = new unsigned short[nr*nc*ncomponents];
	aux_read_file_uint16(nr, nc, ncomponents, argv[1], AA1);

	/* quantized disparity */
	float *quantDM = new float[nr*nc];
	aux_read_file_float(nr, nc, 1, argv[2], quantDM);
	
	
	int *ROWS = new int[nr*nc];
	int *COLS = new int[nr*nc];

	//aux_read_file_int32(nr, nc, 1, argv[7], ROWS); //maybe we use these later ...
	//aux_read_file_int32(nr, nc, 1, argv[8], COLS);


	float *DispTarg = new float[nr*nc];

	float *DispTarg_col = new float[nr*nc];
	float *DispTarg_row = new float[nr*nc];

	unsigned short *Warped = new unsigned short[nr*nc * 3];

	for (int ij = 0; ij < nr*nc; ij++){
		DispTarg[ij] = -1;
		DispTarg_col[ij] = 99999;
		DispTarg_row[ij] = 99999;
	}

	
	for (int ij = 0; ij < nr*nc; ij++)
	{
		float disp0 = quantDM[ij];

		int ix = ij % nr; //row
		int iy = (ij - ix) / nr; //col

		int iynew = iy + (int)round(((float)(ii1 - ii0))*quantDM[ij] + (float)COLS[ij]);
		int ixnew = ix + (int)round(((float)(jj1 - jj0))*quantDM[ij] + (float)ROWS[ij]);

		if (iynew >= 0 && iynew<nc && ixnew >= 0 && ixnew < nr){
			int indnew = ixnew + iynew*nr;
			if (DispTarg[indnew] < disp0){
				DispTarg[indnew] = disp0;
				DispTarg_col[indnew] = iynew - iy;
				DispTarg_row[indnew] = ixnew - ix;
				Warped[indnew] = AA1[ij];
				Warped[indnew + nr*nc] = AA1[ij+nr*nc];
				Warped[indnew + 2*nr*nc] = AA1[ij+2*nr*nc];
			}
		}

	}
	
	aux_write_header_file(nr, nc, 3, 1, argv[9]);

	FILE *f_file = fopen(argv[9], "a+b");

	fwrite(Warped, sizeof(unsigned short), nr*nc*3, f_file);

	fclose(f_file);

	aux_write_header_file(nr, nc, 1, 1, argv[10]);

	f_file = fopen(argv[10], "a+b");

	fwrite(DispTarg, sizeof(float), nr*nc, f_file);

	fclose(f_file);

	aux_write_header_file(nr, nc, 1, 1, argv[11]);

	f_file = fopen(argv[11], "a+b");

	fwrite(DispTarg_col, sizeof(float), nr*nc, f_file);

	fclose(f_file);

	aux_write_header_file(nr, nc, 1, 1, argv[12]);

	f_file = fopen(argv[12], "a+b");

	fwrite(DispTarg_row, sizeof(float), nr*nc, f_file);

	fclose(f_file);
	

	return 0;

}