/* View warping. */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "gen_types.hh"

int main(const int argc, const char** argv) {

	int nr, nc, ncomponents, nvc;

	//float ii0 = atof(argv[3]), jj0 = atof(argv[4]), ii1 = atof(argv[5]), jj1 = atof(argv[6]);

	aux_read_header_file(&nr, &nc, &ncomponents, &nvc, argv[1]);

	/* reference */
	unsigned short *AA1 = new unsigned short[nr*nc*ncomponents];
	aux_read_file_uint16(nr, nc, ncomponents, argv[1], AA1);

	/* disparity */
	float *DM_ROW = new float[nr*nc];
	aux_read_file_float(nr, nc, 1, argv[2], DM_ROW);

	float *DM_COL = new float[nr*nc];
	aux_read_file_float(nr, nc, 1, argv[3], DM_COL);
	
	
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
		float disp0 = abs( DM_COL[ij] ) + abs( DM_ROW[ij] );

		int ix = ij % nr; //row
		int iy = (ij - ix) / nr; //col

		int iynew = iy + (int)round( DM_COL[ij] );
		int ixnew = ix + (int)round( DM_ROW[ij] );

		if (iynew >= 0 && iynew<nc && ixnew >= 0 && ixnew < nr && (int)round(DM_COL[ij])>-9999){
			int indnew = ixnew + iynew*nr;
			if (DispTarg[indnew] < disp0){
				DispTarg[indnew] = disp0;
				DispTarg_col[indnew] = iynew - iy; // warped horizontal disparity
				DispTarg_row[indnew] = ixnew - ix; // warped vertical disparity
				Warped[indnew] = AA1[ij];
				Warped[indnew + nr*nc] = AA1[ij+nr*nc];
				Warped[indnew + 2*nr*nc] = AA1[ij+2*nr*nc];
			}
		}

	}
	
	aux_write_header_file(nr, nc, 3, 1, argv[4]);

	FILE *f_file = fopen(argv[4], "a+b");

	fwrite(Warped, sizeof(unsigned short), nr*nc*3, f_file);

	fclose(f_file);

	aux_write_header_file(nr, nc, 1, 1, argv[5]);

	f_file = fopen(argv[5], "a+b");

	fwrite(DispTarg, sizeof(float), nr*nc, f_file);

	fclose(f_file);

	aux_write_header_file(nr, nc, 1, 1, argv[6]);

	f_file = fopen(argv[6], "a+b");

	fwrite(DispTarg_col, sizeof(float), nr*nc, f_file);

	fclose(f_file);

	aux_write_header_file(nr, nc, 1, 1, argv[7]);

	f_file = fopen(argv[7], "a+b");

	fwrite(DispTarg_row, sizeof(float), nr*nc, f_file);

	fclose(f_file);
	

	return 0;

}