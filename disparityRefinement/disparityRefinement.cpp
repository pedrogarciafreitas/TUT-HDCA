/* Disparity refinement over 11x11 window. */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "gen_types.hh"


int main(const int argc, const char** argv) {

	int nr, nc, ncomponents, nvc;

	int ii0 = atoi(argv[5]), jj0 = atoi(argv[6]), ii1 = atoi(argv[7]), jj1 = atoi(argv[8]);

	aux_read_header_file(&nr, &nc, &ncomponents, &nvc, argv[1]);

	/* reference */
	unsigned short *AA1 = new unsigned short[nr*nc*ncomponents];

	aux_read_file_uint16(nr, nc, ncomponents, argv[1], AA1);

	/* side view */
	unsigned short *AA3 = new unsigned short[nr*nc*ncomponents];

	aux_read_file_uint16(nr, nc, ncomponents, argv[2], AA3);

	/* quantized labels */
	int *St1 = new int[nr*nc];

	aux_read_file_int32(nr, nc, 1, argv[3], St1);

	/* quantized disparity */
	float *quantDM = new float[nr*nc];

	aux_read_file_float(nr, nc, 1, argv[4], quantDM);

	int rowsreg[512];
	int colsreg[512];

	printf("%d %d %d \n", nr, nc, ncomponents);

	for (int is = 1; is < 512; is++)
	{

		std::vector<int> indi;
		std::vector<float> MAE;
		std::vector<int> irows;
		std::vector<int> icols;

		for (int ij = 0; ij < nr*nc; ij++)
		{
			//printf("%f\t", quantDM[ij]);
			if (St1[ij] == is)
				indi.push_back(ij);
		}

		

		std::vector<int> new_indi;
		std::vector<int> iynews, ixnews;

		for (int ij = 0; ij < indi.size(); ij++){
			int ix = indi.at(ij) % nr; //row
			int iy = ( indi.at(ij) - ix )/nr; //col

			int iynew = iy + (int)round(((float)(ii1 - ii0))*quantDM[indi.at(ij)]);
			int ixnew = ix + (int)round(((float)(jj1 - jj0))*quantDM[indi.at(ij)]*(60.0/40.0));

			//int iynew = (int)iynewf;
			//int ixnew = (int)ixnewf;

			if (iynew>=0 && iynew<nc && ixnew>=0 && ixnew < nr){
				new_indi.push_back(indi.at(ij));
				iynews.push_back(iynew);
				ixnews.push_back(ixnew);
			}
		}


		int best_irow = 0, best_icol = 0;
		float best_MAE = 10000000000;

		int NNL = 6;

		for (int icol = -NNL; icol < NNL + 1; icol++){
			for (int irow = -NNL; irow < NNL + 1; irow++){
				float tmp_MAE = 0;
				for (int ij = 0; ij < new_indi.size(); ij++){

					int ixn = ixnews.at(ij) + irow;
					int iyn = iynews.at(ij) + icol;

					if (ixn>=0 && ixn<nr && iyn>=0 && iyn < nc){

						tmp_MAE = tmp_MAE +
							abs((float)AA1[new_indi.at(ij)] - (float)AA3[ixn + nr*iyn])
							+ abs((float)AA1[new_indi.at(ij) + nr*nc] - (float)AA3[ixn + nr*iyn + nr*nc])
							+ abs((float)AA1[new_indi.at(ij) + 2*nr*nc ] - (float)AA3[ixn + nr*iyn + 2*nr*nc]);
					}
					//printf("%f\n", tmp_MAE);
				}
				if ( (tmp_MAE < best_MAE) && tmp_MAE>0 ){
					best_MAE = tmp_MAE;
					best_irow = irow;
					best_icol = icol;
				}
			}
		}

		rowsreg[is] = best_irow;
		colsreg[is] = best_icol;

		printf("%d\t%f\t%d\t%d\t%d\t%d\n", is, best_MAE, best_irow, best_icol, indi.size(), new_indi.size());
	}


	int *COLS = new int[nr*nc];
	int *ROWS = new int[nr*nc];

	for (int is = 1; is < 512; is++)
	{
		for (int ij = 0; ij < nr*nc; ij++)
		{
			if (St1[ij] == is){
				COLS[ij] = colsreg[is];
				ROWS[ij] = rowsreg[is];
			}
		}
	}

	aux_write_header_file(nr, nc, 1, 1, argv[9]);

	FILE *f_file = fopen(argv[9], "a+b");

	fwrite(ROWS, sizeof(int), nr*nc, f_file);

	fclose(f_file);

	aux_write_header_file(nr, nc, 1, 1, argv[10]);

	f_file = fopen(argv[10], "a+b");

	fwrite(COLS, sizeof(int), nr*nc, f_file);

	fclose(f_file);

	return 0;

}