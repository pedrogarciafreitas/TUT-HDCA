/* View warping. */

#include <cmath>


void collectWarpedLS(unsigned short *warpedColorViews[5], float *DispTargs[5], const int nr, const int nc, const int ncomponents, const int n_views, const float* LSw){


	unsigned short *seg_vp = new unsigned short[nr*nc]();

	unsigned short *seg_vp2 = new unsigned short[nr*nc*3]();


	for (int ii = 0; ii < nr*nc; ii++){

		int ci = 0;

		for (int ik = 0; ik < n_views; ik++){
			float *pf = DispTargs[ik];
			if( *(pf+ii)>-1 )
				ci = ci + pow(2, ik);
		}

		seg_vp[ii] = ci;

		seg_vp2[ii] = ci;
		seg_vp2[ii + nr*nc] = ci;
		seg_vp2[ii + 2 * nr*nc] = ci;

	}

	float *AA1 = new float[nr*nc*ncomponents]();
	unsigned short *AA2 = new unsigned short[nr*nc*ncomponents]();

	for (int ii = 0; ii < nr*nc; ii++){

		int ci = seg_vp[ii];

		for (int ik = 0; ik < n_views; ik++){
			unsigned short *ps = warpedColorViews[ik];
			for (int icomp = 0; icomp < 3; icomp++){
				AA1[ii + icomp*nr*nc] = AA1[ii + icomp*nr*nc] + LSw[ci + ik * 32]*((float)(*(ps + ii + icomp*nr*nc)));
			}
		}

		for (int icomp = 0; icomp < 3; icomp++){
			if (AA1[ii + icomp*nr*nc] < 0)
				AA1[ii + icomp*nr*nc] = 0;
			if (AA1[ii + icomp*nr*nc]>(pow(2,16)-1))
				AA1[ii + icomp*nr*nc] = (pow(2, 16) - 1);

			AA2[ii + icomp*nr*nc] = (unsigned short)( floor(AA1[ii + icomp*nr*nc]) );

		}
	}

	memcpy(warpedColorViews[0], AA2, sizeof(unsigned short)*nr*nc*ncomponents);

	delete[](seg_vp);
	delete[](AA1);
	delete[](AA2);
	delete[](seg_vp2);
}


void collectWarped(unsigned short *warpedColorViews[5], float *DispTargs[5], const int nr, const int nc, const int ncomponents, const int n_views){


	//unsigned short *AA1 = new unsigned short[nr*nc*ncomponents];
	//memcpy(AA1, warpedColorViews[0], nr*nc*ncomponents);

	//float *DispTarg1 = new float[nr*nc];
	//memcpy(DispTarg1, DispTargs[0], nr*nc);

	unsigned short *AA1 = warpedColorViews[0];
	float *DispTarg1 = DispTargs[0];

	unsigned short *AA2 = new unsigned short[nr*nc*ncomponents]();
	float *DispTarg2 = new float[nr*nc]();

	for (int ik = 1; ik < n_views; ik++){

		memcpy(AA2, warpedColorViews[ik], sizeof(unsigned short)*nr*nc*ncomponents);
		memcpy(DispTarg2, DispTargs[ik], sizeof(float)*nr*nc);

		for (int ij = 0; ij < nr*nc; ij++){
			if (DispTarg1[ij] < 0 && DispTarg2[ij] >= 0){
				DispTarg1[ij] = DispTarg2[ij];
				AA1[ij] = AA2[ij];
				AA1[ij + nr*nc] = AA2[ij + nr*nc];
				AA1[ij + nr*nc * 2] = AA2[ij + nr*nc * 2];
			}

		}
		//printf("ik:\t%d\n", ik);
	}
	delete[](AA2);
	delete[](DispTarg2);
}

void warpColorView(const unsigned short *AA1, const float *DM_ROW, const float *DM_COL, const int nr, const int nc, unsigned short* Warped, float* DispTarg) {


	int *ROWS = new int[nr*nc]();
	int *COLS = new int[nr*nc]();

	float *DispTarg_col = new float[nr*nc]();
	float *DispTarg_row = new float[nr*nc]();

	for (int ij = 0; ij < nr*nc; ij++){
		DispTarg[ij] = -1;
		DispTarg_col[ij] = 99999;
		DispTarg_row[ij] = 99999;
	}


	for (int ij = 0; ij < nr*nc; ij++)
	{
		float disp0 = abs(DM_COL[ij]) + abs(DM_ROW[ij]);

		int ix = ij % nr; //row
		int iy = (ij - ix) / nr; //col

		int iynew = iy + (int)round(DM_COL[ij]);
		int ixnew = ix + (int)round(DM_ROW[ij]);

		if (iynew >= 0 && iynew < nc && ixnew >= 0 && ixnew < nr){
			int indnew = ixnew + iynew*nr;
			if (DispTarg[indnew] < disp0){
				DispTarg[indnew] = disp0;
				DispTarg_col[indnew] = iynew - iy; // warped horizontal disparity
				DispTarg_row[indnew] = ixnew - ix; // warped vertical disparity
				Warped[indnew] = AA1[ij];
				Warped[indnew + nr*nc] = AA1[ij + nr*nc];
				Warped[indnew + 2 * nr*nc] = AA1[ij + 2 * nr*nc];
			}
		}

	}

	delete[](ROWS);
	delete[](COLS);
	delete[](DispTarg_col);
	delete[](DispTarg_row);
}