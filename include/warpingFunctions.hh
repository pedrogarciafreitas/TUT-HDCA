/* View warping. */

#ifndef WARPING_FUN_HH
#define WARPING_FUN_HH

#include <cmath>

#ifndef BIT_DEPTH
#define BIT_DEPTH 10
#endif

void viewMergingLSWeights(unsigned short *warpedColorViews[5], float *DispTargs[5], const unsigned short* original_intermediate_view, 
	const int nr, const int nc, const bool *bmask, signed short *LScoeffs)
{

	/* This function puts the LS view merging weights into LSw */

	unsigned short *seg_vp = new unsigned short[nr*nc]();

	//unsigned short *seg_vp2 = new unsigned short[nr*nc * 3]();

	int *number_of_pixels_per_region = new int[32]();


	for (int ii = 0; ii < nr*nc; ii++){

		int ci = 0;

		for (int ik = 0; ik < 5; ik++){
			float *pf = DispTargs[ik];
			if (*(pf + ii)>-1)
				ci = ci + pow(2, ik);
		}

		seg_vp[ii] = ci;

		//seg_vp2[ii] = ci;
		//seg_vp2[ii + nr*nc] = ci;
		//seg_vp2[ii + 2 * nr*nc] = ci;

		number_of_pixels_per_region[ci]++;

	}

	/* debug */
	//FILE *filept1;
	//filept1 = fopen("G:/HEVC_HDCA/Berlin_verification_model/Encoder_1st_tests/seg_vp", "wb");
	//fwrite(seg_vp, sizeof(unsigned short), nr*nc, filept1);
	//fclose(filept1);

	
	/* go through all regions, collect regressors from references */
	unsigned short **reference_view_pixels_in_classes = new unsigned short*[32 * 5]();

	/* also collect desired values */
	unsigned short **original_view_in_classes = new unsigned short*[32]();

	for (int ij = 1; ij < 32; ij++){ // region

		int NN = number_of_pixels_per_region[ij];

		original_view_in_classes[ij] = new unsigned short[NN*3]();
		unsigned short *p3s = original_view_in_classes[ij];

		int jj = 0;

		for (int ii = 0; ii < nr*nc; ii++){
			if (seg_vp[ii] == ij){
				for (int icomp = 0; icomp < 3;icomp++)
					*(p3s + jj + icomp*NN) = *(original_intermediate_view+ii+icomp*nr*nc);
				jj++;
			}
		}

		for (int ik = 0; ik < 5; ik++){ // reference view

			if (bmask[ij + ik * 32]){

				/* allocate pixels for region */
				reference_view_pixels_in_classes[ij + ik * 32] = new unsigned short[NN * 3]();

				unsigned short *ps = reference_view_pixels_in_classes[ij + ik * 32];
				unsigned short *pss = warpedColorViews[ik];		

				
				jj = 0;

				for (int ii = 0; ii < nr*nc; ii++){
					if (seg_vp[ii] == ij){
						for (int icomp = 0; icomp < 3; icomp++)
							*(ps + jj + icomp*NN) = *(pss + ii + icomp*nr*nc);
						jj++;
					}
				}

			}

		}

	}

	
	/* run fastOLS on the classes */

	unsigned short *thetas = new unsigned short[32 * 5]();

	for (int ij = 1; ij < 32; ij++){

		/* form A for this class, compute A'*A (phi)
		also compute A'*y (psi), where y is the desired data from the original view */

		int M = 0;

		for( int ik = 0; ik<5; ik++ )
		{
			if (bmask[ij + 32 * ik])
				M++;
		}

		int N = number_of_pixels_per_region[ij] * 3; // number of rows in A

		double *A = new double[N*M]();
		double *Yd = new double[N]();

		unsigned short *ps;

		int ikk = 0;

		for (int ik = 0; ik < 5; ik++){
			if (bmask[ij + ik * 32]){
				ps = reference_view_pixels_in_classes[ij + ik * 32];
				for (int ii = 0; ii < N; ii++){
					*(A + ii + ikk*N) = ((double)*(ps + ii)) / pow(2, BIT_DEPTH);
				}
				ikk++;
			}
		}

		ps = original_view_in_classes[ij];

		for (int ii = 0; ii < N; ii++){
			*(Yd + ii) = ((double)*(ps + ii)) / pow(2, BIT_DEPTH);
		}

		//double *ATA = new double[M*M]();
		//double *ATYd = new double[M]();

		///* make ATA */
		//for (int i1 = 0; i1 < M; i1++){
		//	for (int j1 = 0; j1 < M; j1++){
		//		for (int ii = 0; ii < N; ii++){
		//			*(ATA + i1 + j1*M) += (*(A + ii + i1*N))*(*(A + ii + j1*N));
		//		}
		//	}
		//}

		////for (int i1 = 0; i1 < M; i1++){
		////	for (int j1 = 0; j1 < M; j1++){
		////		std::cout << *(ATA + i1 + j1*M) << "\n";
		////	}
		////}

		////std::cout << "------------------------------------------------\n";

		///* make ATYd */
		//for (int i1 = 0; i1 < M; i1++){
		//	for (int ii = 0; ii < N; ii++){
		//		*(ATYd + i1) += (*(A + ii + i1*N))*(*(Yd + ii));
		//	}
		//}

		///* YdTYd */
		//double YdTYd = 0;
		//for (int ii = 0; ii < N; ii++)
		//	YdTYd += (*(Yd + ii))*(*(Yd + ii));

		/* fastols */

		int *PredRegr0 = new int[M]();
		double *PredTheta0 = new double[M]();

		//int Mtrue = FastOLS(ATA, ATYd, YdTYd, PredRegr0, PredTheta0, M, M, M);

		int Mtrue = FastOLS_new(A, Yd, PredRegr0, PredTheta0, M, M, M, N);

		/* establish the subset of reference views available for class */
		int *iks = new int[M]();
		int ee = 0;
		for (int ik = 0; ik < 5; ik++){
			if (bmask[ij + ik * 32]){
				*(iks + ee) = ik;
				ee++;
			}
		}

		for (int ii = 0; ii < M; ii++){
			thetas[ij + 32 * iks[PredRegr0[ii]]] = (signed short)round(*(PredTheta0 + ii)*pow(2, 14));
		}

		delete[](iks);

		delete[](PredRegr0);
		delete[](PredTheta0);

		//delete[](ATA);
		//delete[](ATYd);
		delete[](A);
		delete[](Yd);

		

	}

	/* columnwise collecting of thetas */
	int in = 0;
	for (int ik = 0; ik < 5; ik++){
		for (int ij = 0; ij < 32; ij++){
			if (bmask[ij + ik * 32]){
				LScoeffs[in] = thetas[ij + 32 * ik];
				in++;
			}
		}
	}

	delete[](thetas);

	for (int ij = 0; ij < 32; ij++){

		delete[](original_view_in_classes[ij]);

		for (int ik = 0; ik < 5; ik++){
			delete[](reference_view_pixels_in_classes[ij+32*ik]);
		}

	}

	delete[](original_view_in_classes);
	delete[](reference_view_pixels_in_classes);
	delete[](seg_vp);
	//delete[](seg_vp2);
	delete[](number_of_pixels_per_region);
}


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
			if (AA1[ii + icomp*nr*nc]>(pow(2, BIT_DEPTH) - 1))
				AA1[ii + icomp*nr*nc] = (pow(2, BIT_DEPTH) - 1);

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

#endif