/* View warping. */

#ifndef WARPING_FUN_HH
#define WARPING_FUN_HH

#include <cmath>

struct view{

	unsigned short *color;
	unsigned short *depth;

	int r, c;

	int nr, nc;

	float y, x;

	int n_references;

	int *references;

	int NCOEFFS; // LS merging number of coefficients

	signed short *merge_weights;
	signed short *sparse_weights;

	float *merge_weights_float;
	float *sparse_weights_float;

	unsigned short *residual_color;
	unsigned short *residual_depth;

	bool *bmask;
	unsigned short *seg_vp;

};

void holefilling(unsigned short *pshort, const int ncomps, const int nr, const int nc, const unsigned short maskval)
{
	std::vector<unsigned short> neighbours;
	int dsz = 1;

	for (int ii = 0; ii < nr*nc; ii++){
		int y, x;
		y = ii%nr;
		x = floor(ii / nr);

		bool is_hole = (pshort[ii] == maskval);
		if (ncomps == 3){
			is_hole = is_hole && (pshort[ii + nr*nc] == maskval) && (pshort[ii + 2 * nr*nc] == maskval);
		}

		if (is_hole){
			for (int icomp = 0; icomp < ncomps; icomp++){
				neighbours.clear();
				for (int dy = -dsz; dy <= dsz; dy++){
					for (int dx = -dsz; dx <= dsz; dx++){
						if (!(dy == 0 && dx == 0)){
							if ((y + dy) >= 0 && (y + dy) < nr && (x + dx) >= 0 && (x + dx) < nc){
								bool is_hole_again = (pshort[y + dy + (x + dx)*nr] == maskval);
								if (ncomps == 3){
									is_hole_again = is_hole_again && (pshort[y + dy + (x + dx)*nr + nr*nc] == maskval) && (pshort[y + dy + (x + dx)*nr + 2 * nr*nc] == maskval);
								}
								if (!is_hole_again){
									neighbours.push_back(pshort[y + dy + (x + dx)*nr + icomp*nr*nc]);
								}
							}
						}
					}
				}
				if (neighbours.size() > 0)
					pshort[ii + icomp*nr*nc] = getMedian(neighbours);
			}
		}
	}
}

void getViewMergingLSWeights_N(view *view0, unsigned short **warpedColorViews, float **DispTargs, const unsigned short *original_intermediate_view)
{

	/* This function puts the LS view merging weights into LSw */

	int n_references = (view0)->n_references;
	int nr = (view0)->nr;
	int nc = (view0)->nc;

	signed short *LScoeffs = (view0)->merge_weights;

	int MMM = pow(2, n_references);
	
	bool *bmask = new bool[MMM * n_references]();
	((view0)->bmask) = bmask;

	for (int ij = 0; ij < MMM; ij++){

		int uu = ij;

		for (int ik = n_references-1; ik >= 0; ik--){

			if (floor(uu / pow(2, ik)) > 0)
			{
				uu = uu - pow(2, ik);
				bmask[ij + ik * MMM] = 1;
			}

		}
	}

	unsigned short *seg_vp = new unsigned short[nr*nc]();
	(view0)->seg_vp = seg_vp;

	int *number_of_pixels_per_region = new int[MMM]();


	for (int ii = 0; ii < nr*nc; ii++){

		unsigned short ci = 0;

		for (int ik = 0; ik < n_references; ik++){
			float *pf = DispTargs[ik];
			if (*(pf + ii)>-1)
				ci = ci + pow(2, ik);
		}

		seg_vp[ii] = ci;



		number_of_pixels_per_region[ci]++;

	}

	//FILE *tmpf;
	//tmpf = fopen("G:/HEVC_HDCA/Berlin_verification_model/TMP_CPP/segvp.short", "wb");
	//fwrite(seg_vp, sizeof(unsigned short), 1080 * 1920, tmpf);
	//fclose(tmpf);


	///* debug */
	////BUGAAA
	//aux_write16PGMPPM("G:/HEVC_HDCA/Berlin_verification_model/TMP_CPP/seg_vp.pgm", nc, nr, 1, seg_vp);

	/* go through all regions, collect regressors from references */
	unsigned short **reference_view_pixels_in_classes = new unsigned short*[MMM * n_references]();

	/* also collect desired values */
	unsigned short **original_view_in_classes = new unsigned short*[MMM]();

	for (int ij = 1; ij < MMM; ij++){ // region

		int NN = number_of_pixels_per_region[ij];

		original_view_in_classes[ij] = new unsigned short[NN * 3]();
		unsigned short *p3s = original_view_in_classes[ij];

		int jj = 0;

		for (int ii = 0; ii < nr*nc; ii++){
			if (seg_vp[ii] == ij){
				for (int icomp = 0; icomp < 3; icomp++)
					*(p3s + jj + icomp*NN) = *(original_intermediate_view + ii + icomp*nr*nc);
				jj++;
			}
		}

		for (int ik = 0; ik < n_references; ik++){ // reference view

			if (bmask[ij + ik * MMM]){

				/* allocate pixels for region */
				reference_view_pixels_in_classes[ij + ik * MMM] = new unsigned short[NN * 3]();

				unsigned short *ps = reference_view_pixels_in_classes[ij + ik * MMM];
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

	unsigned short *thetas = new unsigned short[MMM * n_references]();

	for (int ij = 1; ij < MMM; ij++){

		/* form A for this class, compute A'*A (phi)
		also compute A'*y (psi), where y is the desired data from the original view */

		int M = 0;

		for (int ik = 0; ik<n_references; ik++)
		{
			if (bmask[ij + MMM * ik])
				M++;
		}

		int N = number_of_pixels_per_region[ij] * 3; // number of rows in A

		double *A = new double[N*M]();
		double *Yd = new double[N]();

		unsigned short *ps;

		int ikk = 0;

		for (int ik = 0; ik < n_references; ik++){
			if (bmask[ij + ik * MMM]){
				ps = reference_view_pixels_in_classes[ij + ik * MMM];
				for (int ii = 0; ii < N; ii++){
					*(A + ii + ikk*N) = ((double)*(ps + ii)) / ( pow(2, BIT_DEPTH)-1 );
				}
				ikk++;
			}
		}

		ps = original_view_in_classes[ij];

		for (int ii = 0; ii < N; ii++){
			*(Yd + ii) = ((double)*(ps + ii)) / ( pow(2, BIT_DEPTH)-1);
		}

		/* fastols */

		int *PredRegr0 = new int[M]();
		double *PredTheta0 = new double[M]();

		//int Mtrue = FastOLS(ATA, ATYd, YdTYd, PredRegr0, PredTheta0, M, M, M);

		int Mtrue = FastOLS_new(A, Yd, PredRegr0, PredTheta0, M, M, M, N);

		/* establish the subset of reference views available for class */
		int *iks = new int[M]();
		int ee = 0;
		for (int ik = 0; ik < n_references; ik++){
			if (bmask[ij + ik * MMM]){
				*(iks + ee) = ik;
				ee++;
			}
		}

		for (int ii = 0; ii < M; ii++){
			thetas[ij + MMM * iks[PredRegr0[ii]]] = (signed short)round(*(PredTheta0 + ii)*pow(2, 14));
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
	for (int ik = 0; ik < n_references; ik++){
		for (int ij = 0; ij < MMM; ij++){
			if (bmask[ij + ik * MMM]){
				LScoeffs[in] = thetas[ij + MMM * ik];
				in++;
			}
		}
	}

	delete[](thetas);

	for (int ij = 0; ij < MMM; ij++){

		delete[](original_view_in_classes[ij]);

		for (int ik = 0; ik < n_references; ik++){
			delete[](reference_view_pixels_in_classes[ij + MMM * ik]);
		}

	}

	delete[](original_view_in_classes);
	delete[](reference_view_pixels_in_classes);
	//delete[](seg_vp);
	delete[](number_of_pixels_per_region);

	//delete[](bmask);
}

void mergeWarped_N(unsigned short **warpedColorViews, float **DispTargs, view *view0, const int ncomponents)
{

	int MMM = pow(2, (view0)->n_references);

	bool *bmask = (view0)->bmask;
	
	int uu = 0;

	for (int ii = 0; ii < MMM * (view0)->n_references; ii++){
		if (bmask[ii]){
			(view0)->merge_weights_float[ii] = ((float)(view0)->merge_weights[uu++]) / pow(2.0, 14.0);
			//std::cout << (view0)->merge_weights_float[ii] << "\n";
		}
		else{
			(view0)->merge_weights_float[ii] = 0.0;
		}
	}

	int nr = (view0)->nr;
	int nc = (view0)->nc;
	int n_views = (view0)->n_references;
	float *LSw = (view0)->merge_weights_float;

	unsigned short *seg_vp = (view0)->seg_vp;

	float *AA1 = new float[nr*nc*ncomponents]();
	unsigned short *AA2 = new unsigned short[nr*nc*ncomponents]();

	for (int ii = 0; ii < nr*nc; ii++){

		int ci = seg_vp[ii];

		for (int ik = 0; ik < n_views; ik++){
			unsigned short *ps = warpedColorViews[ik];
			for (int icomp = 0; icomp < 3; icomp++){
				AA1[ii + icomp*nr*nc] = AA1[ii + icomp*nr*nc] + LSw[ci + ik * MMM] * ((float)(*(ps + ii + icomp*nr*nc)));
			}
		}

		for (int icomp = 0; icomp < 3; icomp++){
			if (AA1[ii + icomp*nr*nc] < 0)
				AA1[ii + icomp*nr*nc] = 0;
			if (AA1[ii + icomp*nr*nc]>(pow(2, BIT_DEPTH) - 1))
				AA1[ii + icomp*nr*nc] = (pow(2, BIT_DEPTH) - 1);

			AA2[ii + icomp*nr*nc] = (unsigned short)(floor(AA1[ii + icomp*nr*nc]));

		}
	}

	memcpy((view0)->color, AA2, sizeof(unsigned short)*nr*nc*ncomponents);

	delete[](seg_vp);
	delete[](bmask);
	delete[](AA1);
	delete[](AA2);

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

void warpView0_to_View1(view *view0, view *view1, unsigned short *&warpedColor, unsigned short *&warpedDepth, float *&DispTarg)
{

	float ddy = view0->y - view1->y;
	float ddx = view0->x - view1->x;

	unsigned short *AA1 = view0->color;
	unsigned short *DD1 = view0->depth;

	warpedColor = new unsigned short[view0->nr*view0->nc * 3]();
	warpedDepth = new unsigned short[view0->nr*view0->nc]();
	DispTarg = new float[view0->nr*view0->nc]();

	for (int ij = 0; ij < view0->nr*view0->nc; ij++){
		DispTarg[ij] = -1.0;
	}

	for (int ij = 0; ij < view0->nr*view0->nc; ij++)
	{

		float disp = ((float)DD1[ij]) / pow(2,D_DEPTH);
		float DM_COL = disp*ddx;
		float DM_ROW = -disp*ddy;

		float disp0 = abs(DM_COL) + abs(DM_ROW);

		int ix = ij % view0->nr; //row
		int iy = (ij - ix) / view0->nr; //col

		int iynew = iy + (int)round(DM_COL);
		int ixnew = ix + (int)round(DM_ROW);

		if (iynew >= 0 && iynew < view0->nc && ixnew >= 0 && ixnew < view0->nr){
			int indnew = ixnew + iynew*view0->nr;
			if (DispTarg[indnew] < disp0){
				DispTarg[indnew] = disp0;
				warpedColor[indnew] = AA1[ij];
				warpedColor[indnew + view0->nr*view0->nc] = AA1[ij + view0->nr*view0->nc];
				warpedColor[indnew + 2 * view0->nr*view0->nc] = AA1[ij + 2 * view0->nr*view0->nc];
				warpedDepth[indnew] = DD1[ij];
			}
		}

	}

}

#endif