/* View warping. */

#ifndef WARPING_FUN_HH
#define WARPING_FUN_HH

#include <cmath>
#include "float.h"

struct view{

	unsigned short *color;
	unsigned short *depth;

	unsigned short *segmentation;

	int r, c; // SAI subscript

	int nr, nc; // image height, width

	float y, x; // camera displacement

	int min_inv_d; // needed only if inverse depth has negative values, [0,max]-mind = [-mind,max-mind]

	int n_references, n_depth_references;

	int *references, *depth_references; /* depth references not necessarily the same as color references
													  we can have, for example, depth warping only from the externally obtained depth but we still
													  warp color from neighbors that don't have depth provided. We don't want to propagate depth errors from
													  badly warped depth views, thus we restrict depth warping to some high quality subset (usually meaning the 
													  externally obtained high quality depth maps)*/

	signed short *merge_weights;
	int32_t *sparse_weights;

	unsigned char *sparse_mask;

	float *merge_weights_float;

	int *number_of_pixels_per_region;

	bool *bmask; /* view mask for merging weights */
	unsigned short *seg_vp; /* class segmentation, used for view merging weights */
	int NB;

	float residual_rate_color;
	float residual_rate_depth;

	float stdd;

	int NNt, Ms; //for global sparse, NNt defines the neighborhood size [ -NNt:NNt,-NNt:NNt ], Ms is the filter order

	int has_segmentation;
	int maxL; // number of regions in segmentation

	int ****region_displacements; /* region displacement vectors [iar][iac][iR][xy], e.g., [13][13][25][2], for 13x13 angular views with 25 regions for segmentation */

	char path_input_pgm[1024], path_input_ppm[1024], path_input_seg[1024];
	char path_out_pgm[1024], path_out_ppm[1024];

	float *DM_ROW, *DM_COL; /* for lenslet with region displacement vectors */

	int i_order; /* view position in encoding configuration */

};

void initView(view* view)
{

	view->color = NULL;
	view->depth = NULL;
	view->segmentation = NULL;

	view->r = 0;
	view->c = 0;

	view->min_inv_d = 0;
	view->n_references = 0;

	view->references = NULL;
	view->depth_references = NULL;

	view->merge_weights = NULL;
	view->sparse_weights = NULL;
	view->sparse_mask = NULL;
	view->merge_weights_float = NULL;
	view->number_of_pixels_per_region = NULL;

	view->bmask = NULL;
	view->seg_vp = NULL;

	view->NNt = 0;
	view->Ms = 0;

	view->has_segmentation = 0;
	view->maxL = 0;

	view->DM_ROW = NULL;
	view->DM_COL = NULL;

	view->i_order = 0;


}

float getPSNR(FILE *fileout, const char *path_out_ppm, const char *path_input_ppm, const char *difftest_call)
{

	if (fileout == NULL)
		fileout = stdout;

	/* run psnr here */

	char psnr_call[1024];
	sprintf(psnr_call, "%s%s%s%s", difftest_call, path_out_ppm, " ", path_input_ppm);

	FILE *pfile;
	pfile = _popen(psnr_call, "r");

	char psnr_buffer[1024];
	while (fgets(psnr_buffer, sizeof(psnr_buffer), pfile) != 0) {
		/*...*/
	}
	_pclose(pfile);

	char tmp_char[1024];
	float psnr_value = 0;

	sscanf(psnr_buffer, "%s\t%f", tmp_char, &psnr_value);

	fprintf(fileout,"%s\n", psnr_buffer);

	return psnr_value;
}

void applyGlobalSparseFilter(view *view0){

	unsigned char *Regr0 = view0->sparse_mask;
	int32_t *theta0 = view0->sparse_weights;
	int Ms = view0->Ms;
	int NNt = view0->NNt;
	int nr = view0->nr, nc = view0->nc;


	float *theta = new float[(2 * NNt + 1)*(2 * NNt + 1) + 1]();

	for (int ii = 0; ii < Ms; ii++){
		if (Regr0[ii] > 0){
			theta[Regr0[ii] - 1] = ((float)theta0[ii]) / pow(2, BIT_DEPTH_SPARSE);
			//std::cout << theta[Regr0[ii] - 1] << "\t";
		}
	}

	float *final_view = new float[nr*nc * 3]();

	unsigned short *pshort = view0->color;

	for (int ii = 0; ii < nr*nc * 3; ii++)
		final_view[ii] = pshort[ii];

	for (int rr = NNt; rr < nr - NNt; rr++){
		for (int cc = NNt; cc < nc - NNt; cc++)
		{
			for (int icomp = 0; icomp < 3; icomp++)
				final_view[rr + cc*nr + icomp*nr*nc] = 0;

			int ee = 0;

			for (int dy = -NNt; dy <= NNt; dy++){
				for (int dx = -NNt; dx <= NNt; dx++){
					for (int icomp = 0; icomp < 3; icomp++){
						final_view[rr + cc*nr + icomp*nr*nc] = final_view[rr + cc*nr + icomp*nr*nc] + theta[ee] * ((float)pshort[rr + dy + (cc + dx)*nr + icomp*nr*nc]);
					}
					ee++;
				}
			}

			/* bias term */
			for (int icomp = 0; icomp < 3; icomp++){
				final_view[rr + cc*nr + icomp*nr*nc] = final_view[rr + cc*nr + icomp*nr*nc] + theta[(2 * NNt + 1)*(2 * NNt + 1)];
			}

		}
	}

	//unsigned short *final_view_s = new unsigned short[nr*nc * 3]();

	for (int ii = 0; ii < nr*nc * 3; ii++){
		if (final_view[ii] < 0)
			final_view[ii] = 0;
		if (final_view[ii]> (pow(2, BIT_DEPTH) - 1))
			final_view[ii] = (pow(2, BIT_DEPTH) - 1);

		pshort[ii] = (unsigned short)(final_view[ii]);
	}

	delete[](theta);
	delete[](final_view);

}

void getGlobalSparseFilter(view *view0, unsigned short *original_color_view)
{

	int nr = view0->nr;
	int nc = view0->nc;
	int NNt = view0->NNt;
	int Ms = view0->Ms;

	unsigned char *Regr0 = new unsigned char[Ms]();
	int32_t *theta0 = new int32_t[Ms]();

	view0->sparse_weights = theta0;
	view0->sparse_mask = Regr0;

	int Npp = (nr - NNt * 2)*(nc - NNt * 2) * 3;
	int Npp0 = Npp / 3;

	double *AA = new double[Npp*((NNt * 2 + 1)*(NNt * 2 + 1) + 1)]();
	double *Yd = new double[Npp]();

	for (int ii = 0; ii < Npp; ii++)
		*(AA + ii + (NNt * 2 + 1)*(NNt * 2 + 1)*Npp) = 1.0;

	int iiu = 0;

	unsigned short *pshort = view0->color;

	for (int ir = NNt; ir < nr - NNt; ir++){
		for (int ic = NNt; ic < nc - NNt; ic++){
			int ai = 0;
			for (int dy = -NNt; dy <= NNt; dy++){
				for (int dx = -NNt; dx <= NNt; dx++){
					for (int icomp = 0; icomp < 3; icomp++){

						int offset = ir + dy + nr*(ic + dx) + icomp*nr*nc;

						/* get the desired Yd*/
						if (dy == 0 && dx == 0){
							*(Yd + iiu + icomp*Npp0) = ((double)*(original_color_view + offset)) / (pow(2, BIT_DEPTH) - 1);
						}

						/* get the regressors */
						*(AA + iiu + icomp*Npp0 + ai*Npp) = ((double)*(pshort + offset)) / (pow(2, BIT_DEPTH) - 1);

					}
					ai++;
				}
			}
			iiu++;
		}
	}

	int *PredRegr0 = new int[(2 * NNt + 1)*(2 * NNt + 1) + 1]();
	double *PredTheta0 = new double[(2 * NNt + 1)*(2 * NNt + 1) + 1]();

	int Mtrue = FastOLS_new(AA, Yd, PredRegr0, PredTheta0, Ms, (NNt * 2 + 1)*(NNt * 2 + 1) + 1, (NNt * 2 + 1)*(NNt * 2 + 1) + 1, Npp);

	delete[](AA);
	delete[](Yd);

	for (int ii = 0; ii < Ms; ii++){
		*(Regr0 + ii) = ((unsigned char)*(PredRegr0 + ii) + 1);
		*(theta0 + ii) = (int32_t)round(*(PredTheta0 + ii) * pow(2, BIT_DEPTH_SPARSE));
	}
}

void holefilling(unsigned short *pshort, const int ncomps, const int nr, const int nc, const unsigned short maskval)
{
	std::vector<unsigned short> neighbours;
	int dsz = 1;

	for (int ii = 0; ii < nr*nc; ii++){
		int y, x;
		y = ii%nr;
		x = (ii / nr);

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

void setBMask(view *view0) 
{

	/* sets the binary mask used to derive view availability in each of the MMM classes,
	size of the binary mask is [MMM x n_references] */

	int MMM = pow(2, view0->n_references);

	bool *bmask = new bool[MMM * view0->n_references]();

	(view0)->bmask = bmask;

	for (int ij = 0; ij < MMM; ij++) {

		int uu = ij;

		for (int ik = view0->n_references - 1; ik >= 0; ik--) {

			if (floor(uu / pow(2, ik)) > 0)
			{
				uu = uu - pow(2, ik);
				bmask[ij + ik * MMM] = 1;
			}

		}
	}
}

void getViewMergingLSWeights_N(view *view0, unsigned short **warpedColorViews, float **DispTargs, const unsigned short *original_color_view)
{

	/* This function puts the LS view merging weights into LSw */

	int n_references = (view0)->n_references;
	int nr = (view0)->nr;
	int nc = (view0)->nc;

	signed short *LScoeffs = (view0)->merge_weights;

	int MMM = pow(2, n_references);
	
	bool *bmask = (view0)->bmask;

	unsigned short *seg_vp = (view0)->seg_vp;

	/* go through all regions, collect regressors from references */
	unsigned short **reference_view_pixels_in_classes = new unsigned short*[MMM * n_references]();

	/* also collect desired values */
	unsigned short **original_view_in_classes = new unsigned short*[MMM]();

	int *number_of_pixels_per_region = (view0)->number_of_pixels_per_region;

	for (int ij = 1; ij < MMM; ij++){ // region

		int NN = number_of_pixels_per_region[ij];

		original_view_in_classes[ij] = new unsigned short[NN * 3]();
		unsigned short *p3s = original_view_in_classes[ij];

		int jj = 0;

		for (int ii = 0; ii < nr*nc; ii++){
			if (seg_vp[ii] == ij){
				for (int icomp = 0; icomp < 3; icomp++)
					*(p3s + jj + icomp*NN) = *(original_color_view + ii + icomp*nr*nc);
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

	signed short *thetas = new signed short[MMM * n_references]();

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
			thetas[ij + MMM * iks[PredRegr0[ii]]] = (signed short)floor(*(PredTheta0 + ii)*pow(2, BIT_DEPTH_MERGE)+0.5);
		}

		delete[](iks);

		delete[](PredRegr0);
		delete[](PredTheta0);
		
		

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
	//delete[](number_of_pixels_per_region);

}

void initSegVp(view *view0, float **DispTargs) {

	int nr = view0->nr;
	int nc = view0->nc;
	int n_references = view0->n_references;

	unsigned short *seg_vp = new unsigned short[nr*nc]();

	(view0)->seg_vp = seg_vp;

	int MMM = pow(2, (view0)->n_references);

	int *number_of_pixels_per_region = new int[MMM]();

	for (int ii = 0; ii < nr*nc; ii++) {

		unsigned short ci = 0;

		for (int ik = 0; ik < n_references; ik++) {
			float *pf = DispTargs[ik];
			if (*(pf + ii)>-1)
				ci = ci + pow(2, ik);
		}

		seg_vp[ii] = ci;



		number_of_pixels_per_region[ci]++;

	}

	view0->number_of_pixels_per_region = number_of_pixels_per_region;
}

void initViewW(view *view0, float **DispTargs) {

	/* sets some of the parameters for a view in the light view structure */

	(view0)->NB = (pow(2, (view0)->n_references)*(view0)->n_references);

	signed short *merge_weights = new signed short[(view0)->NB / 2]();

	(view0)->merge_weights = merge_weights;

	float *merge_weights_float = new float[(view0)->NB]();

	(view0)->merge_weights_float = merge_weights_float;

	setBMask(view0);

	initSegVp(view0, DispTargs);

}

void getGeomWeight(view *view0, view *LF, float stdd) {

	int MMM = pow(2, (view0)->n_references);

	float *thetas = new float[view0->NB]();

	bool *bmask = (view0)->bmask;

	for (int ii = 0; ii < MMM; ii++) {
		float sumw = 0;

		for (int ij = 0; ij < (view0)->n_references; ij++) {
			view *view1 = LF + (view0)->references[ij];
			float vdistance = (view0->x - view1->x)*(view0->x - view1->x) + (view0->y - view1->y)*(view0->y - view1->y);
			if (bmask[ii + ij * MMM]) {
				thetas[ii + ij * MMM] = exp(-(vdistance) / (2 * stdd*stdd));
				sumw += thetas[ii + ij * MMM];
			}
		}
		for (int ij = 0; ij < (view0)->n_references; ij++)
			thetas[ii + ij * MMM] = thetas[ii + ij * MMM] / sumw;
	}

	/* columnwise collecting of thetas */
	signed short *LScoeffs = (view0)->merge_weights;
	int in = 0;
	for (int ik = 0; ik < (view0)->n_references; ik++) {
		for (int ij = 0; ij < MMM; ij++) {
			if (bmask[ij + ik * MMM]) {
				LScoeffs[in] = (signed short)floor(thetas[ij + MMM * ik] * pow(2, BIT_DEPTH_MERGE) + 0.5);
				in++;
			}
		}
	}
	delete[](thetas);

	//view0->stdd = stdd;
}

void mergeWarped_N(unsigned short **warpedColorViews, float **DispTargs, view *view0, const int ncomponents)
{

	int MMM = pow(2, (view0)->n_references);

	bool *bmask = (view0)->bmask;
	
	int uu = 0;

	for (int ii = 0; ii < MMM * (view0)->n_references; ii++){
		if (bmask[ii]){
			(view0)->merge_weights_float[ii] = ((float)(view0)->merge_weights[uu++]) / pow(2, BIT_DEPTH_MERGE);
			std::cout << (view0)->merge_weights_float[ii] << "\n";
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

	//delete[](seg_vp);
	//delete[](bmask);
	delete[](AA1);
	delete[](AA2);

}

void warpView0_to_View1(view *view0, view *view1, unsigned short *&warpedColor, unsigned short *&warpedDepth, float *&DispTarg)
{

	/*this function forward warps from view0 to view1 for both color and depth*/

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


		float disp = ((float)DD1[ij] - (float)view0->min_inv_d) / pow(2, D_DEPTH);
		float DM_COL = disp*ddx;
		float DM_ROW = -disp*ddy;
		float disp0 = abs(DM_COL) + abs(DM_ROW);

		//if (view0->DM_ROW != NULL && view0->DM_COL != NULL)
		//{
		//	DM_COL = view0->DM_COL[ij];
		//	DM_ROW = view0->DM_ROW[ij];
		//	disp0 = 1 / (abs(DM_COL) + abs(DM_ROW)); /*for lenslet large disparity means further away ... */
		//}

		int ix = ij % view0->nr; //row
		int iy = (ij - ix) / view0->nr; //col

		int iynew = iy + (int)floor(DM_COL + 0.5);
		int ixnew = ix + (int)floor(DM_ROW + 0.5);

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