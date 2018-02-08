/* View warping. */

#ifndef WARPING_HH
#define WARPING_HH

#include <cmath>
#include <algorithm>

template <class T>
void warpViews(const float d_cen[21][101][2], const int ref_rows[], const int ref_cols[],
	T *input[], float *qDMF[], float *RowDisps[], float *ColDisps[], const int nr, const int nc,
	T *warped[], float *DispTargs[], const int r, const int c, const int nviews, const int ncomponents);

template <class T>
void warpView(const T *input, const float *DM_ROW, const float *DM_COL,
	const int nr, const int nc, T* warped, float* DispTarg, const int ncomponents);

template <class T>
void collectWarped(T *warped[5], float *DispTargs[5], const int nr, const int nc,
	const int ncomponents, const int n_views);

template <class T>
T getMedian(std::vector<T> scores);

template <class T>
void medfilt2D(T* input, T* output, int SZ, int nr, int nc);

template <class T>
T getMedian(std::vector<T> scores)
{
	T median;
	size_t size = scores.size();

	std::sort(scores.begin(), scores.end());

	if (size % 2 == 0)
	{
		median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
	}
	else
	{
		median = scores[size / 2];
	}

	return median;
}

template <class T>
void medfilt2D(T* input, T* output, int SZ, int nr, int nc)
{
	int dsz = floor(SZ / 2);
	std::vector<T> scores;
	for (int y = 0; y < nr; y++){
		for (int x = 0; x < nc; x++){
			scores.clear();
			for (int dy = -dsz; dy < dsz; dy++){
				for (int dx = -dsz; dx < dsz; dx++){
					if ((y + dy) >= 0 && (y + dy) < nr
						&& (x + dx) >= 0 && (x + dx) < nc)
						scores.push_back(input[y + dy + (x + dx)*nr]);
				}
			}
			output[y + x*nr] = getMedian(scores);
		}
	}
}

template<class T>
void fillHoles_T(T *pshort, const int nr, const int nc, const int ncomponents, const int dsz)
{
	std::vector<T> neighbours;

	for (int ii = 0; ii < nr*nc; ii++){
		int y, x;
		y = ii%nr;
		x = floor(ii / nr);
		if (pshort[ii] + pshort[ii + nr*nc] + pshort[ii + 2 * nr*nc] < 1){
			for (int icomp = 0; icomp < ncomponents; icomp++){
				neighbours.clear();
				for (int dy = -dsz; dy < dsz; dy++){
					for (int dx = -dsz; dx < dsz; dx++){
						if ((y + dy) >= 0 && (y + dy) < nr
							&& (x + dx) >= 0 && (x + dx) < nc){
							if (pshort[y + dy + (x + dx)*nr + icomp*nr*nc]>0)
								neighbours.push_back(pshort[y + dy + (x + dx)*nr + icomp*nr*nc]);
						}
					}
				}
				if (neighbours.size()>0)
					pshort[ii + icomp*nr*nc] = getMedian(neighbours);
			}
		}
	}
}

template <class T>
void warpViews(const float d_cen[21][101][2], int *ref_rows, int *ref_cols,
	T *input[], float *qDMF[], float *RowDisps[], float *ColDisps[], const int nr, const int nc,
	T *warped[], float *DispTargs[], const int r, const int c, const int nviews, const int ncomponents)
{
	
	std::vector<pair<float, int>> view_distances(nviews);
	std::vector<pair<float, float>> dydx(nviews);

	for (int ij = 0; ij < nviews; ij++){

		dydx.at(ij).first = d_cen[ref_rows[ij]][ref_cols[ij]][1] - d_cen[r][c][1];
		dydx.at(ij).second = d_cen[ref_rows[ij]][ref_cols[ij]][0] - d_cen[r][c][0];

		view_distances.at(ij).first = dydx.at(ij).first * dydx.at(ij).first + dydx.at(ij).second * dydx.at(ij).second;
		view_distances.at(ij).second = ij;

	}

	sort(view_distances.begin(), view_distances.end());

	for (int uu = 0; uu < nviews; uu++){

		int ij = view_distances.at(uu).second;

		/* to zero for consistency */
		memset(warped[ij], 0x00, sizeof(T)*nr*nc*ncomponents);
		memset(DispTargs[ij], 0x00, sizeof(float)*nr*nc);
		memset(RowDisps[ij], 0x00, sizeof(float)*nr*nc);
		memset(ColDisps[ij], 0x00, sizeof(float)*nr*nc);
		
		float *p = qDMF[ij];
		float *pp = RowDisps[ij];
		float *ppp = ColDisps[ij];
		for (int ii = 0; ii < nr*nc; ii++){
			pp[ii] = -p[ii] * dydx.at(ij).first;
			ppp[ii] = p[ii] * dydx.at(ij).second;
		}

		warpView(input[ij], RowDisps[ij], ColDisps[ij], nr, nc, warped[ij], DispTargs[ij], ncomponents);

	}
	return;
}

template <class T>
void warpView(const T *input, const float *DM_ROW, const float *DM_COL, const int nr, const int nc, T* warped, float* DispTarg, const int ncomponents)
{

	int *ROWS = new int[nr*nc];
	int *COLS = new int[nr*nc];

	float *DispTarg_col = new float[nr*nc];
	float *DispTarg_row = new float[nr*nc];

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
				for (int icomp = 0; icomp < ncomponents; icomp++){
					warped[indnew + icomp*nr*nc] = input[ij + icomp*nr*nc];
				}
			}
		}

	}

	delete(ROWS);
	delete(COLS);
	delete(DispTarg_col);
	delete(DispTarg_row);

	return;
}

template <class T>
void collectWarped(T *warped[5], float *DispTargs[5], const int nr, const int nc, const int ncomponents, const int n_views)
{

	T *AA1 = warped[0];
	float *DispTarg1 = DispTargs[0];

	T *AA2 = new T[nr*nc*ncomponents];
	float *DispTarg2 = new float[nr*nc];

	for (int ik = 1; ik < n_views; ik++){

		memcpy(AA2, warped[ik], sizeof(T)*nr*nc*ncomponents);
		memcpy(DispTarg2, DispTargs[ik], sizeof(float)*nr*nc);

		for (int ij = 0; ij < nr*nc; ij++){
			if (DispTarg1[ij] < 0 && DispTarg2[ij] >= 0){
				DispTarg1[ij] = DispTarg2[ij];
				for (int icomp = 0; icomp < ncomponents; icomp++){
					AA1[ij + icomp*nr*nc] = AA2[ij + icomp*nr*nc];
				}
			}

		}
	}

	delete(AA2);
	delete(DispTarg2);

	return;
}


//void warpInverseDepths(const float d_cen[21][101][2], const int ref_rows[5], const int ref_cols[5],
//	float *qDMF[5], float *RowDisps[5], float *ColDisps[5], const int nr, const int nc,
//	float *warpedInverseDepths[5], float *DispTargs[5], const int r, const int c)
//{
//	//float view_distances[5];
//	std::vector<pair<float, int>> view_distances(5);
//	float ddy[5];
//	float ddx[5];
//
//	for (int ij = 0; ij < 5; ij++){
//
//		ddy[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][1] - d_cen[r][c][1];
//		ddx[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][0] - d_cen[r][c][0];
//
//		view_distances.at(ij).first = ddy[ij] * ddy[ij] + ddx[ij] * ddx[ij];
//		view_distances.at(ij).second = ij;
//
//	}
//
//	sort(view_distances.begin(), view_distances.end());
//
//	for (int uu = 0; uu < 5; uu++){
//
//		int ij = view_distances.at(ij).second;
//
//		float *p = qDMF[ij];
//		float *pp = RowDisps[ij];
//		float *ppp = ColDisps[ij];
//		for (int ii = 0; ii < nr*nc; ii++){
//			pp[ii] = -p[ii] * ddy[ij];
//			ppp[ii] = p[ii] * ddx[ij];
//		}
//
//		/* to zero for consistency */
//		memset(warpedInverseDepths[ij], 0x00, sizeof(unsigned short)*nr*nc * 3);
//		memset(DispTargs[ij], 0x00, sizeof(float)*nr*nc);
//
//		warpInverseDepth(qDMF[ij], RowDisps[ij], ColDisps[ij], nr, nc, warpedInverseDepths[ij], DispTargs[ij]);
//
//	}
//	return;
//}
//
//void warpColor5Ref(float d_cen[21][101][2], int ref_rows[5], int ref_cols[5],
//	float *qDMF[5], float *RowDisps[5], float *ColDisps[5], int nr, int nc,
//	unsigned short *warpedColorViews[5], unsigned short *colorViews[5], float *DispTargs[5],
//	const int r, const int c)
//{
//	float view_distances[5];
//	float ddy[5];
//	float ddx[5];
//
//	for (int ij = 0; ij < 5; ij++){
//
//		ddy[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][1] - d_cen[r][c][1];
//		ddx[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][0] - d_cen[r][c][0];
//
//		view_distances[ij] = ddy[ij] * ddy[ij] + ddx[ij] * ddx[ij];
//
//	}
//
//	for (int uu = 0; uu < 5; uu++){
//
//		/* sorting */
//		float minim = 999999999999;
//		int ij = 0;
//		for (int kk = 0; kk < 5; kk++){
//			if (view_distances[kk] < minim){
//				minim = view_distances[kk];
//				ij = kk;
//			}
//		}
//
//		view_distances[ij] = 9999999999999999999;
//
//		float *p = qDMF[ij];
//		float *pp = RowDisps[ij];
//		float *ppp = ColDisps[ij];
//		for (int ii = 0; ii < nr*nc; ii++){
//			pp[ii] = -p[ii] * ddy[ij];
//			ppp[ii] = p[ii] * ddx[ij];
//		}
//
//		/* to zero for consistency */
//		memset(warpedColorViews[ij], 0x00, sizeof(unsigned short)*nr*nc * 3);
//		memset(DispTargs[ij], 0x00, sizeof(float)*nr*nc);
//
//		warpColorView(colorViews[ij], RowDisps[ij], ColDisps[ij], nr, nc, warpedColorViews[ij], DispTargs[ij]);
//
//	}
//	return;
//}
//
//void collectWarped(unsigned short *warpedColorViews[5], float *DispTargs[5], const int nr, const int nc, const int ncomponents, const int n_views){
//
//
//	//unsigned short *AA1 = new unsigned short[nr*nc*ncomponents];
//	//memcpy(AA1, warpedColorViews[0], nr*nc*ncomponents);
//
//	//float *DispTarg1 = new float[nr*nc];
//	//memcpy(DispTarg1, DispTargs[0], nr*nc);
//
//	unsigned short *AA1 = warpedColorViews[0];
//	float *DispTarg1 = DispTargs[0];
//
//	unsigned short *AA2 = new unsigned short[nr*nc*ncomponents];
//	float *DispTarg2 = new float[nr*nc];
//
//	for (int ik = 1; ik < n_views; ik++){
//
//		memcpy(AA2, warpedColorViews[ik], sizeof(unsigned short)*nr*nc*ncomponents);
//		memcpy(DispTarg2, DispTargs[ik], sizeof(float)*nr*nc);
//
//		for (int ij = 0; ij < nr*nc; ij++){
//			if (DispTarg1[ij] < 0 && DispTarg2[ij] >= 0){
//				DispTarg1[ij] = DispTarg2[ij];
//				AA1[ij] = AA2[ij];
//				AA1[ij + nr*nc] = AA2[ij + nr*nc];
//				AA1[ij + nr*nc * 2] = AA2[ij + nr*nc * 2];
//			}
//
//		}
//		//printf("ik:\t%d\n", ik);
//	}
//	delete(AA2);
//	delete(DispTarg2);
//}
//
//void warpInverseDepth(const float *inverseDepth, const float *DM_ROW, const float *DM_COL, const int nr, const int nc, float* warpedInvDepth, float* DispTarg) {
//
//
//	int *ROWS = new int[nr*nc];
//	int *COLS = new int[nr*nc];
//
//	float *DispTarg_col = new float[nr*nc];
//	float *DispTarg_row = new float[nr*nc];
//
//	for (int ij = 0; ij < nr*nc; ij++){
//		DispTarg[ij] = -1;
//		DispTarg_col[ij] = 99999;
//		DispTarg_row[ij] = 99999;
//	}
//
//
//	for (int ij = 0; ij < nr*nc; ij++)
//	{
//		float disp0 = abs(DM_COL[ij]) + abs(DM_ROW[ij]);
//
//		int ix = ij % nr; //row
//		int iy = (ij - ix) / nr; //col
//
//		int iynew = iy + (int)round(DM_COL[ij]);
//		int ixnew = ix + (int)round(DM_ROW[ij]);
//
//		if (iynew >= 0 && iynew < nc && ixnew >= 0 && ixnew < nr){
//			int indnew = ixnew + iynew*nr;
//			if (DispTarg[indnew] < disp0){
//				DispTarg[indnew] = disp0;
//				DispTarg_col[indnew] = iynew - iy; // warped horizontal disparity
//				DispTarg_row[indnew] = ixnew - ix; // warped vertical disparity
//				warpedInvDepth[indnew] = inverseDepth[ij];
//			}
//		}
//
//	}
//
//	delete(ROWS);
//	delete(COLS);
//	delete(DispTarg_col);
//	delete(DispTarg_row);
//}
//
//void warpColorView(const unsigned short *AA1, const float *DM_ROW, const float *DM_COL, const int nr, const int nc, unsigned short* Warped, float* DispTarg) {
//
//
//	int *ROWS = new int[nr*nc];
//	int *COLS = new int[nr*nc];
//
//	float *DispTarg_col = new float[nr*nc];
//	float *DispTarg_row = new float[nr*nc];
//
//	for (int ij = 0; ij < nr*nc; ij++){
//		DispTarg[ij] = -1;
//		DispTarg_col[ij] = 99999;
//		DispTarg_row[ij] = 99999;
//	}
//
//
//	for (int ij = 0; ij < nr*nc; ij++)
//	{
//		float disp0 = abs(DM_COL[ij]) + abs(DM_ROW[ij]);
//
//		int ix = ij % nr; //row
//		int iy = (ij - ix) / nr; //col
//
//		int iynew = iy + (int)round(DM_COL[ij]);
//		int ixnew = ix + (int)round(DM_ROW[ij]);
//
//		if (iynew >= 0 && iynew < nc && ixnew >= 0 && ixnew < nr){
//			int indnew = ixnew + iynew*nr;
//			if (DispTarg[indnew] < disp0){
//				DispTarg[indnew] = disp0;
//				DispTarg_col[indnew] = iynew - iy; // warped horizontal disparity
//				DispTarg_row[indnew] = ixnew - ix; // warped vertical disparity
//				Warped[indnew] = AA1[ij];
//				Warped[indnew + nr*nc] = AA1[ij + nr*nc];
//				Warped[indnew + 2 * nr*nc] = AA1[ij + 2 * nr*nc];
//			}
//		}
//
//	}
//
//	delete(ROWS);
//	delete(COLS);
//	delete(DispTarg_col);
//	delete(DispTarg_row);
//}



#endif