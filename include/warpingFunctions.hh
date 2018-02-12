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

	delete[](ROWS);
	delete[](COLS);
	delete[](DispTarg_col);
	delete[](DispTarg_row);

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

	delete[](AA2);
	delete[](DispTarg2);

	return;
}



#endif