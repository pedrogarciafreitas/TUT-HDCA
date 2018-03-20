#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include "..\include\gen_types.hh"
#include "..\include\warpingFunctions.hh"

int main(int argc, char** argv) {

	const char* path_to_references = argv[1];
	const char* path_camera_centers = argv[2];
	const char* path_out_dir = argv[3];

	const int nr = 1080;
	const int nc = 1920;

	float *qDMF[5];
	int *qDM[5];

	float *p, *pp, *ppp;

	int ref_cols[] = { 2, 98, 2, 98, 50 };
	int ref_rows[] = { 0, 0, 20, 20, 10 };

	float *d_cen0 = new float[101 * 21 * 2]();

	FILE *filept;

	filept = fopen(path_camera_centers, "rb");
	fread(d_cen0, sizeof(float), 101 * 21 * 2, filept);
	fclose(filept);

	bool *bmask = new bool[32 * 5]();

	for (int ij = 0; ij < 32; ij++){

		int uu = ij;

		for (int ik = 4; ik >= 0; ik--){

			if ( floor(uu /pow(2,ik)) > 0){
				uu = uu - pow(2, ik);
				bmask[ij + ik*32] = 1;
			}

		}
	}


	float d_cen[21][101][2];

	p = &d_cen0[0];

	for (int k = 0; k < 2; k++)
		for (int r = 0; r < 21; r++)
			for (int c = 0; c < 101; c++)
				d_cen[r][c][k] = *(p++);


	bool ref_array[11][33];

	for (int r = 0; r < 11; r++)
		for (int c = 0; c < 33; c++)
			ref_array[r][c] = false;


	char pathLSW[160];
	sprintf(pathLSW, "%s%s", path_to_references, "LS_weights");

	std::cout << pathLSW << "\n";

	FILE *fileptLS;
	fileptLS = fopen(pathLSW, "rb");

	if (1)
	{

		ref_array[0][0] = true;
		ref_array[10][0] = true;
		ref_array[10][32] = true;
		ref_array[0][32] = true;
		ref_array[5][16] = true;


		float *ColDisps[5];
		float *RowDisps[5];

		unsigned short *colorViews[5];
		unsigned short *warpedColorViews[5];
		float *DispTargs[5];

		for (int ij = 0; ij < 5; ij++){
			ColDisps[ij] = new float[nr*nc]();
			RowDisps[ij] = new float[nr*nc]();
			colorViews[ij] = new unsigned short[nr*nc * 3]();
			warpedColorViews[ij] = new unsigned short[nr*nc * 3]();
			DispTargs[ij] = new float[nr*nc]();

			qDM[ij] = new int[nr*nc]();
			qDMF[ij] = new float[nr*nc]();

			char pathtoref[160];
			sprintf(pathtoref, "%s%c%03d_%03d%s", path_to_references, '/', ref_cols[ij], ref_rows[ij], ".ppm");

			filept = fopen(pathtoref, "rb");
			aux_read16ppm(filept, nc, nr, colorViews[ij]);
			fclose(filept);



			char pgm_filename[12];
			sprintf(pgm_filename, "%s%c%03d_%03d%s", path_to_references, '/', ref_cols[ij], ref_rows[ij], "_drf.pgm");
			filept = fopen(pgm_filename, "rb");
			aux_read16pgm(filept, qDM[ij]);
			fclose(filept);

			float *pf = qDMF[ij];
			int *p = qDM[ij];



			for (int ik = 0; ik < nr*nc; ik++){
				int tmpi = *(p + ik);
				*(pf + ik) = ((float)tmpi) / pow(2, 14);
			}

		}

		for (int rr = 0; rr < 11; rr++){
			for (int cc = 0; cc < 33; cc++){

				int r, c;
				r = rr * 2;
				c = cc * 3 + 2;
				char path_out_ppm[160];
				sprintf(path_out_ppm, "%s%c%03d_%03d%s", path_out_dir, '/', c, r, ".ppm");


				if (!ref_array[rr][cc]){

					printf("Decoding to %s\n ", path_out_ppm);
					
					float ddy[5];
					float ddx[5];

					std::vector<std::pair<float, int>> view_distances(5);
					std::vector<std::pair<float, float>> dydx(5);

					std::cout << d_cen[r][c][0] << "\t" << d_cen[r][c][1] << "\t\n\n";

					for (int ij = 0; ij < 5; ij++){

						ddy[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][1] - d_cen[r][c][1];
						ddx[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][0] - d_cen[r][c][0];

						dydx.at(ij).first = ref_rows[ij] - r;
						dydx.at(ij).second = ref_cols[ij] - c;

						view_distances.at(ij).first = sqrt( dydx.at(ij).first * dydx.at(ij).first + dydx.at(ij).second * dydx.at(ij).second );
						view_distances.at(ij).second = ij;
					}

					sort(view_distances.begin(), view_distances.end());


					for (int uu = 0; uu < 5; uu++){

						int ij = view_distances.at(uu).second;

						p = qDMF[ij];
						pp = RowDisps[ij];
						ppp = ColDisps[ij];

						for (int ii = 0; ii < nr*nc; ii++){
							pp[ii] = -p[ii] * ddy[ij];
							ppp[ii] = p[ii] * ddx[ij];
						}

						/* to zero for consistency */
						memset(warpedColorViews[uu], 0x00, sizeof(unsigned short)*nr*nc * 3);
						memset(DispTargs[uu], 0x00, sizeof(float)*nr*nc);

						warpColorView(colorViews[ij], RowDisps[ij], ColDisps[ij], nr, nc, warpedColorViews[uu], DispTargs[uu]);

					}


					signed short *LSw_s = new signed short[80];

					float *LSw = new float[32 * 5]();

					float stdd = 10.0;

					if ( fileptLS==NULL ||  !(fread(LSw_s, sizeof(signed short), 80, fileptLS)==80) )
					{

						std::cout << "Fixed exponential weights, sigma = " << stdd << "\n";

						for (int ii = 0; ii < 32; ii++){
							float sumw = 0;
							for (int ij = 0; ij < 5; ij++){
								if (bmask[ii + ij * 32]){
									LSw[ii + ij * 32] = exp(-(view_distances.at(ij).first * view_distances.at(ij).first) / (2*stdd*stdd));
									sumw = sumw + LSw[ii + ij * 32];
								}
								else{
									LSw[ii + ij * 32] = 0.0;
								}

							}
							for (int ij = 0; ij < 5; ij++)
								if (bmask[ii + ij * 32])
									LSw[ii + ij * 32] = LSw[ii + ij * 32] / sumw;
						}

					}
					else{

						std::cout << "LS weights " << pathLSW << "\n";

						int uu = 0;

						for (int ii = 0; ii < 32 * 5; ii++){
							if (bmask[ii]){
								LSw[ii] = ((float)LSw_s[uu++]) / pow(2.0, 14.0);
							}
							else{
								LSw[ii] = 0.0;
							}
						}
					}


					
					//for (int ii = 0; ii < 32 * 5; ii++)
					//	std::cout << LSw[ii] << "\n";
						
					/* results are collected to warpedColorViews[0] and DispTargs[0] */
					/*here we need the LS weights*/
					collectWarpedLS(warpedColorViews, DispTargs, nr, nc, 3, 5, LSw);
					
					unsigned short *pshort = warpedColorViews[0];

					std::vector<unsigned short> neighbours;
					int dsz = 1;

					for (int ii = 0; ii < nr*nc; ii++){
						int y, x;
						y = ii%nr;
						x = floor(ii / nr);
						if (pshort[ii] + pshort[ii + nr*nc] + pshort[ii + 2 * nr*nc] < 1){
							for (int icomp = 0; icomp < 3; icomp++){
								neighbours.clear();
								for (int dy = -dsz; dy < dsz; dy++){
									for (int dx = -dsz; dx < dsz; dx++){
										if ((y + dy) >= 0 && (y + dy) < nr
											&& (x + dx) >= 0 && (x + dx) < nc){
											if (pshort[y + dy + (x + dx)*nr + icomp*nr*nc] > 0)
												neighbours.push_back(pshort[y + dy + (x + dx)*nr + icomp*nr*nc]);
										}
									}
								}
								if (neighbours.size() > 0)
									pshort[ii + icomp*nr*nc] = getMedian(neighbours);
							}
						}
					}

					filept = fopen(path_out_ppm, "wb");

					aux_write16ppm_16(filept, nc, nr, warpedColorViews[0]);

				}
				else{
					for (int ij = 0; ij < 5; ij++)
					{
						if ((abs(ref_rows[ij] - r) + abs(ref_cols[ij] - c)) < 1)
						{
							printf("Decoding to (this view is a reference) %s\n ", path_out_ppm);
							filept = fopen(path_out_ppm, "wb");
							aux_write16ppm_16(filept, nc, nr, colorViews[ij]);
							break;
						}
					}
				}
				fclose(filept);
			}
		}

		/* clean up ... */
		for (int ij = 0; ij < 5; ij++){
			delete[](ColDisps[ij]);// = new float[nr*nc];
			delete[](RowDisps[ij]);// = new float[nr*nc];
			delete[](colorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete[](warpedColorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete[](DispTargs[ij]);// = new float[nr*nc];
		}/*
		delete[](ccomp);
		delete[](medccomp);*/
	}

	fclose(fileptLS);

	/* clean up ... */
	for (int ij = 0; ij < 5; ij++){
		delete[](qDMF[ij]);// = new float[nr*nc];
		delete[](qDM[ij]);// = new int[nr*nc];
	}



	exit(0);
}