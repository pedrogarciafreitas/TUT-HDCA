#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

//#include "..\CERV\cerv.h"
//#include "..\GolombCoder\golomb_coder.hh"

#include "..\include\gen_types.hh"
#include "..\include\warpingFunctions.hh"

int main(int argc, char** argv) {

	/* parameters for encoder: CB/5 path_to_original_views path_to_unsw path_to_camera_centers bitrate NrQLev output_directory */

	const char* path_to_references = argv[1];
	const char* path_camera_centers = argv[2];
	const char* path_out_dir = argv[3];

	const int nr = 1080;
	const int nc = 1920;

	int *inverse_depths[5]; //inverse depth
	float *qDMF[5]; //quantized inverse depth
	int *qDM[5];
	int *LDM[5]; //labels

	float *p, *pp, *ppp;
	int *ppi; // dummy pointer

	//int ref_cols[] = { 2, 2, 50, 98, 98 };
	//int ref_rows[] = { 0, 20, 10, 0, 20 };

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


	if (0)//(strcmp(MODE, "CB") == 0)
	{
		for (int r = 0; r < 11; r = r + 2)
			for (int c = 0; c < 33; c = c + 2)
				ref_array[r][c] = true;
		for (int r = 1; r < 11; r = r + 2)
			for (int c = 1; c < 33; c = c + 2)
				ref_array[r][c] = true;

	}
	if (1)//(strcmp(MODE, "5Ref") == 0)
	{

		//unsigned short *ccomp = new unsigned short[nr*nc]();
		//unsigned short *medccomp = new unsigned short[nr*nc]();

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
					

					float view_distances[5];
					float ddy[5];
					float ddx[5];

					for (int ij = 0; ij < 5; ij++){

						ddy[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][1] - d_cen[r][c][1];
						ddx[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][0] - d_cen[r][c][0];

						view_distances[ij] = ddy[ij] * ddy[ij] + ddx[ij] * ddx[ij];

					}

					for (int ij = 0; ij < 5; ij++){


						p = qDMF[ij];


						p = qDMF[ij];
						pp = RowDisps[ij];
						ppp = ColDisps[ij];
						for (int ii = 0; ii < nr*nc; ii++){
							pp[ii] = -p[ii] * ddy[ij];
							ppp[ii] = p[ii] * ddx[ij];
						}


						/* to zero for consistency */
						memset(warpedColorViews[ij], 0x00, sizeof(unsigned short)*nr*nc * 3);
						memset(DispTargs[ij], 0x00, sizeof(float)*nr*nc);

						warpColorView(colorViews[ij], RowDisps[ij], ColDisps[ij], nr, nc, warpedColorViews[ij], DispTargs[ij]);

					}


					signed short *LSw_s = new signed short[80];

					char pathLSW[160];
					sprintf(pathLSW, "%s%c%03d_%03d%s", path_to_references, '/', c, r, "_LS_weights");

					std::cout << pathLSW << "\n";


					filept = fopen(pathLSW, "rb");

					if (filept == NULL){
						std::cout << "Unable to open " << pathLSW << "\n";
						return 0;
					}
					
					fread(LSw_s, 80, sizeof(signed short), filept);
					fclose(filept);

					


					float *LSw = new float[32 * 5]();

					int uu = 0;

					for (int ii = 0; ii < 32 * 5; ii++){
						if (bmask[ii]){
							LSw[ii] = ((float)LSw_s[uu++])/pow(2,12);
						}
						else{
							LSw[ii] = 0.0;
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

	/* clean up ... */
	for (int ij = 0; ij < 5; ij++){
		delete[](qDMF[ij]);// = new float[nr*nc];
		delete[](qDM[ij]);// = new int[nr*nc];
	}



	exit(0);
}