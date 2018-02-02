#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <..\CERV\cerv.h>
#include <..\GolombCoder\golomb_coder.hh>

#include "gen_types.hh"
#include "warpingFunctions.hh"

int main(int argc, char** argv) {

	/* parameters for encoder: CB/5 path_to_original_views path_to_unsw path_to_camera_centers bitrate NrQLev output_directory */

	const char* path_to_references = argv[1];
	const char* path_camera_centers = argv[2];
	const char* MODE = argv[3];
	const char* path_out_dir = argv[4];

	// if using hevc, locations for external binaries of x265 encoder and decoder need be defined as well as ffmpeg binary
	const char x265_encoder_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/x265/x265.exe";
	const char x265_decoder_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/x265/TAppDecoder.exe";
	const char ffmpeg_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/ffmpeg/bin/ffmpeg.exe";

	const int nr = 1080;
	const int nc = 1920;

	int *inverse_depths[5]; //inverse depth
	float *qDMF[5]; //quantized inverse depth
	int *qDM[5];
	int *LDM[5]; //labels

	float *p, *pp, *ppp;
	int *ppi; // dummy pointer

	/* for CERV */
	int** SEGM2D = (int**)malloc(nr*sizeof(int*));
	for (int i = 0; i < nr; ++i)
	{
		SEGM2D[i] = (int*)malloc(nc*sizeof(int));
	}
	int* SEGMFINAL = alocaVector((int)nr*nc);

	int ref_cols[] = { 2, 2, 50, 98, 98 };
	int ref_rows[] = { 0, 20, 10, 0, 20 };

	float *d_cen0 = new float[101 * 21 * 2];

	FILE *filept;

	filept = fopen(path_camera_centers, "rb");
	fread(d_cen0, sizeof(float), 101 * 21 * 2, filept);
	fclose(filept);

	float d_cen[21][101][2];

	p = &d_cen0[0];

	for (int k = 0; k < 2; k++)
		for (int r = 0; r < 21; r++)
			for (int c = 0; c < 101; c++)
				d_cen[r][c][k] = *(p++);

	/* Decode inverse depth from .cerv and .gr */
	for (int ij = 0; ij < 5; ij++)
	{

		char cerv_filename[12];
		sprintf(cerv_filename, "%03d_%03d.cerv", ref_cols[ij], ref_rows[ij]);

		cerv_decode(SEGM2D, nr, nc, cerv_filename);

		int number_of_regions = 0;

		for (int i = 0; i < nr; ++i) {
			for (int j = 0; j < nc; ++j) {
				SEGMFINAL[i + j*nr] = SEGM2D[i][j];
				if (SEGMFINAL[i + j*nr] > number_of_regions)
					number_of_regions++;
			}
		}


		number_of_regions = number_of_regions + 1; //account for 0

		printf("NRegions:\t%d\n", number_of_regions);

		std::vector<int> labels_symbols;

		char cerv_labels_filename[12];
		sprintf(cerv_labels_filename, "%03d_%03d.gr", ref_cols[ij], ref_rows[ij]);
		GolombCoder golomb_coder(cerv_labels_filename, 1);
		golomb_coder.decode_symbols(labels_symbols, NBIT_GR);

		printf("\n n symbols = %i", labels_symbols.size());

		qDMF[ij] = new float[nr*nc];
		qDM[ij] = new int[nr*nc];

		p = qDMF[ij];
		ppi = qDM[ij];
		for (int i = 0; i < nr*nc; i++) {
			*(p + i) = ((float)labels_symbols.at(SEGMFINAL[i])) / 16384;
			*(ppi + i) = labels_symbols.at(SEGMFINAL[i]);
		}

		char qDM_filename[12];
		sprintf(qDM_filename, "%03d_%03d_Q.pgm", ref_cols[ij], ref_rows[ij]);
		filept = fopen(qDM_filename, "wb");
		aux_write16pgm(filept, nc, nr, qDM[ij]);
		fclose(filept);
	}

	bool ref_array[11][33];


	if (strcmp(MODE, "CB") == 0)
	{
		for (int r = 0; r < 11; r = r + 2)
			for (int c = 0; c < 33; c = c + 2)
				ref_array[r][c] = true;
		for (int r = 1; r < 11; r = r + 2)
			for (int c = 1; c < 33; c = c + 2)
				ref_array[r][c] = true;

	}
	if (strcmp(MODE, "5Ref") == 0)
	{

		unsigned short *ccomp = new unsigned short[nr*nc];
		unsigned short *medccomp = new unsigned short[nr*nc];

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
			ColDisps[ij] = new float[nr*nc];
			RowDisps[ij] = new float[nr*nc];
			colorViews[ij] = new unsigned short[nr*nc * 3];
			warpedColorViews[ij] = new unsigned short[nr*nc * 3];
			DispTargs[ij] = new float[nr*nc];

			char pathtoref[160];
			sprintf(pathtoref, "%s%c%03d_%03d%s", path_to_references, '/', ref_cols[ij], ref_rows[ij], ".ppm");

			filept = fopen(pathtoref, "rb");
			aux_read16ppm(filept, nc, nr, colorViews[ij]);
			fclose(filept);
		}

		for (int cc = 0; cc < 33; cc++){
			for (int rr = 0; rr < 11; rr++){

				int r, c;
				r = rr * 2;
				c = cc * 3 + 2;
				char path_out_ppm[160];
				sprintf(path_out_ppm, "%s%c%03d_%03d%s", path_out_dir, '/', c, r, ".ppm");


				if (!ref_array[rr][cc]){

					printf("Decoding to %s\n ", path_out_ppm);
					filept = fopen(path_out_ppm, "wb");

					float view_distances[5];
					float ddy[5];
					float ddx[5];

					for (int ij = 0; ij < 5; ij++){

						ddy[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][1] - d_cen[r][c][1];
						ddx[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][0] - d_cen[r][c][0];

						view_distances[ij] = ddy[ij] * ddy[ij] + ddx[ij] * ddx[ij];

					}

					for (int uu = 0; uu < 5; uu++){

						float minim = 999999999999;
						int ij = 0;
						for (int kk = 0; kk < 5; kk++){
							if (view_distances[kk] < minim){
								minim = view_distances[kk];
								ij = kk;
							}
						}

						view_distances[ij] = 9999999999999999999;

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

					/* results are collected to warpedColorViews[0] and DispTargs[0] */
					collectWarped(warpedColorViews, DispTargs, nr, nc, 3, 5);

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

					aux_write16ppm_16(filept, nc, nr, warpedColorViews[0]);
				}
				else{
					for (int ij = 0; ij < 5; ij++)
					{
						if (abs(ref_rows[ij] - rr) + abs(ref_cols[ij] - cc) < 1)
						{
							printf("Decoding to %s\n ", path_out_ppm);
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
			delete(ColDisps[ij]);// = new float[nr*nc];
			delete(RowDisps[ij]);// = new float[nr*nc];
			delete(colorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete(warpedColorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete(DispTargs[ij]);// = new float[nr*nc];
		}
		delete(ccomp);
		delete(medccomp);
	}

	/* clean up ... */
	for (int ij = 0; ij < 5; ij++){
		delete(qDMF[ij]);// = new float[nr*nc];
		delete(qDM[ij]);// = new int[nr*nc];
	}



	exit(0);
}