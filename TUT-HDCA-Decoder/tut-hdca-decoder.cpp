#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include "../CERV/cerv.h"
#include "../GolombCoder/golomb_coder.hh"
#include "../include/gen_types.hh"
#include "../include/warpingFunctions.hh"

int main(int argc, char** argv) {

	/* parameters for encoder: CB/5 path_to_original_views path_to_unsw path_to_camera_centers bitrate NrQLev output_directory */

	const char* path_to_references = argv[1];
	const char* path_camera_centers = argv[2];
	const char* MODE = argv[3];

	const char* input_dir = argv[4];

	const char* path_out_dir = argv[5];

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
		sprintf(cerv_filename, "%s%c%03d_%03d.cerv", input_dir, '/', ref_cols[ij], ref_rows[ij]);

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
		sprintf(cerv_labels_filename, "%s%c%03d_%03d.gr", input_dir, '/', ref_cols[ij], ref_rows[ij]);
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

		char qDM_filename[160];
		sprintf(qDM_filename, "%s%c%03d_%03d_Q.pgm", path_out_dir, '/', ref_cols[ij], ref_rows[ij]);
		filept = fopen(qDM_filename, "wb");
		aux_write16pgm(filept, nc, nr, qDM[ij]);
		fclose(filept);
	}

	bool ref_array[11][33];

	for (int r = 0; r < 11; r++)
		for (int c = 0; c < 33; c++)
			ref_array[r][c] = false;


	if (strcmp(MODE, "CB") == 0)
	{
		/* Read UNSW camera centers */
		int ref_cols[] = { 2, 2, 50, 98, 98 };
		int ref_rows[] = { 0, 20, 10, 0, 20 };

		float *d_cen0 = new float[101 * 21 * 2];

		filept = fopen(path_camera_centers, "rb");
		fread(d_cen0, sizeof(float), 101 * 21 * 2, filept);
		fclose(filept);

		float d_cen[21][101][2];

		float *pff = &d_cen0[0];

		for (int k = 0; k < 2; k++)
			for (int r = 0; r < 21; r++)
				for (int c = 0; c < 101; c++)
					d_cen[r][c][k] = *(p++);

		/* Make checkerboard array */
		bool ref_array[11][33];

		for (int r = 0; r < 11; r++)
			for (int c = 0; c < 33; c++)
				ref_array[r][c] = false;

		for (int r = 0; r < 11; r = r + 2)
			for (int c = 0; c < 33; c = c + 2)
				ref_array[r][c] = true;
		for (int r = 1; r < 11; r = r + 2)
			for (int c = 1; c < 33; c = c + 2)
				ref_array[r][c] = true;

		float *ColDisps[5];
		float *RowDisps[5];

		float *warpedInverseDepths[5];
		float *DispTargs[5];

		unsigned short* colorViews[5];
		unsigned short* warpedColorViews[5];

		for (int ij = 0; ij < 5; ij++){

			ColDisps[ij] = new float[nr*nc];
			RowDisps[ij] = new float[nr*nc];
			warpedInverseDepths[ij] = new float[nr*nc];
			DispTargs[ij] = new float[nr*nc];

			colorViews[ij] = new unsigned short[nr*nc * 3];
			warpedColorViews[ij] = new unsigned short[nr*nc * 3];

		}

		for (int cc = 0; cc < 33; cc++){
			for (int rr = 0; rr < 11; rr++){

				int r, c;
				r = rr * 2;
				c = cc * 3 + 2;
				if (!ref_array[rr][cc]){

					std::vector< int > refscb_row, refscb_col;

					/* collect available neighbourhood */
					for (int dy = -1; dy < 1; dy++){
						for (int dx = -1; dx < 1; dx++){
							if ((dy == 0 || dx == 0) && !(dy == 0 && dx == 0)){
								if (rr + dy < 33 && rr + dy >= 0 && cc + dx < 11 && cc + dx >= 0)
								{
									refscb_row.push_back(rr + dy);
									refscb_col.push_back(cc + dx);
								}
							}
						}
					}

					/* create inverse depths */
					for (int ij = 0; ij < refscb_col.size(); ij++){

						warpViews<float>(d_cen, ref_rows, ref_cols,
							qDMF, qDMF, RowDisps, ColDisps, nr, nc,
							warpedInverseDepths, DispTargs, refscb_row.at(ij), refscb_col.at(ij), 5, 1);

						collectWarped<float>(warpedInverseDepths, DispTargs, nr, nc, 1, 5);

						fillHoles_T<float>(warpedInverseDepths[0], nr, nc, 1, 1);

						/* read color reference */

						char pathtoref[160];
						sprintf(pathtoref, "%s%c%03d_%03d%s", path_to_references, '/', refscb_col[ij], refscb_row[ij], ".ppm");

						filept = fopen(pathtoref, "rb");
						aux_read16ppm(filept, nc, nr, colorViews[ij]);
						fclose(filept);

					}

					/* warp color views */
					warpViews<unsigned short>(d_cen, &refscb_row[0], &refscb_col[0],
						colorViews, warpedInverseDepths, RowDisps, ColDisps, nr, nc,
						warpedColorViews, DispTargs, r, c, refscb_col.size(), 1);

					std::vector<pair<float, int>> mses(refscb_col.size());

					char path_to_vorder[160];
					sprintf(path_to_vorder, "%s%c%03d_%03d%s", path_to_references, '/', cc, rr, ".vorder");
					filept = fopen(path_to_vorder, "rb");
					for (int ij = 0; ij < 5; ij++){
						if (ij >= refscb_col.size()){
						}
						else
						{
							fread(&mses.at(0).second, sizeof(int), 1, filept);
						}
					}
					fclose(filept);
					
					unsigned short *warpedColorViews_sorted[5];

					for (int ij = 0; ij < refscb_col.size(); ij++)
						warpedColorViews_sorted[ij] = warpedColorViews[mses.at(ij).second];	

					collectWarped<unsigned short>(warpedColorViews_sorted, DispTargs, nr, nc, 1, 5);

					fillHoles_T<unsigned short>(warpedColorViews_sorted[0], nr, nc, 1, 1);

					aux_write16ppm_16(filept, nc, nr, warpedColorViews_sorted[0]);
					
				}
			}

		}

		delete[](d_cen0);

		for (int ij = 0; ij < 5; ij++){
			delete[](ColDisps[ij]);// = new float[nr*nc];
			delete[](RowDisps[ij]);// = new float[nr*nc];
			delete[](colorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete[](warpedColorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete[](DispTargs[ij]);// = new float[nr*nc];
			delete[](warpedInverseDepths[ij]);
		}
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
				filept = fopen(path_out_ppm, "wb");

				char path_out_pgm[160];
				sprintf(path_out_pgm, "%s%s%03d_%03d%s", path_out_dir, "/DispTargT_", c, r, ".pgm");

				if (!ref_array[rr][cc]){

					printf("Decoding to %s\n ", path_out_ppm);


					warpViews<unsigned short>(d_cen, &ref_rows[0], &ref_cols[0],
						colorViews, qDMF, RowDisps, ColDisps, nr, nc,
						warpedColorViews, DispTargs, r, c, 5, 3);

					collectWarped<unsigned short>(warpedColorViews, DispTargs, nr, nc, 3, 5);

					/*collectWarped<float>(DispTargs, DispTargs, nr, nc, 1, 5);

					unsigned short *tmp_dt = new unsigned short[nr*nc];
					p = DispTargs[0];
					unsigned short *ps = tmp_dt;
					for (int ijj = 0; ijj < nr*nc; ijj++){
					if (*(p + ijj)>=0)
					tmp_dt[ijj] = (unsigned short)(*(p + ijj));
					}

					FILE *filept1 = fopen(path_out_pgm, "wb");

					aux_write16pgm(filept1, nc, nr, tmp_dt);

					fclose(filept1);*/

					//medfilt2D(warpedColorViews[0], warpedColorViews[0],
					//	nr, nc, 3);

					fillHoles_T<unsigned short>(warpedColorViews[0], nr, nc, 3, 1);

					aux_write16ppm_16(filept, nc, nr, warpedColorViews[0]);

				}
				else{
					for (int ij = 0; ij < 5; ij++)
					{
						if ((abs(ref_rows[ij] - r) + abs(ref_cols[ij] - c)) < 1)
						{
							printf("Decoding to (this view is a reference) %s\n ", path_out_ppm);
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
		}
		delete[](ccomp);
		delete[](medccomp);
	}

	/* clean up ... */
	for (int ij = 0; ij < 5; ij++){
		delete[](qDMF[ij]);// = new float[nr*nc];
		delete[](qDM[ij]);// = new int[nr*nc];
	}



	exit(0);
}