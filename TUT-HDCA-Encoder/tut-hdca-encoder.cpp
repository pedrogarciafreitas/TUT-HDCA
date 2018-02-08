#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include "..\CERV\cerv.h"
#include "..\GolombCoder\golomb_coder.hh"
#include "..\include\gen_types.hh"
#include "..\include\warpingFunctions.hh"

int main(int argc, char** argv) {

	/* parameters for encoder: CB/5 path_to_original_views path_to_unsw path_to_camera_centers bitrate NrQLev output_directory */

	//const char* path_original_views = argv[2];
	const char* path_UNSW = argv[1];
	//const char* path_camera_centers = argv[4];

	//double bitrate = atof(argv[5]);
	int NrQLev = atoi(argv[2]);

	int MED_FILT_SZ = atoi(argv[3]);

	char* MODE = "5Ref";
	char* path_camera_centers = "";
	char * path_to_references = "";

	if (argc > 3){
		MODE = argv[4];
		path_camera_centers = argv[5];
		path_to_references = argv[6];
	}

	const std::string output_path = std::string(argv[4]) + "/";

	// if using hevc, locations for external binaries of x265 encoder and decoder need be defined as well as ffmpeg binary
	/*const char x265_encoder_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/x265/x265.exe";
	const char x265_decoder_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/x265/TAppDecoder.exe";
	const char ffmpeg_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/ffmpeg/bin/ffmpeg.exe";*/

	const int nr = 1080;
	const int nc = 1920;

	int *inverse_depths[5]; //inverse depths
	float *inverse_depths_float[5];
	int *qDM[5]; //quantized inverse depths
	int *LDM[5]; //labelss
	float *qDMF[5]; //quantized floating point

	int *p, *pp; // dummy pointer

	float *ppf;

	/* for CERV */
	int** SEGM2D = (int**)malloc(nr*sizeof(int*));
	for (int i = 0; i < nr; ++i)
	{
		SEGM2D[i] = (int*)malloc(nc*sizeof(int));
	}
	int* SEGMFINAL = alocaVector((int)nr*nc);

	int ref_cols[] = { 2, 2, 50, 98, 98 };
	int ref_rows[] = { 0, 20, 10, 0, 20 };

	std::string filepath_pgm(path_UNSW);
	FILE *filept;

	std::string filename;

	for (int ij = 0; ij < 5; ij++)
	{

		/* Read UNSW inverse depth */
		char buffer[11];
		sprintf(buffer, "%03d_%03d.pgm", ref_cols[ij], ref_rows[ij]);

		filepath_pgm.replace(filepath_pgm.size() - 11, 11, buffer);
		std::cout << filepath_pgm << "\n";

		inverse_depths[ij] = new int[nr * nc];

		inverse_depths_float[ij] = new float[nr * nc];

		filept = fopen(filepath_pgm.c_str(), "rb");
		aux_read16pgm_1080p(filept, inverse_depths[ij]);
		fclose(filept);

		ppf = inverse_depths_float[ij];
		p = inverse_depths[ij];

		for (int ii = 0; ii < nr*nc; ii++)
			*(ppf++) = (static_cast<float>(*(p++))) / 16384;

		/* verify */

		if (0){
			filept = fopen(buffer, "wb");
			int *medianfiltered = new int[nr*nc];
			medfilt2D<int>(inverse_depths[ij], medianfiltered, 3, nr, nc);
			aux_write16pgm(filept, nc, nr, medianfiltered);
			fclose(filept);
		}


		/* Quantize UNSW inverse depth */

		float *medianfiltered = new float[nr*nc];
		medfilt2D<float>(inverse_depths_float[ij], medianfiltered, MED_FILT_SZ, nr, nc);

		delete(inverse_depths[ij]);
		delete(inverse_depths_float[ij]);

		float minDM = *std::min_element(medianfiltered, medianfiltered + nr*nc);
		float maxDM = *std::max_element(medianfiltered, medianfiltered + nr*nc);
		float Delta1 = (maxDM - minDM) / ((float)NrQLev - 1);

		std::cout << "Delta1 =\t" << Delta1 << "\tmaxDM = " << maxDM << "\tminDM = " << minDM << "\n";

		/*QUANTIZE*/

		qDMF[ij] = new float[nr*nc];
		qDM[ij] = new int[nr*nc];
		LDM[ij] = new int[nr*nc];

		ppf = qDMF[ij];
		p = LDM[ij];
		pp = qDM[ij];

		for (int ii = 0; ii < nr*nc; ii++){
			int Label1 = round((medianfiltered[ii] - minDM) / Delta1);

			*(ppf + ii) = ((float)Label1)*Delta1 + minDM; // in float
			*(p + ii) = Label1; // labels
			*(pp + ii) = (int)(*(ppf + ii) * 16384); // in integer (as in original .pgm, but quantized)
		}


		char unsw_quantized[128];
		sprintf(unsw_quantized, "%s%03d_%03d_Q_%d.pgm", output_path.c_str(), ref_cols[ij], ref_rows[ij], NrQLev - 1);

		filept = fopen(unsw_quantized, "wb");
		aux_write16pgm(filept, nc, nr, qDM[ij]);


		delete(medianfiltered);

		/* Encode quantized UNSW inverse depth with CERV */
		/* encode quantized labels with cerv */

		p = LDM[ij];

		for (int i = 0; i < nr; ++i) {
			for (int j = 0; j < nc; ++j) {
				SEGM2D[i][j] = *(p + i + j*nr);
			}
		}

		char cerv_filename[128];
		sprintf(cerv_filename, "%s%03d_%03d.cerv", output_path.c_str(), ref_cols[ij], ref_rows[ij]);

		std::cout << cerv_filename << "\n";

		cerv_encode(SEGM2D, nr, nc, cerv_filename);

		for (int i = 0; i < nr; ++i) {
			for (int j = 0; j < nc; ++j) {
				SEGM2D[i][j] = 0;
			}
		}

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

		std::vector<int> labels_symbols(number_of_regions, 0);

		pp = qDM[ij];
		for (int i = 0; i < nr*nc; ++i)
			labels_symbols.at(SEGMFINAL[i]) = *(pp + i);


		char cerv_labels_filename[128];
		sprintf(cerv_labels_filename, "%s%03d_%03d.gr", output_path.c_str(), ref_cols[ij], ref_rows[ij]);

		std::cout << cerv_labels_filename << "\n";

		GolombCoder golomb_coder(cerv_labels_filename, 0);
		golomb_coder.encode_symbols(labels_symbols, NBIT_GR);

		//delete(qDM[ij]);
		//delete(LDM[ij]);
		//delete(qDMF[ij]);

	}

	/* For checkerboard and other view configurations */
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

					/* reorder based on MSE */

					unsigned short *AA3 = new unsigned short[nr*nc * 3];
					float *SE_i = new float[nr*nc * 3];

					char pathtoref[160];
					sprintf(pathtoref, "%s%c%03d_%03d%s", path_to_references, '/', cc, rr, ".ppm");
					filept = fopen(pathtoref, "rb");
					aux_read16ppm(filept, nc, nr, AA3);
					fclose(filept);

					std::vector<pair<float, int>> mses(refscb_col.size());

					unsigned short *ps;

					for (int ij = 0; ij < refscb_col.size(); ij++){
						memset(SE_i, 0x00, sizeof(float));
						ps = warpedColorViews[ij];
						pff = DispTargs[ij];
						for (int ii = 0; ii < nr*nc * 3; ii++){
							if (*(pff + ii) > 0){
								float E = (float)AA3[ii] - (float)ps[ii];

								SE_i[ii] = E*E;
							}
						}
						for (int ii = 0; ii < nr*nc * 3; ii++)
							mses.at(ij).first = mses.at(ij).first + SE_i[ii];

						mses.at(ij).first = mses.at(ij).first / nr / nc / 3;
						mses.at(ij).second = ij;
					}

					delete(SE_i);
					delete(AA3);

					sort(mses.begin(), mses.end());

					char path_to_vorder[160];
					sprintf(path_to_vorder, "%s%c%03d_%03d%s", path_to_references, '/', cc, rr, ".vorder");
					filept = fopen(path_to_vorder, "wb");
					for (int ij = 0; ij < 5; ij++){
						if (ij >= refscb_col.size())
						{
							fwrite(0, sizeof(int), 1, filept);
						}
						else
						{
							fwrite(&mses.at(0).second, sizeof(int), 1, filept);
						}
					}
					fclose(filept);

					/*
					unsigned short *warpedColorViews_sorted[5];

					for (int ij = 0; ij < refscb_col.size(); ij++){
					warpedColorViews_sorted[ij] = warpedColorViews[mses.at(ij).second];
					}

					collectWarped<unsigned short>(warpedColorViews_sorted, DispTargs, nr, nc, 1, 5);

					fillHoles_T<unsigned short>(warpedColorViews_sorted[0], nr, nc, 1, 1);
					*/
				}
			}
		}

		

		delete(d_cen0);

		for (int ij = 0; ij < 5; ij++){
			delete(ColDisps[ij]);// = new float[nr*nc];
			delete(RowDisps[ij]);// = new float[nr*nc];
			delete(colorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete(warpedColorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete(DispTargs[ij]);// = new float[nr*nc];
			delete(warpedInverseDepths[ij]);
		}

	}

	/* prediction */

	/* All done. */

	for (int ij = 0; ij < 5; ij++){
		delete(qDM[ij]);
		delete(LDM[ij]);
		delete(qDMF[ij]);
	}

	exit(0);
}