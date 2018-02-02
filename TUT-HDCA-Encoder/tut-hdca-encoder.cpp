#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <..\CERV\cerv.h>
#include <..\GolombCoder\golomb_coder.hh>

#include "gen_types.hh"

int main(int argc, char** argv) {

	/* parameters for encoder: CB/5 path_to_original_views path_to_unsw path_to_camera_centers bitrate NrQLev output_directory */

	const char* path_original_views = argv[2];
	const char* path_UNSW = argv[3];
	const char* path_camera_centers = argv[4];

	double bitrate = atof(argv[5]);
	int NrQLev = atoi(argv[6]);

	// if using hevc, locations for external binaries of x265 encoder and decoder need be defined as well as ffmpeg binary
	const char x265_encoder_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/x265/x265.exe";
	const char x265_decoder_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/x265/TAppDecoder.exe";
	const char ffmpeg_path[] = "C:/Local/astolap/Data/JPEG_PLENO/05_Software/ffmpeg/bin/ffmpeg.exe";

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
			*(ppf++) = (static_cast< float >(*(p++))) / 16384;

		/* verify */
		/*
		if (0){
		filept = fopen(buffer, "wb");
		int *medianfiltered = new int[nr*nc];
		medfilt2D(inverse_depths[ij], medianfiltered, 3, nr, nc);
		aux_write16pgm(filept, nc, nr, medianfiltered);
		fclose(filept);
		}
		*/


		/* Quantize UNSW inverse depth */

		/*for (int ij = 0; ij < 5; ij++){*/
		float *medianfiltered = new float[nr*nc];
		medfilt2D(inverse_depths_float[ij], medianfiltered, 3, nr, nc);

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
			*(pp + ii) = (int)( *(ppf + ii) * 16384 ); // in integer (as in original .pgm, but quantized)
		}

		filept = fopen(buffer, "wb");
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

		char cerv_filename[12];
		sprintf(cerv_filename, "%03d_%03d.cerv", ref_cols[ij], ref_rows[ij]);

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

		std::vector<int> labels_symbols(number_of_regions);

		pp = qDM[ij];
		for (int i = 0; i < nr*nc; ++i)
			labels_symbols.at(SEGMFINAL[i]) = *(pp+i);


		char cerv_labels_filename[12];
		sprintf(cerv_labels_filename, "%03d_%03d.gr", ref_cols[ij], ref_rows[ij]);
		GolombCoder golomb_coder(cerv_labels_filename, 0);
		golomb_coder.encode_symbols(labels_symbols, NBIT_GR);

		delete(qDM[ij]);
		delete(LDM[ij]);
		delete(qDMF[ij]);
	}

	/* Encode reference views with X265 */


	/* All done. */


	exit(0);
}