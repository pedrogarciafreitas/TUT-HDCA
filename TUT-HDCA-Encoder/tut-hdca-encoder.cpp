#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <..\CERV\cerv.h>
#include <..\GolombCoder\golomb_coder.hh>

#include "gen_types.hh"

int getMedian(std::vector<int> scores)
{
	int median;
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

void medfilt2D(int* input, int* output, int SZ, int nr, int nc)
{
	int dsz = floor(SZ / 2);
	std::vector<int> scores;
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

	int *inverse_depths[5]; //inverse depth
	int *qDM[5]; //quantized inverse depth
	int *LDM[5]; //labels

	int *p, *pp; // dummy pointer

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


		filept = fopen(filepath_pgm.c_str(), "rb");
		aux_read16pgm_1080p(filept, inverse_depths[ij]);
		fclose(filept);

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
		int *medianfiltered = new int[nr*nc];
		medfilt2D(inverse_depths[ij], medianfiltered, 3, nr, nc);

		int minDM = *std::min_element(medianfiltered, medianfiltered + nr*nc);
		int maxDM = *std::max_element(medianfiltered, medianfiltered + nr*nc);
		int Delta1 = (maxDM - minDM) / (NrQLev - 1);

		if (1){
			qDM[ij] = new int[nr*nc];
			LDM[ij] = new int[nr*nc];
			pp = qDM[ij];
			p = LDM[ij];
			for (int ii = 0; ii < nr*nc; ii++){
				int Label1 = round((medianfiltered[ii] - minDM) / Delta1);

				*(pp + ii) = Label1*Delta1 + minDM;
				*(p + ii) = Label1;
			}

			filept = fopen(buffer, "wb");
			aux_write16pgm(filept, nc, nr, qDM[ij]);
		}


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

		std::vector<int> labels_symbols;

		pp = qDM[ij];
		int prevr = 0;
		for (int i = 0; i < nr*nc; ++i) {
			if (SEGMFINAL[i] == prevr+1){
				prevr = prevr + 1;
				labels_symbols.push_back(floor(*(pp + i) * 2 ^ 10));
			}
		}

		char cerv_labels_filename[12];
		sprintf(cerv_labels_filename, "%03d_%03d.gr", ref_cols[ij], ref_rows[ij]);
		GolombCoder golomb_coder(cerv_labels_filename, 0);
		golomb_coder.encode_symbols(labels_symbols, 10);

		//printf("%d\n", number_of_regions);

	}

	/* Encode reference views with X265 */


	/* All done. */


	exit(0);
}