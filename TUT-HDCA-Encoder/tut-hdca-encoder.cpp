#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include <cstdint>

#include <ctime>

#include "../include/gen_types.hh"
#include "../include/warpingFunctions.hh"


void FiveRefHierarchy_2_disk(const char *hiearchy_file)
{
	const int ref_cols[] = { 2, 98, 2, 98, 50 };
	const int ref_rows[] = { 0, 0, 20, 20, 10 };

	bool vmask[256][256];

	for (int rr = 0; rr < 256; rr++) {
		for (int cc = 0; cc < 256; cc++) {
			vmask[rr][cc] = false;
			//wmask[rr][cc] = false;
		}
	}

	FILE *filept = fopen(hiearchy_file, "wb");

	int val = 33 * 11;

	fwrite(&val, sizeof(int), 1, filept);

	int ee_mask[21][101];
	int ee = 0;

	for (int ik = 0; ik < 5; ik++)
	{
		fwrite(&ref_rows[ik], sizeof(int), 1, filept);

		fwrite(&ref_cols[ik], sizeof(int), 1, filept);

		int rate = 100000 / 2;

		fwrite(&rate, sizeof(int), 1, filept);
		fwrite(&rate, sizeof(int), 1, filept);

		int val = 0;

		fwrite(&val, sizeof(int), 1, filept);

		vmask[ref_rows[ik]][ref_cols[ik]] = true;

		ee_mask[ref_rows[ik]][ref_cols[ik]] = ee;
		ee++;

	}

	for (int rr = 0; rr < 11; rr++) {
		for (int cc = 0; cc < 33; cc++) {

			int r, c;
			r = rr * 2;
			c = cc * 3 + 2;

			if (!vmask[r][c]) {

				fwrite(&r, sizeof(int), 1, filept);

				fwrite(&c, sizeof(int), 1, filept);

				int rate = 100000 / 2;

				fwrite(&rate, sizeof(int), 1, filept);
				fwrite(&rate, sizeof(int), 1, filept);

				int val = 5;

				fwrite(&val, sizeof(int), 1, filept);

				for (int ik = 0; ik < 5; ik++)
				{
					fwrite(&ee_mask[ref_rows[ik]][ref_cols[ik]], sizeof(int), 1, filept);
				}

				vmask[r][c] = 1;

				ee_mask[r][c] = ee;
				ee++;
			}

		}
	}

	fclose(filept);
}

int main(int argc, char** argv) {

	//C:/Local/astolap/Data/JPEG_PLENO/JPEG-PLENO-DATASETS/Fraunhofer_HDCA/ C:/Local/astolap/Data/JPEG_PLENO/TUT-HDCA_tmp_output/ C:/Local/astolap/Data/JPEG_PLENO/Kakadu/ 21 101

	const char* input_dir = argv[1];
	const char* output_dir = argv[2];

	//const char *kdu_compress_path = "\"C:/Program Files (x86)/Kakadu/kdu_compress.exe\"";
	//const char *kdu_expand_path = "\"C:/Program Files (x86)/Kakadu/kdu_expand.exe\"";

	const char *kakadu_dir = argv[3];

	//const float ref_color_rate = atof(argv[5]);
	//const float ref_depth_rate = atof(argv[6]);
	//const float residual_rate = atof(argv[7]);

	const int nar = atoi(argv[4]);
	const int nac = atoi(argv[5]);

	const char *hiearchy_file = argv[6];

	//const int nr = 1080;
	//const int nc = 1920;

	char kdu_compress_path[256];
	char kdu_expand_path[256];

	sprintf(kdu_compress_path, "%s%s", kakadu_dir, "kdu_compress.exe");
	sprintf(kdu_expand_path, "%s%s", kakadu_dir, "kdu_expand.exe");

	//FiveRefHierarchy_2_disk(hiearchy_file);

	const char *difftest_call = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --toycbcr --psnr ";
	const char *difftest_call_pgm = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --psnr ";

	char path_camera_centers[256];
	sprintf(path_camera_centers, "%s%s", input_dir, "/camera_displacement.bin");

	FILE *filept;

	float *d_cen0 = new float[nar * nac * 2]();

	filept = fopen(path_camera_centers, "rb");
	fread(d_cen0, sizeof(float), nar * nac * 2, filept);
	fclose(filept);

	float d_cen[256][256][2];

	float *p = &d_cen0[0];

	for (int k = 0; k < 2; k++)
		for (int r = 0; r < nar; r++)
			for (int c = 0; c < nac; c++)
				d_cen[r][c][k] = *(p++);

	filept = fopen(hiearchy_file, "rb");

	int n_views_total;
	fread(&n_views_total, sizeof(int), 1, filept);

	int Nd; // defines how many of the reference views are used for warping of the depth, for HDCA Nd = 5, for lenslet Nd maybe just 1
	fread(&Nd, sizeof(int), 1, filept);

	view *LF = new view[n_views_total]();

	int ee_mask[256][256];
	int ee = 0;

	for (int ii = 0; ii < n_views_total; ii++) {

		view *SAI = LF+ii;

		fread(&(SAI->r), sizeof(int), 1, filept);
		fread(&(SAI->c), sizeof(int), 1, filept);

		int rate_color, rate_depth, Ms, NNt;

		fread(&rate_color, sizeof(int), 1, filept);
		fread(&rate_depth, sizeof(int), 1, filept);

		fread(&SAI->Ms, sizeof(int), 1, filept);
		fread(&SAI->NNt, sizeof(int), 1, filept);

		SAI->residual_rate_color = ((float)rate_color) / 100000;
		SAI->residual_rate_depth = ((float)rate_depth) / 100000;

		ee_mask[SAI->r][SAI->c] = ee;
		ee++;

		fread(&(SAI->n_references), sizeof(int), 1, filept);

		if (SAI->n_references > 0) {

			SAI->references = new int[SAI->n_references]();

			int *ip = SAI->references;

			for (int uu = 0; uu < SAI->n_references; uu++) {
				fread(ip + uu, sizeof(int), 1, filept);
			}

		}

		SAI->y = d_cen[SAI->r][SAI->c][1];
		SAI->x = d_cen[SAI->r][SAI->c][0];

	}

	fclose(filept);


	for (int ii = 0; ii < n_views_total; ii++) {

		view *SAI = LF+ii;

		char output_results[1024];
		int output_buffer_length = 0;
		output_buffer_length += sprintf(output_results + output_buffer_length, "%03d\t%03d", SAI->r, SAI->c);

		unsigned short *original_color_view = NULL;
		unsigned short *original_depth_view = NULL;

		char path_input_ppm[1024];
		sprintf(path_input_ppm, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, ".ppm");

		int nc1, nr1, ncomp1;
		aux_read16PGMPPM(path_input_ppm, SAI->nc, SAI->nr, ncomp1, original_color_view);

		char path_input_depth_pgm[1024];
		sprintf(path_input_depth_pgm, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, ".pgm");
		bool depth_file_exist = aux_read16PGMPPM(path_input_depth_pgm, nc1, nr1, ncomp1, original_depth_view);

		char path_out_ppm[1024];
		sprintf(path_out_ppm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".ppm");

		char path_out_pgm[1024];
		sprintf(path_out_pgm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".pgm");

		SAI->color = new unsigned short[SAI->nr*SAI->nc * 3]();
		SAI->depth = new unsigned short[SAI->nr*SAI->nc]();

		if (SAI->n_references > 0) {

			/* currently we forward warp the depth from the N (for HDCA N = 5, lenslet maybe 1?) references */
			unsigned short **warped_color_views_0_N = new unsigned short*[SAI->n_references]();
			unsigned short **warped_depth_views_0_N = new unsigned short*[SAI->n_references]();
			float **DispTargs_0_N = new float*[SAI->n_references]();

			for (int ij = 0; ij < Nd; ij++)
			{
				warpView0_to_View1(LF + ij, SAI, warped_color_views_0_N[ij], warped_depth_views_0_N[ij], DispTargs_0_N[ij]);
			}

			/* holds partial warped views for ii */
			unsigned short **warped_color_views = new unsigned short*[SAI->n_references]();
			unsigned short **warped_depth_views = new unsigned short*[SAI->n_references]();
			float **DispTargs = new float*[SAI->n_references]();

			for (int ij = 0; ij < SAI->n_references; ij++)
			{

				int uu = SAI->references[ij];

				/* FORWARD warp color AND depth */
				warpView0_to_View1(LF + uu, SAI, warped_color_views[ij], warped_depth_views[ij], DispTargs[ij]);

				char tmp_str[1024];

				sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (LF + uu)->c, (LF + uu)->r, "_warped_to_", SAI->c, SAI->r, ".ppm");
				aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 3, warped_color_views[ij]);

				sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (LF + uu)->c, (LF + uu)->r, "_warped_to_", SAI->c, SAI->r, ".pgm");
				aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 1, warped_depth_views[ij]);

				//FILE *tmpf;
				//tmpf = fopen("G:/HEVC_HDCA/Berlin_verification_model/TMP_CPP/disptarg.float", "wb");
				//fwrite(DispTargs[ij], sizeof(float), 1080 * 1920, tmpf);
				//fclose(tmpf);

			}

			initViewW(SAI, DispTargs);

			/* get LS weights */
			if (1) {
				getViewMergingLSWeights_N(SAI, warped_color_views, DispTargs, original_color_view);
			}
			else {
				/* we don't use LS weights but something derived on geometric distance in view array*/
				getGeomWeight(SAI, LF, 10.0);
			}

			/* merge color */
			mergeWarped_N(warped_color_views, DispTargs, SAI, 3);

			//aux_write16PGMPPM(path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);

			/* hole filling for color*/
			holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);

			//aux_write16PGMPPM(path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);

			//aux_write16PGMPPM(path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);

			/* merge depth with median*/

			int startt = clock();

			#pragma omp parallel for
			for (int ij = 0; ij < SAI->nr*SAI->nc; ij++) {
				std::vector<unsigned short> depth_values;
				for (int uu = 0; uu < Nd; uu++) {
					//for (int uu = 0; uu < SAI->n_references; uu++) {
					unsigned short *pp = warped_depth_views_0_N[uu];
					float *pf = DispTargs_0_N[uu];
					if (*(pf + ij) > -1) {
						depth_values.push_back(*(pp + ij));
					}
				}
				if (depth_values.size() > 0)
					SAI->depth[ij] = getMedian(depth_values);
			}

			std::cout << "time elapsed in depth merging\t" << (int)clock() - startt << "\n";

			/* hole filling for depth */
			holefilling(SAI->depth, 1, SAI->nr, SAI->nc, 0);

			/* clean */
			for (int ij = 0; ij < SAI->n_references; ij++)
			{
				delete[](warped_color_views[ij]);
				delete[](warped_depth_views[ij]);
				delete[](DispTargs[ij]);
			}
			for (int ij = 0; ij < Nd; ij++)
			{
				delete[](warped_color_views_0_N[ij]);
				delete[](warped_depth_views_0_N[ij]);
				delete[](DispTargs_0_N[ij]);
			}

			delete[](warped_color_views_0_N);
			delete[](warped_depth_views_0_N);
			delete[](DispTargs_0_N);

			delete[](warped_color_views);
			delete[](warped_depth_views);
			delete[](DispTargs);
		}


		float psnr_result;

		if (SAI->n_references > 0) {

			aux_write16PGMPPM(path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
			psnr_result = getPSNR(NULL, path_out_ppm, path_input_ppm, difftest_call);

			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_result);

		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0);
		}

		if (SAI->NNt > 0 && SAI->Ms > 0) 
		{

			int startt = clock();

			getGlobalSparseFilter(SAI, original_color_view);

			std::cout << "time elapsed in getGlobalSparseFilter()\t" << (int)clock() - startt << "\n";

			applyGlobalSparseFilter( SAI );

			aux_write16PGMPPM(path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
			psnr_result = getPSNR(NULL, path_out_ppm, path_input_ppm, difftest_call);
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_result);
		}
		else 
		{
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0);
		}

		/* get residual */
		if (SAI->residual_rate_color > 0)
		{

			/* COLOR residual here */

			FILE *residual_file;

			char ppm_residual_path[512];

			char jp2_residual_path_jp2[512];

			sprintf(ppm_residual_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.ppm");

			sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.jp2");

			encodeResidualJP2(SAI->nr, SAI->nc, original_color_view, SAI->color, ppm_residual_path,
				kdu_compress_path, jp2_residual_path_jp2, SAI->residual_rate_color, 3, pow(2, 10) - 1);

			decodeResidualJP2(SAI->color, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path, ncomp1, pow(2, 10) - 1, pow(2, 10) - 1);

			if (depth_file_exist && SAI->residual_rate_depth > 0) { /* residual depth if needed */

				char pgm_residual_depth_path[512];

				char jp2_residual_depth_path_jp2[512];

				sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.pgm");

				sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.jp2");

				encodeResidualJP2(SAI->nr, SAI->nc, original_depth_view, SAI->depth, pgm_residual_depth_path,
					kdu_compress_path, jp2_residual_depth_path_jp2, SAI->residual_rate_depth, 1, 0);

				decodeResidualJP2(SAI->depth, kdu_expand_path, jp2_residual_depth_path_jp2, pgm_residual_depth_path, ncomp1, 0, pow(2, 16) - 1);

			}
		}

		/* medfilt depth */
		unsigned short *tmp_depth = new unsigned short[SAI->nr*SAI->nc]();
		int startt = clock();
		medfilt2D(SAI->depth, tmp_depth, 3, SAI->nr, SAI->nc);
		std::cout << "time elapsed in depth median filtering\t" << (int)clock() - startt << "\n";
		memcpy(SAI->depth, tmp_depth, sizeof(unsigned short)*SAI->nr*SAI->nc);
		delete[](tmp_depth);

		aux_write16PGMPPM(path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
		aux_write16PGMPPM(path_out_pgm, SAI->nc, SAI->nr, 1, SAI->depth);

		if (depth_file_exist) {
			psnr_result = getPSNR(NULL, path_out_pgm, path_input_depth_pgm, difftest_call_pgm);
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_result);
		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0);
		}


		psnr_result = getPSNR(NULL, path_out_ppm, path_input_ppm, difftest_call);

		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_result);

		delete[](original_color_view);
		delete[](original_depth_view);

		FILE *output_results_file;
		char output_results_filename[512];
		sprintf(output_results_filename, "%s%s", output_dir, "results.txt");
		if (ii < 1) {
			output_results_file = fopen(output_results_filename, "w");
		}
		else {
			output_results_file = fopen(output_results_filename, "a");
		}
		fprintf(output_results_file, "%s\n", output_results);
		fclose(output_results_file);
	}

	for (int ii = 0; ii < n_views_total; ii++) 
	{

		printf("ii=%d\n", ii);

		view *SAI = LF + ii;

		if (SAI->color != NULL)
			delete[](SAI->color);
		if (SAI->depth != NULL)
			delete[](SAI->depth);
		if (SAI->references != NULL)
			delete[](SAI->references);
		if (SAI->depth_references != NULL)
			delete[](SAI->depth_references);
		if (SAI->merge_weights != NULL)
			delete[](SAI->merge_weights);
		if (SAI->sparse_weights != NULL)
			delete[](SAI->sparse_weights);
		if (SAI->bmask != NULL)
			delete[](SAI->bmask);
		if (SAI->seg_vp != NULL) 
			delete[](SAI->seg_vp);
		if (SAI->sparse_mask != NULL)
			delete[](SAI->sparse_mask);

	}

	delete[](LF);

	exit(0);
}