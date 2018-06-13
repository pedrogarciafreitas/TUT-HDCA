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

	const char *input_dir = argv[1];
	const char *output_dir = argv[2];
	const char *kakadu_dir = argv[3];
	const char *hiearchy_file = argv[4];

	char kdu_compress_path[1024];
	char kdu_expand_path[1024];

	sprintf(kdu_compress_path, "%s%s", kakadu_dir, "kdu_compress");
	sprintf(kdu_expand_path, "%s%s", kakadu_dir, "kdu_expand");

	const char *difftest_call = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --toycbcr --psnr ";
	const char *difftest_call_pgm = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --psnr ";

	FILE *filept;

	filept = fopen(hiearchy_file, "rb");

	int n_views_total;
	fread(&n_views_total, sizeof(int), 1, filept);

	view *LF = new view[n_views_total](); /* one dimensional view vector */

	int maxR = 0, maxC = 0;

	for (int ii = 0; ii < n_views_total; ii++) {

		view *SAI = LF + ii;

		fread(&(SAI->r), sizeof(int), 1, filept);
		fread(&(SAI->c), sizeof(int), 1, filept);

		/* find largest row and column */
		maxR = (SAI->r > maxR) ? SAI->r + 1 : maxR;
		maxC = (SAI->c > maxC) ? SAI->c + 1 : maxC;

		int xx = 0, yy = 0;

		fread(&xx, sizeof(int), 1, filept);
		fread(&yy, sizeof(int), 1, filept);

		SAI->x = float(xx) / 100000;
		SAI->y = float(yy) / 100000;

		int rate_color, rate_depth, Ms, NNt;

		fread(&rate_color, sizeof(int), 1, filept);
		fread(&rate_depth, sizeof(int), 1, filept);

		SAI->residual_rate_color = ((float)rate_color) / 100000;
		SAI->residual_rate_depth = ((float)rate_depth) / 100000;

		fread(&SAI->Ms, sizeof(int), 1, filept);
		fread(&SAI->NNt, sizeof(int), 1, filept);

		int stdd = 0;

		fread(&stdd, sizeof(int), 1, filept);
		SAI->stdd = ((float)stdd) / 100000;

		fread(&(SAI->n_references), sizeof(int), 1, filept);

		if (SAI->n_references > 0) {

			SAI->references = new int[SAI->n_references]();

			fread(SAI->references, sizeof(int), SAI->n_references, filept);

		}

		fread(&(SAI->n_depth_references), sizeof(int), 1, filept);
		if (SAI->n_depth_references > 0) {

			SAI->depth_references = new int[SAI->n_depth_references]();

			fread(SAI->depth_references, sizeof(int), SAI->n_depth_references, filept);

		}

		fread(&SAI->has_segmentation, sizeof(int), 1, filept); // new,13.06.18

		sprintf(SAI->path_input_ppm, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, ".ppm");
		sprintf(SAI->path_input_pgm, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, ".pgm");

		sprintf(SAI->path_input_seg, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, "_segmentation.pgm");

		sprintf(SAI->path_out_ppm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".ppm");
		sprintf(SAI->path_out_pgm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".pgm");

	}
	fclose(filept);


	/* 2D array format for views, useful in some cases */
	view ***LF_mat = new view**[maxR]();
	for (int ii = 0; ii < maxC; ii++) {
		LF_mat[ii] = new view*[maxC]();
	}

	for (int r = 0; r < maxR; r++) {
		for (int c = 0; c < maxC; c++) {
			LF_mat[r][c] = NULL;
		}
	}

	for (int ii = 0; ii < n_views_total; ii++) {
		view *SAI = LF + ii;
		LF_mat[SAI->r][SAI->c] = SAI;
	}

	/* get motion vectors/displacements/warping tables for segmentations */
	for (int ii = 0; ii < n_views_total; ii++)
	{
		view *SAI = LF + ii;

		if (SAI->has_segmentation > 0) {

			unsigned short *temp_segmentation;
			int ncomp1; //nc_temp and nr_temp should equal SAI->nr,SAI->nc

			if (!aux_read16PGMPPM(SAI->path_input_seg, SAI->nc, SAI->nr, ncomp1, temp_segmentation)) {
				printf("Error reading segmentation %s\nExiting...\n", SAI->path_input_seg);
				exit(0);
			}

			SAI->segmentation = new int[SAI->nr*SAI->nc]();

			int maxL = 0;

			for (int jj = 0; jj < SAI->nr*SAI->nc; jj++) {
				*(SAI->segmentation + jj) = (int)*(temp_segmentation + jj);

				maxL = *(SAI->segmentation + jj) > maxL ? *(SAI->segmentation + jj) : maxL;

			}

			/* initialize region displacement table */
			SAI->region_displacements = new int***[maxR]();
			for (int iar = 0; iar < maxR; iar++) {
				SAI->region_displacements[iar] = new int**[maxC]();
				for (int iac = 0; iac < maxC; iac++) {
					SAI->region_displacements[iar][iac] = new int*[maxL + 1]();
					for (int iR = 0; iR <= maxL; iR++) {
						SAI->region_displacements[iar][iac][iR] = new int[2]();
					}
				}
			}

			delete[](temp_segmentation);

			/* motion vectors for all regions against all views */

			unsigned short *im0;

			aux_read16PGMPPM(SAI->path_input_ppm, SAI->nc, SAI->nr, ncomp1, im0);

			for (int jj = 0; jj < n_views_total; jj++) {

				view *SAI1 = LF + jj;

				if (jj != ii) {

					unsigned short *im1;
					aux_read16PGMPPM(SAI1->path_input_ppm, SAI1->nc, SAI1->nr, ncomp1, im1);

					int search_radius = 10; // search -search_radius:search_radius in both x,y

					float ***match_score = new float**[maxL + 1]();
					for (int ik = 0; ik <= maxL; ik++) {
						match_score[ik] = new float*[search_radius * 2 + 1]();
						for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
							match_score[ik][isr] = new float[search_radius * 2 + 1]();
						}
					}

					int ***counts = new int**[maxL + 1]();
					for (int ik = 0; ik <= maxL; ik++) {
						counts[ik] = new int*[search_radius * 2 + 1]();
						for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
							counts[ik][isr] = new int[search_radius * 2 + 1]();
						}
					}

					int *segp = SAI->segmentation;

					for (int ijk = 0; ijk < SAI->nr*SAI->nc; ijk++) {
						for (int ik = 0; ik <= maxL; ik++) {
							if (*(segp + ijk) == ik) {
#pragma omp parallel for
								for (int isr = -search_radius; isr <= search_radius; isr++) {
									for (int isc = -search_radius; isc <= search_radius; isc++) {

										int iy = ijk % SAI->nr; //row
										int ix = (ijk - iy) / SAI->nr; //col

										int iy1 = iy + isr;
										int ix1 = ix + isc;

										if (iy1 >= 0 && iy1 < SAI->nr && ix1 >= 0 && ix1 < SAI->nc) {

											int ijk1 = ix1*SAI->nr + iy1;

											for (int ic = 0; ic < 3; ic++) {
												int offc = ic*SAI->nr*SAI->nc;
												match_score[ik][isr + search_radius][isc + search_radius] += abs(((float)*(im0 + ijk + offc) - (float)*(im1 + ijk1 + offc)));
											}

											//printf("%f\n", match_score[ik][isr + search_radius][isc + search_radius]);

											counts[ik][isr + search_radius][isc + search_radius]++;

											//printf("%d\n", counts[ik][isr + search_radius][isc + search_radius]);

										}
									}
								}
							}
						}
					}



					delete[](im1);

					/* find best matching displacement */
					for (int ik = 0; ik <= maxL; ik++) {
						float lowest_score = FLT_MAX;
						for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
							for (int isc = 0; isc < search_radius * 2 + 1; isc++) {

								int nk = counts[ik][isr][isc];

								if (nk > 0) {

									float clow = match_score[ik][isr][isc] / (float)nk;

									//printf("%f\n", clow);

									if (clow < lowest_score) {
										SAI->region_displacements[SAI1->r][SAI1->c][ik][0] = isr - search_radius;
										SAI->region_displacements[SAI1->r][SAI1->c][ik][1] = isc - search_radius;
										lowest_score = clow;
									}

								}
							}
						}
					}

					for (int iR = 0; iR <= maxL; iR++) {
						printf("SAI->region_displacements[%d][%d][%d][0] = %d\t", SAI1->r, SAI1->c, iR, SAI->region_displacements[SAI1->r][SAI1->c][iR][0]);
						printf("SAI->region_displacements[%d][%d][%d][1] = %d\n", SAI1->r, SAI1->c, iR, SAI->region_displacements[SAI1->r][SAI1->c][iR][1]);
					}


					//char dummy;
					//std::cin >> dummy;

					for (int ik = 0; ik <= maxL; ik++) {
						for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
							delete[](match_score[ik][isr]);
							delete[](counts[ik][isr]);
						}
						delete[](match_score[ik]);
						delete[](counts[ik]);
					}
					delete[](match_score);
					delete[](counts);

				}
				else {
					/* nothing to do, region displacements are zero*/
				}

			}

			delete[](im0);

		}
	}

	char path_out_LF_data[1024];
	sprintf(path_out_LF_data, "%s%c%s", output_dir, '/', "output.LF");

	/* our bitstream starts here */
	FILE *output_LF_file;
	output_LF_file = fopen(path_out_LF_data, "wb");
	fwrite(&n_views_total, sizeof(int), 1, output_LF_file);
	fclose(output_LF_file);

	for (int ii = 0; ii < n_views_total; ii++) {

		view *SAI = LF + ii;

		char output_results[1024];
		int output_buffer_length = 0;
		output_buffer_length += sprintf(output_results + output_buffer_length, "%03d\t%03d", SAI->r, SAI->c);

		unsigned short *original_color_view = NULL;
		unsigned short *original_depth_view = NULL;

		int nc1, nr1, ncomp1;
		aux_read16PGMPPM(SAI->path_input_ppm, SAI->nc, SAI->nr, ncomp1, original_color_view);

		bool depth_file_exist = aux_read16PGMPPM(SAI->path_input_pgm, nc1, nr1, ncomp1, original_depth_view);

		SAI->color = new unsigned short[SAI->nr*SAI->nc * 3]();
		SAI->depth = new unsigned short[SAI->nr*SAI->nc]();

		/* forward warp depth */
		if (SAI->n_depth_references > 0) {
			/* currently we forward warp the depth from the N (for HDCA N = 5, lenslet maybe 1?) references */
			unsigned short **warped_color_views_0_N = new unsigned short*[SAI->n_depth_references]();
			unsigned short **warped_depth_views_0_N = new unsigned short*[SAI->n_depth_references]();
			float **DispTargs_0_N = new float*[SAI->n_depth_references]();

			for (int ij = 0; ij < SAI->n_depth_references; ij++)
			{
				view *ref_view = LF + SAI->depth_references[ij];

				int tmp_w, tmp_r, tmp_ncomp;

				aux_read16PGMPPM(ref_view->path_out_pgm, tmp_w, tmp_r, tmp_ncomp, ref_view->depth);
				aux_read16PGMPPM(ref_view->path_out_ppm, tmp_w, tmp_r, tmp_ncomp, ref_view->color);

				warpView0_to_View1(ref_view, SAI, warped_color_views_0_N[ij], warped_depth_views_0_N[ij], DispTargs_0_N[ij]);

				delete[](ref_view->depth);
				delete[](ref_view->color);

				ref_view->depth = NULL;
				ref_view->color = NULL;
			}

			/* merge depth using median*/

			int startt = clock();

#pragma omp parallel for
			for (int ij = 0; ij < SAI->nr*SAI->nc; ij++) {
				std::vector<unsigned short> depth_values;
				for (int uu = 0; uu < SAI->n_depth_references; uu++) {
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

			for (int ij = 0; ij < SAI->n_depth_references; ij++)
			{
				delete[](warped_color_views_0_N[ij]);
				delete[](warped_depth_views_0_N[ij]);
				delete[](DispTargs_0_N[ij]);
			}

			delete[](warped_color_views_0_N);
			delete[](warped_depth_views_0_N);
			delete[](DispTargs_0_N);
		}

		/* forward warp color */
		if (SAI->n_references > 0) {

			/* holds partial warped views for ii */
			unsigned short **warped_color_views = new unsigned short*[SAI->n_references]();
			unsigned short **warped_depth_views = new unsigned short*[SAI->n_references]();
			float **DispTargs = new float*[SAI->n_references]();

			for (int ij = 0; ij < SAI->n_references; ij++)
			{

				view *ref_view = LF + SAI->references[ij];

				int tmp_w, tmp_r, tmp_ncomp;

				aux_read16PGMPPM(ref_view->path_out_pgm, tmp_w, tmp_r, tmp_ncomp, ref_view->depth);
				aux_read16PGMPPM(ref_view->path_out_ppm, tmp_w, tmp_r, tmp_ncomp, ref_view->color);

				/* FORWARD warp color AND depth */
				warpView0_to_View1(ref_view, SAI, warped_color_views[ij], warped_depth_views[ij], DispTargs[ij]);

				delete[](ref_view->depth);
				delete[](ref_view->color);

				ref_view->depth = NULL;
				ref_view->color = NULL;

				char tmp_str[1024];

				sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, ".ppm");
				aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 3, warped_color_views[ij]);

				sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, ".pgm");
				aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 1, warped_depth_views[ij]);

				//FILE *tmpf;
				//tmpf = fopen("G:/HEVC_HDCA/Berlin_verification_model/TMP_CPP/disptarg.float", "wb");
				//fwrite(DispTargs[ij], sizeof(float), 1080 * 1920, tmpf);
				//fclose(tmpf);

			}

			initViewW(SAI, DispTargs);

			/* get LS weights */
			if (SAI->stdd == 0) {
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

			//char tmp_str[1024];
			//sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (LF + 0)->c, (LF + 0)->r, "_warped_to_", SAI->c, SAI->r, "_inpainted.ppm");
			//aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 3, SAI->color);

			//aux_write16PGMPPM(path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);

			//aux_write16PGMPPM(path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);

			/* clean */
			for (int ij = 0; ij < SAI->n_references; ij++)
			{
				delete[](warped_color_views[ij]);
				delete[](warped_depth_views[ij]);
				delete[](DispTargs[ij]);
			}

			delete[](warped_color_views);
			delete[](warped_depth_views);
			delete[](DispTargs);
		}


		float psnr_result;

		if (SAI->n_references > 0) {

			aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
			psnr_result = getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call);

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

			applyGlobalSparseFilter(SAI);

			aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
			psnr_result = getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call);
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_result);

		}
		else
		{
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0);
		}

		char ppm_residual_path[512];

		char jp2_residual_path_jp2[512];

		char pgm_residual_depth_path[512];

		char jp2_residual_depth_path_jp2[512];

		/* get residual */
		if (SAI->residual_rate_color > 0)
		{

			/* COLOR residual here */

			sprintf(ppm_residual_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.ppm");

			sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.jp2");

			encodeResidualJP2(SAI->nr, SAI->nc, original_color_view, SAI->color, ppm_residual_path,
				kdu_compress_path, jp2_residual_path_jp2, SAI->residual_rate_color, 3, pow(2, 10) - 1);

			decodeResidualJP2(SAI->color, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path, ncomp1, pow(2, 10) - 1, pow(2, 10) - 1);

		}

		if (depth_file_exist && SAI->residual_rate_depth > 0) { /* residual depth if needed */

			sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.pgm");

			sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.jp2");

			encodeResidualJP2(SAI->nr, SAI->nc, original_depth_view, SAI->depth, pgm_residual_depth_path,
				kdu_compress_path, jp2_residual_depth_path_jp2, SAI->residual_rate_depth, 1, 0);

			decodeResidualJP2(SAI->depth, kdu_expand_path, jp2_residual_depth_path_jp2, pgm_residual_depth_path, ncomp1, 0, pow(2, 16) - 1);

		}

		/* medfilt depth */
		unsigned short *tmp_depth = new unsigned short[SAI->nr*SAI->nc]();
		int startt = clock();
		medfilt2D(SAI->depth, tmp_depth, 3, SAI->nr, SAI->nc);
		std::cout << "time elapsed in depth median filtering\t" << (int)clock() - startt << "\n";
		memcpy(SAI->depth, tmp_depth, sizeof(unsigned short)*SAI->nr*SAI->nc);
		delete[](tmp_depth);

		aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
		aux_write16PGMPPM(SAI->path_out_pgm, SAI->nc, SAI->nr, 1, SAI->depth);

		if (depth_file_exist) {
			psnr_result = getPSNR(NULL, SAI->path_out_pgm, SAI->path_input_pgm, difftest_call_pgm);
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_result);
		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0);
		}


		psnr_result = getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call);

		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_result);

		delete[](original_color_view);
		delete[](original_depth_view);

		FILE *output_results_file;
		char output_results_filename[1024];
		sprintf(output_results_filename, "%s%s", output_dir, "results.txt");
		if (ii < 1) {
			output_results_file = fopen(output_results_filename, "w");
		}
		else {
			output_results_file = fopen(output_results_filename, "a");
		}
		fprintf(output_results_file, "%s\n", output_results);
		fclose(output_results_file);


		/* write view configuration data to disk */
		output_LF_file = fopen(path_out_LF_data, "ab");

		fwrite(&SAI->r, sizeof(int), 1, output_LF_file);
		fwrite(&SAI->c, sizeof(int), 1, output_LF_file);

		fwrite(&SAI->nr, sizeof(int), 1, output_LF_file);
		fwrite(&SAI->nc, sizeof(int), 1, output_LF_file);

		fwrite(&SAI->x, sizeof(int), 1, output_LF_file);
		fwrite(&SAI->y, sizeof(float), 1, output_LF_file);

		fwrite(&SAI->n_references, sizeof(int), 1, output_LF_file);

		if (SAI->n_references > 0) {
			fwrite(SAI->references, sizeof(int), SAI->n_references, output_LF_file);
		}

		fwrite(&SAI->n_depth_references, sizeof(int), 1, output_LF_file);
		if (SAI->n_depth_references > 0) {
			fwrite(SAI->depth_references, sizeof(int), SAI->n_depth_references, output_LF_file);
		}

		fwrite(&SAI->NNt, sizeof(int), 1, output_LF_file);
		fwrite(&SAI->Ms, sizeof(int), 1, output_LF_file);

		fwrite(&SAI->stdd, sizeof(float), 1, output_LF_file);

		if (SAI->stdd == 0) {

			if (SAI->NB > 0) {

				/* for debug and studying, we write merging weights to standalone files also*/

				FILE *wfile;
				char wfile_name[1024];
				sprintf(wfile_name, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_merging_weights");
				wfile = fopen(wfile_name, "wb");
				fwrite(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, wfile);
				fclose(wfile);


				fwrite(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, output_LF_file);
			}

		}

		if (SAI->Ms > 0 && SAI->NNt > 0) {
			fwrite(SAI->sparse_weights, sizeof(int32_t), SAI->Ms, output_LF_file);
			fwrite(SAI->sparse_mask, sizeof(unsigned char), SAI->Ms, output_LF_file);
		}

		if (SAI->residual_rate_color > 0) {
			int n_bytes_color_residual = aux_GetFileSize(jp2_residual_path_jp2);

			unsigned char *jp2_residual = new unsigned char[n_bytes_color_residual]();
			FILE *jp2_color_residual_file = fopen(jp2_residual_path_jp2, "rb");
			fread(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, jp2_color_residual_file);
			fclose(jp2_color_residual_file);

			fwrite(&n_bytes_color_residual, sizeof(int), 1, output_LF_file);
			fwrite(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, output_LF_file);

			delete[](jp2_residual);
		}
		else {
			int n_bytes_color_residual = 0;
			fwrite(&n_bytes_color_residual, sizeof(int), 1, output_LF_file);
		}

		if (depth_file_exist && SAI->residual_rate_depth > 0) {
			int n_bytes_depth_residual = aux_GetFileSize(jp2_residual_depth_path_jp2);

			unsigned char *jp2_depth_residual = new unsigned char[n_bytes_depth_residual]();
			FILE *jp2_depth_residual_file = fopen(jp2_residual_depth_path_jp2, "rb");
			fread(jp2_depth_residual, sizeof(unsigned char), n_bytes_depth_residual, jp2_depth_residual_file);
			fclose(jp2_depth_residual_file);

			fwrite(&n_bytes_depth_residual, sizeof(int), 1, output_LF_file);
			fwrite(jp2_depth_residual, sizeof(unsigned char), n_bytes_depth_residual, output_LF_file);

			delete[](jp2_depth_residual);
		}
		else {
			int n_bytes_depth_residual = 0;
			fwrite(&n_bytes_depth_residual, sizeof(int), 1, output_LF_file);
		}

		fclose(output_LF_file);

		/* to reduce memory usage */
		if (SAI->color != NULL) {
			delete[](SAI->color);
			SAI->color = NULL;
		}

		if (SAI->depth != NULL) {
			delete[](SAI->depth);
			SAI->depth = NULL;
		}

		if (SAI->seg_vp != NULL) {
			delete[](SAI->seg_vp);
			SAI->seg_vp = NULL;
		}

		if (SAI->segmentation != NULL) {
			delete[](SAI->segmentation);
			SAI->segmentation = NULL;
		}

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
		if (SAI->segmentation != NULL)
			delete[](SAI->segmentation);

	}

	delete[](LF);

	for (int ii = 0; ii < maxC; ii++) {
		delete[](LF_mat[ii]);
	}
	delete[](LF_mat);

	exit(0);
}