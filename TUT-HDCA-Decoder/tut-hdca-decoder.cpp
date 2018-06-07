#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include <cstdint>

#include <ctime>

#include "../include/gen_types.hh"
#include "../include/warpingFunctions.hh"


int main(int argc, char** argv) {

	const char *input_file = argv[1];
	const char *output_dir = argv[2];
	const char *kakadu_dir = argv[3];

	char kdu_expand_path[256];

	sprintf(kdu_expand_path, "%s%s", kakadu_dir, "kdu_expand");

	FILE *input_LF;
	input_LF = fopen(input_file, "rb");

	int n_views_total, Nd;

	int f_status;

	f_status = fread(&n_views_total, sizeof(int), 1,input_LF);
	//f_status = fread(&Nd, sizeof(int), 1,input_LF);

	view *LF = new view[n_views_total]();

	int ii = 0; /*view index*/

	while (f_status>0) {

		view *SAI = LF+ii;

		ii++;

		f_status = fread(&SAI->r, sizeof(int), 1, input_LF);
		f_status = fread(&SAI->c, sizeof(int), 1, input_LF);

		f_status = fread(&SAI->nr, sizeof(int), 1, input_LF);
		f_status = fread(&SAI->nc, sizeof(int), 1, input_LF);

		f_status = fread(&SAI->x, sizeof(int), 1, input_LF);
		f_status = fread(&SAI->y, sizeof(float), 1, input_LF);

		f_status = fread(&SAI->n_references, sizeof(int), 1, input_LF);

		if (SAI->n_references > 0) {
			SAI->references = new int[SAI->n_references]();
			f_status = fread(SAI->references, sizeof(int), SAI->n_references, input_LF);
		}

		f_status = fread(&SAI->n_depth_references, sizeof(int), 1, input_LF);
		if (SAI->n_depth_references > 0) {
			SAI->depth_references = new int[SAI->n_depth_references]();
			f_status = fread(SAI->depth_references, sizeof(int), SAI->n_depth_references, input_LF);
		}

		f_status = fread(&SAI->NNt, sizeof(int), 1, input_LF);
		f_status = fread(&SAI->Ms, sizeof(int), 1, input_LF);

		f_status = fread(&SAI->stdd, sizeof(float), 1, input_LF);

		f_status = fread(&SAI->stdd, sizeof(float), 1, input_LF);

		if (!(f_status > 0)) {
			break;
		}

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

			/* merge depth with median*/

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

				//int uu = SAI->references[ij];

				/* FORWARD warp color AND depth */
				warpView0_to_View1(ref_view, SAI, warped_color_views[ij], warped_depth_views[ij], DispTargs[ij]);

				delete[](ref_view->depth);
				delete[](ref_view->color);

				ref_view->depth = NULL;
				ref_view->color = NULL;

				//char tmp_str[1024];

				//sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (LF + uu)->c, (LF + uu)->r, "_warped_to_", SAI->c, SAI->r, ".ppm");
				//aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 3, warped_color_views[ij]);

				//sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (LF + uu)->c, (LF + uu)->r, "_warped_to_", SAI->c, SAI->r, ".pgm");
				//aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 1, warped_depth_views[ij]);

			}

			initViewW(SAI, DispTargs);

			if (!(SAI->stdd == 0)) {
				if (SAI->NB > 0) {
					f_status = fread(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, input_LF);
				}
			}
			else{
				/* we don't use LS weights but something derived on geometric distance in view array*/
				getGeomWeight(SAI, LF, SAI->stdd);
			}

			/* merge color */
			mergeWarped_N(warped_color_views, DispTargs, SAI, 3);

			/* hole filling for color*/
			holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
			
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

		if (SAI->NNt > 0 && SAI->Ms > 0)
		{

			SAI->sparse_weights = new int32_t[SAI->Ms]();
			SAI->sparse_mask = new unsigned char[SAI->Ms]();

			f_status = fread(SAI->sparse_weights, sizeof(int32_t), SAI->Ms, input_LF);
			f_status = fread(SAI->sparse_mask, sizeof(unsigned char), SAI->Ms, input_LF);

			applyGlobalSparseFilter(SAI);
		}


		int n_bytes_color_residual = 0, n_bytes_depth_residual = 0;

		f_status = fread(&n_bytes_color_residual, sizeof(int), 1,  input_LF);

		/* get residual */
		if (n_bytes_color_residual > 0)
		{

			int ncomp1 = 0; /* temporary to hold the number of components */

			char ppm_residual_path[512];

			char jp2_residual_path_jp2[512];

			sprintf(ppm_residual_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.ppm");

			sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.jp2");

			unsigned char *jp2_residual = new unsigned char[n_bytes_color_residual]();
			f_status = fread(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, input_LF);

			FILE *jp2_res_file;
			jp2_res_file = fopen(jp2_residual_path_jp2, "wb");
			fwrite(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, jp2_res_file);
			fclose(jp2_res_file);

			delete[](jp2_residual);

			decodeResidualJP2(SAI->color, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path, ncomp1, pow(2, 10) - 1, pow(2, 10) - 1);

			
		}

		f_status = fread(&n_bytes_depth_residual, sizeof(int), 1, input_LF);

		if (n_bytes_depth_residual>0) { /* residual depth if needed */

			int ncomp1 = 0; /* temporary to hold the number of components */

			char pgm_residual_depth_path[512];

			char jp2_residual_depth_path_jp2[512];

			sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.pgm");

			sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.jp2");

			unsigned char *jp2_depth_residual = new unsigned char[n_bytes_depth_residual]();
			f_status = fread(jp2_depth_residual, sizeof(unsigned char), n_bytes_depth_residual, input_LF);

			FILE *jp2_depth_res_file;
			jp2_depth_res_file = fopen(jp2_residual_depth_path_jp2, "wb");
			fwrite(jp2_depth_residual, sizeof(unsigned char), n_bytes_depth_residual, jp2_depth_res_file);
			fclose(jp2_depth_res_file);

			delete[](jp2_depth_residual);

			decodeResidualJP2(SAI->depth, kdu_expand_path, jp2_residual_depth_path_jp2, pgm_residual_depth_path, ncomp1, 0, pow(2, 16) - 1);

		}

		/* medfilt depth */
		unsigned short *tmp_depth = new unsigned short[SAI->nr*SAI->nc]();
		int startt = clock();
		medfilt2D(SAI->depth, tmp_depth, 3, SAI->nr, SAI->nc);
		std::cout << "time elapsed in depth median filtering\t" << (int)clock() - startt << "\n";
		memcpy(SAI->depth, tmp_depth, sizeof(unsigned short)*SAI->nr*SAI->nc);
		delete[](tmp_depth);

		sprintf(SAI->path_out_ppm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".ppm");
		sprintf(SAI->path_out_pgm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".pgm");

		aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
		aux_write16PGMPPM(SAI->path_out_pgm, SAI->nc, SAI->nr, 1, SAI->depth);

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

	}

	fclose(input_LF);

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