#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include <cstdint>

#include "../include/gen_types.hh"
#include "../include/warpingFunctions.hh"

int main(int argc, char** argv) {

	const char* input_dir = argv[1];
	const char* output_dir = argv[2];

	//const char *kdu_compress_path = "\"C:/Program Files (x86)/Kakadu/kdu_compress.exe\"";
	//const char *kdu_expand_path = "\"C:/Program Files (x86)/Kakadu/kdu_expand.exe\"";

	const char *kdu_compress_path = argv[3];
	const char *kdu_expand_path = argv[4];

	const float ref_color_rate = atof(argv[5]);
	const float ref_depth_rate = atof(argv[6]);
	const float residual_rate = atof(argv[7]);

	const int nar = atoi(argv[8]);
	const int nac = atoi(argv[9]);

	const int nr = 1080;
	const int nc = 1920;

	const char *difftest_call = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --toycbcr --psnr ";

	char path_camera_centers[256];
	sprintf(path_camera_centers, "%s%s", input_dir, "/camera_displacement.bin");

	float d_cen[256][256][2];

	float *d_cen0 = new float[nar * nac * 2]();
	float *p = &d_cen0[0];

	for (int k = 0; k < 2; k++)
		for (int r = 0; r < nar; r++)
			for (int c = 0; c < nac; c++)
				d_cen[r][c][k] = *(p++);

	view *LF = new view[11 * 33]();

	int ikiv = 0;

	for (int rr = 0; rr < 11; rr++){
		for (int cc = 0; cc < 33; cc++){

			int r, c;
			r = rr * 2;
			c = cc * 3 + 2;

			(LF + ikiv)->r = r;
			(LF + ikiv)->c = c;

			(LF + ikiv)->n_references = 0;

			(LF + ikiv)->y = d_cen[r][c][1];
			(LF + ikiv)->x = d_cen[r][c][0];

			(LF + ikiv)->nr = nr;
			(LF + ikiv)->nc = nc;

			ikiv++;

		}
	}


	for (int ii = 0; ii < 11 * 33; ii++){
	
		(LF + ii)->color = new unsigned short[nr*nc * 3]();
		(LF + ii)->depth = new unsigned short[nr*nc]();

		unsigned short *original_intermediate_view = new unsigned short[(LF + ii)->nr*(LF + ii)->nc * 3]();

		char path_input_ppm[160];
		sprintf(path_input_ppm, "%s%c%03d_%03d%s", input_dir, '/', (LF + ii)->c, (LF + ii)->r, ".ppm");

		int nc1, nr1;
		aux_read16ppm(path_input_ppm, nc1, nr1, original_intermediate_view);

		if ((LF + ii)->n_references > 0){

			/* holds partial warped views for ii */

			unsigned short **warped_color_views = new unsigned short*[(LF + ii)->n_references]();
			unsigned short **warped_depth_views = new unsigned short*[(LF + ii)->n_references]();
			float **DispTargs = new float*[(LF + ii)->n_references]();

			for (int ij = 0; ij < (LF + ii)->n_references; ij++)
			{

				warped_color_views[ij] = new unsigned short[nr*nc * 3]();
				warped_depth_views[ij] = new unsigned short[nr*nc]();
				DispTargs[ij] = new float[nr*nc]();

				int uu = (LF + ii)->references[ij];

				/* warp color AND warp depth */
				warpView0_2_View1(LF + ii, LF + uu, warped_color_views[ij], warped_depth_views[ij], DispTargs[ij]);

			}

			int NCOEFFS = (pow(2, (LF + ii)->n_references)*(LF + ii)->n_references) / 2;

			(LF + ii)->merge_weights = new signed short[NCOEFFS]();
			(LF + ii)->merge_weights_float = new float[NCOEFFS]();

			/* get LS weights */

			getViewMergingLSWeights_N((LF + ii)->n_references, warped_color_views, DispTargs, original_intermediate_view,
				(LF + ii)->nr, (LF + ii)->nc, (LF + ii)->merge_weights);

			/* merge color */

			for (int ij = 0; ij < NCOEFFS; ij++)
				(LF + ii)->merge_weights_float[ij] = ((float)(LF + ii)->merge_weights[ij]) / pow(2, 14);

			mergeWarpedLS_N(warped_color_views, DispTargs, (LF + ii)->nr, (LF + ii)->nc, 3, (LF + ii)->n_references, (LF + ii)->merge_weights_float);

			/* hole filling */
			holefilling((LF + ii)->color, 3, (LF + ii)->nr, (LF + ii)->nc);

			/* merge depth */

			for (int ij = 0; ij < nr*nc; ij++)
				(LF + ii)->depth[ij] = 9999;

			for (int ij = 0; ij < nr*nc; ij++){
				for (int uu = 0; uu < (LF + ii)->n_references; uu++){
					unsigned short *pp = warped_depth_views[uu];
					if (*(pp + ij) < (LF + ii)->depth[ij])
						(LF + ii)->depth[ij] = *(pp + ij);
				}
			}

			for (int ij = 0; ij < nr*nc; ij++)
				if ((LF + ii)->depth[ij] == 9999)
					(LF + ii)->depth[ij] = 0;

			/* hole filling */
			holefilling((LF + ii)->depth, 1, (LF + ii)->nr, (LF + ii)->nc);

			/* clean */
			for (int ij = 0; ij < (LF + ii)->n_references; ij++)
			{
				delete[](warped_color_views[ij]);
				delete[](warped_depth_views[ij]);
			}
			delete[](warped_color_views);
			delete[](warped_depth_views);

		}

		/* get sparse weights */
		/* get residual */

		char path_out_ppm[160];
		sprintf(path_out_ppm, "%s%c%03d_%03d%s", output_dir, '/', (LF+ii)->c, (LF+ii)->r, ".ppm");

		if (residual_rate > 0.0001)
		{

			/* residual here */

			FILE *residual_file;

			char ppm_residual_path[256];

			char jp2_residual_path_jp2[256];

			sprintf(ppm_residual_path, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_residual.ppm");

			sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_residual.jp2");

			/*establish residual*/
			unsigned short *residual_image = new unsigned short[(LF + ii)->nr*(LF + ii)->nc * 3]();

			unsigned short *ps = (LF + ii)->color;

			for (int iir = 0; iir < ((LF + ii)->nr) * ((LF + ii)->nc) * 3; iir++){
				unsigned short res_val = (unsigned short)(((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + (pow(2, BIT_DEPTH) - 1));
				if (res_val>pow(2, BIT_DEPTH_RESIDUAL) - 1)
					res_val = pow(2, BIT_DEPTH_RESIDUAL) - 1;
				if (res_val < 0)
					res_val = 0;
				*(residual_image + iir) = res_val;
			}

			aux_write16ppm(ppm_residual_path, (LF + ii)->nc, (LF + ii)->nr, residual_image);

			delete[](residual_image);

			/* here encode residual with kakadu */

			char kdu_compress_s[256];
			sprintf(kdu_compress_s, "\"%s\"%s%s%s%s%s%f", kdu_compress_path, " -i ", ppm_residual_path, " -o ", jp2_residual_path_jp2, " -no_weights -precise -full -rate ", residual_rate);

			//std::cout << kdu_compress_s << "\n";

			int status = system_1(kdu_compress_s);

			/* decode residual with kakadu */
			char kdu_expand_s[256];
			sprintf(kdu_expand_s, "\"%s\"%s%s%s%s", kdu_expand_path, " -i ", jp2_residual_path_jp2, " -o ", ppm_residual_path);

			//std::cout << kdu_expand_s << "\n";

			status = system_1(kdu_expand_s);

			/* apply residual */

			unsigned short* jp2_residual;

			if (aux_read16ppm(ppm_residual_path, nc1, nr1, jp2_residual))
			{

				for (int iir = 0; iir < nr*nc * 3; iir++)
				{
					signed int val = ((signed int)*(ps + iir)) + ((signed int)jp2_residual[iir]) - (pow(2, BIT_DEPTH) - 1);
					if (val < 0)
						val = 0;
					if (val >(pow(2, BIT_DEPTH) - 1))
						val = pow(2, BIT_DEPTH) - 1;
					*(ps + iir) = (unsigned short)(val);
				}

				delete[](jp2_residual);
			}
		}

		aux_write16ppm(path_out_ppm, (LF + ii)->nc, (LF + ii)->nr, (LF + ii)->color);

		if (1){

			/* run psnr here */

			char tmp_original_intermediate_view[256];
			sprintf(tmp_original_intermediate_view, "%s%s", output_dir, "tmp_iv.ppm");
			aux_write16ppm(tmp_original_intermediate_view, (LF + ii)->nc, (LF + ii)->nr, original_intermediate_view);

			char psnr_call[256];
			sprintf(psnr_call, "%s%s%s%s", difftest_call, path_out_ppm, " ", tmp_original_intermediate_view);

			FILE *pfile;
			pfile = _popen(psnr_call, "r");

			char psnr_buffer[256];
			while (fgets(psnr_buffer, sizeof(psnr_buffer), pfile) != 0) {
				/*...*/
			}
			_pclose(pfile);

			printf("PSNR %s\n", psnr_buffer);

		}
	
		delete[](original_intermediate_view);

	}

	for (int ii = 0; ii < 11 * 33; ii++){

		delete[]((LF + ii)->color);
		delete[]((LF + ii)->depth);
		delete[]((LF + ii)->merge_weights);
		delete[]((LF + ii)->sparse_weights);
	}


	exit(0);
}