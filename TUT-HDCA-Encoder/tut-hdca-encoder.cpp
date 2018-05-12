#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include <cstdint>

#include "../include/gen_types.hh"
#include "../include/warpingFunctions.hh"

void FiveRefHierarchy_2_disk(const char *input_dir)
{
	const int ref_cols[] = { 2, 98, 2, 98, 50 };
	const int ref_rows[] = { 0, 0, 20, 20, 10 };

	bool vmask[256][256];
	//bool wmask[256][256];

	for (int rr = 0; rr < 256; rr++){
		for (int cc = 0; cc < 256; cc++){
			vmask[rr][cc] = false;
			//wmask[rr][cc] = false;
		}
	}

	//for (int rr = 0; rr < 11; rr++){
	//	for (int cc = 0; cc < 33; cc++){

	//		int r, c;
	//		r = rr * 2;
	//		c = cc * 3 + 2;

	//		wmask[r][c] = true;

	//	}
	//}

	char hname[256];
	sprintf(hname, "%s%s", input_dir, "five.h");

	FILE *filept = fopen(hname, "wb");

	int val = 33*11;

	fwrite(&val, sizeof(int), 1, filept);


	int ee_mask[21][101];
	int ee = 0;


	for (int ik = 0; ik < 5; ik++)
	{
		fwrite(&ref_rows[ik], sizeof(int), 1, filept);
		
		fwrite(&ref_cols[ik], sizeof(int), 1, filept);

		int val = 0;

		fwrite(&val, sizeof(int), 1, filept);

		vmask[ref_rows[ik]][ref_cols[ik]] = true;

		ee_mask[ref_rows[ik]][ref_cols[ik]] = ee;
		ee++;

	}

	for (int rr = 0; rr < 11; rr++){
		for (int cc = 0; cc < 33; cc++){

			int r, c;
			r = rr * 2;
			c = cc * 3 + 2;

			if (!vmask[r][c]){

				fwrite(&r, sizeof(int), 1, filept);

				fwrite(&c, sizeof(int), 1, filept);

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

	FiveRefHierarchy_2_disk(input_dir);

	const char *difftest_call = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --toycbcr --psnr ";

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


	char hname[256];
	sprintf(hname, "%s%s", input_dir, "five.h");

	filept = fopen(hname, "rb");

	view *LF = new view[11 * 33]();

	int n_views_total;
	fread(&n_views_total, sizeof(int), 1, filept);

	int ee_mask[21][101];
	int ee = 0;

	for (int ikiv = 0; ikiv < n_views_total; ikiv++){

		fread(&( (LF + ikiv)->r), sizeof(int), 1, filept);
		fread(&( (LF + ikiv)->c), sizeof(int), 1, filept);

		ee_mask[(LF + ikiv)->r][(LF + ikiv)->c] = ee;
		ee++;

		fread(&((LF + ikiv)->n_references), sizeof(int), 1, filept);

		if ((LF + ikiv)->n_references > 0){

			(LF + ikiv)->references = new int[(LF + ikiv)->n_references]();

			int *ip = (LF + ikiv)->references;

			for (int uu = 0; uu < (LF + ikiv)->n_references; uu++){
				fread(ip + uu, sizeof(int), 1, filept);
			}

		}

		(LF + ikiv)->y = d_cen[(LF + ikiv)->r][(LF + ikiv)->c][1];
		(LF + ikiv)->x = d_cen[(LF + ikiv)->r][(LF + ikiv)->c][0];

		(LF + ikiv)->nr = nr;
		(LF + ikiv)->nc = nc;

	}


	for (int ii = 0; ii < 11 * 33; ii++){
	
		(LF + ii)->color = new unsigned short[(LF+ii)->nr*(LF+ii)->nc * 3]();
		(LF + ii)->depth = new unsigned short[(LF + ii)->nr*(LF + ii)->nc]();

		unsigned short *original_intermediate_view = NULL;
		unsigned short *original_depth_view = NULL;

		char path_input_ppm[160];
		sprintf(path_input_ppm, "%s%c%03d_%03d%s", input_dir, '/', (LF + ii)->c, (LF + ii)->r, ".ppm");

		int nc1, nr1, ncomp1;
		aux_read16PGMPPM(path_input_ppm, nc1, nr1, ncomp1, original_intermediate_view);

		char path_input_depth_pgm[160];
		sprintf(path_input_depth_pgm, "%s%c%03d_%03d%s", input_dir, '/', (LF + ii)->c, (LF + ii)->r, ".pgm");
		bool depth_file_exist = aux_read16PGMPPM(path_input_depth_pgm, nc1, nr1, ncomp1, original_depth_view);

		char path_out_ppm[160];
		sprintf(path_out_ppm, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, ".ppm");

		char path_out_pgm[160];
		sprintf(path_out_pgm, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, ".pgm");

		if ((LF + ii)->n_references > 0){

			/* holds partial warped views for ii */

			unsigned short **warped_color_views = new unsigned short*[(LF + ii)->n_references]();
			unsigned short **warped_depth_views = new unsigned short*[(LF + ii)->n_references]();
			float **DispTargs = new float*[(LF + ii)->n_references]();

			for (int ij = 0; ij < (LF + ii)->n_references; ij++)
			{

				int uu = (LF + ii)->references[ij];

				/* warp color AND warp depth */
				warpView0_to_View1(LF + uu, LF+ii, warped_color_views[ij], warped_depth_views[ij], DispTargs[ij]);

				char tmp_str[256];
				sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (LF + uu)->c, (LF + uu)->r, "_warped_to_", (LF + ii)->c, (LF + ii)->r,".ppm");
				aux_write16PGMPPM(tmp_str, nc, nr, 3, warped_color_views[ij]);

				//FILE *tmpf;
				//tmpf = fopen("G:/HEVC_HDCA/Berlin_verification_model/TMP_CPP/disptarg.float", "wb");
				//fwrite(DispTargs[ij], sizeof(float), 1080 * 1920, tmpf);
				//fclose(tmpf);

			}

			(LF + ii)->NCOEFFS = (pow(2, (LF + ii)->n_references)*(LF + ii)->n_references);
			(LF + ii)->merge_weights = new signed short[(LF + ii)->NCOEFFS]();
			(LF + ii)->merge_weights_float = new float[(LF + ii)->NCOEFFS]();

			/* get LS weights */

			getViewMergingLSWeights_N( (LF + ii), warped_color_views, DispTargs, original_intermediate_view);

			/* merge color */

			mergeWarped_N(warped_color_views, DispTargs, LF+ii, 3);

			aux_write16PGMPPM(path_out_ppm, (LF + ii)->nc, (LF + ii)->nr, 3, (LF + ii)->color);

			/* hole filling */
			holefilling((LF + ii)->color, 3, (LF + ii)->nr, (LF + ii)->nc,0);

			aux_write16PGMPPM(path_out_ppm, (LF + ii)->nc, (LF + ii)->nr, 3, (LF + ii)->color);

			/* merge depth */

			unsigned short *pp = warped_depth_views[0];
			float *pf = DispTargs[0];
			for (int ij = 0; ij < nr*nc; ij++){
				if (*(pf + ij) > -1){
					(LF + ii)->depth[ij] = *(pp + ij);
				}
				else{
					(LF + ii)->depth[ij] = 0;
				}
			}

			for (int ij = 0; ij < nr*nc; ij++){
				for (int uu = 1; uu < (LF + ii)->n_references; uu++){
					unsigned short *pp = warped_depth_views[uu];
					float *pf = DispTargs[uu];
					if ((*(pf + ij)>-1) && (*(pp + ij) > (LF + ii)->depth[ij]))
						(LF + ii)->depth[ij] = *(pp + ij);
				}
			}

			/* hole filling */
			holefilling((LF + ii)->depth, 1, (LF + ii)->nr, (LF + ii)->nc, 0);

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
		if (residual_rate > 0.0001)
		{

			/* COLOR residual here */

			FILE *residual_file;

			char ppm_residual_path[256];

			char jp2_residual_path_jp2[256];

			sprintf(ppm_residual_path, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_residual.ppm");

			sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_residual.jp2");

			encodeResidualJP2((LF + ii)->nr, (LF + ii)->nc, original_intermediate_view, (LF + ii)->color, ppm_residual_path,
				kdu_compress_path, jp2_residual_path_jp2, residual_rate,3,pow(2,10)-1);

			decodeResidualJP2((LF + ii)->color, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path,ncomp1,pow(2,10)-1,pow(2,10)-1);

			if (depth_file_exist){

				char pgm_residual_depth_path[256];

				char jp2_residual_depth_path_jp2[256];

				sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_depth_residual.pgm");

				sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_depth_residual.jp2");

				encodeResidualJP2((LF + ii)->nr, (LF + ii)->nc, original_depth_view, (LF + ii)->depth, pgm_residual_depth_path,
					kdu_compress_path, jp2_residual_depth_path_jp2, residual_rate, 1, 0);

				decodeResidualJP2((LF + ii)->depth, kdu_expand_path, jp2_residual_depth_path_jp2, pgm_residual_depth_path, ncomp1,0,pow(2,16)-1);

			}
		}

		aux_write16PGMPPM(path_out_ppm, (LF + ii)->nc, (LF + ii)->nr, 3, (LF + ii)->color);
		aux_write16PGMPPM(path_out_pgm, (LF + ii)->nc, (LF + ii)->nr, 1, (LF + ii)->depth);

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

			printf("%s\n", psnr_buffer);

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