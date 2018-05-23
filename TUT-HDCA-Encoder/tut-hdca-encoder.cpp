#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include <cstdint>

#include <chrono> //C++11
#include <ctime>

#include "../include/gen_types.hh"
#include "../include/warpingFunctions.hh"


void FiveRefHierarchy_2_disk(const char *hiearchy_file)
{
	const int ref_cols[] = { 2, 98, 2, 98, 50 };
	const int ref_rows[] = { 0, 0, 20, 20, 10 };

	bool vmask[256][256];

	for (int rr = 0; rr < 256; rr++){
		for (int cc = 0; cc < 256; cc++){
			vmask[rr][cc] = false;
			//wmask[rr][cc] = false;
		}
	}

	FILE *filept = fopen(hiearchy_file, "wb");

	int val = 33*11;

	fwrite(&val, sizeof(int), 1, filept);

	int ee_mask[21][101];
	int ee = 0;

	for (int ik = 0; ik < 5; ik++)
	{
		fwrite(&ref_rows[ik], sizeof(int), 1, filept);
		
		fwrite(&ref_cols[ik], sizeof(int), 1, filept);

		int rate = 100000/2;

		fwrite(&rate, sizeof(int), 1, filept);
		fwrite(&rate, sizeof(int), 1, filept);

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

	view *LF = new view[n_views_total]();

	int ee_mask[256][256];
	int ee = 0;

	for (int ikiv = 0; ikiv < n_views_total; ikiv++){

		fread(&( (LF + ikiv)->r), sizeof(int), 1, filept);
		fread(&( (LF + ikiv)->c), sizeof(int), 1, filept);

		int rate_color, rate_depth;

		fread(&rate_color, sizeof(int), 1, filept);
		fread(&rate_depth, sizeof(int), 1, filept);

		(LF + ikiv)->residual_rate_color = ((float)rate_color) / 100000;
		(LF + ikiv)->residual_rate_depth = ((float)rate_depth) / 100000;

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

	}


	for (int ii = 0; ii < n_views_total; ii++){

		char output_results[1024];
		int output_buffer_length = 0;
		output_buffer_length += sprintf(output_results+output_buffer_length, "%03d\t%03d", (LF + ii)->r, (LF + ii)->c);
	
		unsigned short *original_color_view = NULL;
		unsigned short *original_depth_view = NULL;

		char path_input_ppm[160];
		sprintf(path_input_ppm, "%s%c%03d_%03d%s", input_dir, '/', (LF + ii)->c, (LF + ii)->r, ".ppm");

		int nc1, nr1, ncomp1;
		aux_read16PGMPPM(path_input_ppm, (LF + ii)->nc, (LF + ii)->nr, ncomp1, original_color_view);

		char path_input_depth_pgm[160];
		sprintf(path_input_depth_pgm, "%s%c%03d_%03d%s", input_dir, '/', (LF + ii)->c, (LF + ii)->r, ".pgm");
		bool depth_file_exist = aux_read16PGMPPM(path_input_depth_pgm, nc1, nr1, ncomp1, original_depth_view);

		char path_out_ppm[160];
		sprintf(path_out_ppm, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, ".ppm");

		char path_out_pgm[160];
		sprintf(path_out_pgm, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, ".pgm");

		(LF + ii)->color = new unsigned short[(LF + ii)->nr*(LF + ii)->nc * 3]();
		(LF + ii)->depth = new unsigned short[(LF + ii)->nr*(LF + ii)->nc]();

		if ((LF + ii)->n_references > 0){

			/* holds partial warped views for ii */

			unsigned short **warped_color_views = new unsigned short*[(LF + ii)->n_references]();
			unsigned short **warped_depth_views = new unsigned short*[(LF + ii)->n_references]();
			float **DispTargs = new float*[(LF + ii)->n_references]();

			/* currently we forward warp the depth from the 5 references*/
			unsigned short **warped_color_views_0_5 = new unsigned short*[(LF + ii)->n_references]();
			unsigned short **warped_depth_views_0_5 = new unsigned short*[(LF + ii)->n_references]();
			float **DispTargs_0_5 = new float*[(LF + ii)->n_references]();

			for (int ij = 0; ij < 5; ij++)
			{
				warpView0_to_View1(LF + ij, LF + ii, warped_color_views_0_5[ij], warped_depth_views_0_5[ij], DispTargs_0_5[ij]);
			}

			for (int ij = 0; ij < (LF + ii)->n_references; ij++)
			{

				int uu = (LF + ii)->references[ij];

				/* FORWARD warp color AND depth */
				warpView0_to_View1(LF + uu, LF+ii, warped_color_views[ij], warped_depth_views[ij], DispTargs[ij]);

				char tmp_str[256];

				sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (LF + uu)->c, (LF + uu)->r, "_warped_to_", (LF + ii)->c, (LF + ii)->r,".ppm");
				aux_write16PGMPPM(tmp_str, (LF + ii)->nc, (LF + ii)->nr, 3, warped_color_views[ij]);

				sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (LF + uu)->c, (LF + uu)->r, "_warped_to_", (LF + ii)->c, (LF + ii)->r, ".pgm");
				aux_write16PGMPPM(tmp_str, (LF + ii)->nc, (LF + ii)->nr, 1, warped_depth_views[ij]);

				//FILE *tmpf;
				//tmpf = fopen("G:/HEVC_HDCA/Berlin_verification_model/TMP_CPP/disptarg.float", "wb");
				//fwrite(DispTargs[ij], sizeof(float), 1080 * 1920, tmpf);
				//fclose(tmpf);

			}

			/* filtering step to remove erroneus pixels from merging */
			/*for (int ij = 0; ij < (LF + ii)->n_references; ij++)
			{

				unsigned short *tmpd = new unsigned short[(LF + ii)->nr*(LF + ii)->nc];
				medfilt2D(warped_depth_views[ij], tmpd, 3, (LF + ii)->nr, (LF + ii)->nc);

				unsigned short *pp = warped_depth_views[ij];
				float *pf = DispTargs[ij];

				for (int jj = 0; jj < (LF + ii)->nr*(LF + ii)->nc; jj++) {
					if (*(pf + jj) > -1 && abs((float)(*(tmpd + jj)) - (float)(*(pp + jj))) >(float)(*(pp + jj))*0.75)
					{
						*(pf + jj) = -1;
					}

				}

				delete[](tmpd);

			}*/


			/* get LS weights */

			getViewMergingLSWeights_N( (LF + ii), warped_color_views, DispTargs, original_color_view);

			/* merge color */

			mergeWarped_N(warped_color_views, DispTargs, LF+ii, 3);

			//aux_write16PGMPPM(path_out_ppm, (LF + ii)->nc, (LF + ii)->nr, 3, (LF + ii)->color);

			/* hole filling for color*/
			holefilling((LF + ii)->color, 3, (LF + ii)->nr, (LF + ii)->nc,0);

			//aux_write16PGMPPM(path_out_ppm, (LF + ii)->nc, (LF + ii)->nr, 3, (LF + ii)->color);

			/* merge depth with minimum winning*/
			//unsigned short *pp = warped_depth_views[0];
			//float *pf = DispTargs[0];
			//for (int ij = 0; ij < (LF + ii)->nr*(LF + ii)->nc; ij++){
			//	if (*(pf + ij) > -1){
			//		(LF + ii)->depth[ij] = *(pp + ij);
			//	}
			//	else{
			//		(LF + ii)->depth[ij] = 0;
			//	}
			//}
			//for (int ij = 0; ij < (LF + ii)->nr*(LF + ii)->nc; ij++) {
			//	for (int uu = 1; uu < (LF + ii)->n_references; uu++) {
			//		unsigned short *pp = warped_depth_views[uu];
			//		float *pf = DispTargs[uu];
			//		if ((*(pf + ij)>-1) && (*(pp + ij) >(LF + ii)->depth[ij]))
			//			(LF + ii)->depth[ij] = *(pp + ij);
			//	}
			//}

			/* merge depth with median*/

			int startt = clock();

			#pragma omp parallel for
			for (int ij = 0; ij < (LF + ii)->nr*(LF + ii)->nc; ij++) {
				std::vector<unsigned short> depth_values;
				for (int uu = 0; uu < 5; uu++) {
				//for (int uu = 0; uu < (LF + ii)->n_references; uu++) {
					unsigned short *pp = warped_depth_views_0_5[uu];
					float *pf = DispTargs_0_5[uu];
					if (*(pf + ij) > -1) {
						depth_values.push_back(*(pp + ij));
					}
				}
				if( depth_values.size()>0)
					(LF + ii)->depth[ij] = getMedian(depth_values);
			}

			std::cout << "time elapsed in depth merging\t" << (int)clock() - startt << "\n";

			/* hole filling for depth */
			holefilling((LF + ii)->depth, 1, (LF + ii)->nr, (LF + ii)->nc, 0);

			/* clean */
			for (int ij = 0; ij < (LF + ii)->n_references; ij++)
			{
				delete[](warped_color_views[ij]);
				delete[](warped_depth_views[ij]);
				delete[](DispTargs[ij]);
			}
			delete[](warped_color_views_0_5);
			delete[](warped_depth_views_0_5);
			delete[](DispTargs_0_5);

			delete[](warped_color_views);
			delete[](warped_depth_views);
			delete[](DispTargs);
		}


		float psnr_result;

		if ((LF + ii)->n_references > 0) {

			aux_write16PGMPPM(path_out_ppm, (LF + ii)->nc, (LF + ii)->nr, 3, (LF + ii)->color);
			psnr_result = getPSNR(NULL, LF + ii, path_out_ppm, path_input_ppm, difftest_call);

			output_buffer_length += sprintf(output_results+ output_buffer_length, "\t%f", psnr_result);

		}
		else {
			output_buffer_length += sprintf(output_results+ output_buffer_length, "\t%f", 0);
		}

		if ((LF + ii)->NNt > 0 && (LF + ii)->Ms > 0 && ii>4){ /*we need to put Ms and NNt in to the input config and remove ii>5*/

			int startt = clock();

			getGlobalSparseFilter(LF + ii, original_color_view);

			std::cout <<  "time elapsed in getGlobalSparseFilter()\t" << (int)clock()-startt << "\n";

			applyGlobalSparseFilter(LF + ii);

			aux_write16PGMPPM(path_out_ppm, (LF + ii)->nc, (LF + ii)->nr, 3, (LF + ii)->color);
			psnr_result = getPSNR(NULL, LF + ii, path_out_ppm, path_input_ppm, difftest_call);
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_result);
		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0);
		}
		
		/* get residual */
		if ((LF + ii)->residual_rate_color > 0)
		{

			/* COLOR residual here */

			FILE *residual_file;

			char ppm_residual_path[512];

			char jp2_residual_path_jp2[512];

			sprintf(ppm_residual_path, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_residual.ppm");

			sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_residual.jp2");

			encodeResidualJP2((LF + ii)->nr, (LF + ii)->nc, original_color_view, (LF + ii)->color, ppm_residual_path,
				kdu_compress_path, jp2_residual_path_jp2, (LF + ii)->residual_rate_color, 3, pow(2, 10) - 1);

			decodeResidualJP2((LF + ii)->color, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path,ncomp1,pow(2,10)-1,pow(2,10)-1);

			if (depth_file_exist && (LF + ii)->residual_rate_depth>0){

				char pgm_residual_depth_path[512];

				char jp2_residual_depth_path_jp2[512];

				sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_depth_residual.pgm");

				sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', (LF + ii)->c, (LF + ii)->r, "_depth_residual.jp2");

				encodeResidualJP2((LF + ii)->nr, (LF + ii)->nc, original_depth_view, (LF + ii)->depth, pgm_residual_depth_path,
					kdu_compress_path, jp2_residual_depth_path_jp2, (LF + ii)->residual_rate_depth, 1, 0);

				decodeResidualJP2((LF + ii)->depth, kdu_expand_path, jp2_residual_depth_path_jp2, pgm_residual_depth_path, ncomp1,0,pow(2,16)-1);

			}
		}

		/* medfilt depth */
		unsigned short *tmp_depth = new unsigned short[(LF + ii)->nr*(LF + ii)->nc];
		int startt = clock();
		medfilt2D((LF + ii)->depth, tmp_depth, 3, (LF + ii)->nr, (LF + ii)->nc);
		std::cout << "time elapsed in depth median filtering\t" << (int)clock() - startt << "\n";
		memcpy((LF + ii)->depth, tmp_depth, sizeof(unsigned short)*(LF + ii)->nr*(LF + ii)->nc);

		aux_write16PGMPPM(path_out_ppm, (LF + ii)->nc, (LF + ii)->nr, 3, (LF + ii)->color);
		aux_write16PGMPPM(path_out_pgm, (LF + ii)->nc, (LF + ii)->nr, 1, (LF + ii)->depth);

		if (depth_file_exist) {
			psnr_result = getPSNR(NULL, LF + ii, path_out_pgm, path_input_depth_pgm, difftest_call_pgm);
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_result);
		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0);
		}


		psnr_result = getPSNR(NULL,LF + ii, path_out_ppm, path_input_ppm, difftest_call);

		output_buffer_length += sprintf(output_results+ output_buffer_length, "\t%f", psnr_result);

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
		fprintf(output_results_file, "%s\n",output_results);
		fclose(output_results_file);
	}

	for (int ii = 0; ii < 11 * 33; ii++){

		delete[]((LF + ii)->color);
		delete[]((LF + ii)->depth);
		delete[]((LF + ii)->merge_weights);
		delete[]((LF + ii)->sparse_weights);
		delete[]((LF + ii)->sparse_mask);
	}


	exit(0);
}