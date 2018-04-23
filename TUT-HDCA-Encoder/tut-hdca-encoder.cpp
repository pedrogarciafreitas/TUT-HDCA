#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include <cstdint>

#include "../include/gen_types.hh"
#include "../include/warpingFunctions.hh"

int main(int argc, char** argv) {

	const char* path_to_references = argv[1];
	const char* path_camera_centers = argv[2];
	const char* path_out_dir = argv[3];
	const char* path_to_original_views = argv[4];

	const char *kdu_compress_path = "\"C:/Program Files (x86)/Kakadu/kdu_compress.exe\"";
	const char *kdu_expand_path = "\"C:/Program Files (x86)/Kakadu/kdu_expand.exe\"";

	const int nr = 1080;
	const int nc = 1920;
	const int r_crop = 209 + 540;
	const int c_crop = 85 + 960;

	int nr0 = 0, nc0 = 0;

	float *qDMF[5];
	int *qDM[5];

	float *p, *pp, *ppp;

	const int ref_cols[] = { 2, 98, 2, 98, 50 };
	const int ref_rows[] = { 0, 0, 20, 20, 10 };

	float *d_cen0 = new float[101 * 21 * 2]();

	FILE *filept;

	filept = fopen(path_camera_centers, "rb");
	fread(d_cen0, sizeof(float), 101 * 21 * 2, filept);
	fclose(filept);

	bool *bmask = new bool[32 * 5]();

	for (int ij = 0; ij < 32; ij++){

		int uu = ij;

		for (int ik = 4; ik >= 0; ik--){

			if (floor(uu / pow(2, ik)) > 0)
			{
				uu = uu - pow(2, ik);
				bmask[ij + ik * 32] = 1;
			}

		}
	}


	float d_cen[21][101][2];

	p = &d_cen0[0];

	for (int k = 0; k < 2; k++)
		for (int r = 0; r < 21; r++)
			for (int c = 0; c < 101; c++)
				d_cen[r][c][k] = *(p++);


	bool ref_array[11][33];

	for (int r = 0; r < 11; r++)
		for (int c = 0; c < 33; c++)
			ref_array[r][c] = false;


	char pathLSW[160];
	sprintf(pathLSW, "%s%s", path_to_references, "LS_weights");

	char pathSPARSEGLOBAL[256];
	sprintf(pathSPARSEGLOBAL, "%s%s", path_to_references, "/sparse_global_weights");

	std::cout << pathLSW << "\n";
	std::cout << pathSPARSEGLOBAL << "\n";

	unsigned short NNt, Ms;


	FILE *fileptLS;
	fileptLS = fopen(pathLSW, "rb");

	FILE *fileptSPARSEGLOBAL;
	fileptSPARSEGLOBAL = fopen(pathSPARSEGLOBAL, "rb");

	if (!(fileptSPARSEGLOBAL == NULL)){
		fread(&NNt, sizeof(unsigned short), 1, fileptSPARSEGLOBAL);
		fread(&Ms, sizeof(unsigned short), 1, fileptSPARSEGLOBAL);
	}

	std::cout << NNt << "\t" << Ms << "\n";

	if (1)
	{

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

			ColDisps[ij] = new float[nr*nc]();
			RowDisps[ij] = new float[nr*nc]();
			colorViews[ij] = new unsigned short[nr*nc * 3]();
			warpedColorViews[ij] = new unsigned short[nr*nc * 3]();
			DispTargs[ij] = new float[nr*nc]();

			qDM[ij] = new int[nr*nc]();
			qDMF[ij] = new float[nr*nc]();

			char pathtoref[160];
			sprintf(pathtoref, "%s%c%03d_%03d%s", path_to_references, '/', ref_cols[ij], ref_rows[ij], ".ppm");

			filept = fopen(pathtoref, "rb");
			aux_read16ppm(filept, nc0, nr0, colorViews[ij]);
			fclose(filept);

			char pgm_filename[12];
			sprintf(pgm_filename, "%s%c%03d_%03d%s", path_to_references, '/', ref_cols[ij], ref_rows[ij], "_drf.pgm");
			filept = fopen(pgm_filename, "rb");
			aux_read16pgm(filept, qDM[ij]);
			fclose(filept);

			float *pf = qDMF[ij];
			int *p = qDM[ij];

			for (int ik = 0; ik < nr*nc; ik++){
				int tmpi = *(p + ik);
				*(pf + ik) = ((float)tmpi) / pow(2, 14);
			}

		}

		/* Encoding of color and disparity here or before */

		unsigned short *original_intermediate_view = new unsigned short[nr*nc * 3]();

		for (int rr = 0; rr < 11; rr++){
			for (int cc = 0; cc < 33; cc++){

				int r, c;
				r = rr * 2;
				c = cc * 3 + 2;

				char path_out_ppm[160];
				sprintf(path_out_ppm, "%s%c%03d_%03d%s", path_out_dir, '/', c, r, ".ppm");

				char path_input_ppm[160];
				sprintf(path_input_ppm, "%s%c%03d_%03d%s", path_to_original_views, '/', c, r, ".ppm");

				std::cout << "READING " << path_input_ppm << "\n";

				/* Read the original intermediate view to be encoded */
				int nc1, nr1;
				unsigned short *im_read;
				filept = fopen(path_input_ppm, "rb");
				aux_read16ppm_1(filept, nc1, nr1, im_read);
				fclose(filept);

				memset(original_intermediate_view, 0x00, sizeof(unsigned short)*nr*nc * 3);

				/* crop the image to [nr,nc] with offset r_crop,c_crop (y,x) */
				for (int x = 0; x < nc; x++){
					for (int y = 0; y < nr; y++){
						for (int icomp = 0; icomp < 3; icomp++){
							*(original_intermediate_view + y + x*nr + icomp*nr*nc) = *(im_read + y + r_crop + (x + c_crop)*nr1 + icomp*nr1*nc1) * 64; // *64 to make range 16bit
						}

					}
				}

				delete[](im_read);

				if (!ref_array[rr][cc]){

					printf("Decoding to %s\n ", path_out_ppm);

					float ddy[5];
					float ddx[5];

					std::vector<std::pair<float, int>> view_distances(5);
					std::vector<std::pair<float, float>> dydx(5);

					//std::cout << d_cen[r][c][0] << "\t" << d_cen[r][c][1] << "\t\n\n";

					for (int ij = 0; ij < 5; ij++){

						ddy[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][1] - d_cen[r][c][1];
						ddx[ij] = d_cen[ref_rows[ij]][ref_cols[ij]][0] - d_cen[r][c][0];

						dydx.at(ij).first = ref_rows[ij] - r;
						dydx.at(ij).second = ref_cols[ij] - c;

						view_distances.at(ij).first = sqrt(dydx.at(ij).first * dydx.at(ij).first + dydx.at(ij).second * dydx.at(ij).second);
						view_distances.at(ij).second = ij;
					}

					sort(view_distances.begin(), view_distances.end());


					for (int uu = 0; uu < 5; uu++){

						int ij = view_distances.at(uu).second;

						p = qDMF[ij];
						pp = RowDisps[ij];
						ppp = ColDisps[ij];

						for (int ii = 0; ii < nr*nc; ii++){
							pp[ii] = -p[ii] * ddy[ij];
							ppp[ii] = p[ii] * ddx[ij];
						}

						/* to zero for consistency */
						memset(warpedColorViews[uu], 0x00, sizeof(unsigned short)*nr*nc * 3);
						memset(DispTargs[uu], 0x00, sizeof(float)*nr*nc);

						warpColorView(colorViews[ij], RowDisps[ij], ColDisps[ij], nr, nc, warpedColorViews[uu], DispTargs[uu]);

					}

					signed short *LSw_s = new signed short[80]();

					float *LSw = new float[32 * 5]();

					float stdd = 10.0;

					/* LS merging here, we need fastOLS for solving the LS weights */

					//float *LScoeffs;

					viewMergingLSWeights(warpedColorViews, DispTargs, original_intermediate_view,
						nr, nc, bmask, LSw_s);

					/*delete[](LScoeffs);*/


					/*if (fileptLS == NULL || !(fread(LSw_s, sizeof(signed short), 80, fileptLS) == 80))
					{

					std::cout << "Fixed exponential weights, sigma = " << stdd << "\n";

					for (int ii = 0; ii < 32; ii++){
					float sumw = 0;
					for (int ij = 0; ij < 5; ij++){
					if (bmask[ii + ij * 32]){
					LSw[ii + ij * 32] = exp(-(view_distances.at(ij).first * view_distances.at(ij).first) / (2 * stdd*stdd));
					sumw = sumw + LSw[ii + ij * 32];
					}
					else{
					LSw[ii + ij * 32] = 0.0;
					}

					}
					for (int ij = 0; ij < 5; ij++)
					if (bmask[ii + ij * 32])
					LSw[ii + ij * 32] = LSw[ii + ij * 32] / sumw;
					}

					}
					else{*/

					std::cout << "LS weights " << pathLSW << "\n";

					int uu = 0;

					for (int ii = 0; ii < 32 * 5; ii++){
						if (bmask[ii]){
							LSw[ii] = ((float)LSw_s[uu++]) / pow(2.0, 14.0);
							//std::cout << LSw[ii] << "\n";
						}
						else{
							LSw[ii] = 0.0;
						}
					}
					//}


					//for (int ii = 0; ii < 32 * 5; ii++)
					//	std::cout << LSw[ii] << "\n";

					/* results are collected to warpedColorViews[0] and DispTargs[0] */
					/*here we need the LS weights*/
					collectWarpedLS(warpedColorViews, DispTargs, nr, nc, 3, 5, LSw);

					unsigned short *pshort = warpedColorViews[0];

					std::vector<unsigned short> neighbours;
					int dsz = 1;

					for (int ii = 0; ii < nr*nc; ii++){
						int y, x;
						y = ii%nr;
						x = floor(ii / nr);
						if ((pshort[ii] < 1) && (pshort[ii + nr*nc] < 1) && (pshort[ii + 2 * nr*nc] < 1)){
							for (int icomp = 0; icomp < 3; icomp++){
								neighbours.clear();
								for (int dy = -dsz; dy <= dsz; dy++){
									for (int dx = -dsz; dx <= dsz; dx++){
										if (!(dy == 0 && dx == 0)){
											if ((y + dy) >= 0 && (y + dy) < nr && (x + dx) >= 0 && (x + dx) < nc){
												//if (pshort[y + dy + (x + dx)*nr + icomp*nr*nc] > 0)
												if (!((pshort[y + dy + (x + dx)*nr] < 1) && (pshort[y + dy + (x + dx)*nr + nr*nc] < 1) && (pshort[y + dy + (x + dx)*nr + 2 * nr*nc] < 1)))
													neighbours.push_back(pshort[y + dy + (x + dx)*nr + icomp*nr*nc]);
											}
										}
									}
								}
								if (neighbours.size() > 0)
									pshort[ii + icomp*nr*nc] = getMedian(neighbours);
							}
						}
					}


					//filept = fopen(path_out_ppm, "wb");
					//aux_write16ppm_16(filept, nc, nr, warpedColorViews[0]);
					//fclose(filept);

					//exit(0);

					/* global sparse here */
					int NNt = 3; // window -NNt:NNt
					int Npp = (nr - NNt * 2)*(nc - NNt * 2) * 3;
					int Npp0 = Npp / 3;

					double *AA = new double[Npp*((NNt * 2 + 1)*(NNt * 2 + 1) + 1)]();
					double *Yd = new double[Npp]();

					for (int ii = 0; ii < Npp; ii++)
						*(AA + ii + (NNt * 2 + 1)*(NNt * 2 + 1)*Npp) = 1.0;

					int iiu = 0;

					for (int ir = NNt; ir < nr - NNt; ir++){
						for (int ic = NNt; ic < nc - NNt; ic++){
							int ai = 0;
							for (int dy = -NNt; dy <= NNt; dy++){
								for (int dx = -NNt; dx <= NNt; dx++){
									for (int icomp = 0; icomp < 3; icomp++){

										int offset = ir + dy + nr*(ic + dx) + icomp*nr*nc;

										/* get the desired Yd*/
										if (dy == 0 && dx == 0){
											*(Yd + iiu + icomp*Npp0) = ((double)*(original_intermediate_view + offset)) / 65535;
										}

										/* get the regressors */
										*(AA + iiu + icomp*Npp0 + ai*Npp) = ((double)*(pshort + offset)) / 65535;

									}
									ai++;
								}
							}
							iiu++;
						}
					}

					int *PredRegr0 = new int[(2 * NNt + 1)*(2 * NNt + 1) + 1]();
					double *PredTheta0 = new double[(2 * NNt + 1)*(2 * NNt + 1) + 1]();

					int Mtrue = FastOLS_new(AA, Yd, PredRegr0, PredTheta0, 25, (NNt * 2 + 1)*(NNt * 2 + 1) + 1, (NNt * 2 + 1)*(NNt * 2 + 1) + 1, Npp);

					delete[](AA);
					delete[](Yd);

					unsigned char *Regr0 = new unsigned char[Ms]();
					int32_t *theta0 = new int32_t[Ms]();

					double *theta = new double[(2 * NNt + 1)*(2 * NNt + 1) + 1]();

					//int n_regr = fread(Regr0, sizeof(unsigned char), Ms, fileptSPARSEGLOBAL);
					//int n_theta = fread(theta0, sizeof(int32_t), Ms, fileptSPARSEGLOBAL);

					for (int ii = 0; ii < Ms; ii++){
						*(Regr0 + ii) = ((unsigned char)*(PredRegr0 + ii) + 1);
						*(theta0 + ii) = (int32_t)round(*(PredTheta0 + ii) * pow(2, 20));
						/* now these should be added to bitstream */
					}

					//if (fileptSPARSEGLOBAL != NULL && n_regr == Ms && n_theta == Ms)
					//{

					for (int ii = 0; ii < Ms; ii++){
						if (Regr0[ii] > 0){
							theta[Regr0[ii] - 1] = ((double)theta0[ii]) / pow(2, 20);
							//std::cout << theta[Regr0[ii] - 1] << "\t";
						}
					}

					double *final_view = new double[nr*nc * 3]();

					pshort = warpedColorViews[0];

					for (int ii = 0; ii < nr*nc * 3; ii++)
						final_view[ii] = pshort[ii];

					for (int rr = NNt; rr < nr - NNt; rr++){
						for (int cc = NNt; cc < nc - NNt; cc++)
						{
							for (int icomp = 0; icomp < 3; icomp++)
								final_view[rr + cc*nr + icomp*nr*nc] = 0;

							int ee = 0;

							for (int dy = -NNt; dy <= NNt; dy++){
								for (int dx = -NNt; dx <= NNt; dx++){
									for (int icomp = 0; icomp < 3; icomp++){
										final_view[rr + cc*nr + icomp*nr*nc] = final_view[rr + cc*nr + icomp*nr*nc] + theta[ee] * ((double)pshort[rr + dy + (cc + dx)*nr + icomp*nr*nc]);
									}
									ee++;
								}
							}

							/* bias term */
							for (int icomp = 0; icomp < 3; icomp++){
								final_view[rr + cc*nr + icomp*nr*nc] = final_view[rr + cc*nr + icomp*nr*nc] + theta[(2 * NNt + 1)*(2 * NNt + 1)];
							}

						}
					}

					unsigned short *final_view_s = new unsigned short[nr*nc * 3]();

					for (int ii = 0; ii < nr*nc * 3; ii++){
						if (final_view[ii] < 0)
							final_view[ii] = 0;
						if (final_view[ii]> (pow(2, 16) - 1))
							final_view[ii] = (pow(2, 16) - 1);

						final_view_s[ii] = (unsigned short)(final_view[ii]);
					}

					memcpy(pshort, final_view_s, sizeof(unsigned short)*nr*nc * 3);

					delete[] final_view;
					delete[] final_view_s;

					//}

					delete[] theta;
					delete[] Regr0;
					delete[] theta0;

					if (1){

						/* residual here */

						FILE *residual_file;

						char ppm_residual_path[256];

						char jp2_residual_path_jp2[256];

						sprintf(ppm_residual_path, "%s%c%03d_%03d%s", path_out_dir, '/', c, r, "_residual.ppm");

						sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", path_out_dir, '/', c, r, "_residual.jp2");

						/*establish residual*/
						unsigned short *residual_image = new unsigned short[nr*nc * 3]();
						for (int ii = 0; ii < nr*nc*3; ii++){
							unsigned short res_val = (unsigned short)( ( (signed int)*(original_intermediate_view + ii) ) - ( (signed int)*(pshort + ii) ) + pow(2,15) );
							if (res_val>pow(2, 16)-1)
								res_val = pow(2, 16)-1;
							if (res_val < 0)
								res_val = 0;
							*(residual_image + ii) = res_val;
						}

						residual_file = fopen(ppm_residual_path, "wb");
						aux_write16ppm_16(residual_file, nc, nr, residual_image);
						fclose(residual_file);

						delete[](residual_image);

						/* here encode residual with kakadu */

						char kdu_compress_s[256];
						sprintf(kdu_compress_s, "%s%s%s%s%s%s%f", kdu_compress_path, " -i ", ppm_residual_path, " -o ", jp2_residual_path_jp2, " -no_weights -precise -full -rate ", 0.02);

						std::cout << kdu_compress_s << "\n";

						int status = system(kdu_compress_s);

						/* decode residual with kakadu */
						char kdu_expand_s[256];
						sprintf(kdu_expand_s, "%s%s%s%s%s", kdu_expand_path, " -i ", jp2_residual_path_jp2, " -o ", ppm_residual_path);

						std::cout << kdu_expand_s << "\n";

						status = system(kdu_expand_s);

						/* apply residual */

						residual_file = fopen(ppm_residual_path, "rb");

						if (!(residual_file == NULL)){

							unsigned short* jp2_residual = new unsigned short[nr*nc * 3]();

							aux_read16ppm(residual_file, nc, nr, jp2_residual);
							fclose(residual_file);

							for (int ii = 0; ii < nr*nc * 3; ii++)
							{
								signed int val = ((signed int)pshort[ii]) + ((signed int)jp2_residual[ii]) - pow(2,15);
								if (val < 0)
									val = 0;
								if (val >(pow(2, 16) - 1))
									val = pow(2, 16) - 1;
								pshort[ii] = (unsigned short)(val);
							}

							delete[] jp2_residual;
						}
					}

					filept = fopen(path_out_ppm, "wb");

					aux_write16ppm_16(filept, nc, nr, warpedColorViews[0]);

					/* DEBUG */
					//exit(0);

				}
				else{
					for (int ij = 0; ij < 5; ij++)
					{
						if ((abs(ref_rows[ij] - r) + abs(ref_cols[ij] - c)) < 1)
						{
							printf("Decoding to (this view is a reference) %s\n ", path_out_ppm);
							filept = fopen(path_out_ppm, "wb");
							aux_write16ppm_16(filept, nc, nr, colorViews[ij]);
							break;
						}
					}
				}
				fclose(filept);
			}
		}

		delete[](original_intermediate_view);

		/* clean up ... */
		for (int ij = 0; ij < 5; ij++){
			delete[](ColDisps[ij]);// = new float[nr*nc];
			delete[](RowDisps[ij]);// = new float[nr*nc];
			delete[](colorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete[](warpedColorViews[ij]);// = new unsigned short[nr*nc * 3];
			delete[](DispTargs[ij]);// = new float[nr*nc];
		}/*
		 delete[](ccomp);
		 delete[](medccomp);*/

		if (fileptLS != NULL)
			fclose(fileptLS);
		if (fileptSPARSEGLOBAL != NULL)
			fclose(fileptSPARSEGLOBAL);
	}



	/* clean up ... */
	for (int ij = 0; ij < 5; ij++){
		delete[](qDMF[ij]);// = new float[nr*nc];
		delete[](qDM[ij]);// = new int[nr*nc];
	}



	exit(0);
}