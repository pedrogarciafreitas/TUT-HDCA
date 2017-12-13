/* Reads Ostendo disparity (.dsp) and RGB image (.ppm), and creates quantized disparity of 512 levels, based on aggregation of disparity+color into consistent regions. */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "gen_types.hh"
#include "msImageProcessor.h"

#ifdef _MSC_VER
typedef __int8 int8_t;
typedef __int32 int32_t;
typedef unsigned __int8 uint8_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

const float max_dsp = 13.0, min_dsp = 4.0;

float inline dsp_to_float(const float dspvalue){
	return dspvalue*(max_dsp - min_dsp) / 65535 - max_dsp;
}

float getMedian(std::vector<float> scores)
{
	float median;
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

int main(const int argc, const char** argv) {

	// exe rgb disparity sigmaS sigmaR minRegion labels quantDisparity

	// Read in the input variables
	int sigmaS = atoi(argv[3]);
	float sigmaR = atof(argv[4]);
	int minRegion = atoi(argv[5]);

	msImageProcessor im_proc;

	int nr, nc, ncomponents, nvc;

	aux_read_header_file(&nr, &nc, &ncomponents, &nvc, argv[1]);

	printf("nr %d \t nc %d \t ncomponents %d \n", nr, nc, ncomponents);

	/* READ RGB UINT8, here add .ppm */

	unsigned char *rgb_data_uint8 = new unsigned char[nr*nc*ncomponents];

	aux_read_file_uint8(nr, nc, ncomponents, argv[1], rgb_data_uint8);

	/* READ DISPARITY UINT16 */

	unsigned short *disparity_data_uint16 = new unsigned short[nr*nc];

	aux_readDSP_uint16(nr, nc, 1, argv[2], disparity_data_uint16);

	/* DISPARITY INTO FLOATING POINT */

	float *ref_dsp = new float[nr*nc];
	float *BB = ref_dsp;
	for (int row = 0; row < nr; row++){
		for (int col = 0; col < nc; col++){
			*BB++ = -dsp_to_float((float)disparity_data_uint16[col + row*nc]);
		}
	}

	//printf("debug at 5678 \t %d\n", rgb_data_uint8[5678 - 1]);
	//printf("debug at 5678 \t %f\n", ref_dsp[5678 - 1]);

	//return 0;

	/* REPLACE R COMPONENT WITH DISPARITY */

	unsigned char *rgb_image = new unsigned char[nr*nc*ncomponents];

	unsigned char *B = rgb_image;

	for (int row = 0; row < nr; row++){
		for (int col = 0; col < nr*nc; col += nr){
			//*B++ = rgb_data_uint8[row + col];
			*B++ = (unsigned char)( ref_dsp[row*nc + col / nr]*25 + 0.5); //  +0.5 because of truncation by casting
			*B++ = rgb_data_uint8[row + col + nr*nc];
			*B++ = rgb_data_uint8[row + col + 2 * nr*nc];
		}
	}


	///* debug dsp uint8 */
	//unsigned char *debug_dsp = new unsigned char[nr*nc];
	//for (int row = 0; row < nr; row++){
	//	for (int col = 0; col < nr*nc; col += nr){
	//		debug_dsp[row + col] = (unsigned char)(ref_dsp[row*nc + col / nr] * 25 + 0.5); //  +0.5 because of truncation by casting
	//	}
	//}
	//aux_write_header_file(nr, nc, 1, 1, "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/HDCA-warping/debug_dsp.uint8");

	//FILE* f_file = fopen("C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/HDCA-warping/debug_dsp.uint8", "a+b");

	//fwrite(debug_dsp, sizeof(unsigned char), nr*nc, f_file);

	//fclose(f_file);

	//float *debug_dspf = new float[nr*nc];
	//for (int row = 0; row < nr; row++){
	//	for (int col = 0; col < nr*nc; col += nr){
	//		debug_dspf[row + col] = ref_dsp[row*nc + col / nr];
	//	}
	//}
	//aux_write_header_file(nr, nc, 1, 1, "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/HDCA-warping/debug_dsp.float");

	//f_file = fopen("C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/HDCA-warping/debug_dsp.float", "a+b");

	//fwrite(debug_dspf, sizeof(float), nr*nc, f_file);

	//fclose(f_file);

	//return 0;

	/* RUN MEANSHIFT */
	im_proc.DefineImage(rgb_image, COLOR, nr, nc);

	delete[] rgb_image;

	printf("meanshift ...\n");
	im_proc.Segment(sigmaS, sigmaR, minRegion, HIGH_SPEEDUP);

	int *labels = im_proc.GetLabels();

	uint32_t *labelIm = new uint32_t[nr*nc];

	std::vector< std::vector<float> > region_values;

	for (int row = 0; row < nr; row++){
		for (int col = 0; col < nr*nc; col += nr){
			labelIm[row + col] = (uint32_t)(*labels++) + 1;
			
			if (labelIm[row + col] > region_values.size()){

				std::vector<float> tmp;

				region_values.push_back(tmp);
			}
			else{
				region_values.at(labelIm[row + col] - 1).push_back(ref_dsp[row*nc + col / nr]);
			}
			
		}
	}

	/* WRITE LABELS TO DISK */
	aux_write_header_file(nr, nc, 1, 1, argv[6]);

	FILE* f_file = fopen(argv[6], "a+b");

	fwrite(labelIm, sizeof(uint32_t), nr*nc, f_file);

	fclose(f_file);

	//return 0;

	/* MEDIAN FILTER DISPARITY USING REGIONS */

	printf("median filter, nreg = %d ...\n", region_values.size());

	std::vector< float > dsp_medians(region_values.size());

	for (int ij = 0; ij < region_values.size(); ij++){
		dsp_medians.at(ij) = getMedian(region_values.at(ij));
	}

	printf("median assign ...\n");

	float *DMS1t = new float[nr*nc];
	for (int ij = 0; ij < nr*nc; ij++){
		DMS1t[ij] = dsp_medians.at(labelIm[ij]-1);
	}

	/* QUANTIZE DISPARITY */

	printf("quantize disparity ...\n");

	float minDM = *std::min_element(dsp_medians.begin(),dsp_medians.end());
	float maxDM = *std::max_element(dsp_medians.begin(), dsp_medians.end());

	float Delta1 = (maxDM - minDM) / 510;

	printf("%f\t%f\t%f\n", minDM, maxDM, Delta1);

	float *vals0 = new float[511];
	vals0[0] = minDM;
	//for (int ij = 1; ij < 511; ij++)
	int ijj = 1;
	while (vals0[ijj-1]<maxDM)
	{
		vals0[ijj] = vals0[ijj - 1] + Delta1;
		ijj++;
	}

	float *D = new float[nr*nc];

	int *Labels = new int[nr*nc];

	float *quantDM = new float[nr*nc];

	for (int ij = 0; ij < nr*nc; ij++){
		Labels[ij] = 1;
		quantDM[ij] = minDM;
		D[ij] = 10 ^ 30;
	}

	for (int ij = 0; ij < nr*nc; ij++){
		//int min_ind = 0;
		//int min_eucl = 10000000000;
		for (int ik = 0; ik < 511; ik++){
			float distance = (DMS1t[ij] - vals0[ik])*(DMS1t[ij] - vals0[ik]);
			//for (int ijj = 0; ijj < nr*nc; ij++){
				if (D[ij]>distance){
					D[ij] = distance;
					Labels[ij] = ik+1;
					quantDM[ij] = vals0[ik];
				}
			//}
			//if (distance < min_eucl){
			//	min_eucl = distance;
			//	min_ind = ik;
			//}
		}
		//quantDM[ij] = vals0[min_ind];
		//Labels[ij] = min_ind + 1;
	}

	/* MAKE COLUMN WISE (transpose) AND WRITE TO DISK */

	//printf("transpose for matlab ...\n");

	//float *quantDM_todisk = new float[nr*nc];
	//int *Labels_todisk = new int[nr*nc];

	////for (int row = 0; row < nr; row++){
	////	for (int col = 0; col < nc; col++){
	////		quantDM_todisk[row + col*nr] = quantDM[col + row*nc];
	////		Labels_todisk[row + col*nr] = Labels[col + row*nc];
	////	}
	////}

	//int *AA = Labels;
	//float *AAA = quantDM;

	//for (int row = 0; row < nr; row++){
	//	for (int col = 0; col < nr*nc; col += nr){
	//		quantDM_todisk[row + col] = (*AAA++);
	//		Labels_todisk[row + col] = (*AA++);
	//	}
	//}

	/* WRITE QUANTIZED DISPARITY AND LABELS TO DISK */
	aux_write_header_file(nr, nc, 1, 1, argv[7]);

	f_file = fopen(argv[7], "a+b");

	fwrite(Labels, sizeof(int), nr*nc, f_file);

	fclose(f_file);

	aux_write_header_file(nr, nc, 1, 1, argv[8]);

	f_file = fopen(argv[8], "a+b");

	fwrite(quantDM, sizeof(float), nr*nc, f_file);

	fclose(f_file);

	return 0;

}