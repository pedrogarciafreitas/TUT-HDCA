#ifndef GEN_TYPES_HH
#define GEN_TYPES_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sys/stat.h>

//#include <emmintrin.h>	// Need this for SSE compiler intrinsics

#define IO_V false
#define NBIT_GR 32
#define SYSTEM_VERBOSE 0

#define BIT_DEPTH 10
#define D_DEPTH 14
#define BIT_DEPTH_RESIDUAL 16
#define BIT_DEPTH_MERGE 14
#define BIT_DEPTH_SPARSE 20

#define MEDFILT_DEPTH false /* 3x3 median filter depth after decoding  */
#define SAVE_PARTIAL_WARPED_VIEWS false /* save partial warped views (with missing regions) to disk, good for debug */
#define YUV_422 0 /* otherwise YUV 444, has effect only if if YUV_TRANSFORM true. */

#define RESIDUAl_16BIT true

#ifdef __unix__
#define _popen popen
#define _pclose pclose
#endif

const bool verbose = false;

// number of pixels in the original lenslet image.
const int num_pixels_in_lenslet = 41483904;

int system_1(char *str) {

	std::string sys_call_str(str);

#ifdef SYSTEM_VERBOSE

#ifdef _WIN32
	sys_call_str.append("> nul");
#endif
#ifdef __unix__
	sys_call_str.append("> /dev/null");
#endif

#endif

	return system(sys_call_str.c_str());

}

float getPSNR(FILE *fileout, const char *path_out_ppm, const char *path_input_ppm, const char *difftest_call)
{

	if (fileout == NULL)
		fileout = stdout;

	/* run psnr here */

	char psnr_call[1024];
	sprintf(psnr_call, "%s%s%s%s", difftest_call, path_out_ppm, " ", path_input_ppm);

	FILE *pfile;
	pfile = _popen(psnr_call, "r");

	char psnr_buffer[1024];
	while (fgets(psnr_buffer, sizeof(psnr_buffer), pfile) != 0) {
		/*...*/
	}
	_pclose(pfile);

	char tmp_char[1024];
	float psnr_value = 0;

	sscanf(psnr_buffer, "%s\t%f", tmp_char, &psnr_value);

	fprintf(fileout, "%s\n", psnr_buffer);

	return psnr_value;
}

template <class T>
T clip(T in, const T min, const T max) {

	if (in > max) {
		return max;
	}
	if (in < min) {
		return min;
	}
	return in;
}

template <class T>
T getMedian(std::vector<T> scores)
{
	T median;
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

template <class T>
void medfilt2D(T* input, T* output, int SZ, int nr, int nc)
{
	int dsz = floor(SZ / 2);
	std::vector<T> scores;

	for (int y = 0; y < nr; y++) {
		for (int x = 0; x < nc; x++) {
			scores.clear();
			for (int dy = -dsz; dy < dsz; dy++) {
				for (int dx = -dsz; dx < dsz; dx++) {
					if ((y + dy) >= 0 && (y + dy) < nr
						&& (x + dx) >= 0 && (x + dx) < nc)
						scores.push_back(input[y + dy + (x + dx)*nr]);
				}
			}
			output[y + x*nr] = getMedian(scores);
		}
	}
}

bool aux_read16PGMPPM(const char* filename, int &width, int &height, int &ncomp, unsigned short *&img)
{
	if (IO_V)
		printf("Reading %s\n", filename);

	int  max, x, y;
	int red, green, blue;
	char dummy[100];

	unsigned short *Image16bit = NULL;

	FILE *filept;

	filept = fopen(filename, "rb");

	if (filept == NULL) {
		printf("%s does not exist\n", filename);
		return false;
	}


	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);

	if (!strncmp(dummy, "P6", 2)) {
		ncomp = 3;
	}
	else if (!strncmp(dummy, "P5", 2)) {
		ncomp = 1;
	}
	else { printf("ERROR NOT PGM OR PPM\n"); return false; }


	//std::cout << width << "\t" << height << "\n";

	img = new unsigned short[width*height * ncomp]();

	Image16bit = new unsigned short[width*height * ncomp]();

	/*--< Read 16bit ppm image from filept >--*/
	int nread = fread(Image16bit, sizeof(unsigned short), width*height * ncomp, filept);

	if (nread != width*height * ncomp)
	{
		fprintf(stderr, "READ ERROR aux_read16ppm() %s\n", filename);
		delete[](img);
		delete[](Image16bit);
		return false;
	}

	fclose(filept);

	int i = 0;

	for (x = 0; x < width; x++) {
		for (y = 0; y < height; y++) {

			red = Image16bit[(x + y*width) * ncomp];
			if (ncomp == 3) {
				green = Image16bit[(x + y*width) * ncomp + 1];
				blue = Image16bit[(x + y*width) * ncomp + 2];
			}

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
			if (ncomp == 3) {
				green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
				blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);
			}

			img[i] = red;
			if (ncomp == 3) {
				img[i + height*width] = green;
				img[i + 2 * height*width] = blue;
			}

			i++;

		}
	}

	delete[](Image16bit);

	return true;

}

bool aux_write16PGMPPM(const char* filename, const int width, const int height, const int ncomp, unsigned short *img)
{

	if (IO_V)
		printf("Writing %s\n", filename);

	unsigned char *p;
	int i, tmp, j;

	unsigned short maxi = 0;

	FILE *filept;

	filept = fopen(filename, "wb");

	if (filept == NULL) {
		printf("Cannot open %s\n", filename);
		return false;
	}

	unsigned short* img16bit = new unsigned short[height*width * ncomp]();

	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = img[j + i*height];
			if (ncomp == 3) {
				img16bit[lin_ind + 1] = img[j + i*height + width*height];
				img16bit[lin_ind + 2] = img[j + i*height + 2 * width*height];
				lin_ind = lin_ind + 3;
			}
			else {
				lin_ind = lin_ind + 1;
			}

		}
	}

	for (i = 0; i < width*height * ncomp; i++)
	{
		if (*(img16bit + i) > maxi)
			maxi = *(img16bit + i);
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++) {
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		if (ncomp == 3) {
			tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
			tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		}
	}

	if (maxi > 1023)
		maxi = 65535;
	else
		maxi = 1023;

	if (ncomp == 3) {
		fprintf(filept, "P6\n%d %d\n%d\n", width, height, maxi);
	}
	else {
		fprintf(filept, "P5\n%d %d\n%d\n", width, height, maxi);
	}

	fwrite(img16bit, sizeof(unsigned short), width*height * ncomp, filept);

	fclose(filept);

	delete[](img16bit);

	return true;
}

void decodeResidualJP2(unsigned short *ps, const char *kdu_expand_path, const char *jp2_residual_path_jp2, const char *ppm_residual_path, int ncomp, const int offset, const int maxvali)
{
	/* decode residual with kakadu */
	char kdu_expand_s[1024];
	sprintf(kdu_expand_s, "\"%s\"%s%s%s%s", kdu_expand_path, " -i ", jp2_residual_path_jp2, " -o ", ppm_residual_path);

	//std::cout << kdu_expand_s << "\n";

	int status = system_1(kdu_expand_s);

	/* apply residual */

	unsigned short* jp2_residual;

	int nc1, nr1;

	if (aux_read16PGMPPM(ppm_residual_path, nc1, nr1, ncomp, jp2_residual))
	{

		for (int iir = 0; iir < nc1*nr1 * ncomp; iir++)
		{
			signed int val = ((signed int)*(ps + iir)) + ((signed int)jp2_residual[iir]) - offset;
			if (val < 0)
				val = 0;
			if (val > maxvali)
				val = maxvali;
			*(ps + iir) = (unsigned short)(val);
		}

		delete[](jp2_residual);
	}
}

void RGB2YCbCr(unsigned short *rgb, unsigned short *ycbcr, const int nr, const int nc,const int N) 
{

	/* N-bit RGB 444 -> YCbCr 444 conversion */

	double M[] = { 0.212600000000000, -0.114572000000000,   0.500000000000000,
		0.715200000000000, -0.385428000000000, -0.454153000000000,
		0.072200000000000,   0.500000000000000, -0.045847000000000, };

	double *rgbD = new double[nr*nc * 3]();
	double *ycbcrD = new double[nr*nc * 3]();

	double nd = pow(2, (double)N - 8);

	double clipval = pow(2, N) - 1;

	for (int ii = 0; ii < nr*nc*3; ii++) {

		*(rgbD + ii) = (double) *(rgb + ii);
		*(rgbD + ii) = *(rgbD + ii) / clipval;

	}

	for (int ii = 0; ii < nr*nc; ii++) {
		for (int icomp = 0; icomp < 3; icomp++) {

			*(ycbcrD + ii + icomp*nr*nc) = *(rgbD + ii)*M[icomp + 0] + 
				*(rgbD + ii + 1*nr*nc)*M[icomp + 3] + 
				*(rgbD + ii + 2*nr*nc)*M[icomp + 6];

			//printf("rgb\t%f\tycbcr\t%f\n", *(rgbD + ii + icomp*nr*nc), *(ycbcrD + ii + icomp*nr*nc));

			if (icomp < 1) {
				*(ycbcrD + ii + icomp*nr*nc) = (219 * (*(ycbcrD + ii + icomp*nr*nc)) + 16)*nd;
			}
			else {
				*(ycbcrD + ii + icomp*nr*nc) = (224 * (*(ycbcrD + ii + icomp*nr*nc)) + 128)*nd;
			}

			*(ycbcr + ii + icomp*nr*nc) = (unsigned short)*(ycbcrD + ii + icomp*nr*nc);
		}

	}

	delete[](rgbD);
	delete[](ycbcrD);

}

void YCbCr2RGB(unsigned short *ycbcr, unsigned short *rgb, const int nr, const int nc, const int N)
{

	double M[] = {		1.000000000000000,		1.000000000000000,				1.000000000000000,
						0,						-0.187330000000000,				1.855630000000000,
						1.574800000000000,		-0.468130000000000,                   0 };

	double *rgbD = new double[nr*nc * 3]();
	double *ycbcrD = new double[nr*nc * 3]();

	double nd = pow(2, (double)N - 8);

	unsigned short clipval = pow(2, N) - 1;

	double sval1 = 16 * nd;
	double sval2 = 219 * nd;
	double sval3 = 128 * nd;
	double sval4 = 224 * nd;

	for (int ii = 0; ii < nr*nc; ii++) {

		for (int icomp = 0; icomp < 3; icomp++) {

			*(ycbcrD + ii + icomp*nr*nc) = (double) *(ycbcr + ii + icomp*nr*nc);
			
			if (icomp < 1) {
				*(ycbcrD + ii + icomp*nr*nc) = clip( (*(ycbcrD + ii + icomp*nr*nc) - sval1) / sval2, 0.0, 1.0);
			}
			else {
				*(ycbcrD + ii + icomp*nr*nc) = clip( (*(ycbcrD + ii + icomp*nr*nc) - sval3) / sval4, -0.5, 0.5);
			}

		}

		for (int icomp = 0; icomp < 3; icomp++) {

			*(rgbD + ii + icomp*nr*nc) = *(ycbcrD + ii)*M[icomp + 0] +
				*(ycbcrD + ii + 1 * nr*nc)*M[icomp + 3] +
				*(ycbcrD + ii + 2 * nr*nc)*M[icomp + 6];

				*(rgb + ii + icomp*nr*nc) = (unsigned short)clip(( *(rgbD + ii + icomp*nr*nc)*clipval ), 0.0, (double)clipval);
		}

	}

	delete[](rgbD);
	delete[](ycbcrD);

}

double PSNR(unsigned short *im0, unsigned short* im1, const int NR, const int NC, const int NCOMP, double maxval)
{

	double se = 0;

	for (int ii = 0; ii < NR*NC*NCOMP; ii++)

	{
		double dx = (double)(*(im0 + ii)) - (double)(*(im1 + ii));

		se += dx*dx;
	}

	double mse = se / NR / NC / NCOMP;

	return 10 * log10((maxval*maxval) / mse);

}

double PSNR(unsigned short *im0, unsigned short* im1, const int NR, const int NC, const int NCOMP)
{

	double se = 0;

	double maxval = 0;

	for (int ii = 0; ii < NR*NC*NCOMP; ii++)

	{
		double dx = (double)(*(im0 + ii)) - (double)(*(im1 + ii));

		maxval = *(im1 + ii) > maxval ? *(im1 + ii) : maxval;

		se += dx*dx;
	}

	double mse = se / NR / NC / NCOMP;

	return 10 * log10((maxval*maxval) / mse);

}

void RGB2YUV422(unsigned short *rgb, unsigned short **yy, unsigned short **cbb, unsigned short **crr,
	const int NR, const int NC, const int NCOMP, const int N) 
{

	unsigned short *ycbcr;

	ycbcr = new unsigned short[NR*NC*NCOMP]();

	RGB2YCbCr(rgb, ycbcr, NR, NC, N);

	unsigned short *y, *cb, *cr;

	*yy = new unsigned short[NR*NC*NCOMP]();
	*cbb = new unsigned short[NR*NC / 2 * NCOMP]();
	*crr = new unsigned short[NR*NC / 2 * NCOMP]();

	y = *yy;
	cb = *cbb;
	cr = *crr;

	memcpy(y, ycbcr, sizeof(unsigned short)*NR*NC);

	for (int cc = 0; cc < NC; cc += 2) {
		memcpy(cb + (cc / 2)*NR, ycbcr + cc*NR + NR*NC, sizeof(unsigned short)*NR);
		memcpy(cr + (cc / 2)*NR, ycbcr + cc*NR + NR*NC * 2, sizeof(unsigned short)*NR);
	}

	delete[](ycbcr);
}

double getYCbCr_422_PSNR(unsigned short *im0, unsigned short* im1, const int NR, const int NC, const int NCOMP, const int N)
{

	unsigned short 
		*im0_y, *im1_y, 
		*im0_cb, *im1_cb, 
		*im0_cr, *im1_cr;

	RGB2YUV422(im0, &im0_y, &im0_cb, &im0_cr, NR, NC, NCOMP, N);
	RGB2YUV422(im1, &im1_y, &im1_cb, &im1_cr, NR, NC, NCOMP, N);

	double nd = pow(2, N) - 1;

	double PSNR_Y = PSNR(im0_y, im1_y, NR, NC, 1, nd);
	double PSNR_Cb = PSNR(im0_cb, im1_cb, NR, NC/2, 1, nd);
	double PSNR_Cr = PSNR(im0_cb, im1_cb, NR, NC/2, 1, nd);

	delete[](im0_y);
	delete[](im1_y);
	delete[](im0_cb);
	delete[](im1_cb);
	delete[](im0_cr);
	delete[](im1_cr);

	return (6 * PSNR_Y + PSNR_Cb + PSNR_Cr) / 8;

}

void decodeResidualJP2_YUV(unsigned short *ps, const char *kdu_expand_path, char *ycbcr_jp2_names[], char *ycbcr_pgm_names[], const int ncomp, const int offset, const int maxvali)
{
	/* decode residual with kakadu */
	char kdu_expand_s[1024];

	for (int icomp = 0; icomp < ncomp; icomp++) {
		sprintf(kdu_expand_s, "\"%s\"%s%s%s%s", kdu_expand_path, " -i ", ycbcr_jp2_names[icomp], " -o ", ycbcr_pgm_names[icomp]);
		int status = system_1(kdu_expand_s);
		if (status < 0) {
			printf("KAKADU ERROR\nTERMINATING ... \n");
			exit(0);
		}
	}

	unsigned short *ycbcr;

	unsigned short *jp2_residual;

	int nc1, nr1,ncomp1;

	for (int icomp = 0; icomp < ncomp; icomp++) {
		if (aux_read16PGMPPM(ycbcr_pgm_names[icomp], nc1, nr1, ncomp1, jp2_residual))
		{
			if (icomp < 1) {
				ycbcr = new unsigned short[nc1*nr1*ncomp]();
			}

			if (YUV_422) {
				if (icomp < 1) {
					memcpy(ycbcr + icomp*nr1*nc1, jp2_residual, sizeof(unsigned short)*nr1*nc1);
				}
				else {
					for (int cc = 0; cc < nc1; cc++) {
						memcpy(ycbcr + cc * 2 * nr1 + icomp*nr1*nc1*2, jp2_residual + cc*nr1, sizeof(unsigned short)*nr1);
						memcpy(ycbcr + (cc*2+1)*nr1 + icomp*nr1*nc1*2, jp2_residual + cc*nr1, sizeof(unsigned short)*nr1);
					}
					nc1 = nc1 * 2; // since we keep using this variable...
				}
			}
			else {
				memcpy(ycbcr + icomp*nr1*nc1, jp2_residual, sizeof(unsigned short)*nr1*nc1);
			}

			

			delete[](jp2_residual);
		}
	}

	signed int dv = RESIDUAl_16BIT ? 1 : 2;
	signed int BP = RESIDUAl_16BIT ? 16 : 10;
	signed int maxval = pow(2, BP)-1;

	for (int ii = 0; ii < nr1*nc1*ncomp; ii++) {
		*(ycbcr + ii) = clip(*(ycbcr + ii), (unsigned short)0, (unsigned short)maxval);
	}

	unsigned short *rgb = new unsigned short[nr1*nc1*ncomp]();


	if (RESIDUAl_16BIT) {
		YCbCr2RGB(ycbcr, rgb, nr1, nc1, 16);
	}
	else {
		YCbCr2RGB(ycbcr, rgb, nr1, nc1, 10);
	}

	/* apply residual */

	for (int iir = 0; iir < nc1*nr1 * ncomp; iir++)
	{
		signed int val = (signed int)*(ps + iir) + (signed int)(rgb[iir]*dv) - offset; // we assume that for 10bit case we have offset as 2^10-1, so go from 2^11 range to 2^10 and lose 1 bit of precision
		val = clip(val, 0, maxvali);
		*(ps + iir) = (unsigned short)(val);
	}

	delete[](ycbcr);
	delete[](rgb);
}

void encodeResidualJP2_YUV(const int nr, const int nc, unsigned short *original_intermediate_view, unsigned short *ps, char *ycbcr_pgm_names[],
	const char *kdu_compress_path, char *ycbcr_jp2_names[], const float residual_rate, const int ncomp, const int offset, float rate_a)
{
	/*establish residual*/
	unsigned short *residual_image = new unsigned short[nr*nc * ncomp]();

	signed int dv = RESIDUAl_16BIT ? 1 : 2;
	signed int BP = RESIDUAl_16BIT ? 16 : 10;
	signed int maxval = pow(2, BP)-1;

	for (int iir = 0; iir < nr*nc*ncomp; iir++) {
		signed int res_val = ( (((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + offset) )/dv;
		res_val = clip(res_val, 0, maxval);
		*(residual_image + iir) = (unsigned short)(res_val);
	}

	unsigned short *ycbcr = new unsigned short[nr*nc*ncomp]();

	if (RESIDUAl_16BIT) {
		RGB2YCbCr(residual_image, ycbcr, nr, nc, 16);
	}
	else {
		RGB2YCbCr(residual_image, ycbcr, nr, nc, 10);
	}

	unsigned short *tmp_im = new unsigned short[nr*nc]();

	char kdu_compress_s[1024];

	for (int icomp = 0; icomp < ncomp; icomp++) {

		float rateR = residual_rate;

		if (icomp < 1) {
			rateR = rate_a*rateR;
		}
		else {
			rateR = (1-rate_a)*rateR/2;
		}

		if (YUV_422) {
			if (icomp < 1) {
				memcpy(tmp_im, ycbcr + icomp*nr*nc, sizeof(unsigned short)*nr*nc);
				aux_write16PGMPPM(ycbcr_pgm_names[icomp], nc, nr, 1, tmp_im);
			}
			else {
				for (int cc = 0; cc < nc; cc += 2) {
					memcpy(tmp_im + (cc / 2)*nr, ycbcr + cc*nr + nr*nc*icomp, sizeof(unsigned short)*nr);
				}
				aux_write16PGMPPM(ycbcr_pgm_names[icomp], nc/2, nr, 1, tmp_im);
			}
		}
		else {
			memcpy(tmp_im, ycbcr + icomp*nr*nc, sizeof(unsigned short)*nr*nc);
			aux_write16PGMPPM(ycbcr_pgm_names[icomp], nc, nr, 1, tmp_im);
		}

		sprintf(kdu_compress_s, "\"%s\"%s%s%s%s%s%f", kdu_compress_path, " -i ", ycbcr_pgm_names[icomp], " -o ", ycbcr_jp2_names[icomp], " -no_weights -no_info -precise -full -rate ", rateR);

		int status = system_1(kdu_compress_s);

		if (status < 0) {
			printf("KAKADU ERROR\nTERMINATING ... \n");
			exit(0);
		}


	}

	delete[](tmp_im);
	delete[](ycbcr);
	delete[](residual_image);


}

void encodeResidualJP2(const int nr, const int nc, unsigned short *original_intermediate_view, unsigned short *ps, const char *ppm_residual_path,
	const char *kdu_compress_path, const char *jp2_residual_path_jp2, const float residual_rate, const int ncomp, const int offset)
{
	/*establish residual*/
	unsigned short *residual_image = new unsigned short[nr*nc * ncomp]();

	for (int iir = 0; iir < nr*nc*ncomp; iir++) {
		signed int res_val = (((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + offset);
		if (res_val > pow(2, 16) - 1)
			res_val = pow(2, 16) - 1;
		if (res_val < 0)
			res_val = 0;
		*(residual_image + iir) = (unsigned short)res_val;
	}

	aux_write16PGMPPM(ppm_residual_path, nc, nr, ncomp, residual_image);

	delete[](residual_image);

	/* here encode residual with kakadu */

	char kdu_compress_s[1024]; // tolerance 0 ?
	sprintf(kdu_compress_s, "\"%s\"%s%s%s%s%s%f", kdu_compress_path, " -i ", ppm_residual_path, " -o ", jp2_residual_path_jp2, " -no_weights -no_info -precise -full -rate ", residual_rate);

	std::cout << kdu_compress_s << "\n";

	int status = system_1(kdu_compress_s);
}


inline long aux_GetFileSize(char* filename)
{
	struct stat stat_buf;
	int rc = stat(filename, &stat_buf);
	return rc == 0 ? stat_buf.st_size : -1;
}

int FastOLS_new(double **AAA, double **Ydd, int *PredRegr0, double *PredTheta0, const int Ms, const int MT, const int MPHI, const int N)
{
	int mTheta, M, iM, iM1;
	double *B, *C, sigerr, *Ag, *g;
	double C1, valm1, temp, crit, sabsval;
	int p, j_p, i, j, k, itemp;

	double *AA = *AAA;
	double *Yd = *Ydd;

	double *PHI = new double[MT*MT]();
	double *PSI = new double[MT]();

	//int startt = clock();

	/* make ATA. this is slow. */
	for (int i1 = 0; i1 < MT; i1++) {
#pragma omp parallel for shared(i1)
		for (int j1 = 0; j1 < MT; j1++) {
			for (int ii = 0; ii < N; ii++) {
				*(PHI + i1 + j1*MT) += (*(AA + ii + i1*N))*(*(AA + ii + j1*N));
			}
		}
	}

	///* try ATA with simple SSE optimization, down to around 9sec for 1080p*/
	//double* result = (double*)_aligned_malloc(2 * sizeof(double), 16);

	//__m128d x, y, z;

	//for (int i1 = 0; i1 < MT; i1++) {
	//	for (int j1 = 0; j1 < MT; j1++) {
	//		int ii = 0;
	//		while (ii + 2 < N)
	//		{

	//			x = _mm_set_pd(*(AA + ii + i1*N), *(AA + ii + i1*N + 1));
	//			y = _mm_set_pd(*(AA + ii + j1*N), *(AA + ii + j1*N + 1));
	//			z = _mm_mul_pd(x, y);

	//			_mm_store_pd(result, z);

	//			*(PHI + i1 + j1*MT) += result[0] + result[1];

	//			ii += 2;
	//		}
	//		for (int ee = ii; ee < N; ee++) {
	//			*(PHI + i1 + j1*MT) += (*(AA + ee + i1*N))*(*(AA + ee + j1*N));
	//		}
	//	}
	//}
	//_aligned_free(result);

	//std::cout << "time elapsed in ATA\t" << (int)clock() - startt << "\n";

	//for (int i1 = 0; i1 < MT; i1++){
	//	for (int j1 = 0; j1 < MT; j1++){
	//		std::cout << *(ATA + i1 + j1*MT) << "\n";
	//	}
	//}

	//std::cout << "------------------------------------------------\n";

	/* make ATYd */
#pragma omp parallel for
	for (int i1 = 0; i1 < MT; i1++) {
		for (int ii = 0; ii < N; ii++) {
			*(PSI + i1) += (*(AA + ii + i1*N))*(*(Yd + ii));
		}
	}

	delete[](*AAA);
	*AAA = NULL;

	/* YdTYd */
	double yd2 = 0;
#pragma omp parallel for
	for (int ii = 0; ii < N; ii++)
		yd2 += (*(Yd + ii))*(*(Yd + ii));

	delete[](*Ydd);
	*Ydd = NULL;


	// Usage example: Ms= 3 says the sparsity (length of final predictor) and MT =42 tells how many regressors are available
	// Finally, MPHI = 63 tells the dimensions of the matrices, for getting linear indices in PHI
	M = MT + 1;
	//B = alocadoubleVector(M*M);// ((M+2)*(M+2));
	//C = alocadoubleVector(M*M);// ((M+2)*(M+2));
	//Ag = alocadoubleVector(M*M);// ((M+2)*(M+2));
	//g = alocadoubleVector(M);

	B = new double[M*M]();
	C = new double[M*M]();
	Ag = new double[M*M]();
	g = new double[M]();

	// Inputs: PHI is MTxMT, PSI is MTx1;
	// Outputs: PredRegr and PredTheta are also MTx1 but only first Ms entries are needed
	// Internal variables: B and C are (MT+1)x(MT+1) i.e. MxM

	B[MT + MT*M] = yd2; //B[MT,MT] = yd2; // we start from B[0,0]
	for (iM = 0; iM < MT; iM++)
	{
		PredRegr0[iM] = iM;
		B[iM + MT*M] = PSI[iM];// B[iM,MT] = PSI[iM];
		B[MT + iM*M] = PSI[iM];//B[MT,iM] = PSI[iM];
		for (iM1 = 0; iM1 < MT; iM1++)
		{
			B[iM + iM1*M] = PHI[iM + iM1*MPHI];//B[iM,iM1]=PHI[iM,iM1];
		}
	}
	for (iM = 0; iM < M; iM++)
		for (iM1 = 0; iM1 < M; iM1++)
			C[iM + iM1*M] = 0;//C[iM,iM1] = 0
	for (iM = 0; iM < MT; iM++)
		C[iM + iM*M] = 1;// C[iM,iM] = 1;
	crit = B[MT + MT*M];//crit = B[MT,MT];
	if (crit < 0.0000001)
	{
		//printf("crir, yd2 [%f] [%f] ", crit, yd2);
		i = 0;
		return i;
	}


	for (p = 0; p < Ms; p++)
	{
		valm1 = 0; j_p = 0; // pick the max value in next loop
		for (j = p; j < MT; j++)
		{
			//if(B[j+j*M] > 0.00000000000000001)
			sigerr = B[j + MT*M] * B[j + MT*M] / B[j + j*M];//sigerr  = B[j,MT]*B[j,MT]/B[j,j];
															//else
															//	sigerr = 0;
			if (sigerr > valm1)
			{
				valm1 = sigerr;
				j_p = j;
			}
		} // j_p is the index of maximum
		crit = crit - valm1;
		itemp = PredRegr0[j_p]; PredRegr0[j_p] = PredRegr0[p]; PredRegr0[p] = itemp;
		for (j = p; j < M; j++)
		{
			//% interchange B(p:end,j_p) with B(p:end,p)
			temp = B[j + j_p*M]; //temp = B[j,j_p];
			B[j + j_p*M] = B[j + p*M]; //B[j,j_p] = B[j,p];
			B[j + p*M] = temp;// B[j,p] = temp;
		}
		for (j = p; j < M; j++)
		{
			//% interchange B(j_p,p:end) with B(p,p:end)
			temp = B[j_p + j*M]; // temp = B[j_p,j];
			B[j_p + j*M] = B[p + j*M]; //B[j_p,j] = B[p,j];
			B[p + j*M] = temp;//B[p,j] = temp;
		}

		//% Fast
		for (j = 0; j <= p - 1; j++)
		{
			// % interchange C(1:p-1,j_p) with C(1:p-1,p)
			temp = C[j + j_p*M]; // temp = C[j,j_p];
			C[j + j_p*M] = C[j + p*M]; //C[j,j_p] = C[j,p];
			C[j + p*M] = temp; //C[j,p] = temp;
		}
		for (j = (p + 1); j < M; j++)
		{
			//if(B[p+p*M] > 0.00000000000000000000001)
			C[p + j*M] = B[p + j*M] / B[p + p*M];//C[p,j] = B[p,j]/B[p,p];
												 //else
												 //C[p+j*M] = 0;
		}

		for (j = (p + 1); j < MT; j++)
			for (k = j; k <= MT; k++)
			{
				B[j + k*M] = B[j + k*M] - C[p + j*M] * C[p + k*M] * B[p + p*M];//B[j,k] = B[j,k]-C[p,j]*C[p,k]*B[p,p];
			}
		for (j = (p + 1); j < MT; j++)
			for (k = j; k <= MT; k++)
			{
				B[k + j*M] = B[j + k*M];//B[k,j] = B[j,k];
			}
		//for j = (p+1):M
		//    for k = j:(M+1)
		//        B(j,k) = B(j,k)-C(p,j)*C(p,k)*B(p,p);
		//        %B(j,k) = B(j,k)-C(p,j)*B(p,k);
		//        B(k,j) = B(j,k);
		//    end
		//end
		//        for( iM=0; iM<Ms; iM++)
		//        {
		//        for( iM1=0; iM1<Ms; iM1++)
		//        	printf("C[%f] ",C[iM+iM1*M]);
		//        printf(" \n" );
		//        }
		// scanf("%d",&i);
	} //% for( p=0; p<M; p++ )


	  // final triangular backsolving
	for (i = 0; i < Ms; i++)
	{
		g[i] = C[i + MT*M];//g[i] = C[i,MT];
		for (j = 0; j < Ms; j++)
			Ag[i + j*M] = C[i + j*M];//Ag[i,j] = C[i,j];
	}
	PredTheta0[Ms - 1] = g[Ms - 1];
	for (i = Ms - 2; i >= 0; i--)
	{
		PredTheta0[i] = g[i];
		for (j = i + 1; j < Ms; j++)
			PredTheta0[i] = PredTheta0[i] - Ag[i + j*M] * PredTheta0[j];//PredTheta[i] = PredTheta[i]- Ag[i,j]*PredTheta[j];
	}
	//printf("pred FASTOLS [%f][%f][%f][%f][%f][%f]\n",PredTheta0[0], PredTheta0[1], PredTheta0[2], PredTheta0[3], PredTheta0[4], PredTheta0[5]);
	// printf("pred [%d][%d][%d][%d][%d]\n",PredRegr0[0], PredRegr0[1], PredRegr0[2], PredRegr0[3], PredRegr0[4] );

	if (PredTheta0[0] != PredTheta0[0])
	{// if is nan
		//printf("PredTheta0[0]  is NaN\n");
		PredTheta0[0] = 1.0;
		for (i = 1; i < Ms; i++)
		{
			PredTheta0[i] = 0.0;
		}
	}

	sabsval = 0;
	for (i = 0; i < Ms; i++)
	{

		if (PredTheta0[i] != PredTheta0[i])
		{// if is nan
			PredTheta0[i] = 0.0;
			//printf("PredTheta0[%i]  is NaN\n", i);
		}

		if (PredTheta0[i] > 0)
			sabsval = sabsval + PredTheta0[i];
		else
			sabsval = sabsval - PredTheta0[i];
	}
	//printf("%f\n", sabsval);
	if (sabsval > 2 * Ms) // if average coefficients are too high forget about intrpolation
	{
		PredTheta0[0] = C[0 + MT*M];//g[0] = C[0,MT];

									//PredTheta0[0] = 1;

		if (PredTheta0[0] != PredTheta0[0])
		{// if is nan
			//printf("C[0+MT*M]  is NaN\n", i);
			PredTheta0[0] = 1.0;
		}

		// fix ???
		//if (abs(PredTheta0[0]) > 100)
		//PredTheta0[0] = PredTheta0[0] / abs(PredTheta0[0]);

		for (i = 1; i < Ms; i++) {
			PredTheta0[i] = 0.0;
		}

		i = 1;
	}
	else {
		i = Ms;
	}

	delete[](B);
	delete[](C);
	delete[](Ag);
	delete[](g);

	delete[](PSI);
	delete[](PHI);

	//printf("pred Ms [%d]",Ms);
	return i;

}

int FastOLS_new(double *AA, double *Yd, int *PredRegr0, double *PredTheta0, const int Ms, const int MT, const int MPHI, const int N,
	double *PHI, double *PSI)
{
	int mTheta, M, iM, iM1;
	double *B, *C, sigerr, *Ag, *g;
	double C1, valm1, temp, crit, sabsval;
	int p, j_p, i, j, k, itemp;


	//double *PHI = new double[MT*MT]();
	//double *PSI = new double[MT]();

	//int startt = clock();

//	/* make ATA. this is slow. */
//	for (int i1 = 0; i1 < MT; i1++) {
//#pragma omp parallel for shared(i1)
//		for (int j1 = 0; j1 < MT; j1++) {
//			for (int ii = 0; ii < N; ii++) {
//				*(PHI + i1 + j1*MT) += (*(AA + ii + i1*N))*(*(AA + ii + j1*N));
//			}
//		}
//	}

	///* try ATA with simple SSE optimization, down to around 9sec for 1080p*/
	//double* result = (double*)_aligned_malloc(2 * sizeof(double), 16);

	//__m128d x, y, z;

	//for (int i1 = 0; i1 < MT; i1++) {
	//	for (int j1 = 0; j1 < MT; j1++) {
	//		int ii = 0;
	//		while (ii + 2 < N)
	//		{

	//			x = _mm_set_pd(*(AA + ii + i1*N), *(AA + ii + i1*N + 1));
	//			y = _mm_set_pd(*(AA + ii + j1*N), *(AA + ii + j1*N + 1));
	//			z = _mm_mul_pd(x, y);

	//			_mm_store_pd(result, z);

	//			*(PHI + i1 + j1*MT) += result[0] + result[1];

	//			ii += 2;
	//		}
	//		for (int ee = ii; ee < N; ee++) {
	//			*(PHI + i1 + j1*MT) += (*(AA + ee + i1*N))*(*(AA + ee + j1*N));
	//		}
	//	}
	//}
	//_aligned_free(result);

	//std::cout << "time elapsed in ATA\t" << (int)clock() - startt << "\n";

	//for (int i1 = 0; i1 < MT; i1++){
	//	for (int j1 = 0; j1 < MT; j1++){
	//		std::cout << *(ATA + i1 + j1*MT) << "\n";
	//	}
	//}

	//std::cout << "------------------------------------------------\n";

//	/* make ATYd */
//#pragma omp parallel for
//	for (int i1 = 0; i1 < MT; i1++) {
//		for (int ii = 0; ii < N; ii++) {
//			*(PSI + i1) += (*(AA + ii + i1*N))*(*(Yd + ii));
//		}
//	}
//
//	delete[](AA);
//	AA = NULL;

	/* YdTYd */
	double yd2 = 0;
#pragma omp parallel for
	for (int ii = 0; ii < N; ii++)
		yd2 += (*(Yd + ii))*(*(Yd + ii));

	delete[](Yd);
	Yd = NULL;


	// Usage example: Ms= 3 says the sparsity (length of final predictor) and MT =42 tells how many regressors are available
	// Finally, MPHI = 63 tells the dimensions of the matrices, for getting linear indices in PHI
	M = MT + 1;
	//B = alocadoubleVector(M*M);// ((M+2)*(M+2));
	//C = alocadoubleVector(M*M);// ((M+2)*(M+2));
	//Ag = alocadoubleVector(M*M);// ((M+2)*(M+2));
	//g = alocadoubleVector(M);

	B = new double[M*M]();
	C = new double[M*M]();
	Ag = new double[M*M]();
	g = new double[M]();

	// Inputs: PHI is MTxMT, PSI is MTx1;
	// Outputs: PredRegr and PredTheta are also MTx1 but only first Ms entries are needed
	// Internal variables: B and C are (MT+1)x(MT+1) i.e. MxM

	B[MT + MT*M] = yd2; //B[MT,MT] = yd2; // we start from B[0,0]
	for (iM = 0; iM < MT; iM++)
	{
		PredRegr0[iM] = iM;
		B[iM + MT*M] = PSI[iM];// B[iM,MT] = PSI[iM];
		B[MT + iM*M] = PSI[iM];//B[MT,iM] = PSI[iM];
		for (iM1 = 0; iM1 < MT; iM1++)
		{
			B[iM + iM1*M] = PHI[iM + iM1*MPHI];//B[iM,iM1]=PHI[iM,iM1];
		}
	}
	for (iM = 0; iM < M; iM++)
		for (iM1 = 0; iM1 < M; iM1++)
			C[iM + iM1*M] = 0;//C[iM,iM1] = 0
	for (iM = 0; iM < MT; iM++)
		C[iM + iM*M] = 1;// C[iM,iM] = 1;
	crit = B[MT + MT*M];//crit = B[MT,MT];
	if (crit < 0.0000001)
	{
		//printf("crir, yd2 [%f] [%f] ", crit, yd2);
		i = 0;
		return i;
	}


	for (p = 0; p < Ms; p++)
	{
		valm1 = 0; j_p = 0; // pick the max value in next loop
		for (j = p; j < MT; j++)
		{
			//if(B[j+j*M] > 0.00000000000000001)
			sigerr = B[j + MT*M] * B[j + MT*M] / B[j + j*M];//sigerr  = B[j,MT]*B[j,MT]/B[j,j];
															//else
															//	sigerr = 0;
			if (sigerr > valm1)
			{
				valm1 = sigerr;
				j_p = j;
			}
		} // j_p is the index of maximum
		crit = crit - valm1;
		itemp = PredRegr0[j_p]; PredRegr0[j_p] = PredRegr0[p]; PredRegr0[p] = itemp;
		for (j = p; j < M; j++)
		{
			//% interchange B(p:end,j_p) with B(p:end,p)
			temp = B[j + j_p*M]; //temp = B[j,j_p];
			B[j + j_p*M] = B[j + p*M]; //B[j,j_p] = B[j,p];
			B[j + p*M] = temp;// B[j,p] = temp;
		}
		for (j = p; j < M; j++)
		{
			//% interchange B(j_p,p:end) with B(p,p:end)
			temp = B[j_p + j*M]; // temp = B[j_p,j];
			B[j_p + j*M] = B[p + j*M]; //B[j_p,j] = B[p,j];
			B[p + j*M] = temp;//B[p,j] = temp;
		}

		//% Fast
		for (j = 0; j <= p - 1; j++)
		{
			// % interchange C(1:p-1,j_p) with C(1:p-1,p)
			temp = C[j + j_p*M]; // temp = C[j,j_p];
			C[j + j_p*M] = C[j + p*M]; //C[j,j_p] = C[j,p];
			C[j + p*M] = temp; //C[j,p] = temp;
		}
		for (j = (p + 1); j < M; j++)
		{
			//if(B[p+p*M] > 0.00000000000000000000001)
			C[p + j*M] = B[p + j*M] / B[p + p*M];//C[p,j] = B[p,j]/B[p,p];
												 //else
												 //C[p+j*M] = 0;
		}

		for (j = (p + 1); j < MT; j++)
			for (k = j; k <= MT; k++)
			{
				B[j + k*M] = B[j + k*M] - C[p + j*M] * C[p + k*M] * B[p + p*M];//B[j,k] = B[j,k]-C[p,j]*C[p,k]*B[p,p];
			}
		for (j = (p + 1); j < MT; j++)
			for (k = j; k <= MT; k++)
			{
				B[k + j*M] = B[j + k*M];//B[k,j] = B[j,k];
			}
		//for j = (p+1):M
		//    for k = j:(M+1)
		//        B(j,k) = B(j,k)-C(p,j)*C(p,k)*B(p,p);
		//        %B(j,k) = B(j,k)-C(p,j)*B(p,k);
		//        B(k,j) = B(j,k);
		//    end
		//end
		//        for( iM=0; iM<Ms; iM++)
		//        {
		//        for( iM1=0; iM1<Ms; iM1++)
		//        	printf("C[%f] ",C[iM+iM1*M]);
		//        printf(" \n" );
		//        }
		// scanf("%d",&i);
	} //% for( p=0; p<M; p++ )


	  // final triangular backsolving
	for (i = 0; i < Ms; i++)
	{
		g[i] = C[i + MT*M];//g[i] = C[i,MT];
		for (j = 0; j < Ms; j++)
			Ag[i + j*M] = C[i + j*M];//Ag[i,j] = C[i,j];
	}
	PredTheta0[Ms - 1] = g[Ms - 1];
	for (i = Ms - 2; i >= 0; i--)
	{
		PredTheta0[i] = g[i];
		for (j = i + 1; j < Ms; j++)
			PredTheta0[i] = PredTheta0[i] - Ag[i + j*M] * PredTheta0[j];//PredTheta[i] = PredTheta[i]- Ag[i,j]*PredTheta[j];
	}
	//printf("pred FASTOLS [%f][%f][%f][%f][%f][%f]\n",PredTheta0[0], PredTheta0[1], PredTheta0[2], PredTheta0[3], PredTheta0[4], PredTheta0[5]);
	// printf("pred [%d][%d][%d][%d][%d]\n",PredRegr0[0], PredRegr0[1], PredRegr0[2], PredRegr0[3], PredRegr0[4] );

	if (PredTheta0[0] != PredTheta0[0])
	{// if is nan
		//printf("PredTheta0[0]  is NaN\n");
		PredTheta0[0] = 1.0;
		for (i = 1; i < Ms; i++)
		{
			PredTheta0[i] = 0.0;
		}
	}

	sabsval = 0;
	for (i = 0; i < Ms; i++)
	{

		if (PredTheta0[i] != PredTheta0[i])
		{// if is nan
			PredTheta0[i] = 0.0;
			//printf("PredTheta0[%i]  is NaN\n", i);
		}

		if (PredTheta0[i] > 0)
			sabsval = sabsval + PredTheta0[i];
		else
			sabsval = sabsval - PredTheta0[i];
	}
	//printf("%f\n", sabsval);
	if (sabsval > 2 * Ms) // if average coefficients are too high forget about intrpolation
	{
		PredTheta0[0] = C[0 + MT*M];//g[0] = C[0,MT];

									//PredTheta0[0] = 1;

		if (PredTheta0[0] != PredTheta0[0])
		{// if is nan
			//printf("C[0+MT*M]  is NaN\n", i);
			PredTheta0[0] = 1.0;
		}

		// fix ???
		//if (abs(PredTheta0[0]) > 100)
		//PredTheta0[0] = PredTheta0[0] / abs(PredTheta0[0]);

		for (i = 1; i < Ms; i++) {
			PredTheta0[i] = 0.0;
		}

		i = 1;
	}
	else {
		i = Ms;
	}

	delete[](B);
	delete[](C);
	delete[](Ag);
	delete[](g);

	delete[](PSI);
	delete[](PHI);

	//printf("pred Ms [%d]",Ms);
	return i;

}

#endif
