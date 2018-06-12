#ifndef CERV_H
#define CERV_H

#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

#include "openjpeg.h"
#include "mqc.h"
#include "opj_malloc.h"

	//void cerv_encode(int** img, int NR, int NC, char* filename);
	//void cerv_decode(int** img, int NR, int NC, char* filename);
#ifdef __cplusplus
}
#endif

#endif

static const opj_mqc_state_t mqc_states[47 * 2] = {
	{ 0x5601, 0, &mqc_states[2], &mqc_states[3] },
	{ 0x5601, 1, &mqc_states[3], &mqc_states[2] },
	{ 0x3401, 0, &mqc_states[4], &mqc_states[12] },
	{ 0x3401, 1, &mqc_states[5], &mqc_states[13] },
	{ 0x1801, 0, &mqc_states[6], &mqc_states[18] },
	{ 0x1801, 1, &mqc_states[7], &mqc_states[19] },
	{ 0x0ac1, 0, &mqc_states[8], &mqc_states[24] },
	{ 0x0ac1, 1, &mqc_states[9], &mqc_states[25] },
	{ 0x0521, 0, &mqc_states[10], &mqc_states[58] },
	{ 0x0521, 1, &mqc_states[11], &mqc_states[59] },
	{ 0x0221, 0, &mqc_states[76], &mqc_states[66] },
	{ 0x0221, 1, &mqc_states[77], &mqc_states[67] },
	{ 0x5601, 0, &mqc_states[14], &mqc_states[13] },
	{ 0x5601, 1, &mqc_states[15], &mqc_states[12] },
	{ 0x5401, 0, &mqc_states[16], &mqc_states[28] },
	{ 0x5401, 1, &mqc_states[17], &mqc_states[29] },
	{ 0x4801, 0, &mqc_states[18], &mqc_states[28] },
	{ 0x4801, 1, &mqc_states[19], &mqc_states[29] },
	{ 0x3801, 0, &mqc_states[20], &mqc_states[28] },
	{ 0x3801, 1, &mqc_states[21], &mqc_states[29] },
	{ 0x3001, 0, &mqc_states[22], &mqc_states[34] },
	{ 0x3001, 1, &mqc_states[23], &mqc_states[35] },
	{ 0x2401, 0, &mqc_states[24], &mqc_states[36] },
	{ 0x2401, 1, &mqc_states[25], &mqc_states[37] },
	{ 0x1c01, 0, &mqc_states[26], &mqc_states[40] },
	{ 0x1c01, 1, &mqc_states[27], &mqc_states[41] },
	{ 0x1601, 0, &mqc_states[58], &mqc_states[42] },
	{ 0x1601, 1, &mqc_states[59], &mqc_states[43] },
	{ 0x5601, 0, &mqc_states[30], &mqc_states[29] },
	{ 0x5601, 1, &mqc_states[31], &mqc_states[28] },
	{ 0x5401, 0, &mqc_states[32], &mqc_states[28] },
	{ 0x5401, 1, &mqc_states[33], &mqc_states[29] },
	{ 0x5101, 0, &mqc_states[34], &mqc_states[30] },
	{ 0x5101, 1, &mqc_states[35], &mqc_states[31] },
	{ 0x4801, 0, &mqc_states[36], &mqc_states[32] },
	{ 0x4801, 1, &mqc_states[37], &mqc_states[33] },
	{ 0x3801, 0, &mqc_states[38], &mqc_states[34] },
	{ 0x3801, 1, &mqc_states[39], &mqc_states[35] },
	{ 0x3401, 0, &mqc_states[40], &mqc_states[36] },
	{ 0x3401, 1, &mqc_states[41], &mqc_states[37] },
	{ 0x3001, 0, &mqc_states[42], &mqc_states[38] },
	{ 0x3001, 1, &mqc_states[43], &mqc_states[39] },
	{ 0x2801, 0, &mqc_states[44], &mqc_states[38] },
	{ 0x2801, 1, &mqc_states[45], &mqc_states[39] },
	{ 0x2401, 0, &mqc_states[46], &mqc_states[40] },
	{ 0x2401, 1, &mqc_states[47], &mqc_states[41] },
	{ 0x2201, 0, &mqc_states[48], &mqc_states[42] },
	{ 0x2201, 1, &mqc_states[49], &mqc_states[43] },
	{ 0x1c01, 0, &mqc_states[50], &mqc_states[44] },
	{ 0x1c01, 1, &mqc_states[51], &mqc_states[45] },
	{ 0x1801, 0, &mqc_states[52], &mqc_states[46] },
	{ 0x1801, 1, &mqc_states[53], &mqc_states[47] },
	{ 0x1601, 0, &mqc_states[54], &mqc_states[48] },
	{ 0x1601, 1, &mqc_states[55], &mqc_states[49] },
	{ 0x1401, 0, &mqc_states[56], &mqc_states[50] },
	{ 0x1401, 1, &mqc_states[57], &mqc_states[51] },
	{ 0x1201, 0, &mqc_states[58], &mqc_states[52] },
	{ 0x1201, 1, &mqc_states[59], &mqc_states[53] },
	{ 0x1101, 0, &mqc_states[60], &mqc_states[54] },
	{ 0x1101, 1, &mqc_states[61], &mqc_states[55] },
	{ 0x0ac1, 0, &mqc_states[62], &mqc_states[56] },
	{ 0x0ac1, 1, &mqc_states[63], &mqc_states[57] },
	{ 0x09c1, 0, &mqc_states[64], &mqc_states[58] },
	{ 0x09c1, 1, &mqc_states[65], &mqc_states[59] },
	{ 0x08a1, 0, &mqc_states[66], &mqc_states[60] },
	{ 0x08a1, 1, &mqc_states[67], &mqc_states[61] },
	{ 0x0521, 0, &mqc_states[68], &mqc_states[62] },
	{ 0x0521, 1, &mqc_states[69], &mqc_states[63] },
	{ 0x0441, 0, &mqc_states[70], &mqc_states[64] },
	{ 0x0441, 1, &mqc_states[71], &mqc_states[65] },
	{ 0x02a1, 0, &mqc_states[72], &mqc_states[66] },
	{ 0x02a1, 1, &mqc_states[73], &mqc_states[67] },
	{ 0x0221, 0, &mqc_states[74], &mqc_states[68] },
	{ 0x0221, 1, &mqc_states[75], &mqc_states[69] },
	{ 0x0141, 0, &mqc_states[76], &mqc_states[70] },
	{ 0x0141, 1, &mqc_states[77], &mqc_states[71] },
	{ 0x0111, 0, &mqc_states[78], &mqc_states[72] },
	{ 0x0111, 1, &mqc_states[79], &mqc_states[73] },
	{ 0x0085, 0, &mqc_states[80], &mqc_states[74] },
	{ 0x0085, 1, &mqc_states[81], &mqc_states[75] },
	{ 0x0049, 0, &mqc_states[82], &mqc_states[76] },
	{ 0x0049, 1, &mqc_states[83], &mqc_states[77] },
	{ 0x0025, 0, &mqc_states[84], &mqc_states[78] },
	{ 0x0025, 1, &mqc_states[85], &mqc_states[79] },
	{ 0x0015, 0, &mqc_states[86], &mqc_states[80] },
	{ 0x0015, 1, &mqc_states[87], &mqc_states[81] },
	{ 0x0009, 0, &mqc_states[88], &mqc_states[82] },
	{ 0x0009, 1, &mqc_states[89], &mqc_states[83] },
	{ 0x0005, 0, &mqc_states[90], &mqc_states[84] },
	{ 0x0005, 1, &mqc_states[91], &mqc_states[85] },
	{ 0x0001, 0, &mqc_states[90], &mqc_states[86] },
	{ 0x0001, 1, &mqc_states[91], &mqc_states[87] },
	{ 0x5601, 0, &mqc_states[92], &mqc_states[92] },
	{ 0x5601, 1, &mqc_states[93], &mqc_states[93] },
};

class CERV
{

	std::vector<unsigned char> mqc_bytes;

	int NC, NR; //width, height

	int NRgr, NCgr;

	const int NS_NCON = 1024 * 256, NS = 2, NCON = 18;
	const int LargestCount = 1080 * 1920 * 2 * 2;

	int global_memory_counter = 0;

	int *segmentation = NULL; // array of dimensions height x width

	int drc[17][4] = { { -1, -1, -1, -1 },{ -1, 1, -1, 1 },{ 0, -2, -2, 0 },{ 0, -4, 0, -2 },{ -2, 2, -3, -1 },{ -2, -2, -2, -2 },{ -2, 0, -3, 1 },{ -3, -1, -2, 2 },{ -3, 1, -1, 3 },{ -1, 3, -1, -3 },{ -1, -3, -4, 0 },{ -2, -4, -3, -3 },{ -2, 4, -3, 3 },{ -3, -3, -4, 2 },{ -3, 3, -4, -2 },{ -5, -1, -1, -5 },{ -5, 1, -1, 5 } };
	int two_pow[17] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536 };

	int **img = NULL, **gr = NULL, **gSR = NULL;
	int **Cimg = NULL, **Himg = NULL;

	int newregindex;

	int *Labels = NULL;

	

	OPJ_BYTE *bp;

	void encodeContour();
	void decodeContour();

	int** alocaMatrice(const int m, const int n);
	double** DoubleAlocaMatrice(const int m, const int n);
	int*** TripleAlocaMatrice(const int m, const int n, const int ell);
	int* alocaVector(const int m);

	int MarkRegions();

	int RegTreatSegment(int istart, int iend,
		int newregindex, int* sameregs,
		int nsame, int ir, int* minlabels, int* usamereg);

	void get_gSR_matrix();
	void get_gSR_from_gr_matrix();

	void clear_allocations();

public:

	int n_bytes[2] = { 0,0 };
	int *depth_values = NULL;

	void encode(int* Labels, int NRi, int NCi);
	void decode(unsigned char* bitstream, int* n_bytes_in, int *depth_values, int NRi, int NCi);

	unsigned char* getBitstream();
	int getBitstreamLength();

	void getSegmentation(int *SEGMFINAL);

	int getNR();
	int getNC();

	~CERV()
	{

		if (segmentation != NULL) {
			delete[](segmentation);
		}

		if (depth_values != NULL) {
			free(depth_values);
		}

		clear_allocations();
	}

};

void CERV::encode(int* Labels, int NRi, int NCi)
{
	int i, j;
	int DATA[10];
	int nrBM, nrBN, biti;
	float sec;

	NR = NRi;
	NC = NCi;

	gSR = alocaMatrice(NR, NC);
	img = alocaMatrice(NR, NC);
	Cimg = alocaMatrice(NR, NC);

	int minL = 65536; int maxL = 0;
	for (int i = 0; i < NR; ++i) {
		for (int j = 0; j < NC; ++j) {
			img[i][j] = Labels[i + j*NR];
			if (minL > Labels[i + j*NR])
				minL = Labels[i + j*NR];
			if (maxL < Labels[i + j*NR])
				maxL = Labels[i + j*NR];
		}
	}
	printf("minL maxL %d %d\n", minL, maxL);

	NRgr = 2 * NR + 4;
	NCgr = 2 * NC + 4;

	gr = alocaMatrice(NRgr, NCgr);

	for (i = 0; i < NR; i++)
		for (j = 0; j < NC; j++)
			gr[3 + 2 * i][3 + 2 * j] = img[i][j];

	// get gSR image
	get_gSR_matrix();
	// Find the regions
	newregindex = MarkRegions();

	printf("newregindex =%d\n", newregindex);

	encodeContour();

	printf("done encodeContour2\n");

	depth_values = (int*)malloc(sizeof(int)*(newregindex + 1));

	for (i = 0; i < NR; i++)
		for (j = 0; j < NC; j++)
			depth_values[Cimg[i][j]] = img[i][j];

	//fwrite(&newregindex, sizeof(int), 1, outputFile);
	//fwrite(depth_values, sizeof(int), newregindex, outputFile);

	/*free(depth_values);*/


	//if (Himg != NULL)
	//{
	//	for (i = 0; i < NR; i++)
	//		free(Himg[i]);
	//	free(Himg);
	//}

	//if (Cimg != NULL)
	//{
	//	for (i = 0; i < NR; i++)
	//		free(Cimg[i]);
	//	free(Cimg);
	//}

	//if (gSR != NULL)
	//{
	//	for (i = 0; i < NR; i++)
	//		free(gSR[i]);
	//	free(gSR);
	//}

	//if (img != NULL)
	//{
	//	for (i = 0; i < NR; i++)
	//		free(img[i]);
	//	free(img);
	//}

}

void CERV::decode(unsigned char* bitstream, int* n_bytes_in, int *depth_values_in, int NRi, int NCi)
{

	bp = static_cast<OPJ_BYTE*>( bitstream);

	n_bytes[0] = n_bytes_in[0];
	n_bytes[1] = n_bytes_in[1];

	NR = NRi;
	NC = NCi;

	gSR = alocaMatrice(NR, NC);
	img = NULL;

	Himg = alocaMatrice(NR, NC);
	Cimg = alocaMatrice(NR, NC);

	decodeContour();

	newregindex = MarkRegions();
	depth_values = depth_values_in;

	img = alocaMatrice(NR, NC);

	for (int i = 0; i < NR; i++)
		for (int j = 0; j < NC; j++)
			img[i][j] = depth_values[Cimg[i][j]];

}

void CERV::getSegmentation(int *SEGMFINAL) {

	for (int i = 0; i < NR; ++i) {
		for (int j = 0; j < NC; ++j) {
			SEGMFINAL[i + j*NR] = img[i][j];
		}
	}

}

unsigned char* CERV::getBitstream()
{
	return &mqc_bytes.at(0);
}

int CERV::getNR() {
	return NR;
}

int CERV::getNC() {
	return NC;
}

int CERV::getBitstreamLength()
{
	return mqc_bytes.size();
}

int** CERV::alocaMatrice(const int m, const int n)
{
	int **matrix, i, j;
	matrix = (int**)malloc(m * sizeof(int*));
	for (i = 0; i < m; i++)
	{
		matrix[i] = (int*)malloc(n * sizeof(int));
		for (j = 0; j < n; j++)
			matrix[i][j] = 0;
	}
	//printf("alocamatrice,  matrix=%d, m=%d, n=%d\n", matrix, m, n);

	global_memory_counter = global_memory_counter + m*n * sizeof(int);

	return matrix;
}

double** CERV::DoubleAlocaMatrice(const int m, const int n)
{
	double **matrix;
	int i, j;
	matrix = (double**)malloc(m * sizeof(double*));
	for (i = 0; i < m; i++)
	{
		matrix[i] = (double*)malloc(n * sizeof(double));
		for (j = 0; j < n; j++)
			matrix[i][j] = 0;
	}

	global_memory_counter = global_memory_counter + m*n * sizeof(double);
	return matrix;
}

int*** CERV::TripleAlocaMatrice(const int m, const int n, const int ell)
{
	int ***matrix, i, j, k;
	matrix = (int***)malloc(m * sizeof(int**));
	for (i = 0; i < m; i++)
	{
		matrix[i] = (int**)malloc(n * sizeof(int*));
		for (k = 0; k < n; k++)
		{
			matrix[i][k] = (int*)malloc(ell * sizeof(int));
			for (j = 0; j < ell; j++)
				matrix[i][k][j] = 0;
		}
	}
	global_memory_counter = global_memory_counter + m*n * ell * sizeof(int);
	return matrix;
}

int* CERV::alocaVector(const int m)
{
	int *vector, i;
	vector = (int*)malloc(m * sizeof(int));
	for (i = 0; i < m; i++)
	{
		vector[i] = 0;
	}
	global_memory_counter = global_memory_counter + m * sizeof(int);
	return vector;
}

void CERV::encodeContour()
{
	int i, j, ell, k, ic, ir, ns, context, context0, symbol, s;
	int* putere;

	double delta, CLchildren, CLTot, AddCost;
	int contexti[20];  // for safety reason
	int nBITSTREAM, nsb; //  BITSTREAM[NS_NCON],
	int* BITSTREAM = (int*)calloc(NS_NCON, sizeof(int));


	double* slog2n = (double*)calloc(LargestCount, sizeof(double));
	double* slog2deln = (double*)calloc(LargestCount, sizeof(double));
	int nerr;
	char stopit, dummy;
	int Cgr, ver, ver2, ir2, ic2, icon, Cgr3, nCgr3, nnCgr3;

	int **Ncounts;
	int **SplitT, **SplitT1, **ToknowT, **ToknowT1;
	double **CL, **CLprop;
	int ***Dist;

	opj_mqc_t *mqc;
	//OPJ_BYTE *bp;
	OPJ_UINT32 iopj;
	OPJ_UINT32 dopj;
	OPJ_BYTE *bp_dec;

	FILE *filept1;

	/*these take a lot of memory*/
	Ncounts = alocaMatrice(NS_NCON, NS);
	SplitT = alocaMatrice(NS_NCON, NCON);
	SplitT1 = alocaMatrice(NS_NCON, NCON);
	ToknowT = alocaMatrice(NS_NCON, NCON);
	ToknowT1 = alocaMatrice(NS_NCON, NCON);
	CL = DoubleAlocaMatrice(NS_NCON, NCON);
	CLprop = DoubleAlocaMatrice(NS_NCON, NCON);
	Dist = TripleAlocaMatrice(NS_NCON, NS, NCON);

	nsb = 2;
	ns = 2;
	delta = 0.5;

	// construct auxiliary vectors putere, slog2n, slog2deln (for speed only)
	putere = alocaVector(NCON + 2);
	putere[0] = ns;
	for (i = 1; i <= NCON + 1; i++)
		putere[i] = putere[i - 1] * ns;
	//  set the values sum_(i=0)^(n-1)(log(i+ns*delta)) and sum_(i=0)^(n-1)(log(i+ns*delta)), for n=0,...
	slog2n[0] = 0;
	slog2deln[0] = 0;
	for (i = 1; i < LargestCount; i++)
	{
		slog2n[i] = slog2n[i - 1] + log(delta*ns + i - 1) / log(2.);
		slog2deln[i] = slog2deln[i - 1] + log(delta + i - 1) / log(2.);
	}


	for (i = 0; i < NS_NCON; i++)
	{
		for (ell = 0; ell < ns; ell++)
			Ncounts[i][ell] = 0;
	}

	// Colect the counts for encoding the whole gr image
	nnCgr3 = 0;
	for (ir2 = 3; ir2 <= (2 * NR + 1); ir2++)
	{

		if (ir2 % 2 == 1)
		{
			ver = 1;
			ver2 = 2;
		}
		else
		{
			ver = 0;
			ver2 = 0;
		}
		for (ic2 = ver + 3; ic2 <= (2 * NC + 1); ic2 = ic2 + 2)
		{
			symbol = gr[ir2][ic2];
			Cgr = ver;
			Cgr3 = 0;
			nCgr3 = 0;
			if (ir2 == 3)
			{
				context = NS_NCON + 2;
			}
			else
			{
				for (icon = 1; icon<NCON; icon++)
				{
					Cgr = 2 * Cgr;
					if ((ir2 + drc[icon - 1][ver2] >= 0) && (ic2 + drc[icon - 1][ver2 + 1] >= 0) && (ic2 + drc[icon - 1][ver2 + 1] <= 2 * NC + 1))
					{
						if ((icon <= 3) && (ver == 1))
						{
							nCgr3 = nCgr3 + 1;
							Cgr3 = Cgr3 + gr[ir2 + drc[icon - 1][ver2]][ic2 + drc[icon - 1][ver2 + 1]];
						}
						if (gr[ir2 + drc[icon - 1][ver2]][ic2 + drc[icon - 1][ver2 + 1]] == 1)
							Cgr = Cgr + 1;
					}
				}
				if (Cgr>NS_NCON - 1)
					printf("Error context:  [%d] \n", Cgr);
				context = Cgr;
				if ((ver == 1) && (nCgr3 == 3) && (Cgr3 < 2) && (ir2>3))
				{	// Deterministic cases
					//if(Cgr3==0)
					//	gr[ir2][ic2] = 0;
					//else
					//	gr[ir2][ic2] = 1;
					nnCgr3 = nnCgr3 + 1;
				}
				else
				{
					Ncounts[context][symbol] = Ncounts[context][symbol] + 1;
				}
			}
		}
	}
	printf("encoder1 NR,NC,nnCgr3 [%d][%d][%d]\n ", NR, NC, nnCgr3);
	//	Tree = SPL_OptimizeTree(NC,permi,nsy,ncon,log2n,log2deln)

	// Assign the counts to the last depth of distributions
	for (i = 0; i < putere[NCON - 1]; i++)
	{
		for (j = 0; j < ns; j++)
		{
			for (ell = 0; ell < NCON - 1; ell++)
				Dist[i][j][ell] = 0;
			Dist[i][j][NCON - 1] = Ncounts[i][j];
		}
	}


	// Propagate the best costs from deepest level to root
	// the cost for specifying the type of nodes in last depth level is zero

	for (ell = 0; ell < NCON; ell++)   // current depth level in the tree
		for (i = 0; i < putere[ell]; i++) // current context
		{
			CL[i][ell] = 0.;
			CLprop[i][ell] = 0.;
			SplitT[i][ell] = 0;
		}


	AddCost = 0;
	for (ell = NCON - 1; ell >= 0; ell--)   // current depth level in the tree
	{
		for (i = 0; i < putere[ell]; i++) // current context
		{
			if (ell < NCON - 1)  // construct the distribution at context i and depth ell
			{
				for (k = 0; k < ns; k++)  // path in tree
					for (j = 0; j < ns; j++)  // symbol for the distribution
						Dist[i][j][ell] = Dist[i][j][ell] + Dist[i*ns + k][j][ell + 1];
				// indexes at the next depth level in the tree are i*ns +k
			}
			// compute the codelength at context i and depth ell
			CL[i][ell] = 0;
			s = 0;
			for (j = 0; j < ns; j++)
			{
				if ((Dist[i][j][ell] >= 0) &   (Dist[i][j][ell] < LargestCount))
					CL[i][ell] = CL[i][ell] - slog2deln[Dist[i][j][ell]];
				else
					printf("Too large index:  [%d] \n", Dist[i][j][ell]);

				s = s + Dist[i][j][ell];
			}
			if ((s >= 0) &   (s < LargestCount))
			{
				CL[i][ell] = CL[i][ell] + slog2n[s];
			}
			else
				printf("Too large index:  [%d] \n", s);
			CLprop[i][ell] = CL[i][ell];
			// compare the codelength at the father with CL at  children
			if (ell < NCON - 1)
			{
				CLchildren = 0;
				for (k = 0; k < ns; k++)  // path in tree
					CLchildren = CLchildren + CLprop[i*ns + k][ell + 1];
				if (CLchildren < CL[i][ell])
				{
					CLprop[i][ell] = CLchildren + AddCost;
					SplitT[i][ell] = 1;
				}
				else
				{
					CLprop[i][ell] = CL[i][ell] + AddCost;
					SplitT[i][ell] = 0;
				}
			}
		}
		AddCost = 1; // the cost for specifying the type of nodes is one
	}

	CLTot = 0;
	for (j = 0; j < ns; j++)
		CLTot = CLTot + CLprop[j][0];
	printf("Code length OPT:  [%.3f] \n", CLTot);
	CLTot = 0;
	for (i = 0; i < putere[NCON - 1]; i++) // current context
		CLTot = CLTot + CL[i][NCON - 1];
	printf("Code length:  [%.3f] \n", CLTot);



	mqc = (opj_mqc_t*)opj_malloc(sizeof(opj_mqc_t));

	bp = (OPJ_BYTE*)opj_calloc(100000, sizeof(OPJ_BYTE));


	for (iopj = 0; iopj < MQC_NUMCTXS; iopj++) {
		mqc->ctxs[iopj] = mqc_states;
		//printf("%d\n",iopj);
	}



	opj_mqc_init_enc(mqc, bp);

	// Initialize ToknowT
	for (ell = 0; ell < NCON; ell++)   // current depth level in the tree
	{
		for (i = 0; i < putere[ell]; i++) // current context
		{
			ToknowT[i][ell] = 0;
		}
	}
	for (j = 0; j < ns; j++)
		ToknowT[j][0] = 1;

	// opj_mqc_setstate(mqc,0,0,46);

	nBITSTREAM = 0;
	symbol = 1;
	for (ell = 0; ell < NCON - 1; ell++)   // current depth level in the tree
	{
		for (i = 0; i < putere[ell]; i++) { // current context
			if (ToknowT[i][ell] == 1)
			{
				if (SplitT[i][ell] == 1)
				{
					BITSTREAM[nBITSTREAM] = 1;
					context = symbol - 1;
					symbol = 2;

					opj_mqc_setcurctx(mqc, context);
					opj_mqc_encode(mqc, symbol - 1);


					for (k = 0; k < ns; k++)
						ToknowT[i*ns + k][ell + 1] = 1;
				}
				else
				{
					BITSTREAM[nBITSTREAM] = 0;
					context = symbol - 1;
					symbol = 1;

					opj_mqc_setcurctx(mqc, context);
					opj_mqc_encode(mqc, symbol - 1);


				}
				nBITSTREAM = nBITSTREAM + 1;
			}
		}
	}

	opj_mqc_flush(mqc);

	dopj = opj_mqc_numbytes(mqc);

	printf("numbytes %d nbitstream %d\n BITSTREAM \n", dopj, nBITSTREAM);
	for (i = 0; i < nBITSTREAM; i++)
		printf("%d", BITSTREAM[i]);
	printf(" \n");

	//fwrite(&dopj, sizeof(OPJ_UINT32), 1, outputFile);
	//fwrite(bp, sizeof(OPJ_BYTE), dopj, outputFile);

	for (int ikj = 0; ikj < dopj; ikj++) {
		mqc_bytes.push_back(bp[ikj]);
	}

	n_bytes[0] = dopj;

	free(mqc);
	free(bp);

	// Decode from BITSTREAM the structure of the tree
	// Initialize ToknowT and SplitT1
	for (ell = 0; ell < NCON; ell++)   // current depth level in the tree
	{
		for (i = 0; i < putere[ell]; i++) // current context
		{
			SplitT1[i][ell] = 0;
			ToknowT1[i][ell] = 0;
		}
	}
	for (j = 0; j < ns; j++)
		ToknowT1[j][0] = 1;

	nBITSTREAM = 0;
	for (ell = 0; ell < NCON - 1; ell++)   // current depth level in the tree
	{
		for (i = 0; i < putere[ell]; i++) // current context
			if (ToknowT1[i][ell] == 1)
			{
				symbol = BITSTREAM[nBITSTREAM];



				nBITSTREAM = nBITSTREAM + 1;
				if (symbol == 1)
				{
					SplitT1[i][ell] = 1;
					for (k = 0; k < ns; k++)
						ToknowT1[i*ns + k][ell + 1] = 1;
				}
				else
				{
					SplitT1[i][ell] = 0;
				}
			}
	}

	mqc = (opj_mqc_t*)opj_malloc(sizeof(opj_mqc_t));

	bp = (OPJ_BYTE*)opj_calloc(10000000, sizeof(OPJ_BYTE));


	for (iopj = 0; iopj < MQC_NUMCTXS; iopj++) {
		mqc->ctxs[iopj] = mqc_states;
		//printf("%d\n",i);
	}

	opj_mqc_init_enc(mqc, bp);

	nnCgr3 = 0;
	CLTot = 0.;
	for (ir2 = 3; ir2 <= (2 * NR + 1); ir2++)
	{
		if (ir2 % 2 == 1)
		{
			ver = 1;
			ver2 = 2;
		}
		else
		{
			ver = 0;
			ver2 = 0;
		}
		for (ic2 = ver + 3; ic2 <= (2 * NC + 1); ic2 = ic2 + 2)
		{
			symbol = gr[ir2][ic2] + 1;
			Cgr = ver;
			Cgr3 = 0;
			nCgr3 = 0;
			if (ir2 == 3)
			{
				context = NS_NCON + 2 + 2 * gr[ir2][ic2 - 4] + gr[ir2][ic2 - 2];
			}
			else
			{
				stopit = 0;
				for (icon = 1; icon<NCON; icon++)
				{
					if (SplitT1[Cgr][icon - 1] == 0)
					{
						stopit = 1;
					}
					Cgr = 2 * Cgr;
					if ((ir2 + drc[icon - 1][ver2] >= 0) && (ic2 + drc[icon - 1][ver2 + 1] >= 0) && (ic2 + drc[icon - 1][ver2 + 1] <= 2 * NC + 1))
					{
						if (icon <= 3)
						{
							nCgr3 = nCgr3 + 1;
							Cgr3 = Cgr3 + gr[ir2 + drc[icon - 1][ver2]][ic2 + drc[icon - 1][ver2 + 1]];
						}
						if (stopit == 0)
							if (gr[ir2 + drc[icon - 1][ver2]][ic2 + drc[icon - 1][ver2 + 1]] == 1)
								Cgr = Cgr + 1;
					}
				}
				if (Cgr>NS_NCON - 1)
					printf("Error context:  [%d] \n", Cgr);
				context = Cgr;
			}

			if ((ver == 1) && (nCgr3 == 3) && (Cgr3 < 2) && (ir2>3))
			{	// Deterministic cases
				nnCgr3 = nnCgr3 + 1;
			}
			else
			{
				opj_mqc_setcurctx(mqc, context);
				opj_mqc_encode(mqc, symbol - 1);
			}
		}
	}

	printf("encoder nnCgr3 [%d] \n ", nnCgr3);

	opj_mqc_flush(mqc);

	dopj = opj_mqc_numbytes(mqc);

	//fwrite(&dopj, sizeof(OPJ_UINT32), 1, outputFile);
	//fwrite(bp, sizeof(OPJ_BYTE), dopj, outputFile);

	for (int ikj = 0; ikj < dopj; ikj++) {
		mqc_bytes.push_back(bp[ikj]);
	}

	n_bytes[1] = dopj;

	free(mqc);
	free(bp);

	free(BITSTREAM);
	free(slog2deln);
	free(slog2n);
	free(putere);
	for (i = 0; i < NS_NCON; i++)
	{
		free(Ncounts[i]);
		free(SplitT[i]);
		free(SplitT1[i]);
		free(ToknowT[i]);
		free(ToknowT1[i]);
		free(CL[i]);
		free(CLprop[i]);
		for (j = 0; j < NS; j++)
			free(Dist[i][j]);
		free(Dist[i]);
	}
	free(Ncounts);
	free(SplitT);
	free(SplitT1);
	free(ToknowT);
	free(ToknowT1);
	free(CL);
	free(CLprop);
	free(Dist);

}

void CERV::decodeContour()
{
int i, ic, ir, k, ell, j, ns, context, context0, symbol;
int* putere;
// char SplitT1[NS_NCON][NCON],ToknowT[NS_NCON][NCON];
int **SplitT1, **ToknowT;
int nBITSTREAM; // BITSTREAM[NS_NCON];
int* BITSTREAM = (int*)calloc(NS_NCON, sizeof(int));
int contexti[20];  // for safety reason
char stopit;
int Cgr, ver, ver2, ir2, ic2, icon, Cgr3, nCgr3, nnCgr3;

OPJ_UINT32 n_bytes_bs = 0;

opj_mqc_t *mqc;
OPJ_BYTE bpp;
//OPJ_BYTE *bp;
OPJ_UINT32 iopj;
OPJ_UINT32 dopj;
//OPJ_BYTE *bp_dec;

FILE *filept1;

ns = 2;

SplitT1 = alocaMatrice(NS_NCON, NCON);
ToknowT = alocaMatrice(NS_NCON, NCON);

putere = alocaVector(NCON + 1);
putere[0] = ns;
for (i = 1; i <= NCON; i++)
putere[i] = putere[i - 1] * ns;


mqc = (opj_mqc_t*)opj_malloc(sizeof(opj_mqc_t));

/*opj_mqc_resetstates(mqc);*/

for (iopj = 0; iopj < MQC_NUMCTXS; iopj++) {
	mqc->ctxs[iopj] = mqc_states;

}

opj_mqc_init_dec(mqc, bp, n_bytes[0], 200);

// Decode from BITSTREAM the structure of the tree
// Initialize ToknowT and SplitT
for (ell = 0; ell < NCON; ell++)   // current depth level in the tree
{
	for (i = 0; i < putere[ell]; i++) // current context
	{
		SplitT1[i][ell] = 0;
		ToknowT[i][ell] = 0;
	}
}
nBITSTREAM = 0;
symbol = 0;
for (j = 0; j < 4; j++)
	ToknowT[j][0] = 1;
for (ell = 0; ell < (NCON - 1); ell++)   // current depth level in the tree
{
	for (i = 0; i < putere[ell]; i++) // current context
		if (ToknowT[i][ell] == 1)
		{
			context = symbol;


			opj_mqc_setcurctx(mqc, context);
			opj_mqc_decode(dopj, mqc);

			symbol = dopj;



			BITSTREAM[nBITSTREAM] = symbol;
			nBITSTREAM = nBITSTREAM + 1;
			if (symbol == 1)
			{
				SplitT1[i][ell] = 1;
				for (k = 0; k < ns; k++)
					ToknowT[i*ns + k][ell + 1] = 1;
			}
			else
			{
				SplitT1[i][ell] = 0;
			}
		}
}
printf(" nBITSTREAM [%d]  \n", nBITSTREAM);

opq_mqc_finish_dec(mqc);
free(mqc);

gr = alocaMatrice(2 * NR + 4, 2 * NC + 4);

for (ir2 = 0; ir2 <= (2 * NR + 1); ir2++)
for (ic2 = 0; ic2 <= (2 * NC + 1); ic2++)
gr[ir2][ic2] = 0;

mqc = (opj_mqc_t*)opj_malloc(sizeof(opj_mqc_t));

for (iopj = 0; iopj < MQC_NUMCTXS; iopj++) {
	mqc->ctxs[iopj] = mqc_states;
	//printf("%d\n",i);
}

opj_mqc_init_dec(mqc, bp+n_bytes[0], n_bytes[1], 200);


nnCgr3 = 0;
for (ir2 = 3; ir2 <= (2 * NR + 1); ir2++)
{
	if (ir2 % 2 == 1)
	{
		ver = 1;
		ver2 = 2;
	}
	else
	{
		ver = 0;
		ver2 = 0;
	}
	for (ic2 = ver + 3; ic2 <= (2 * NC + 1); ic2 = ic2 + 2)
	{
		Cgr = ver;
		Cgr3 = 0;
		nCgr3 = 0;
		if (ir2 == 3)
		{
			// context = NS_NCON +2;
			// context = NS_NCON +2 + gr[ir2][ic2-2];
			context = NS_NCON + 2 + 2 * gr[ir2][ic2 - 4] + gr[ir2][ic2 - 2];
		}
		else
		{
			stopit = 0;
			for (icon = 1; icon<NCON; icon++)
			{
				if (SplitT1[Cgr][icon - 1] == 0)
				{
					stopit = 1;
				}
				Cgr = 2 * Cgr;
				if ((ir2 + drc[icon - 1][ver2] >= 0) && (ic2 + drc[icon - 1][ver2 + 1] >= 0) && (ic2 + drc[icon - 1][ver2 + 1] <= 2 * NC + 1))
				{
					if (icon <= 3)
					{
						nCgr3 = nCgr3 + 1;
						Cgr3 = Cgr3 + gr[ir2 + drc[icon - 1][ver2]][ic2 + drc[icon - 1][ver2 + 1]];
					}
					if (stopit == 0)
						if (gr[ir2 + drc[icon - 1][ver2]][ic2 + drc[icon - 1][ver2 + 1]] == 1)
							Cgr = Cgr + 1;
				}
			}
			if (Cgr>NS_NCON - 1)
				printf("Error context:  [%d] \n", Cgr);
			context = Cgr;
		}

		if ((ver == 1) && (nCgr3 == 3) && (Cgr3 < 2) && (ir2>3))
		{	// Deterministic cases
			if (Cgr3 == 0)
				gr[ir2][ic2] = 0;
			else
				gr[ir2][ic2] = 1;
			nnCgr3 = nnCgr3 + 1;
		}
		else
		{
			//symbol = decode_symbol(cum_freq2[context], outputFile);
			//update_model2(symbol, context);

			opj_mqc_setcurctx(mqc, context);
			opj_mqc_decode(dopj, mqc);

			//gr[ir2][ic2] = symbol - 1;
			gr[ir2][ic2] = dopj;
			// printf(" ir2 ic2  [%d][%d] \n ", ir2,ic2);
		}
	}
}

opq_mqc_finish_dec(mqc);
free(mqc);

//free(bp);

printf(" nnCgr3 [%d] \n ", nnCgr3);

// Now decode gSR
get_gSR_from_gr_matrix();

free(putere);
free(BITSTREAM);


for (i = 0; i < NS_NCON; i++)
	free(SplitT1[i]);
free(SplitT1);

for (i = 0; i < NS_NCON; i++)
	free(ToknowT[i]);
free(ToknowT);
}

int CERV::MarkRegions()
{
	int i, ic, ir, j, ns, context, symbol;
	int* minlabels;
	int* finallabels;
	int nbrutelab, nreg, nsame, newregindex;
	int* sameregs;
	int* usamereg;
	int istart, iend;
	int change, labeli, labeli1;


	minlabels = alocaVector(NR*NC + 10);
	for (i = 0; i < NR*NC + 10; i++)
		minlabels[i] = i;
	finallabels = alocaVector(NR*NC + 10);
	sameregs = alocaVector(NC + 10);
	usamereg = alocaVector(NC + 10);   // unique values in sameregs
	newregindex = 0;
	for (ir = 0; ir < NR; ir++)
	{
		// Inspect the current row only once
		// Prepare for the first segment
		istart = 0;
		for (i = 0; i < NC + 2; i++)
			sameregs[i] = 0;
		nsame = 0;
		for (ic = 0; ic < NC; ic++)   // ic=1; ... ???????????
		{
			if (gSR[ir][ic] % 2 == 1)
			{
				// There is a vertical edge at the left of ic
				// stop here for a while, treat the segment istart:(ic-1)
				iend = ic - 1;
				newregindex = RegTreatSegment(istart, iend, newregindex, sameregs, nsame, ir, minlabels, usamereg);
				// now start a new segment
				istart = ic;
				for (i = 1; i < NC + 2; i++)
					sameregs[i] = 0;
				nsame = 0;
			}
			// if there is no hor. edge above a pixel, the pixel will have the same region index as the above pixel
			if (gSR[ir][ic] < 2)
			{
				if (ir > 0)
				{

					sameregs[nsame] = Cimg[ir - 1][ic];
					nsame = nsame + 1;
				}
			}

			if (ic == (NC - 1))  // end of row
			{
				iend = (NC - 1);
				newregindex = RegTreatSegment(istart, iend, newregindex, sameregs, nsame, ir, minlabels, usamereg);
			}
		}
	}

	// Find the minimum label representative for each label

	change = 1;
	while (change == 1)
	{
		change = 0;
		for (i = 0; i < newregindex; i++)
		{
			labeli = minlabels[i];
			labeli1 = labeli;
			while (1)
			{
				if (minlabels[labeli] < labeli)
				{
					labeli = minlabels[labeli];
				}
				else
				{
					break;
				}
			}
			if (labeli1 != labeli)
				change = 1;
			minlabels[i] = labeli;
		}
	}
	nreg = 0;
	for (i = 1; i <= newregindex; i++)
		if (minlabels[i] == i)
		{
			nreg = nreg + 1;
			finallabels[i] = nreg;
		}

	for (ir = 0; ir < NR; ir++)
		for (ic = 0; ic < NC; ic++)
		{
			labeli = minlabels[Cimg[ir][ic]];
			Cimg[ir][ic] = finallabels[labeli];
		}
	//for (j=0;j<NC;j++)
	//	for (i=0;i<NR;i++)
	//	{
	//		fwrite(&Cimg[i][j], sizeof(int), 1, C6Matlab);
	//	}
	//printf("nreg: [%d] \n",nreg);
	free(minlabels);
	free(finallabels);
	free(sameregs);
	free(usamereg);
	//return newregindex;
	return nreg;
}

int CERV::RegTreatSegment(int istart, int iend, int newregindex, int* sameregs, int nsame, int ir, int* minlabels, int* usamereg)
{
	int i, j, unique;
	int nusr, minusame;

	if (nsame == 0)
	{	// no need to update anything, just mark a new region
		newregindex = newregindex + 1;
		minlabels[newregindex] = newregindex;
		for (i = istart; i <= iend; i++)
			Cimg[ir][i] = newregindex;
		return newregindex;
	}

	// We have some neighbors: Find the distinct neighbours in sameregs
	nusr = 1;
	usamereg[0] = sameregs[0];
	for (i = 1; i < nsame; i++)
	{
		unique = 1;
		for (j = 0; j < nusr; j++)
			if (usamereg[j] == sameregs[i])
			{
				unique = 0;
				break;
			}
		if (unique == 1)
		{
			nusr = nusr + 1;
			usamereg[nusr - 1] = sameregs[i];
		}
	}

	// find the minimum of the representatives of sameregs
	minusame = minlabels[usamereg[0]];
	for (j = 0; j < nusr; j++)
		if (minlabels[usamereg[j]] < minusame)
		{
			minusame = minlabels[usamereg[j]];
		}
	// Update the minimum label for all
	for (j = 0; j < nusr; j++)
		minlabels[usamereg[j]] = minusame;

	for (i = istart; i <= iend; i++)
		Cimg[ir][i] = minusame;
	return newregindex;
}

void CERV::get_gSR_matrix() // get gSR and gr(odd,even) and gr(even,odd) from gr(odd,odd)\equiv Initial image
{
	int ic, ir;
	int i, j;

	for (ir = 0; ir < NR; ir++)
		for (ic = 0; ic < NC; ic++)
		{
			if (ic>0)
				if (gr[2 * ir + 3][2 * (ic - 1) + 3] != gr[2 * ir + 3][2 * ic + 3])
				{
					gr[2 * ir + 3][2 * ic + 2] = 1;
					gSR[ir][ic] = 1;
				}

			if (ir>0)
				if (gr[2 * (ir - 1) + 3][2 * ic + 3] != gr[2 * ir + 3][2 * ic + 3])
				{
					gr[2 * ir + 2][2 * ic + 3] = 1;
					gSR[ir][ic] += 2;
				}
		}

	for (i = 0; i < NR; i++)
		for (j = 0; j < NC; j++)
			gr[3 + 2 * i][3 + 2 * j] = 0;
}

void CERV::get_gSR_from_gr_matrix() // get gSR FROM gr(odd,even) and gr(even,odd); Reset gr(odd,odd)
{
	int ic, ir;
	int i, j;

	for (ir = 0; ir < NR; ir++)
		for (ic = 0; ic < NC; ic++)
		{
			if (ic>0)
				if (gr[2 * ir + 3][2 * ic + 2] == 1)
				{
					gSR[ir][ic] = 1;
				}

			if (ir>0)
				if (gr[2 * ir + 2][2 * ic + 3] == 1)
				{
					gSR[ir][ic] += 2;
				}
		}
	for (i = 0; i < NR; i++)
		for (j = 0; j < NC; j++)
			gr[3 + 2 * i][3 + 2 * j] = 0;
}

void CERV::clear_allocations() {

	int i;

	if (gr != NULL)
	{
		for (i = 0; i < NRgr; i++)
			free(gr[i]);
		free(gr);
	}

	if (img != NULL)
	{
		for (i = 0; i < NR; i++)
			free(img[i]);
		free(img);
	}

	if (gSR != NULL)
	{
		for (i = 0; i < NR; i++)
			free(gSR[i]);
		free(gSR);
	}

	if (Himg != NULL)
	{
		for (i = 0; i < NR; i++)
			free(Himg[i]);
		free(Himg);
	}

	if (Cimg != NULL)
	{
		for (i = 0; i < NR; i++)
			free(Cimg[i]);
		free(Cimg);
	}

}



























