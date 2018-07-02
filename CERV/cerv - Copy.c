#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h> 
//#include <conio.h>
#include <math.h>
#include <assert.h>

#include "model.h"

#include "openjpeg.h"

#include "mqc.h"
#include "opj_malloc.h"


//file
FILE *Matlab2C, *C2Matlab, *C3Matlab, *C4Matlab, *C5Matlab, *C6Matlab, *outputFile, *resultFile;
char *compImgName = NULL;

long long unsigned int global_memory_counter = 0;

//time
clock_t Begin, End, Begin1, End1;

//general 
int action = 0, NR, NC, NRgr, NCgr;
int **img, **gr, **gSR;
int **Cimg, **Himg;
int NRREG;
int *ordine;
double  CLval[15], CLind[15], CLsw[15];
int drc[17][4] = { { -1, -1, -1, -1 }, { -1, 1, -1, 1 }, { 0, -2, -2, 0 }, { 0, -4, 0, -2 }, { -2, 2, -3, -1 }, { -2, -2, -2, -2 }, { -2, 0, -3, 1 }, { -3, -1, -2, 2 }, { -3, 1, -1, 3 }, { -1, 3, -1, -3 }, { -1, -3, -4, 0 }, { -2, -4, -3, -3 }, { -2, 4, -3, 3 }, { -3, -3, -4, 2 }, { -3, 3, -4, -2 }, { -5, -1, -1, -5 }, { -5, 1, -1, 5 } };
int two_pow[17] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536 };

// int drc[NCON-1][4] = { {-1,-1,-1, -1},{-1,1,-1,1},{0,-2,-2,0},{0,-3,0,-2},{-2,2,-3,-1},{-2,-2,-2,-2},{-2,0,-3,1},{-3,-1,-2,2},{-3,1,-1,3},{-1,3,-1,-3} }; //,{-1,-3,-4,0},{-2,-4,-3,-3},{-2,4,-3,3},{-3,-3,-4,2},{-3,3,-4,-2}  ,{-5,-1,-1,-5}   ,{-5,1,-1,5} };
// int two_pow[NCON-1] = { 1,2,4,8,16,32,64,128,256,512} ; //,1024 ,2048,4096 ,8192,16384,32768  ,65536 };

// int drc[NCON-1][4] = { {-1,-1,-1, -1},{-1,1,-1,1},{0,-2,-2,0},{0,-3,0,-2},{-2,2,-3,-1},{-2,-2,-2,-2},{-2,0,-3,1},{-3,-1,-2,2},{-3,1,-1,3},{-1,3,-1,-3} };
// int two_pow[NCON-1] = { 1,2,4,8,16,32,64,128,256,512 };

static const opj_mqc_state_t mqc_states[47 * 2] = {
	{0x5601, 0, &mqc_states[2], &mqc_states[3]},
	{0x5601, 1, &mqc_states[3], &mqc_states[2]},
	{0x3401, 0, &mqc_states[4], &mqc_states[12]},
	{0x3401, 1, &mqc_states[5], &mqc_states[13]},
	{0x1801, 0, &mqc_states[6], &mqc_states[18]},
	{0x1801, 1, &mqc_states[7], &mqc_states[19]},
	{0x0ac1, 0, &mqc_states[8], &mqc_states[24]},
	{0x0ac1, 1, &mqc_states[9], &mqc_states[25]},
	{0x0521, 0, &mqc_states[10], &mqc_states[58]},
	{0x0521, 1, &mqc_states[11], &mqc_states[59]},
	{0x0221, 0, &mqc_states[76], &mqc_states[66]},
	{0x0221, 1, &mqc_states[77], &mqc_states[67]},
	{0x5601, 0, &mqc_states[14], &mqc_states[13]},
	{0x5601, 1, &mqc_states[15], &mqc_states[12]},
	{0x5401, 0, &mqc_states[16], &mqc_states[28]},
	{0x5401, 1, &mqc_states[17], &mqc_states[29]},
	{0x4801, 0, &mqc_states[18], &mqc_states[28]},
	{0x4801, 1, &mqc_states[19], &mqc_states[29]},
	{0x3801, 0, &mqc_states[20], &mqc_states[28]},
	{0x3801, 1, &mqc_states[21], &mqc_states[29]},
	{0x3001, 0, &mqc_states[22], &mqc_states[34]},
	{0x3001, 1, &mqc_states[23], &mqc_states[35]},
	{0x2401, 0, &mqc_states[24], &mqc_states[36]},
	{0x2401, 1, &mqc_states[25], &mqc_states[37]},
	{0x1c01, 0, &mqc_states[26], &mqc_states[40]},
	{0x1c01, 1, &mqc_states[27], &mqc_states[41]},
	{0x1601, 0, &mqc_states[58], &mqc_states[42]},
	{0x1601, 1, &mqc_states[59], &mqc_states[43]},
	{0x5601, 0, &mqc_states[30], &mqc_states[29]},
	{0x5601, 1, &mqc_states[31], &mqc_states[28]},
	{0x5401, 0, &mqc_states[32], &mqc_states[28]},
	{0x5401, 1, &mqc_states[33], &mqc_states[29]},
	{0x5101, 0, &mqc_states[34], &mqc_states[30]},
	{0x5101, 1, &mqc_states[35], &mqc_states[31]},
	{0x4801, 0, &mqc_states[36], &mqc_states[32]},
	{0x4801, 1, &mqc_states[37], &mqc_states[33]},
	{0x3801, 0, &mqc_states[38], &mqc_states[34]},
	{0x3801, 1, &mqc_states[39], &mqc_states[35]},
	{0x3401, 0, &mqc_states[40], &mqc_states[36]},
	{0x3401, 1, &mqc_states[41], &mqc_states[37]},
	{0x3001, 0, &mqc_states[42], &mqc_states[38]},
	{0x3001, 1, &mqc_states[43], &mqc_states[39]},
	{0x2801, 0, &mqc_states[44], &mqc_states[38]},
	{0x2801, 1, &mqc_states[45], &mqc_states[39]},
	{0x2401, 0, &mqc_states[46], &mqc_states[40]},
	{0x2401, 1, &mqc_states[47], &mqc_states[41]},
	{0x2201, 0, &mqc_states[48], &mqc_states[42]},
	{0x2201, 1, &mqc_states[49], &mqc_states[43]},
	{0x1c01, 0, &mqc_states[50], &mqc_states[44]},
	{0x1c01, 1, &mqc_states[51], &mqc_states[45]},
	{0x1801, 0, &mqc_states[52], &mqc_states[46]},
	{0x1801, 1, &mqc_states[53], &mqc_states[47]},
	{0x1601, 0, &mqc_states[54], &mqc_states[48]},
	{0x1601, 1, &mqc_states[55], &mqc_states[49]},
	{0x1401, 0, &mqc_states[56], &mqc_states[50]},
	{0x1401, 1, &mqc_states[57], &mqc_states[51]},
	{0x1201, 0, &mqc_states[58], &mqc_states[52]},
	{0x1201, 1, &mqc_states[59], &mqc_states[53]},
	{0x1101, 0, &mqc_states[60], &mqc_states[54]},
	{0x1101, 1, &mqc_states[61], &mqc_states[55]},
	{0x0ac1, 0, &mqc_states[62], &mqc_states[56]},
	{0x0ac1, 1, &mqc_states[63], &mqc_states[57]},
	{0x09c1, 0, &mqc_states[64], &mqc_states[58]},
	{0x09c1, 1, &mqc_states[65], &mqc_states[59]},
	{0x08a1, 0, &mqc_states[66], &mqc_states[60]},
	{0x08a1, 1, &mqc_states[67], &mqc_states[61]},
	{0x0521, 0, &mqc_states[68], &mqc_states[62]},
	{0x0521, 1, &mqc_states[69], &mqc_states[63]},
	{0x0441, 0, &mqc_states[70], &mqc_states[64]},
	{0x0441, 1, &mqc_states[71], &mqc_states[65]},
	{0x02a1, 0, &mqc_states[72], &mqc_states[66]},
	{0x02a1, 1, &mqc_states[73], &mqc_states[67]},
	{0x0221, 0, &mqc_states[74], &mqc_states[68]},
	{0x0221, 1, &mqc_states[75], &mqc_states[69]},
	{0x0141, 0, &mqc_states[76], &mqc_states[70]},
	{0x0141, 1, &mqc_states[77], &mqc_states[71]},
	{0x0111, 0, &mqc_states[78], &mqc_states[72]},
	{0x0111, 1, &mqc_states[79], &mqc_states[73]},
	{0x0085, 0, &mqc_states[80], &mqc_states[74]},
	{0x0085, 1, &mqc_states[81], &mqc_states[75]},
	{0x0049, 0, &mqc_states[82], &mqc_states[76]},
	{0x0049, 1, &mqc_states[83], &mqc_states[77]},
	{0x0025, 0, &mqc_states[84], &mqc_states[78]},
	{0x0025, 1, &mqc_states[85], &mqc_states[79]},
	{0x0015, 0, &mqc_states[86], &mqc_states[80]},
	{0x0015, 1, &mqc_states[87], &mqc_states[81]},
	{0x0009, 0, &mqc_states[88], &mqc_states[82]},
	{0x0009, 1, &mqc_states[89], &mqc_states[83]},
	{0x0005, 0, &mqc_states[90], &mqc_states[84]},
	{0x0005, 1, &mqc_states[91], &mqc_states[85]},
	{0x0001, 0, &mqc_states[90], &mqc_states[86]},
	{0x0001, 1, &mqc_states[91], &mqc_states[87]},
	{0x5601, 0, &mqc_states[92], &mqc_states[92]},
	{0x5601, 1, &mqc_states[93], &mqc_states[93]},
};


typedef struct region
{
	int nrPixels;
	int depth;
	int nrV;
	int  *Vecin;
	int *dVecin;
	int ok;
}REG;

REG** listaREG;

//===================	MY FUNCTIONS	===================
void usage(void)
{
	printf("\n======== PROGRAM DESCRIPTION ========\n");
	printf(" Please run the program as follows:\n\n");
	printf(" 1) DepthImgComp.exe -encode -output out.txt\n");
	printf("   To encode the data from Matlab2C.txt to out.txt file\n\n");
	printf(" 2) DepthImgComp.exe -decode -output out.txt\n");
	printf("   To decode the data from out.txt\n\n to C2Matlba.txt");
}

//=================	MALLOC FUNCTIONS	==================
void closeFiles(void)
{
	//fclose(resultFile);
	//	if (Matlab2C != NULL)
	//		fclose(Matlab2C);
	if (outputFile != NULL)
		fclose(outputFile);
	//	if (C2Matlab != NULL)
	//		fclose(C2Matlab);
	//	if (C3Matlab != NULL)
	//		close(C3Matlab);
	//	if (C4Matlab != NULL)
	//		fclose(C4Matlab);
	//	if (C5Matlab != NULL)
	//		fclose(C5Matlab);
	//	if (C6Matlab != NULL)
	//		fclose(C6Matlab);
}

int** alocaMatrice(const int m, const int n)
{
	int **matrix, i, j;
	matrix = (int**)malloc(m*sizeof(int*));
	for (i = 0; i < m; i++)
	{
		matrix[i] = (int*)malloc(n*sizeof(int));
		for (j = 0; j < n; j++)
			matrix[i][j] = 0;
	}
	//printf("alocamatrice,  matrix=%d, m=%d, n=%d\n", matrix, m, n);

	global_memory_counter = global_memory_counter + m*n * sizeof(int);

	return matrix;
}
int*** TripleAlocaMatrice(const int m, const int n, const int ell)
{
	int ***matrix, i, j, k;
	matrix = (int***)malloc(m*sizeof(int**));
	for (i = 0; i < m; i++)
	{
		matrix[i] = (int**)malloc(n*sizeof(int*));
		for (k = 0; k < n; k++)
		{
			matrix[i][k] = (int*)malloc(ell*sizeof(int));
			for (j = 0; j < ell; j++)
				matrix[i][k][j] = 0;
		}
	}
	global_memory_counter = global_memory_counter + m*n * ell* sizeof(int);
	return matrix;
}
double** DoubleAlocaMatrice(const int m, const int n)
{
	double **matrix;
	int i, j;
	matrix = (double**)malloc(m*sizeof(double*));
	for (i = 0; i < m; i++)
	{
		matrix[i] = (double*)malloc(n*sizeof(double));
		for (j = 0; j < n; j++)
			matrix[i][j] = 0;
	}
	global_memory_counter = global_memory_counter + m*n * sizeof(double);
	return matrix;
}

int* alocaVector(const int m)
{
	int *vector, i;
	vector = (int*)malloc(m*sizeof(int));
	for (i = 0; i < m; i++)
	{
		vector[i] = 0;
	}
	global_memory_counter = global_memory_counter + m* sizeof(int);
	return vector;
}


REG* newREG(int nrPixels, int depth)
{
	int k;
	REG *p;

	p = (REG*)malloc(sizeof(REG));
	assert(p != NULL);

	p->nrPixels = nrPixels;
	p->depth = depth;
	p->nrV = 0;
	k = 2 * nrPixels + 2;
	p->Vecin = alocaVector(k);
	p->dVecin = alocaVector(k);
	p->ok = 0;
	return p;
}


void freee(void)
{
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
			free(Himg[i]); //KAATUU
		free(Himg);
	}
	if (listaREG != NULL)
	{
		for (i = 1; i < NRREG; i++)
		{
			free(listaREG[i]->Vecin);
			free(listaREG[i]->dVecin);
			free(listaREG[i]);
		}
		free(listaREG);
	}

	if (ordine != NULL)
	{
		free(ordine);
	}


}

//===================	FILE FUNCTIONS	===================
void openFiles(int mode)
{
	if (mode == 1)
	{
		//	resultFile = fopen ( "result.txt" , "wb" );
		//	if (resultFile == NULL) 
		//	{
		//		fputs("File error",stderr); 
		//		exit(1);
		//	}
	}
	else if (mode == 2) // encode mode
	{
		//	C3Matlab = fopen ( "C3Matlab.txt" , "wb" );
		//	C4Matlab = fopen ( "C4Matlab.txt" , "wb" );
		//	C5Matlab = fopen ( "C5Matlab.txt" , "wb" );
		//	C6Matlab = fopen ( "C6Matlab.txt" , "wb" );
		//Matlab2C = fopen ( "Matlab2C.txt" , "rt" );
		//if (Matlab2C == NULL) 
		//{
		//	fprintf(resultFile,"1\nFile error [Matlab2C.txt]\0");
		//	fclose(resultFile);
		//	exit(1);
		//}
		outputFile = fopen(compImgName, "wb");
		if (outputFile == NULL)
		{
			//fprintf(resultFile,"1\nFile error [out.txt]");
			//fclose(resultFile);
			exit(1);
		}
	}
	else if (mode == 3)// decode mode
	{
		//		C2Matlab = fopen ( "C2Matlab.txt" , "wb" );
		//		C5Matlab = fopen ( "C5Matlab.txt" , "wb" );
		//		C6Matlab = fopen ( "C6Matlab.txt" , "wb" );
		//		if (C2Matlab == NULL) 
		//		{
		//			fprintf(resultFile,"1\nFile error [C2Matlab.txt]\0");
		//			fclose(resultFile);
		//			exit(1);
		//		}
		outputFile = fopen(compImgName, "rb");
		if (outputFile == NULL)
		{
			//fprintf(resultFile,"1\nFile error [out.txt]");
			//fclose(resultFile);
			exit(1);
		}
	}
}

void readData(void)
{
	int i, j;

	//fscanf(Matlab2C, "%d", &NR);
	//fscanf(Matlab2C, "%d", &NC);
	gSR = alocaMatrice(NR, NC);
	img = alocaMatrice(NR, NC);
	Cimg = alocaMatrice(NR, NC);
	NRgr = 2 * NR + 4;
	NCgr = 2 * NC + 4;
	gr = alocaMatrice(NRgr, NCgr);


}

void get_gSR_matrix(void) // get gSR and gr(odd,even) and gr(even,odd) from gr(odd,odd)\equiv Initial image
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
		// put to one the horiz. edges above the first line
		//		for (j=0; j<NC; j++)
		//			gr[2][2+2*j] = 1;
		for (i = 0; i < NR; i++)
			for (j = 0; j < NC; j++)
				gr[3 + 2 * i][3 + 2 * j] = 0;
}
void get_gSR_from_gr_matrix(void) // get gSR FROM gr(odd,even) and gr(even,odd); Reset gr(odd,odd)
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
void encode(int size, int nrS, int vec[], int add)
{
	int i, symbol;

	My_no_of_symbols = nrS;
	start_model();
	for (i = 0; i < size; i++)
	{
		symbol = vec[i] + add;
		encode_symbol(symbol, cum_freq, outputFile);
		update_model(symbol);
	}
}

void decode(int size, int nrS, int* vec, int add)
{
	int i, symbol;

	My_no_of_symbols = nrS;
	start_model();
	for (i = 0; i < size; i++)
	{
		// printf("decodefunc\n");
		symbol = decode_symbol(cum_freq, outputFile);
		update_model(symbol);
		vec[i] = symbol - add;
	}
}



int lengthBin(int v)
{
	int i = 0;

	if (v == 0)
		i = 1;
	while (v != 0)
	{
		i++;
		v = (v >> 1);
	}
	return i;
}



int int_cmp(const void *a, const void *b)
{
	const int *ia = (const int *)a;
	const int *ib = (const int *)b;
	return -*ia + *ib;
}


void decodeContour2(int method)
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
	OPJ_BYTE *bp;
	OPJ_UINT32 iopj;
	OPJ_UINT32 dopj;
	//OPJ_BYTE *bp_dec;

	FILE *filept1;	

	ns = 2;
	My_no_of_symbols = ns;
	My_no_of_rows = NS_NCON * 2; // 1 + 3*(1 + 4 + 4^2 + 4^3 + 4^4)
	MODEL = method;




	SplitT1 = alocaMatrice(NS_NCON, NCON);
	ToknowT = alocaMatrice(NS_NCON, NCON);

	putere = alocaVector(NCON + 1);
	putere[0] = ns;
	for (i = 1; i <= NCON; i++)
		putere[i] = putere[i - 1] * ns;

	//start_model();
	//cum_freq[0] = 1001;
	//cum_freq[1] = 501;
	//cum_freq[2] = 1;

	//start_model2();

	//filept1 = fopen("temp.mqc","rb");

	//bp = (OPJ_BYTE*)opj_calloc(100000,sizeof(OPJ_BYTE));

	////scanf(&stopit);

	//nBITSTREAM = 0;
	//while( (fread((bp+nBITSTREAM),sizeof(OPJ_BYTE),1,filept1))>0){
	//	//bpp = 0;
	//	//bp[0] = bpp;
	//	//*(bp+nBITSTREAM) = bpp;

	//	//if(nBITSTREAM<10)
	//	//printf("bp[%d] = %d\n",nBITSTREAM,bp[nBITSTREAM]);

	//	nBITSTREAM++;
	//	
	//}
	//fclose(filept1);

	fread(&n_bytes_bs,sizeof(OPJ_UINT32),1,outputFile);

	bp = (OPJ_BYTE*)opj_calloc(n_bytes_bs,sizeof(OPJ_BYTE));

	fread(bp,sizeof(OPJ_BYTE),n_bytes_bs,outputFile);

	printf("nBITSTREAM %d\n",nBITSTREAM);

	mqc = (opj_mqc_t*) opj_malloc(sizeof(opj_mqc_t));

	/*opj_mqc_resetstates(mqc);*/

	for (iopj = 0; iopj < MQC_NUMCTXS; iopj++) {
		mqc->ctxs[iopj] = mqc_states;
		//printf("%d\n",i);
	}

	opj_mqc_init_dec( mqc, bp, n_bytes_bs, 200 );

	// opj_mqc_setstate(mqc,0,0,46);

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
				//// symbol = decode_symbol(cum_freq, outputFile)-1;
				//symbol = decode_symbol(cum_freq2[context], outputFile);
				//update_model2(symbol, context);
				//symbol = symbol - 1;

				opj_mqc_setcurctx(mqc, context);
				opj_mqc_decode(dopj,mqc);

				symbol = dopj;

				//if( symbol == dopj ){
				//	printf("symbol = %d dopj = %d , nBITSTREAM=%dERROR\n",symbol,dopj,nBITSTREAM);
				//	//scanf(&stopit);
				//}

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
	printf(" nBITSTREAM [%d]  \n",nBITSTREAM);

	opq_mqc_finish_dec(mqc);
	free(mqc);
	free(bp);

	// exit(0);
	// for(i=0;i<nBITSTREAM;i++)
	// printf("%d \n",BITSTREAM[i]);



	// Now decode gr
	// start_model2();

	// initializeaza
	//for (i = 0; i < My_no_of_rows; i++)
	//{
	//	for (j = 0; j <= My_no_of_symbols; j++)
	//	{
	//		freq2[i][j] = 1;
	//		cum_freq2[i][j] = NS - j;  // cum_freq2[i][j] = My_no_of_symbols-j;
	//	}
	//	freq2[i][0] = 0;
	//}

	gr = alocaMatrice(2 * NR + 4, 2 * NC + 4);

	for (ir2 = 0; ir2 <= (2 * NR + 1); ir2++)
		for (ic2 = 0; ic2 <= (2 * NC + 1); ic2++)
			gr[ir2][ic2] = 0;

	// put to zero the horiz. edges above the first line

	//	for (j=0; j<NC; j++)
	//		gr[2][2+2*j] = 1;



	/*filept1 = fopen("temp1.mqc","rb");*/

	//bp = (OPJ_BYTE*)opj_calloc(10000000,sizeof(OPJ_BYTE));

	//nBITSTREAM = 0;
	//while( (fread((bp+nBITSTREAM),sizeof(OPJ_BYTE),1,filept1))>0){

	//	nBITSTREAM++;
	//	
	//}
	//fclose(filept1);



	fread(&n_bytes_bs,sizeof(OPJ_UINT32),1,outputFile);

	bp = (OPJ_BYTE*)opj_calloc(n_bytes_bs,sizeof(OPJ_BYTE));

	fread(bp,sizeof(OPJ_BYTE),n_bytes_bs,outputFile);

	printf("nBITSTREAM EDGES %d\n",nBITSTREAM);



	mqc = (opj_mqc_t*) opj_malloc(sizeof(opj_mqc_t));


	for (iopj = 0; iopj < MQC_NUMCTXS; iopj++) {
		mqc->ctxs[iopj] = mqc_states;
		//printf("%d\n",i);
	}

	opj_mqc_init_dec( mqc, bp, n_bytes_bs, 200 );


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
				opj_mqc_decode(dopj,mqc);

				//gr[ir2][ic2] = symbol - 1;
				gr[ir2][ic2] = dopj;
				// printf(" ir2 ic2  [%d][%d] \n ", ir2,ic2);
			}
		}
	}

	opq_mqc_finish_dec(mqc);
	free(mqc);
	free(bp);

	printf(" nnCgr3 [%d] \n ", nnCgr3);

	// Now decode gSR
	get_gSR_from_gr_matrix();
	free_model();
	free(putere);
	free(BITSTREAM);


	for (i = 0; i < NS_NCON; i++)
		free(SplitT1[i]);
	free(SplitT1);

	for (i = 0; i < NS_NCON; i++)
		free(ToknowT[i]);
	free(ToknowT);
}


void encodeContour2(int method)
{
	int i, j, ell, k, ic, ir, ns, context, context0, symbol, s;
	int* putere;
	//	int Ncounts[NS_NCON][NS];   // [index][crt_symb]
	//	int Dist[NS_NCON][NS][NCON];   // [index][crt_symb][depth]
	//	char SplitT[NS_NCON][NCON],SplitT1[NS_NCON][NCON];   // [index][depth]
	//	char ToknowT[NS_NCON][NCON],ToknowT1[NS_NCON][NCON];   // [index][depth]
	//	double CL[NS_NCON][NCON], CLprop[NS_NCON][NCON];  // [index][depth]
	//	double slog2n[100000],slog2deln[100000];
	double delta, CLchildren, CLTot, AddCost;
	int contexti[20];  // for safety reason
	int nBITSTREAM, nsb; //  BITSTREAM[NS_NCON],
	int* BITSTREAM = (int*)calloc(NS_NCON, sizeof(int));
	//
	// LargetsCount = number_rows*number_columns +10= NR*NC+10 = 8294410 for 4K images
	//

	double* slog2n = (double*)calloc(LargetsCount, sizeof(double));
	double* slog2deln = (double*)calloc(LargetsCount, sizeof(double));
	int nerr;
	char stopit, dummy;
	int Cgr, ver, ver2, ir2, ic2, icon, Cgr3, nCgr3, nnCgr3;

	int **Ncounts;
	int **SplitT, **SplitT1, **ToknowT, **ToknowT1;
	double **CL, **CLprop;
	int ***Dist;

	opj_mqc_t *mqc;
	OPJ_BYTE *bp;
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
	for (i = 1; i < LargetsCount; i++)
	{
		slog2n[i] = slog2n[i - 1] + log(delta*ns + i - 1) / log(2.);
		slog2deln[i] = slog2deln[i - 1] + log(delta + i - 1) / log(2.);
	}
	// put to zero the horiz. edges above the first line
	//	for (j=0; j<NC; j++)
	//		gr[2][2+2*j] = 1;

	for (i = 0; i < NS_NCON; i++)
	{
		for (ell = 0; ell < ns; ell++)
			Ncounts[i][ell] = 0;
	}

	My_no_of_symbols = ns;
	My_no_of_rows = NS_NCON * 2; // 1 + 3*(1 + 4 + 4^2 + 4^3 + 4^4)
	MODEL = method;

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
	printf("encoder1 NR,NC,nnCgr3 [%d][%d][%d]\n ", NR,NC,nnCgr3);
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
					if ((Dist[i][j][ell] >= 0) &   (Dist[i][j][ell] < LargetsCount))
						CL[i][ell] = CL[i][ell] - slog2deln[Dist[i][j][ell]];
					else
						printf("Too large index:  [%d] \n", Dist[i][j][ell]);

					s = s + Dist[i][j][ell];
				}
				if ((s >= 0) &   (s < LargetsCount))
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

		//for (j=0; j<NCON; j++)
		//{
		//	for (i=0; i<putere[NCON-1]; i++)
		//		fwrite(&CL[i][j], sizeof(double), 1, C4Matlab);
		//}

		CLTot = 0;
		for (j = 0; j < ns; j++)
			CLTot = CLTot + CLprop[j][0];
		printf("Code length OPT:  [%.3f] \n",CLTot);
		CLTot = 0;
		for (i = 0; i < putere[NCON - 1]; i++) // current context
			CLTot = CLTot + CL[i][NCON - 1];
		printf("Code length:  [%.3f] \n",CLTot);



		mqc = (opj_mqc_t*) opj_malloc(sizeof(opj_mqc_t));

		bp = (OPJ_BYTE*)opj_calloc(100000,sizeof(OPJ_BYTE));


		for (iopj = 0; iopj < MQC_NUMCTXS; iopj++) {
			mqc->ctxs[iopj] = mqc_states;
			//printf("%d\n",iopj);
		}

		//exit(0);

		//opj_mqc_resetstates(mqc);



		opj_mqc_init_enc( mqc, bp );



		/*for(iopj = 0; iopj < 100; iopj++){
		dopj = iopj%2;
		opj_mqc_setcurctx(mqc, 0);
		opj_mqc_encode(mqc, dopj);
		}

		opj_mqc_flush(mqc);

		for(iopj=0;iopj<50;iopj++)
		printf("ii %d = %d\n",iopj,*(bp+iopj+1));

		free(mqc);

		mqc = (opj_mqc_t*) opj_malloc(sizeof(opj_mqc_t));

		bp_dec = (OPJ_BYTE*)opj_calloc(1000,sizeof(OPJ_BYTE));

		opj_mqc_resetstates(mqc);

		opj_mqc_init_dec( mqc, bp, 50, 20 );

		for(iopj=0;iopj<50;iopj++){
		opj_mqc_setcurctx(mqc, 0);
		opj_mqc_decode(dopj,mqc);
		printf("ii %d = %d\n",iopj,dopj);
		}



		exit(0);*/

		// Encode in BITSTREAM the structure of the tree

		//start_model();
		//cum_freq[0] = 1001;
		//cum_freq[1] = 501;
		//cum_freq[2] = 1;

		//start_model2();

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
			for (i = 0; i < putere[ell]; i++){ // current context
				if (ToknowT[i][ell] == 1)
				{
					if (SplitT[i][ell] == 1)
					{
						BITSTREAM[nBITSTREAM] = 1;
						context = symbol - 1;
						symbol = 2;

						opj_mqc_setcurctx(mqc, context);
						opj_mqc_encode(mqc, symbol-1);

						// encode_symbol(symbol, cum_freq, outputFile);
						//encode_symbol(symbol, cum_freq2[context], outputFile);
						//update_model2(symbol, context);
						for (k = 0; k < ns; k++)
							ToknowT[i*ns + k][ell + 1] = 1;
					}
					else
					{
						BITSTREAM[nBITSTREAM] = 0;
						context = symbol - 1;
						symbol = 1;
						// encode_symbol(symbol, cum_freq, outputFile);

						opj_mqc_setcurctx(mqc, context);
						opj_mqc_encode(mqc, symbol-1);

						//encode_symbol(symbol, cum_freq2[context], outputFile);
						//update_model2(symbol, context);
					}
					nBITSTREAM = nBITSTREAM + 1;
				}
			}
		}

		opj_mqc_flush(mqc);

		dopj = opj_mqc_numbytes(mqc);

		printf("numbytes %d nbitstream %d\n BITSTREAM \n",dopj,nBITSTREAM);
		for (i = 0; i < nBITSTREAM; i++)
			printf("%d",BITSTREAM[i]);
		printf(" \n" );


		//for(i=0;i<dopj;i++)
		//printf("bp[%d]=%d\n",i,bp[i]);

		//filept1 = fopen("temp.mqc","wb");
		/*fwrite(bp,sizeof(OPJ_BYTE),dopj,filept1);*/
		//fclose(filept1);

		fwrite(&dopj,sizeof(OPJ_UINT32),1,outputFile);
		fwrite(bp,sizeof(OPJ_BYTE),dopj,outputFile);

		free(mqc);
		free(bp);

		//exit(0);

		/*opj_mqc_resetstates(mqc);

		opj_mqc_init_dec( mqc, bp, 50, 20 );

		for(iopj=0;iopj<50;iopj++){
		opj_mqc_setcurctx(mqc, 0);
		opj_mqc_decode(dopj,mqc);
		printf("ii %d = %d\n",iopj,dopj);
		}*/


		//printf("nBITSTREAM [%d]  \n",nBITSTREAM);

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



		//nerr= 0;
		//for(ell=0; ell<4; ell++)   // current depth level in the tree
		//	for (i=0; i<putere[ell]; i++) // current context
		//		if(SplitT1[i][ell] != SplitT[i][ell])
		//		{
		//			nerr = nerr +1;
		//			printf("context:  [%d] depth [%d] True  [%d]  Reconstr  [%d]  \n",i,ell,SplitT[i][ell],SplitT1[i][ell] );
		//		}
		//printf("Number of errors:  [%d] nBITSTREAM [%d]  \n",nerr,nBITSTREAM);


		// for(i=0;i<nBITSTREAM;i++)
		//	printf("%d \n",BITSTREAM[i]);

		mqc = (opj_mqc_t*) opj_malloc(sizeof(opj_mqc_t));

		bp = (OPJ_BYTE*)opj_calloc(10000000,sizeof(OPJ_BYTE));


		for (iopj = 0; iopj < MQC_NUMCTXS; iopj++) {
			mqc->ctxs[iopj] = mqc_states;
			//printf("%d\n",i);
		}

		opj_mqc_init_enc( mqc, bp );

		// start_model2();
		// initializeaza
		//for (i = 0; i < My_no_of_rows; i++)
		//{
		//	for (j = 0; j <= My_no_of_symbols; j++)
		//	{
		//		freq2[i][j] = 1;
		//		cum_freq2[i][j] = NS - j;  // cum_freq2[i][j] = My_no_of_symbols-j;
		//	}
		//	freq2[i][0] = 0;
		//}

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
					//if(Cgr3==0)
					//	gr[ir2][ic2] = 0;
					//else
					//	gr[ir2][ic2] = 1;
					nnCgr3 = nnCgr3 + 1;
				}
				else
				{
					//encode_symbol(symbol, cum_freq2[context], outputFile);
					//CLTot = CLTot - log((double)(cum_freq2[context][symbol - 1] - cum_freq2[context][symbol]) / cum_freq2[context][0]) / log(2);
					//update_model2(symbol, context);

					opj_mqc_setcurctx(mqc, context);
					opj_mqc_encode(mqc,symbol-1);
				}
			}
		}

		printf("encoder nnCgr3 [%d] \n ", nnCgr3);

		opj_mqc_flush(mqc);

		dopj = opj_mqc_numbytes(mqc);

		//filept1 = fopen("temp1.mqc","wb");
		//fwrite(bp,sizeof(OPJ_BYTE),dopj,filept1);
		//fclose(filept1);

		fwrite(&dopj,sizeof(OPJ_UINT32),1,outputFile);
		fwrite(bp,sizeof(OPJ_BYTE),dopj,outputFile);

		free(mqc);
		free(bp);


		//printf("nnCgr3 [%d] \n",nnCgr3);
		//printf("Code length Real: [%.3f] \n",CLTot);

		// scanf(&dummy);

		//free_model();

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

//
//double TreatSegment(int istart, int iend, int known, int val, int* excluded, int nex, int ir, double CL, int ENCODEIT, int* valregions, int* CLval, int** cum_freq30)
//{
//	double CL1;
//	int i, j, k, context, ref;
//	int* uex;
//	int* of_interest;
//	int nuex;
//	double center[40000];
//	int icenter[40000], card_center[40000], ncen;
//	int isemptyindex, possible, possible_value, index_poss, i1, i2, unique;
//	int classified, n_groups;
//	double CL10;
//	int ip;
//	int contextp;
//	int nbest, ibest, ncou;
//	int index_poss_cr, valr, symbol;
//	int CODING_ACTION;
//
//	CL1 = CL;
//	if (known > 0)
//	{	// no need to encode anything
//		for (i = istart; i <= iend; i++)
//		{
//			Himg[ir][i] = known;
//			valregions[Cimg[ir][i]] = Himg[ir][i];
//		}
//		return CL1;
//	}
//
//	
//	if (valregions[Cimg[ir][istart]] > 0)
//	{	// no need to encode anything
//		if ((valregions[Cimg[ir][istart]] != val) && (ENCODEIT == 1))
//			printf(" Errors  [%d] [%d]  \n", valregions[Cimg[ir][istart]], val);
//		for (i = istart; i <= iend; i++)
//		{
//			Himg[ir][i] = valregions[Cimg[ir][istart]];
//			valregions[Cimg[ir][i]] = Himg[ir][i];
//		}
//		return CL1;
//	}
//
//	uex = alocaVector(NC + 1);   // unique values in excluded
//	of_interest = alocaVector(2 * MAX_NICE + 1);  // displacements of interests
//	of_interest[0] = 0;
//	for (i = 1; i < MAX_NICE; i++)
//	{
//		of_interest[2 * i] = -i;
//		of_interest[2 * i + 1] = i;
//	}
//
//
//	//  We need to encode: Here the simplest case: no neighbors
//
//	//	if(ir==4)
//	//		scanf(&i);
//	if (nex == 0)  // 0 known neighbors: Encode and exit
//	{
//		if (ENCODEIT == 1)
//		{
//			//printf("case1\n");
//			context = 0;
//			CODING_ACTION = ENCODEIT;
//
//			for (i = istart; i <= iend; i++)
//			{
//				Himg[ir][i] = val;
//				valregions[Cimg[ir][i]] = Himg[ir][i];
//			}	
//			
//			// encode_symbol(val, cum_freq30[context], outputFile);   // val is larger than 0 already
//			// update_model3(val, context);
//			CL1 = CL1 - log((double)(cum_freq30[context][val - 1] - cum_freq30[context][val]) / cum_freq30[context][0]) / log(2);
//			CL10 = -log((double)(cum_freq30[context][val - 1] - cum_freq30[context][val]) / cum_freq30[context][0]) / log(2);
//			cum_freq30[context][val] = cum_freq30[context][val]+1;
//			CLval[context] = CLval[context] + CL10;
//			
//			free(uex);
//			free(of_interest);
//			return CL1;
//		}
//		else  // DECODE = 1
//		{
//			context = 0;
//			valr = decode_symbol(cum_freq30[context], outputFile);
//			update_model3(valr, context);
//			for (i = istart; i <= iend; i++)
//			{
//				Himg[ir][i] = valr;
//				valregions[Cimg[ir][i]] = Himg[ir][i];
//			}
//			free(uex);
//			free(of_interest);
//			return CL1;
//		}
//	}
//
//
//	// We have some neighbors: Find the distinct neighbours and put them in uex
//	nuex = 1;
//	uex[0] = excluded[0];
//	for (i = 1; i < nex; i++)
//	{
//		unique = 1;
//		for (j = 0; j < nuex; j++)
//			if (uex[j] == excluded[i])
//			{
//			unique = 0;
//			break;
//			}
//		if (unique == 1)
//		{
//			nuex = nuex + 1;
//			uex[nuex - 1] = excluded[i];
//		}
//	}
//
//	// the main separation to cases depending on number of distinct neighbors
//	switch (nuex)	// depending on how many distinct values
//	{
//	case 1:
//	{
//		if (ENCODEIT == 1)
//		{
//			//printf("case2\n");
//			// regions with 1 distinct neighbour
//			context = 0;
//			isemptyindex = 1;
//			for (j = 1; j < 2 * MAX_NICE + 1; j++)
//			{
//				possible_value = uex[0] + of_interest[j];
//				if (possible_value == val)
//				{
//					index_poss = j - 1;
//					isemptyindex = 0;
//				}
//			}
//		}
//		else // if(DECODE == 1)
//		{
//			context = 0;
//			// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
//			symbol = decode_symbol(cum_freq30[context + 10], outputFile);
//			update_model3(symbol, context + 10);
//			isemptyindex = symbol - 1;
//			// decode first isemptyindex
//			if (isemptyindex == 0)
//			{
//				symbol = decode_symbol(cum_freq30[context + 5], outputFile);
//				update_model3(symbol, context + 5);
//				index_poss = symbol - 1;
//				valr = uex[0] + of_interest[index_poss + 1];
//			}
//			else
//			{
//				// decode val
//				context = 0;
//				valr = decode_symbol(cum_freq30[context], outputFile);
//				update_model3(valr, context);
//			}
//		}
//	}
//		break;
//	case 2:
//	{  // regions with 2 distinct neighbours close to each other
//		if (abs(uex[0] - uex[1]) < MAX_NICE)
//		{
//			if (ENCODEIT == 1)
//			{
//				//printf("case3\n");
//				context = 1;
//				ref = (uex[0] + uex[1]) / 2;
//				isemptyindex = 1;
//				index_poss = 0;
//				for (j = 0; j < 2 * MAX_NICE + 1; j++)
//				{
//					possible_value = ref + of_interest[j];
//					if ((possible_value != uex[0]) & (possible_value != uex[1]))
//						if (possible_value == val)
//						{
//						isemptyindex = 0;
//						break;
//						}
//						else
//						{
//							index_poss = index_poss + 1;
//						}
//				}
//			}
//			else  //  if(DECODE == 1)
//			{
//				context = 1;
//				// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
//				symbol = decode_symbol(cum_freq30[context + 10], outputFile);
//				update_model3(symbol, context + 10);
//				isemptyindex = symbol - 1;
//				// decode first isemptyindex
//				if (isemptyindex == 0)
//				{
//					symbol = decode_symbol(cum_freq30[context + 5], outputFile);
//					update_model3(symbol, context + 5);
//					index_poss = symbol - 1;
//					ref = (uex[0] + uex[1]) / 2;
//					index_poss_cr = 0;
//					for (j = 0; j < 2 * MAX_NICE + 1; j++)
//					{
//						possible_value = ref + of_interest[j];
//						if ((possible_value != uex[0]) & (possible_value != uex[1]))
//							if (index_poss_cr == index_poss)
//							{
//							valr = possible_value;
//							break;
//							}
//							else
//							{
//								index_poss_cr = index_poss_cr + 1;
//							}
//					}
//				}
//				else
//				{ // Decode val
//					context = 0;
//					valr = decode_symbol(cum_freq30[context], outputFile);
//					update_model3(valr, context);
//				}
//			}
//		}
//		else
//		{	// regions with 2 distinct neighbours far apart
//			if (ENCODEIT == 1)
//			{
//				//printf("case4\n");
//				context = 2;
//				isemptyindex = 1;
//				for (j = 1; j < MAX_NICE; j++)
//				{
//					possible_value = uex[0] + of_interest[j];
//					if (possible_value == val)
//					{
//						index_poss = 2 * (j - 1);
//						isemptyindex = 0;
//						break;
//					}
//					possible_value = uex[1] + of_interest[j];
//					if (possible_value == val)
//					{
//						index_poss = 2 * (j - 1) + 1;
//						isemptyindex = 0;
//						break;
//					}
//				}
//			}
//			else    // if(DECODE == 1)
//			{
//				context = 2;
//				// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
//				symbol = decode_symbol(cum_freq30[context + 10], outputFile);
//				update_model3(symbol, context + 10);
//				isemptyindex = symbol - 1;
//				// decode first isemptyindex
//				if (isemptyindex == 0)
//				{
//					symbol = decode_symbol(cum_freq30[context + 5], outputFile);
//					update_model3(symbol, context + 5);
//					index_poss = symbol - 1;
//					for (j = 1; j < MAX_NICE; j++)
//					{
//						possible_value = uex[0] + of_interest[j];
//						if (index_poss == 2 * (j - 1))
//						{
//							valr = possible_value;
//							break;
//						}
//						possible_value = uex[1] + of_interest[j];
//						if (index_poss == 2 * (j - 1) + 1)
//						{
//							valr = possible_value;
//							break;
//						}
//					}
//				}
//				else
//				{ // decode val
//					context = 0;
//					valr = decode_symbol(cum_freq30[context], outputFile);
//					update_model3(valr, context);
//				}
//			}
//		}
//	}
//		break;
//	default:
//	{
//		if (0)
//		{
//			// find  the median  of uex
//			nbest = nuex; ibest = 0;
//			for (i = 0; i < nuex; i++)
//			{ // candidate median is uex[i]
//				ncou = 0;
//				for (j = 0; j < nuex; j++)
//					if (uex[j] < uex[i])
//						ncou = ncou + 1;
//				if (abs(ncou - (nuex / 2)) < nbest)
//				{
//					nbest = abs(ncou - (nuex / 2));
//					ibest = i;
//				}
//			}
//			ref = uex[ibest];
//			//for (i=0;i<nex;i++) 
//			//   printf(" excluded [%d]   [%d]    \n",i,excluded[i]);
//			//for (i=0;i<nuex;i++) 
//			//   printf(" uex [%d]   [%d]    \n",i,uex[i]);
//			//printf(" ibest [%d]  uexbest  [%d]    \n",ibest,uex[ibest]);
//			n_groups = 1;
//		}
//		else
//		{
//
//			// more than 2 neighbors; Partition them into distinct groups, around the points in "centers"
//			icenter[0] = uex[0];
//			center[0] = (double)uex[0];
//			card_center[0] = 1;
//			ncen = 1;
//			for (i = 1; i < nuex; i++)
//			{
//				classified = 0;
//				for (j = 0; j < ncen; j++)
//				{
//					if (abs(uex[i] - icenter[j]) < MAX_NICE)
//					{
//						card_center[j] = card_center[j] + 1;
//						center[j] = (uex[i] + (card_center[j] - 1)*center[j]) / card_center[j];
//						icenter[j] = (int)floor(center[j] + 0.5);
//						classified = 1;
//						break;
//					}
//				}   // end of for(i=1; i<nuex; i++)
//				if (classified == 0)
//				{
//					ncen = ncen + 1;
//					center[ncen - 1] = (double)uex[i];
//					icenter[ncen - 1] = uex[i];
//					card_center[ncen - 1] = 1;
//				}
//			}   // for(i=1; i<nuex; i++)
//			n_groups = 0;
//			// select the two most populated regions
//			if (ncen == 1)
//			{
//				n_groups = 1;
//				i1 = 0;    // index of (the only one and best) center
//				ref = icenter[i1];
//			}
//			else
//			{
//				if (ncen == 2)
//				{
//					n_groups = 2;
//					i1 = 0; i2 = 1;
//				}
//				else   // else of if(ncen== 2
//				{
//					n_groups = 2;
//					i1 = 0;
//					for (i = 1; i < ncen; i++)
//						if (card_center[i1] < card_center[i])
//							i1 = i;
//					if (i1 != 0)
//						i2 = 0;
//					else
//						i2 = 1;
//
//					for (i = 0; i < ncen; i++)
//						if ((card_center[i2] < card_center[i]) & (i != i1))
//							i2 = i;
//				}  // end of if(ncen == 2
//			}   //  end of  if(ncen == 1)
//		}
//		if (n_groups == 2)
//		{
//			if (abs(icenter[i1] - icenter[i2]) > MAX_NICE)
//			{	// regions with 2 neighbours far apart
//				context = 4;
//				if (ENCODEIT == 1)
//				{
//					isemptyindex = 1;
//					index_poss = 0;
//					for (j = 0; j < MAX_NICE; j++)
//					{
//						possible_value = icenter[i1] + of_interest[j];
//						possible = 1;
//						for (k = 0; k < nuex; k++)
//							if (possible_value == uex[k])
//								possible = 0;
//						if (possible)
//						{
//							if (possible_value == val)
//							{
//								isemptyindex = 0;
//								break;
//							}
//							else
//							{
//								index_poss = index_poss + 1;
//							}
//						}
//						possible_value = icenter[i2] + of_interest[j];
//						possible = 1;
//						for (k = 0; k < nuex; k++)
//							if (possible_value == uex[k])
//								possible = 0;
//						if (possible)
//						{
//							if (possible_value == val)
//							{
//								isemptyindex = 0;
//								break;
//							}
//							else
//							{
//								index_poss = index_poss + 1;
//							}
//						}
//					}   // end of for(j=1; j<MAX_NICE; j++)
//				}
//				else    ///if(DECODE == 1)
//				{	// decode first isemptyindex
//					// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
//					symbol = decode_symbol(cum_freq30[context + 10], outputFile);
//					update_model3(symbol, context + 10);
//					isemptyindex = symbol - 1;
//					if (isemptyindex == 0)
//					{
//						symbol = decode_symbol(cum_freq30[context + 5], outputFile);
//						update_model3(symbol, context + 5);
//						index_poss = symbol - 1;
//						index_poss_cr = 0;
//						for (j = 0; j < MAX_NICE; j++)
//						{
//							possible_value = icenter[i1] + of_interest[j];
//							possible = 1;
//							for (k = 0; k < nuex; k++)
//								if (possible_value == uex[k])
//									possible = 0;
//							if (possible)
//							{
//								if (index_poss_cr == index_poss)
//								{
//									valr = possible_value;
//									isemptyindex = 0;
//									break;
//								}
//								else
//								{
//									index_poss_cr = index_poss_cr + 1;
//								}
//							}
//							possible_value = icenter[i2] + of_interest[j];
//							possible = 1;
//							for (k = 0; k < nuex; k++)
//								if (possible_value == uex[k])
//									possible = 0;
//							if (possible)
//							{
//								if (index_poss_cr == index_poss)
//								{
//									valr = possible_value;
//									isemptyindex = 0;
//									break;
//								}
//								else
//								{
//									index_poss_cr = index_poss_cr + 1;
//								}
//							}
//						}   // end of for(j=1; j<MAX_NICE; j++)
//					}
//					else
//					{	// decode val
//						context = 0;
//						valr = decode_symbol(cum_freq30[context], outputFile);
//						update_model3(valr, context);
//					}
//				}
//			}
//			else // if(abs(icenter[i1]-icenter[i2]) > MAX_NICE)
//			{   // two regions close to each other  
//				n_groups = 1;
//				ref = (icenter[i1] + icenter[i2]) / 2;
//			}
//		}   // end of if(n_groups == 2)
//		if (n_groups == 1)
//		{
//			context = 3;
//			if (ENCODEIT == 1)
//			{
//				isemptyindex = 1;
//				index_poss = 0;
//				for (j = 0; j < 2 * MAX_NICE + 1; j++)
//				{
//					possible_value = ref + of_interest[j];
//					possible = 1;
//					for (k = 0; k < nuex; k++)
//						if (possible_value == uex[k])
//							possible = 0;
//					if (possible == 1)
//						if (possible_value == val)
//						{
//						isemptyindex = 0;
//						break;
//						}
//						else
//						{
//							index_poss = index_poss + 1;
//						}
//				}  // end of for(j=1; j<2*MAX_NICE+1; j++)
//			}
//			else    //  if(DECODE == 1)
//			{	// first decode isemptyindex
//				// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
//				symbol = decode_symbol(cum_freq30[context + 10], outputFile);
//				update_model3(symbol, context + 10);
//				isemptyindex = symbol - 1;
//				if (isemptyindex == 0)
//				{  // decode index_poss
//					symbol = decode_symbol(cum_freq30[context + 5], outputFile);
//					update_model3(symbol, context + 5);
//					index_poss = symbol - 1;
//					index_poss_cr = 0;
//					for (j = 0; j < 2 * MAX_NICE + 1; j++)
//					{
//						possible_value = ref + of_interest[j];
//						possible = 1;
//						for (k = 0; k < nuex; k++)
//							if (possible_value == uex[k])
//								possible = 0;
//						if (possible == 1)
//							if (index_poss == index_poss_cr)
//							{
//							valr = possible_value;
//							break;
//							}
//							else
//							{
//								index_poss_cr = index_poss_cr + 1;
//							}
//					}  // end of for(j=1; j<2*MAX_NICE+1; j++)
//				}
//				else
//				{	//decode val
//					context = 0;
//					valr = decode_symbol(cum_freq30[context], outputFile);
//					update_model3(valr, context);
//				}
//			}
//		}  // end of if(n_groups == 1)
//	}  // end of default
//		break;
//	}		// end switch
//
//	// Encode the decision
//	//if(ir==273)
//	//	printf(" istart [%d] iend [%d]  isemptyindex  [%d] \n",istart,iend,isemptyindex);
//
//	if (ENCODEIT == 1)
//	{
//		encode_symbol(isemptyindex + 1, cum_freq30[context + 10], outputFile);
//		CL1 = CL - log((double)(cum_freq30[context + 10][isemptyindex] - cum_freq30[context + 10][isemptyindex + 1]) / cum_freq30[context + 10][0]) / log(2);
//		CL10 = -log((double)(cum_freq30[context + 10][isemptyindex] - cum_freq30[context + 10][isemptyindex + 1]) / cum_freq30[context + 10][0]) / log(2);
//		//printf("CL10 %f\n", CL10);
//		CLsw[context] = CLsw[context] + CL10;
//		update_model3(isemptyindex + 1, context + 10);
//		//printf("update model3\n");
//		//if(ir==273)
//		//	printf(" ista [%d] iend [%d]  CL1:  [%f] \n",istart,iend,CL1);
//		if (isemptyindex == 1)
//		{	// encode the value
//			//printf("start if0\n");
//			context = 0;
//			//printf("val %d context %d\n", val, context);
//			//printf("cum_freq30[context][val - 1] %d\n cum_freq30[context][val] %d \n", cum_freq30[context][val - 1],
//				//cum_freq30[context][val]);
//			encode_symbol(val, cum_freq30[context], outputFile);
//			CL1 = CL1 - log((double)(cum_freq30[context][val - 1] - cum_freq30[context][val]) / cum_freq30[context][0]) / log(2); //goes to -Inf
//			//fwrite(&ir, sizeof(int), 1, C4Matlab);
//			//fwrite(&istart, sizeof(int), 1, C4Matlab);
//			CL10 = -log((double)(cum_freq30[context][val - 1] - cum_freq30[context][val]) / cum_freq30[context][0]) / log(2);
//			CLval[context] = CLval[context] + CL10;
//			update_model3(val, context);
//			//printf("end if0\n");
//		}
//		else
//		{	// encode the index
//			//printf("start if1\n");
//			encode_symbol(index_poss + 1, cum_freq30[context + 5], outputFile);
//			CL1 = CL1 - log((double)(cum_freq30[context + 5][index_poss] - cum_freq30[context + 5][index_poss + 1]) / cum_freq30[context + 5][0]) / log(2);
//			CL10 = -log((double)(cum_freq30[context + 5][index_poss] - cum_freq30[context + 5][index_poss + 1]) / cum_freq30[context + 5][0]) / log(2);
//			CLind[context] = CLind[context] + CL10;
//			update_model3(index_poss + 1, context + 5);
//			//printf("end if1\n");
//		}
//		for (i = istart; i <= iend; i++)
//		{
//			Himg[ir][i] = val;
//			valregions[Cimg[ir][i]] = Himg[ir][i];
//		}
//		//printf("end if2\n");
//	}
//	else  // DECODE =1
//	{
//		for (i = istart; i <= iend; i++)
//		{
//			Himg[ir][i] = valr;
//			valregions[Cimg[ir][i]] = Himg[ir][i];
//		}
//	}
//
//	free(uex);
//	free(of_interest);
//
//	return CL1;
//}
//




double TreatSegment(int istart, int iend, int known, int val, int* excluded, int nex, int ir, double CL, int ENCODEIT, int* valregions)
{
	double CL1;
	int i, j, k, context, ref;
	int* uex;
	int* of_interest;
	int nuex;
	double center[40000];
	int icenter[40000], card_center[40000], ncen;
	int isemptyindex, possible, possible_value, index_poss, i1, i2, unique;
	int classified, n_groups;
	double CL10;
	int ip;
	int contextp;
	int nbest, ibest, ncou;
	int index_poss_cr, valr, symbol;

	CL1 = CL;
	if (known > 0)
	{	// no need to encode anything
		for (i = istart; i <= iend; i++)
		{
			Himg[ir][i] = known;
			valregions[Cimg[ir][i]] = Himg[ir][i];
		}
		return CL1;
	}


	if (valregions[Cimg[ir][istart]] > 0)
	{	// no need to encode anything
		if ((valregions[Cimg[ir][istart]] != val) && (ENCODEIT == 1))
			printf(" Errors  [%d] [%d]  \n", valregions[Cimg[ir][istart]], val);
		for (i = istart; i <= iend; i++)
		{
			Himg[ir][i] = valregions[Cimg[ir][istart]];
			valregions[Cimg[ir][i]] = Himg[ir][i];
		}
		return CL1;
	}

	uex = alocaVector(NC + 1);   // unique values in excluded
	of_interest = alocaVector(2 * MAX_NICE + 1);  // displacements of interests
	of_interest[0] = 0;
	for (i = 1; i < MAX_NICE; i++)
	{
		of_interest[2 * i] = -i;
		of_interest[2 * i + 1] = i;
	}


	//  We need to encode: Here the simplest case: no neighbors

	//	if(ir==4)
	//		scanf(&i);
	if (nex == 0)  // 0 known neighbors: Encode and exit
	{
		if (ENCODEIT == 1)
		{
			//printf("case1\n");
			context = 0;
			//encode_symbol(val, cum_freq30[context], outputFile);   // val is larger than 0 already
			//CL1 = CL1 - log((double)(cum_freq30[context][val - 1] - cum_freq30[context][val]) / cum_freq30[context][0]) / log(2);
			//CL10 = -log((double)(cum_freq30[context][val - 1] - cum_freq30[context][val]) / cum_freq30[context][0]) / log(2);
			//CLval[context] = CLval[context] + CL10;
			//update_model3(val, context);
			for (i = istart; i <= iend; i++)
			{
				Himg[ir][i] = val;
				valregions[Cimg[ir][i]] = Himg[ir][i];
			}
			free(uex);
			free(of_interest);
			return CL1;
		}
		else  // DECODE = 1
		{
			context = 0;
			//valr = decode_symbol(cum_freq30[context], outputFile);
			//update_model3(valr, context);
			for (i = istart; i <= iend; i++)
			{
				Himg[ir][i] = valr;
				valregions[Cimg[ir][i]] = Himg[ir][i];
			}
			free(uex);
			free(of_interest);
			return CL1;
		}
	}


	// We have some neighbors: Find the distinct neighbours and put them in uex
	nuex = 1;
	uex[0] = excluded[0];
	for (i = 1; i < nex; i++)
	{
		unique = 1;
		for (j = 0; j < nuex; j++)
			if (uex[j] == excluded[i])
			{
				unique = 0;
				break;
			}
			if (unique == 1)
			{
				nuex = nuex + 1;
				uex[nuex - 1] = excluded[i];
			}
	}

	// the main separation to cases depending on number of distinct neighbors
	switch (nuex)	// depending on how many distinct values
	{
	case 1:
		{
			if (ENCODEIT == 1)
			{
				//printf("case2\n");
				// regions with 1 distinct neighbour
				context = 0;
				isemptyindex = 1;
				for (j = 1; j < 2 * MAX_NICE + 1; j++)
				{
					possible_value = uex[0] + of_interest[j];
					if (possible_value == val)
					{
						index_poss = j - 1;
						isemptyindex = 0;
					}
				}
			}
			else // if(DECODE == 1)
			{
				context = 0;
				// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
				//symbol = decode_symbol(cum_freq30[context + 10], outputFile);
				//update_model3(symbol, context + 10);
				isemptyindex = symbol - 1;
				// decode first isemptyindex
				if (isemptyindex == 0)
				{
					//symbol = decode_symbol(cum_freq30[context + 5], outputFile);
					//update_model3(symbol, context + 5);
					index_poss = symbol - 1;
					valr = uex[0] + of_interest[index_poss + 1];
				}
				else
				{
					// decode val
					context = 0;
					//valr = decode_symbol(cum_freq30[context], outputFile);
					//update_model3(valr, context);
				}
			}
		}
		break;
	case 2:
		{  // regions with 2 distinct neighbours close to each other
			if (abs(uex[0] - uex[1]) < MAX_NICE)
			{
				if (ENCODEIT == 1)
				{
					//printf("case3\n");
					context = 1;
					ref = (uex[0] + uex[1]) / 2;
					isemptyindex = 1;
					index_poss = 0;
					for (j = 0; j < 2 * MAX_NICE + 1; j++)
					{
						possible_value = ref + of_interest[j];
						if ((possible_value != uex[0]) & (possible_value != uex[1]))
							if (possible_value == val)
							{
								isemptyindex = 0;
								break;
							}
							else
							{
								index_poss = index_poss + 1;
							}
					}
				}
				else  //  if(DECODE == 1)
				{
					context = 1;
					// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
					//symbol = decode_symbol(cum_freq30[context + 10], outputFile);
					//update_model3(symbol, context + 10);
					isemptyindex = symbol - 1;
					// decode first isemptyindex
					if (isemptyindex == 0)
					{
						//symbol = decode_symbol(cum_freq30[context + 5], outputFile);
						//update_model3(symbol, context + 5);
						index_poss = symbol - 1;
						ref = (uex[0] + uex[1]) / 2;
						index_poss_cr = 0;
						for (j = 0; j < 2 * MAX_NICE + 1; j++)
						{
							possible_value = ref + of_interest[j];
							if ((possible_value != uex[0]) & (possible_value != uex[1]))
								if (index_poss_cr == index_poss)
								{
									valr = possible_value;
									break;
								}
								else
								{
									index_poss_cr = index_poss_cr + 1;
								}
						}
					}
					else
					{ // Decode val
						context = 0;
						//valr = decode_symbol(cum_freq30[context], outputFile);
						//update_model3(valr, context);
					}
				}
			}
			else
			{	// regions with 2 distinct neighbours far apart
				if (ENCODEIT == 1)
				{
					//printf("case4\n");
					context = 2;
					isemptyindex = 1;
					for (j = 1; j < MAX_NICE; j++)
					{
						possible_value = uex[0] + of_interest[j];
						if (possible_value == val)
						{
							index_poss = 2 * (j - 1);
							isemptyindex = 0;
							break;
						}
						possible_value = uex[1] + of_interest[j];
						if (possible_value == val)
						{
							index_poss = 2 * (j - 1) + 1;
							isemptyindex = 0;
							break;
						}
					}
				}
				else    // if(DECODE == 1)
				{
					context = 2;
					// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
					//symbol = decode_symbol(cum_freq30[context + 10], outputFile);
					//update_model3(symbol, context + 10);
					isemptyindex = symbol - 1;
					// decode first isemptyindex
					if (isemptyindex == 0)
					{
						//symbol = decode_symbol(cum_freq30[context + 5], outputFile);
						//update_model3(symbol, context + 5);
						index_poss = symbol - 1;
						for (j = 1; j < MAX_NICE; j++)
						{
							possible_value = uex[0] + of_interest[j];
							if (index_poss == 2 * (j - 1))
							{
								valr = possible_value;
								break;
							}
							possible_value = uex[1] + of_interest[j];
							if (index_poss == 2 * (j - 1) + 1)
							{
								valr = possible_value;
								break;
							}
						}
					}
					else
					{ // decode val
						context = 0;
						//valr = decode_symbol(cum_freq30[context], outputFile);
						//update_model3(valr, context);
					}
				}
			}
		}
		break;
	default:
		{
			if (0)
			{
				// find  the median  of uex
				nbest = nuex; ibest = 0;
				for (i = 0; i < nuex; i++)
				{ // candidate median is uex[i]
					ncou = 0;
					for (j = 0; j < nuex; j++)
						if (uex[j] < uex[i])
							ncou = ncou + 1;
					if (abs(ncou - (nuex / 2)) < nbest)
					{
						nbest = abs(ncou - (nuex / 2));
						ibest = i;
					}
				}
				ref = uex[ibest];
				//for (i=0;i<nex;i++) 
				//   printf(" excluded [%d]   [%d]    \n",i,excluded[i]);
				//for (i=0;i<nuex;i++) 
				//   printf(" uex [%d]   [%d]    \n",i,uex[i]);
				//printf(" ibest [%d]  uexbest  [%d]    \n",ibest,uex[ibest]);
				n_groups = 1;
			}
			else
			{

				// more than 2 neighbors; Partition them into distinct groups, around the points in "centers"
				icenter[0] = uex[0];
				center[0] = (double)uex[0];
				card_center[0] = 1;
				ncen = 1;
				for (i = 1; i < nuex; i++)
				{
					classified = 0;
					for (j = 0; j < ncen; j++)
					{
						if (abs(uex[i] - icenter[j]) < MAX_NICE)
						{
							card_center[j] = card_center[j] + 1;
							center[j] = (uex[i] + (card_center[j] - 1)*center[j]) / card_center[j];
							icenter[j] = (int)floor(center[j] + 0.5);
							classified = 1;
							break;
						}
					}   // end of for(i=1; i<nuex; i++)
					if (classified == 0)
					{
						ncen = ncen + 1;
						center[ncen - 1] = (double)uex[i];
						icenter[ncen - 1] = uex[i];
						card_center[ncen - 1] = 1;
					}
				}   // for(i=1; i<nuex; i++)
				n_groups = 0;
				// select the two most populated regions
				if (ncen == 1)
				{
					n_groups = 1;
					i1 = 0;    // index of (the only one and best) center
					ref = icenter[i1];
				}
				else
				{
					if (ncen == 2)
					{
						n_groups = 2;
						i1 = 0; i2 = 1;
					}
					else   // else of if(ncen== 2
					{
						n_groups = 2;
						i1 = 0;
						for (i = 1; i < ncen; i++)
							if (card_center[i1] < card_center[i])
								i1 = i;
						if (i1 != 0)
							i2 = 0;
						else
							i2 = 1;

						for (i = 0; i < ncen; i++)
							if ((card_center[i2] < card_center[i]) & (i != i1))
								i2 = i;
					}  // end of if(ncen == 2
				}   //  end of  if(ncen == 1)
			}
			if (n_groups == 2)
			{
				if (abs(icenter[i1] - icenter[i2]) > MAX_NICE)
				{	// regions with 2 neighbours far apart
					context = 4;
					if (ENCODEIT == 1)
					{
						isemptyindex = 1;
						index_poss = 0;
						for (j = 0; j < MAX_NICE; j++)
						{
							possible_value = icenter[i1] + of_interest[j];
							possible = 1;
							for (k = 0; k < nuex; k++)
								if (possible_value == uex[k])
									possible = 0;
							if (possible)
							{
								if (possible_value == val)
								{
									isemptyindex = 0;
									break;
								}
								else
								{
									index_poss = index_poss + 1;
								}
							}
							possible_value = icenter[i2] + of_interest[j];
							possible = 1;
							for (k = 0; k < nuex; k++)
								if (possible_value == uex[k])
									possible = 0;
							if (possible)
							{
								if (possible_value == val)
								{
									isemptyindex = 0;
									break;
								}
								else
								{
									index_poss = index_poss + 1;
								}
							}
						}   // end of for(j=1; j<MAX_NICE; j++)
					}
					else    ///if(DECODE == 1)
					{	// decode first isemptyindex
						// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
						//symbol = decode_symbol(cum_freq30[context + 10], outputFile);
						//update_model3(symbol, context + 10);
						isemptyindex = symbol - 1;
						if (isemptyindex == 0)
						{
							/*	symbol = decode_symbol(cum_freq30[context + 5], outputFile);
							update_model3(symbol, context + 5);
							*/			index_poss = symbol - 1;
						index_poss_cr = 0;
						for (j = 0; j < MAX_NICE; j++)
						{
							possible_value = icenter[i1] + of_interest[j];
							possible = 1;
							for (k = 0; k < nuex; k++)
								if (possible_value == uex[k])
									possible = 0;
							if (possible)
							{
								if (index_poss_cr == index_poss)
								{
									valr = possible_value;
									isemptyindex = 0;
									break;
								}
								else
								{
									index_poss_cr = index_poss_cr + 1;
								}
							}
							possible_value = icenter[i2] + of_interest[j];
							possible = 1;
							for (k = 0; k < nuex; k++)
								if (possible_value == uex[k])
									possible = 0;
							if (possible)
							{
								if (index_poss_cr == index_poss)
								{
									valr = possible_value;
									isemptyindex = 0;
									break;
								}
								else
								{
									index_poss_cr = index_poss_cr + 1;
								}
							}
						}   // end of for(j=1; j<MAX_NICE; j++)
						}
						else
						{	// decode val
							context = 0;
							//valr = decode_symbol(cum_freq30[context], outputFile);
							//update_model3(valr, context);
						}
					}
				}
				else // if(abs(icenter[i1]-icenter[i2]) > MAX_NICE)
				{   // two regions close to each other  
					n_groups = 1;
					ref = (icenter[i1] + icenter[i2]) / 2;
				}
			}   // end of if(n_groups == 2)
			if (n_groups == 1)
			{
				context = 3;
				if (ENCODEIT == 1)
				{
					isemptyindex = 1;
					index_poss = 0;
					for (j = 0; j < 2 * MAX_NICE + 1; j++)
					{
						possible_value = ref + of_interest[j];
						possible = 1;
						for (k = 0; k < nuex; k++)
							if (possible_value == uex[k])
								possible = 0;
						if (possible == 1)
							if (possible_value == val)
							{
								isemptyindex = 0;
								break;
							}
							else
							{
								index_poss = index_poss + 1;
							}
					}  // end of for(j=1; j<2*MAX_NICE+1; j++)
				}
				else    //  if(DECODE == 1)
				{	// first decode isemptyindex
					// Decode first isemptyindex, and if isemptyindex==0 decode also index_poss
					//symbol = decode_symbol(cum_freq30[context + 10], outputFile);
					//update_model3(symbol, context + 10);
					isemptyindex = symbol - 1;
					if (isemptyindex == 0)
					{  // decode index_poss
						//symbol = decode_symbol(cum_freq30[context + 5], outputFile);
						//update_model3(symbol, context + 5);
						index_poss = symbol - 1;
						index_poss_cr = 0;
						for (j = 0; j < 2 * MAX_NICE + 1; j++)
						{
							possible_value = ref + of_interest[j];
							possible = 1;
							for (k = 0; k < nuex; k++)
								if (possible_value == uex[k])
									possible = 0;
							if (possible == 1)
								if (index_poss == index_poss_cr)
								{
									valr = possible_value;
									break;
								}
								else
								{
									index_poss_cr = index_poss_cr + 1;
								}
						}  // end of for(j=1; j<2*MAX_NICE+1; j++)
					}
					else
					{	//decode val
						context = 0;
						//valr = decode_symbol(cum_freq30[context], outputFile);
						//update_model3(valr, context);
					}
				}
			}  // end of if(n_groups == 1)
		}  // end of default
		break;
	}		// end switch

	// Encode the decision
	//if(ir==273)
	//	printf(" istart [%d] iend [%d]  isemptyindex  [%d] \n",istart,iend,isemptyindex);

	if (ENCODEIT == 1)
	{
		/*encode_symbol(isemptyindex + 1, cum_freq30[context + 10], outputFile);
		CL1 = CL - log((double)(cum_freq30[context + 10][isemptyindex] - cum_freq30[context + 10][isemptyindex + 1]) / cum_freq30[context + 10][0]) / log(2);
		CL10 = -log((double)(cum_freq30[context + 10][isemptyindex] - cum_freq30[context + 10][isemptyindex + 1]) / cum_freq30[context + 10][0]) / log(2);*/
		//printf("CL10 %f\n", CL10);
		//CLsw[context] = CLsw[context] + CL10;
		//update_model3(isemptyindex + 1, context + 10);
		//printf("update model3\n");
		//if(ir==273)
		//	printf(" ista [%d] iend [%d]  CL1:  [%f] \n",istart,iend,CL1);
		if (isemptyindex == 1)
		{	// encode the value
			//printf("start if0\n");
			context = 0;
			//printf("val %d context %d\n", val, context);
			//printf("cum_freq30[context][val - 1] %d\n cum_freq30[context][val] %d \n", cum_freq30[context][val - 1],
			//cum_freq30[context][val]);
			//encode_symbol(val, cum_freq30[context], outputFile);
			//CL1 = CL1 - log((double)(cum_freq30[context][val - 1] - cum_freq30[context][val]) / cum_freq30[context][0]) / log(2); //goes to -Inf
			////fwrite(&ir, sizeof(int), 1, C4Matlab);
			////fwrite(&istart, sizeof(int), 1, C4Matlab);
			//CL10 = -log((double)(cum_freq30[context][val - 1] - cum_freq30[context][val]) / cum_freq30[context][0]) / log(2);
			//CLval[context] = CLval[context] + CL10;
			//update_model3(val, context);
			//printf("end if0\n");
		}
		else
		{	// encode the index
			//printf("start if1\n");
			/*encode_symbol(index_poss + 1, cum_freq30[context + 5], outputFile);
			CL1 = CL1 - log((double)(cum_freq30[context + 5][index_poss] - cum_freq30[context + 5][index_poss + 1]) / cum_freq30[context + 5][0]) / log(2);
			CL10 = -log((double)(cum_freq30[context + 5][index_poss] - cum_freq30[context + 5][index_poss + 1]) / cum_freq30[context + 5][0]) / log(2);
			CLind[context] = CLind[context] + CL10;
			update_model3(index_poss + 1, context + 5);*/
			//printf("end if1\n");
		}
		for (i = istart; i <= iend; i++)
		{
			Himg[ir][i] = val;
			valregions[Cimg[ir][i]] = Himg[ir][i];
		}
		//printf("end if2\n");
	}
	else  // DECODE =1
	{
		for (i = istart; i <= iend; i++)
		{
			Himg[ir][i] = valr;
			valregions[Cimg[ir][i]] = Himg[ir][i];
		}
	}

	free(uex);
	free(of_interest);

	return CL1;
}


void encodeDepth(int method, int newregindex)
{
	int i, ic, ir, ns, context, symbol;
	int* putere;
	int* of_interest;
	int* excluded;
	double CL;
	int istart, iend, val, known, nex;
	int ENCODEIT = 1;
	int* valregions;

	ns = 4;
	My_no_of_symbols = ns;
	My_no_of_rows = NS_NCON; // 1 + 3*(1 + 4 + 4^2 + 4^3 + 4^4)
	MODEL = method;

	for (context = 0; context < 5; context++)
	{
		CLval[context] = 0;
		CLind[context] = 0;
		CLsw[context] = 0;
	}
	start_model3();

	putere = alocaVector(ns);
	putere[0] = ns;
	for (i = 1; i < ns; i++)
		putere[i] = putere[i - 1] * ns;

	of_interest = alocaVector(2 * MAX_NICE + 1);
	of_interest[0] = 0;
	for (i = 1; i < MAX_NICE; i++)
	{
		of_interest[2 * i] = -i;
		of_interest[2 * i + 1] = i;
	}
	excluded = alocaVector(NC + 2);

	valregions = alocaVector(newregindex + 2);

	CL = 0.;
	for (ir = 0; ir < NR; ir++)
	{
		// printf("%d\n", ir);
		// Inspect the current line only once
		// Prepare for the first segment
		istart = 0;
		known = 0;
		for (i = 0; i < NC + 2; i++)
			excluded[i] = 0;
		nex = 0;
		for (ic = 0; ic < NC; ic++)
		{
			if (gSR[ir][ic] % 2 == 1)
			{
				// There is a vertical edge at the left of ic
				// stop here for a while, treat the segment istart:(ic-1)
				iend = ic - 1;
				val = img[ir][iend] + 1;
				CL = TreatSegment(istart, iend, known, val, excluded, nex, ir, CL, ENCODEIT, valregions);
				// printf(" istart[%d]   iend [%d]  known [%d] val [%d]  nex [%d]  \n",istart,iend,known,val,nex);
				// now start a new segment
				istart = ic;
				known = 0;
				for (i = 1; i < NC + 2; i++)
					excluded[i] = 0;
				excluded[0] = Himg[ir][ic - 1]; nex = 1;
			}
			// if there is no hor. edge above a pixel, the pixel will take the value from the above pixel
			if (gSR[ir][ic] < 2)
			{
				if (ir > 0)
				{
					Himg[ir][ic] = Himg[ir - 1][ic];
					known = Himg[ir - 1][ic];
					valregions[Cimg[ir][ic]] = Himg[ir][ic];
				}
			}
			else
				// There is horizontal edge above gray2(ir,ic)
			{
				if (ir > 0)
				{
					nex = nex + 1;
					excluded[nex - 1] = Himg[ir - 1][ic];
				}
			}
			if (ic == (NC - 1))  // end of row
			{
				iend = (NC - 1);
				val = img[ir][iend] + 1;
				CL = TreatSegment(istart, iend, known, val, excluded, nex, ir, CL, ENCODEIT, valregions);
			}

		}
	}

	for (context = 0; context < 5; context++)
	{
		printf("cont [%d] CLval,CLind ,CLsw [%.3f] [%.3f] [%.3f] \n",context, CLval[context] ,CLind[context] ,CLsw[context]);
	}

	//printf("Code length: [%.3f] \n",CL);
	free_model3();
	free(putere);
	free(of_interest);
	free(valregions);
}
int RegTreatSegment(int istart, int iend, int newregindex, int* sameregs, int nsame, int ir, int* minlabels, int* usamereg)
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


int MarkRegions()
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
			return newregindex;
}

void decodeDepth(int method, int newregindex)
{
	int i, ic, ir, ns, context, symbol;
	int* putere;
	int* of_interest;
	int* excluded;
	double CL;
	int istart, iend, val = 0, known, nex;
	int ENCODEIT = 0;
	int* valregions;

	valregions = alocaVector(newregindex + 2);

	ns = 4;
	My_no_of_symbols = ns;
	My_no_of_rows = NS_NCON; // 1 + 3*(1 + 4 + 4^2 + 4^3 + 4^4)
	MODEL = method;

	for (context = 0; context < 5; context++)
	{
		CLval[context] = 0;
		CLind[context] = 0;
		CLsw[context] = 0;
	}
	start_model3();

	putere = alocaVector(ns);
	putere[0] = ns;
	for (i = 1; i < ns; i++)
		putere[i] = putere[i - 1] * ns;

	of_interest = alocaVector(2 * MAX_NICE + 1);
	of_interest[0] = 0;
	for (i = 1; i < MAX_NICE; i++)
	{
		of_interest[2 * i] = -i;
		of_interest[2 * i + 1] = i;
	}
	excluded = alocaVector(NC + 2);
	CL = 0.;

	//printf("%d\n", global_memory_counter);

	for (ir = 0; ir < NR; ir++)
	{
		// Inspect the current line only once
		// Prepare for the first segment
		istart = 0;
		known = 0;
		for (i = 0; i < NC + 2; i++)
			excluded[i] = 0;
		nex = 0;
		for (ic = 0; ic < NC; ic++)
		{
			if (gSR[ir][ic] % 2 == 1)
			{
				// There is a vertical edge at the left of ic
				// stop here for a while, treat the segment istart:(ic-1)
				iend = ic - 1;
				// val = img[ir][iend]+1;
				CL = TreatSegment(istart, iend, known, val, excluded, nex, ir, CL, ENCODEIT, valregions);
				//   printf("Himg [%d] [%d] [%d] [%d]\n",Himg[0][0],Himg[0][1],Himg[0][2],Himg[0][3]);
				// now start a new segment
				istart = ic;
				known = 0;
				for (i = 1; i < NC + 2; i++)
					excluded[i] = 0;
				excluded[0] = Himg[ir][ic - 1]; nex = 1;
			}
			// if there is no hor. edge above a pixel, the pixel will take the value from the above pixel
			if (gSR[ir][ic] < 2)
			{
				if (ir > 0)
				{
					Himg[ir][ic] = Himg[ir - 1][ic];
					known = Himg[ir - 1][ic];
					valregions[Cimg[ir][ic]] = Himg[ir][ic];
				}
			}
			else
				// There is horizontal edge above gray2(ir,ic)
			{
				if (ir > 0)
				{
					nex = nex + 1;
					excluded[nex - 1] = Himg[ir - 1][ic];
				}
			}
			if (ic == (NC - 1))  // end of row
			{
				iend = NC - 1;
				// val = img[ir][iend] + 1;
				CL = TreatSegment(istart, iend, known, val, excluded, nex, ir, CL, ENCODEIT, valregions);
			}

		}
		//printf("%d\n", global_memory_counter);
	}
	for (context = 0; context < 5; context++)
	{
		//printf("cont [%d] CLval,CLind ,CLsw [%.3f] [%.3f] [%.3f] \n",context, CLval[context] ,CLind[context] ,CLsw[context]);
	}

	//printf("Himg [%d] [%d] [%d] [%d]\n",Himg[765][84],Himg[765][85],Himg[765][86],Himg[765][87]);
	//printf("Himg [%d] [%d] [%d] [%d]\n",Himg[766][84],Himg[766][85],Himg[766][86],Himg[766][87]);


	// scanf(&i);
	//printf("Code length: [%.3f] \n",CL);
	free_model3();
	free(putere);
	free(of_interest);
	free(valregions);
}


void cerv_encode(int** iimg, int NRi, int NCi, char* compImageNamei) {
	int i, j;
	int DATA[10];
	int nrBM, nrBN, biti;
	float sec;
	int newregindex;

	int *depth_values;

	NR = NRi;
	NC = NCi;

	printf("%d\n", global_memory_counter);

	//Begin = clock() * CLK_TCK;
	openFiles(1);

	compImgName = compImageNamei;

	openFiles(2); // encode	

	// read input data
	readData();
	//	for (i=0; i<NR; i++)
	//		for (j=0; j<NC; j++)
	//			img[i][j] = iimg[i][j];
	//
	//	// gr(4:2:(2*nr+2), 4:2:(2*nc+2)) = gray2;
	//	for (i=0; i<NR; i++)
	//		for (j=0; j<NC; j++)
	//			gr[3+2*i][3+2*j] = iimg[i][j];


	//fscanf(Matlab2C, "%d", &NR);
	//fscanf(Matlab2C, "%d", &NC);

	gSR = alocaMatrice(NR, NC);
	img = alocaMatrice(NR, NC);
	Cimg = alocaMatrice(NR, NC);

	NRgr = 2 * NR + 4;
	NCgr = 2 * NC + 4;

	gr = alocaMatrice(NRgr, NCgr);

	for (i = 0; i < NR; i++)
		for (j = 0; j < NC; j++)
			img[i][j] = iimg[i][j];

	// gr(4:2:(2*nr+2), 4:2:(2*nc+2)) = gray2;
	for (i = 0; i < NR; i++)
		for (j = 0; j < NC; j++)
			gr[3 + 2 * i][3 + 2 * j] = img[i][j];

	// get gSR image
	get_gSR_matrix();
	// Find the regions
	Begin1 = clock() * CLK_TCK;
	newregindex = MarkRegions();
	printf("newregindex=%d\n",newregindex);
	End1 = clock() * CLK_TCK;
	sec =(float)(End1-Begin1)/1000000;
	printf("Timp de executie MarkRegions: [%.3f] secunde\n",sec);
	//Begin1 = clock() * CLK_TCK;
	//Himg = getRegions(img, NR, NC);
	//End1 = clock() * CLK_TCK;
	//sec =(float)(End1-Begin1)/1000000;
	//printf("Timp de executie GetRegions: [%.3f] secunde\n",sec);
	// start encoding
	//start_outputing_bits();
	//start_encoding();

	// send length and width
	//nrBM = lengthBin(NR);
	//nrBN = lengthBin(NC);
	//biti = (nrBM > nrBN) ? nrBM : nrBN;
	//DATA[0] = biti;
	//encode(1, 16, DATA, NULL);
	//DATA[0] = NR;
	//DATA[1] = NC;
	//encode(2, ((1 << biti) - 1), DATA, NULL);

	fwrite(&NR,sizeof(int),1,outputFile);
	fwrite(&NC,sizeof(int),1,outputFile);

	/* cimg is available */

	printf("newregindex=%d\n",newregindex);

	//depth_values = malloc( (newregindex+1)*sizeof(int));

	//for(i = 0; i < newregindex+1; i++){
	//	depth_values[i] = 0;
	//}

	//for (i = 0; i < NR; i++){
	//	for (j = 0; j < NC; j++){
	//		depth_values[ Cimg[i][j] ] = iimg[i][j];
	//	}
	//}

	//for(i=0;i<newregindex+1;i++){
	//	printf("depth_values[%d]=%d\n",i,depth_values[i]);
	//}

	//free(depth_values);

	// encodeContour1(2);// 1 - Laplace ; 2 - KT
	// encodeContour2(2);
	encodeContour2(1);

	//fwrite(&newregindex,sizeof(int),1,outputFile);
	//fwrite(depth_values,sizeof(int),newregindex+1,outputFile);

	printf("done encodeContour2\n");

	/* NO LABELS ENCODED */
	if (0){
		Himg = alocaMatrice(NR, NC);
		encodeDepth(1, newregindex);
		//printf("done\n");
		// stop encoding
	}

	//done_encoding(outputFile);
	//done_outputing_bits(outputFile);
	// scanf(&i);

	//*/
	//fprintf(resultFile, "0");
	closeFiles();


	if (Himg != NULL)
	{
		for (i = 0; i < NR; i++)
			free(Himg[i]); //KAATUU?
		free(Himg);
	}

	if (Cimg != NULL)
	{
		for (i = 0; i < NR; i++)
			free(Cimg[i]); //KAATUU?
		free(Cimg);
	}

	if (gSR != NULL)
	{
		for (i = 0; i < NR; i++)
			free(gSR[i]); //KAATUU?
		free(gSR);
	}

	printf("%d\n", global_memory_counter);

	//freee();
}

void cerv_decode(int** iimg, int NRi, int NCi, char* compImageNamei) {

	int i, j;
	int DATA[10];
	int nrBM, nrBN, biti;
	float sec;
	int newregindex;

	int *depth_values;

	NR = NRi;
	NC = NCi;

	compImgName = compImageNamei;

	C3Matlab = NULL;
	C4Matlab = NULL;

	printf("%d\n", global_memory_counter);

	openFiles(1);
	openFiles(3); // decode
	//printf("decoding\n");
	//start_inputing_bits();
	//start_decoding(outputFile);
	//printf("decoding\n");


	// read m n
	//decode(1, 16, DATA, 0);
	//biti = DATA[0];
	//decode(2, ((1 << biti) - 1), DATA, 0);

	fread(&NR,sizeof(int),1,outputFile);
	fread(&NC,sizeof(int),1,outputFile);

	// NR = DATA[0];
	// NC = DATA[1];

	printf("NR NC %d %d\n", NR, NC);


	gSR = alocaMatrice(NR, NC);
	img = NULL;

	Himg = alocaMatrice(NR, NC);
	Cimg = alocaMatrice(NR, NC);


	decodeContour2(1);// 1 - Laplace ; 2 - KT
	printf("NR  %d\n", NR );
	newregindex = MarkRegions();
	printf("  NC %d %d\n", NC);

	/*fread(&newregindex,sizeof(int),1,outputFile);

	printf("newregindex=%d\n",newregindex);

	depth_values = malloc((newregindex+1)*sizeof(int));*/

	//for(i = 0; i < newregindex+1; i++){
	//	depth_values[i] = 0;
	//}

	//fread(depth_values,sizeof(int),newregindex+1,outputFile);

	//for(i=0;i<newregindex+1;i++){
	//	printf("depth_values[%d]=%d\n",i,depth_values[i]);
	//}

	/* Cimg contains labels  */
	for (j = 0; j < NC; j++){
		printf("Cimg[i][j]=%d\n",Cimg[i][j]);
		for (i = 0; i < NR; i++)
		{
			iimg[i][j] = Cimg[i][j];
			//iimg[i][j] = depth_values[ Cimg[i][j] ];
		}
	}

	

	//free(depth_values);

	printf("%d\n", global_memory_counter);

	/* NO LABELS DECODED */
	if (0){

		decodeDepth(1, newregindex);

		for (j = 0; j < NC; j++){
			for (i = 0; i < NR; i++){
				//Himg[i][j] = Himg[i][j] - 1;
				iimg[i][j] = Himg[i][j] - 1;
			}
		}

		// clean up
		//printf("done\n");
		//fprintf(resultFile, "0");
	}
	closeFiles();

	if (Himg != NULL)
	{
		for (i = 0; i < NR; i++)
			free(Himg[i]); //KAATUU?
		free(Himg);
	}



	if (gSR != NULL)
	{
		for (i = 0; i < NR; i++)
			free(gSR[i]); //KAATUU?
		free(gSR);
	}

	printf("%d\n", global_memory_counter);


	//freee();
}

