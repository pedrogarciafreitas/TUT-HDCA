#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>
//#include <ctime>

#include "cerv.h"

int main(const int argc, char** argv) {

	int nr, nc;

	int minL,maxL;

	FILE *filept;
	filept = fopen(argv[1], "rb");

	fread(&nr, sizeof(int), 1, filept);
	fread(&nc, sizeof(int), 1, filept);

	printf("NR NC %d %d\n", nr, nc);

	int *Labels = new int[nr*nc]();
	fread(Labels, sizeof(int), nr*nc, filept);
	fclose(filept);

	/* encode quantized labels with cerv */
	int** SEGM2D = (int**)malloc(nr*sizeof(int*));
	for (int i = 0; i < nr; ++i) { SEGM2D[i] = (int*)malloc(nc*sizeof(int)); }
	minL = 65536; maxL = 0;
	for (int i = 0; i < nr; ++i) {
		for (int j = 0; j < nc; ++j) {
			SEGM2D[i][j] = Labels[i + j*nr];
			if(minL > Labels[i + j*nr])
				minL =  Labels[i + j*nr];
			if(maxL < Labels[i + j*nr])
				maxL =  Labels[i + j*nr];
		}
	}
	printf("minL maxL %d %d\n", minL, maxL);
	cerv_encode(SEGM2D, nr, nc, argv[2]);

	printf("Finished encode %d %d\n", minL, maxL);

	int *SEGMFINAL = new int[nr*nc]();

	for (int i = 0; i < nr; ++i) {
		for (int j = 0; j < nc; ++j) {
			SEGM2D[i][j] = 0;
		}
	}

	if(1)
	{
		cerv_decode(SEGM2D, nr, nc, argv[2]);

		int number_of_regions = 0;

		for (int i = 0; i < nr; ++i) {
			for (int j = 0; j < nc; ++j) {
				SEGMFINAL[i + j*nr] = SEGM2D[i][j];
				if (SEGMFINAL[i + j*nr] > number_of_regions)
					number_of_regions++;
			}
		}

		number_of_regions = number_of_regions + 1; //account for 0



		printf("%d\n", number_of_regions);

		filept = fopen(argv[3], "wb");
		fwrite(&nr, sizeof(int), 1, filept);
		fwrite(&nc, sizeof(int), 1, filept);
		fwrite(SEGMFINAL, sizeof(int), nr*nc, filept);
		fclose(filept);
	}

	delete[](Labels);
	delete[](SEGMFINAL);

	return 0;

}