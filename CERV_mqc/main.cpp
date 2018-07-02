#include <iostream>
#include <stdio.h>

#include "cerv.h"

int main(const int argc, char** argv) {

	int nr, nc;

	FILE *filept;
	filept = fopen(argv[1], "rb");

	fread(&nr, sizeof(int), 1, filept);
	fread(&nc, sizeof(int), 1, filept);

	printf("NR NC %d %d\n", nr, nc);

	int *Labels = new int[nr*nc]();
	fread(Labels, sizeof(int), nr*nc, filept);
	fclose(filept);

	CERV *cerv_coder = new CERV;

	cerv_coder->encode(Labels, nr, nc);

	printf("numbytes=%d\n", cerv_coder->getBitstreamLength());

	CERV *cerv_decoder = new CERV;

	cerv_decoder->decode(cerv_coder->getBitstream(), &cerv_coder->n_bytes[0], cerv_coder->depth_values, 
		cerv_coder->getNR(), cerv_coder->getNC());

	int *Labels_out = new int[cerv_coder->getNR()*cerv_coder->getNC()]();

	cerv_decoder->getSegmentation(Labels_out);

	for (int ij = 0; ij < nr*nc; ij++)
		if (Labels_out[ij] != Labels[ij])
			printf("ERROR at %d\n", ij);


	filept = fopen(argv[3], "wb");
	fwrite(&nr, sizeof(int), 1, filept);
	fwrite(&nc, sizeof(int), 1, filept);
	fwrite(Labels_out, sizeof(int), nr*nc, filept);
	fclose(filept);

	delete[](Labels);
	delete[](Labels_out);

	delete cerv_decoder;
	delete cerv_coder;

	return 0;

}