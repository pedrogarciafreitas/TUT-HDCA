/* Collect all warped color views. */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "gen_types.hh"

int main(const int argc, const char** argv) {

	int nr, nc, ncomponents, nvc;

	aux_read_header_file(&nr, &nc, &ncomponents, &nvc, argv[1]);

	int n_views = (argc - 2) / 2;

	printf("n_views:\t%d\n", n_views);

	unsigned short *AA1 = new unsigned short[nr*nc*ncomponents];
	aux_read_file_uint16(nr, nc, ncomponents, argv[1], AA1);

	float *DispTarg1 = new float[nr*nc];
	aux_read_file_float(nr, nc, 1, argv[2], DispTarg1);


	if (n_views > 1){
		for (int ik = 3; ik < n_views*2; ik += 2){
			unsigned short *AA2 = new unsigned short[nr*nc*ncomponents];
			aux_read_file_uint16(nr, nc, ncomponents, argv[ik], AA2);

			float *DispTarg2 = new float[nr*nc];
			aux_read_file_float(nr, nc, 1, argv[ik + 1], DispTarg2);

			for (int ij = 0; ij < nr*nc; ij++){
				if (DispTarg1[ij] < 0 && DispTarg2[ij]>=0){
					DispTarg1[ij] = DispTarg2[ij];
					AA1[ij] = AA2[ij];
					AA1[ij + nr*nc] = AA2[ij + nr*nc];
					AA1[ij + nr*nc * 2] = AA2[ij + nr*nc * 2];
				}

			}
			printf("ik:\t%d\n", ik);
		}
	}

	aux_write_header_file(nr, nc, 3, 1, argv[argc - 1]);

	FILE *f_file = fopen(argv[argc - 1], "a+b");

	fwrite(AA1, sizeof(unsigned short), nr*nc * 3, f_file);

	fclose(f_file);

	return 0;

}