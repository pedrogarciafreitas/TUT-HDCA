/* Disparity refinement over 11x11 window. */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "gen_types.hh"

int main(const int argc, const char** argv) {

	int nr, nc, ncomponents, nvc;

	aux_read_header_file(&nr, &nc, &ncomponents, &nvc, argv[1]);

	int n_views = (argc - 2) / 2;

	unsigned char *AA1 = new unsigned char[nr*nc*ncomponents];
	aux_read_file_uint8(nr, nc, ncomponents, argv[1], AA1);

	float *DispTarg1 = new float[nr*nc];
	aux_read_file_float(nr, nc, ncomponents, argv[2], DispTarg1);


	if (n_views > 1){
		for (int ik = 3; ik < n_views; ik += 2){
			unsigned char *AA2 = new unsigned char[nr*nc*ncomponents];
			aux_read_file_uint8(nr, nc, ncomponents, argv[ik], AA1);

			float *DispTarg2 = new float[nr*nc];
			aux_read_file_float(nr, nc, ncomponents, argv[ik + 1], DispTarg1);

			for (int ij = 0; ij < nr*nc; ij++){


				if (DispTarg1[ij] < 0 & DispTarg2[ij] >= 0){
					DispTarg1[ij] = 1;
					AA1[ij] = AA2[ij];
					AA1[ij + nr*nc] = AA2[ij + nr*nc];
					AA1[ij + nr*nc * 2] = AA2[ij + nr*nc * 2];
				}

			}
		}
	}

	aux_write_header_file(nr, nc, 3, 1, argv[argc - 1]);

	FILE *f_file = fopen(argv[argc - 1], "a+b");

	fwrite(AA1, sizeof(unsigned char), nr*nc * 3, f_file);

	fclose(f_file);

	return 0;

}