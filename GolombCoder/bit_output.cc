// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Simple bit-output.


#include <iostream>
#include "bit_output.hh"

BitOutput::BitOutput(std::string filename): buffer(0), bits_to_go(8), n_bits_(0)
{   
   outputFile = fopen(filename.c_str(), "wb");
}


/* OUTPUT A BIT. */
void BitOutput::output_bit(int bit)
{   buffer >>= 1; if (bit) buffer |= 0x80;	/* Put bit in top of buffer.*/
    bits_to_go -= 1;
    if (bits_to_go==0) {			/* Output buffer if it is   */
        //putc(buffer,stdout);			/* now full.                */
        fwrite(&buffer,1,1,outputFile);
        n_bits_ += 8;

        bits_to_go = 8;
    }
}


/* FLUSH OUT THE LAST BITS. */
BitOutput::~BitOutput()
{   
	buffer = (buffer>>bits_to_go);
	fwrite(&buffer,1,1,outputFile);
   fclose(outputFile);
}
