// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Simple bit-input.

#include <iostream>
#include "bit_input.hh"

BitInput::BitInput(std::string filename): buffer(0), bits_to_go(0)
{   
   outputFile = fopen(filename.c_str(), "rb");
}


BitInput::~BitInput()
{   
}

int BitInput::input_bit()
{   int t;
    bits_to_go -= 1;
    if (bits_to_go<0) {				/* Read the next byte if no */
        if(fread(&buffer,1,1,outputFile) != 1) {
          //printf("Can not read!\n");
        }
		//buffer = getc(stdin);			/* bits are left in the     */
        bits_to_go = 7;				/* buffer. Return anything  */
    }			   			/* after end-of-file.       */
    t = buffer&1;
    buffer >>= 1;				/* Return the next bit from */
    return t;					/* the bottom of the byte.  */
}
