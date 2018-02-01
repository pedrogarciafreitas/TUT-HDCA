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

#ifndef BIT_OUTPUT_HH
#define BIT_OUTPUT_HH

#include <vector>
#include <string>
#include <cstdio>

using namespace std;

class BitOutput {

   public:
      BitOutput(std::string filename);
      ~BitOutput();
      void output_bit(int bit);

   private:

      int buffer;		/* Bits buffered for output                 */
      int bits_to_go;		/* Number of bits free in buffer            */
      FILE* outputFile;
      int n_bits_;
};

#endif
