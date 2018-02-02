// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
// Golomb-rice coder.

#ifndef GOLOMB_CODER_HH
#define GOLOMB_CODER_HH 

#include <vector>
#include <cstdint>
#include <string>
#include "bit_output.hh"
#include "bit_input.hh"

#define P_BITS 5

using namespace std;

class GolombCoder {

   public:
      GolombCoder(std::string filename, bool decode);
      ~GolombCoder();

      void encode_symbols(vector< int >& symbols, int N_bits);
      void decode_symbols(vector< int >& symbols, int N_bits);
   private:

      bool decode_;
      BitOutput* bitoutputter;
      BitInput* bitinputter;
      int bit_counter;


};

#endif