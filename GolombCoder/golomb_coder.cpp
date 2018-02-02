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

#include "golomb_coder.hh"
#include <cstdint>
#include <iostream>
#include <climits>

GolombCoder::GolombCoder(std::string filename, bool decode): decode_(decode), bit_counter(0)
{   
  if(decode) {
    bitinputter = new BitInput(filename);
  }
  else {
    bitoutputter = new BitOutput(filename);
  }
}

void GolombCoder::decode_symbols(vector< int >& symbols, int N_bits)
{
  unsigned int N = 0;
  for(int i = 0; i < N_bits; ++i ) {
    unsigned int bit = bitinputter->input_bit();
    N += bit << i;
  }

  unsigned int p = 0;
  for (int i = 0; i < P_BITS; ++i) {
    unsigned int bit = bitinputter->input_bit();
    p += bit << i;
  }

  for( int i = 0; i < N; ++i ) {
    unsigned int sign = bitinputter->input_bit();
    unsigned int q = 0;
    while( bitinputter->input_bit() == 0 ) {
      ++q;
    }
    unsigned int r = 0;
    for(int i = 0; i < p; ++i ) {
      unsigned int bit = bitinputter->input_bit();
      r += bit << i;
    }
    int theta = (q << p) + r;
    if( sign == 1 ) {
      theta = -theta;
      //cout << "golomb decoder: " << theta << endl;
    }
    symbols.push_back(theta);
  }
//  return symbol;
}


void GolombCoder::encode_symbols(vector< int >& symbols, int N_bits)
{
  //bit_thetawv = [];
  vector<int> p_opts;
  vector<int> p_crit;
  for( int p = 0; p < 20; ++p) {
    p_opts.push_back(p);
    p_crit.push_back(0); 
  }
  for( int pi = 0; pi < p_opts.size(); ++pi) {
    int p = p_opts[pi];
    for( int i = 0; i < symbols.size(); ++i ) {
      int theta = symbols.at(i);
      unsigned int q = abs(theta) >> p;
      unsigned int r = abs(theta)-(q << p);
      p_crit[pi] += q + p + 2;
    }
	p_crit[pi] += N_bits + P_BITS;
  }

  unsigned int p = 0;
  unsigned int crit_best = UINT_MAX;
  
  for( int pi = 0; pi < p_opts.size(); ++pi) {
    if( p_crit[pi] < crit_best ) {
      crit_best = p_crit[pi];
      p = p_opts[pi];
    }
  }

  bit_counter += crit_best;

  unsigned int N = symbols.size();
  for(int i = 0; i < N_bits; ++i ) {
    bitoutputter->output_bit((N >> i) & 1);
  }

  for(int i = 0; i < P_BITS; ++i ) {
    bitoutputter->output_bit((p >> i) & 1);
  }

  for( int i = 0; i < symbols.size(); ++i ) {
    int theta = symbols.at(i);
    unsigned int q = abs(theta) >> p;
    unsigned int r = abs(theta)-(q << p);
    unsigned int sign = 0;
    if( theta < 0) {
      sign = 1;
    }
    bitoutputter->output_bit(sign);
    for( int i = 0; i < q; ++i ) {
      bitoutputter->output_bit(0);
    }
    bitoutputter->output_bit(1);
    for(int i = 0; i < p; ++i ) {
      bitoutputter->output_bit((r >> i) & 1);
    }
  }
  //bit_thetawv = [bit_thetawv sign zeros(1,q) 1];
  //bit_thetawv = [bit_thetawv de2bi(r,p,'left-msb')];

}

GolombCoder::~GolombCoder() {
  cout << "golomb bytes: " << bit_counter/8 << endl;
  if( decode_ ) {
    delete bitinputter;
  } else {
    delete bitoutputter;
  }
}