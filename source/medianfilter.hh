#ifndef MEDIANFILTER_HH
#define MEDIANFILTER_HH

#include <vector>

template <class T>
T getMedian(std::vector<T> scores);

template <class T>
void medfilt2D(T* input, T* output, int SZ, int nr, int nc);

#endif