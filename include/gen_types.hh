// Authors: Ioan Tabus, Petri Helin, Pekka Astola
// Copyright(c) Statistical Inference and Signal Compression (SISC),
//             Tampere University of Technology (TUT)
//             http://www.cs.tut.fi/~tabus/
// All rights reserved.
// Using of this module for developing the JPEG-PLENO standard is allowed
// Using of this module in producing any publication is allowed only if reference to this original publication is inserted:
// Ioan Tabus, Petri Helin, Pekka Astola, "Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000", ICIP 2017, Beijing, China, 17-20 September 2017.17
//
//
// This header contains constant definitions and auxiliary functions used in many places in the lightfield compressor.

#ifndef GEN_TYPES_HH
#define GEN_TYPES_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <math.h>


#define NBIT_GR 32

const bool verbose = false;

// Replacing the variables irMat01p_const, icMat01p_const and Neigh9_165p_const changes the processing order of the views.
// The directory aux_matlab contains a script "create_spiral.m" for this purpose.
const int irMat01p_const[] = { 0, 0, -1, -1, -1, 0, 1, 1, 1, 1, 0, -1, -2, -2, -2, -2, -2, -1, 0, 1, 2, 2, 2, 2, 2, 2, 1, 0, -1, -2, -3, -3, -3, -3, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 0, -1, -2, -3, -4, -4, -4, -4, -4, -4, -4, -4, -4, -3, -2, -1, 0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 };
const int icMat01p_const[] = { 0, -1, -1, 0, 1, 1, 1, 0, -1, -2, -2, -2, -2, -1, 0, 1, 2, 2, 2, 2, 2, 1, 0, -1, -2, -3, -3, -3, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 1, 0, -1, -2, -3, -4, -4, -4, -4, -4, -4, -4, -4, -3, -2, -1, 0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6 };
const int Neigh9_165p_const[] = { 3, 12, 13, 14, 15, 4, 1, 2, 11, 28, 29, 30, 31, 32, 33, 34, 35, 16, 5, 6, 7, 8, 9, 10, 27, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 36, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 51, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 64, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 83, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 100, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 123, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 133, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 144, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 122, 4, 3, 14, 15, 16, 5, 6, 1, 2, 11, 12, 13, 32, 33, 34, 35, 36, 17, 18, 19, 20, 7, 8, 9, 10, 27, 28, 29, 30, 31, 58, 59, 60, 61, 62, 63, 64, 37, 38, 39, 40, 41, 42, 21, 22, 23, 24, 25, 26, 51, 52, 53, 54, 55, 56, 57, 92, 93, 94, 95, 96, 97, 98, 99, 100, 65, 66, 67, 68, 69, 70, 71, 72, 43, 44, 45, 46, 47, 48, 49, 50, 83, 84, 85, 86, 87, 88, 89, 90, 91, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 134, 134, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 121, 5, 4, 15, 16, 17, 18, 19, 6, 1, 2, 3, 14, 33, 34, 35, 36, 37, 38, 39, 40, 41, 20, 7, 8, 9, 10, 11, 12, 13, 32, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 42, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 58, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 72, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 92, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 110, 73, 74, 75, 76, 77, 78, 79, 80, 81, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 132, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 101, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 158, 156, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 169, 2, 11, 12, 3, 4, 1, 8, 9, 10, 27, 28, 29, 30, 13, 14, 15, 16, 5, 6, 7, 22, 23, 24, 25, 26, 51, 52, 53, 54, 55, 56, 31, 32, 33, 34, 35, 36, 17, 18, 19, 20, 21, 44, 45, 46, 47, 48, 49, 50, 83, 84, 85, 86, 87, 88, 89, 90, 57, 58, 59, 60, 61, 62, 63, 64, 37, 38, 39, 40, 41, 42, 43, 74, 75, 76, 77, 78, 79, 80, 81, 82, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 65, 66, 67, 68, 69, 70, 71, 72, 73, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 121, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 91, 132, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 146, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 168, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 169, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 0, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 0, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 158, 0, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 0, 6, 1, 4, 5, 18, 19, 20, 7, 8, 9, 2, 3, 14, 15, 16, 17, 38, 39, 40, 41, 42, 21, 22, 23, 24, 25, 10, 11, 12, 13, 32, 33, 34, 35, 36, 37, 66, 67, 68, 69, 70, 71, 72, 43, 44, 45, 46, 47, 48, 49, 26, 27, 28, 29, 30, 31, 58, 59, 60, 61, 62, 63, 64, 65, 102, 103, 104, 105, 106, 107, 108, 109, 110, 73, 74, 75, 76, 77, 78, 79, 80, 81, 50, 51, 52, 53, 54, 55, 56, 57, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 168, 121, 82, 83, 84, 85, 86, 87, 88, 89, 90, 0, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 146, 0, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 0, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 0, 9, 10, 11, 2, 1, 8, 23, 24, 25, 26, 27, 28, 29, 12, 3, 4, 5, 6, 7, 22, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 30, 13, 14, 15, 16, 17, 18, 19, 20, 21, 44, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 56, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 74, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 90, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 112, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 9, 2, 1, 6, 7, 22, 23, 24, 25, 10, 11, 12, 3, 4, 5, 18, 19, 20, 21, 44, 45, 46, 47, 48, 49, 26, 27, 28, 29, 30, 13, 14, 15, 16, 17, 38, 39, 40, 41, 42, 43, 74, 75, 76, 77, 78, 79, 80, 81, 50, 51, 52, 53, 54, 55, 56, 31, 32, 33, 34, 35, 36, 37, 66, 67, 68, 69, 70, 71, 72, 73, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 82, 83, 84, 85, 86, 87, 88, 89, 90, 57, 58, 59, 60, 61, 62, 63, 64, 65, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 8, 1, 6, 19, 20, 21, 22, 23, 24, 9, 2, 3, 4, 5, 18, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 25, 10, 11, 12, 13, 14, 15, 16, 17, 38, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 49, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 66, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 81, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 102, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

const int hexag_even_C_const[] = { -1, 0, 1, 2, -4, -3, -2, -1, 0, 1, 2, 3, 4, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, -4, -3, -2, -1, 0, 1, 2, 3, 4, -1, 0, 1, 2 };
const int hexag_odd_C_const[] = { -2, -1, 0, 1, -4, -3, -2, -1, 0, 1, 2, 3, 4, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, -4, -3, -2, -1, 0, 1, 2, 3, 4, -2, -1, 0, 1 };
const int hexag_even_R_const[] = { -9, -9, -9, -9, -8, -8, -8, -8, -8, -8, -8, -8, -8, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9 };

// number of pixels in the original lenslet image.
const int num_pixels_in_lenslet = 41483904;

//unsigned short getMedian(std::vector<unsigned short> scores)
//{
//	unsigned short median;
//	size_t size = scores.size();
//
//	std::sort(scores.begin(), scores.end());
//
//	if (size % 2 == 0)
//	{
//		median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
//	}
//	else
//	{
//		median = scores[size / 2];
//	}
//
//	return median;
//}
//
//
//float getMedian(std::vector<float> scores)
//{
//	float median;
//	size_t size = scores.size();
//
//	std::sort(scores.begin(), scores.end());
//
//	if (size % 2 == 0)
//	{
//		median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
//	}
//	else
//	{
//		median = scores[size / 2];
//	}
//
//	return median;
//}

//void medfilt2D(int* input, int* output, int SZ, int nr, int nc)
//{
//	int dsz = floor(SZ / 2);
//	std::vector<int> scores;
//	for (int y = 0; y < nr; y++){
//		for (int x = 0; x < nc; x++){
//			scores.clear();
//			for (int dy = -dsz; dy < dsz; dy++){
//				for (int dx = -dsz; dx < dsz; dx++){
//					if ((y + dy) >= 0 && (y + dy) < nr
//						&& (x + dx) >= 0 && (x + dx) < nc)
//						scores.push_back(input[y + dy + (x + dx)*nr]);
//				}
//			}
//			output[y + x*nr] = getMedian(scores);
//		}
//	}
//}
//
//void medfilt2D(unsigned short* input, unsigned short* output, int SZ, int nr, int nc)
//{
//	int dsz = floor(SZ / 2);
//	std::vector<unsigned short> scores;
//	for (int y = 0; y < nr; y++){
//		for (int x = 0; x < nc; x++){
//			scores.clear();
//			for (int dy = -dsz; dy < dsz; dy++){
//				for (int dx = -dsz; dx < dsz; dx++){
//					if ((y + dy) >= 0 && (y + dy) < nr
//						&& (x + dx) >= 0 && (x + dx) < nc)
//						scores.push_back(input[y + dy + (x + dx)*nr]);
//				}
//			}
//			output[y + x*nr] = getMedian(scores);
//		}
//	}
//}
//
//
//void medfilt2D(float* input, float* output, int SZ, int nr, int nc)
//{
//	int dsz = floor(SZ / 2);
//	std::vector<float> scores;
//	for (int y = 0; y < nr; y++){
//		for (int x = 0; x < nc; x++){
//			scores.clear();
//			for (int dy = -dsz; dy < dsz; dy++){
//				for (int dx = -dsz; dx < dsz; dx++){
//					if ((y + dy) >= 0 && (y + dy) < nr
//						&& (x + dx) >= 0 && (x + dx) < nc)
//						scores.push_back(input[y + dy + (x + dx)*nr]);
//				}
//			}
//			output[y + x*nr] = getMedian(scores);
//		}
//	}
//}

// write a metadata file about sizes of images
inline void aux_write_header_file(int nr, int nc, int nvr, int nvc, const char* filename) {
	FILE* f_main_header = fopen(filename, "wb");
	fwrite(&nr, sizeof(int), 1, f_main_header);
	fwrite(&nc, sizeof(int), 1, f_main_header);
	fwrite(&nvr, sizeof(int), 1, f_main_header);
	fwrite(&nvc, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

inline void aux_read_header_file(int* nr, int* nc, int* nvr, int* nvc, const char* filename) {
	FILE* f_main_header = fopen(filename, "rb");
	fread(nr, sizeof(int), 1, f_main_header);
	fread(nc, sizeof(int), 1, f_main_header);
	fread(nvr, sizeof(int), 1, f_main_header);
	fread(nvc, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

inline void aux_read_file_float(const int nr, const int nc, const int ncomponents, const char* filename, float *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(float), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_read_file_uint32(const int nr, const int nc, const int ncomponents, const char* filename, unsigned int *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(unsigned int), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_read_file_int32(const int nr, const int nc, const int ncomponents, const char* filename, int *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(int), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_read_file_uint8(const int nr, const int nc, const int ncomponents, const char* filename, unsigned char *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(unsigned char), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_read_file_uint16(const int nr, const int nc, const int ncomponents, const char* filename, unsigned short *data) {

	FILE* f_file = fopen(filename, "rb");

	fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(unsigned short), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

inline void aux_readDSP_uint16(const int nr, const int nc, const int ncomponents, const char* filename, unsigned short *data) {

	FILE* f_file = fopen(filename, "rb");

	//fseek(f_file, 4 * sizeof(int), SEEK_SET); // skip header

	fread(data, sizeof(unsigned short), nr*nc*ncomponents, f_file);

	fclose(f_file);

}

// write a metadata file about sizes of images
inline void aux_write_header(int nr, int nc, int nvr, int nvc) {
	FILE* f_main_header = fopen("HDR", "wb");
	fwrite(&nr, sizeof(int), 1, f_main_header);
	fwrite(&nc, sizeof(int), 1, f_main_header);
	fwrite(&nvr, sizeof(int), 1, f_main_header);
	fwrite(&nvc, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

// read the metadata file
inline void aux_read_header(int* nr, int* nc, int* nvr, int* nvc) {
	FILE* f_main_header = fopen("HDR", "rb");
	fread(nr, sizeof(int), 1, f_main_header);
	fread(nc, sizeof(int), 1, f_main_header);
	fread(nvr, sizeof(int), 1, f_main_header);
	fread(nvc, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

// write a header about prediction order and maximum region number
inline void aux_write_pred_header(int Ms, int maxiS) {
	FILE* f_main_header = fopen("HDRP", "wb");
	fwrite(&Ms, sizeof(int), 1, f_main_header);
	fwrite(&maxiS, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}

// read the header about prediction order and maximum region number
inline void aux_read_pred_header(int* Ms, int* maxiS) {
	FILE* f_main_header = fopen("HDRP", "rb");
	fread(Ms, sizeof(int), 1, f_main_header);
	fread(maxiS, sizeof(int), 1, f_main_header);
	fclose(f_main_header);
}


// read an 8 bit ppm file (here only one channel)
inline void aux_read8ppm(FILE *filept, int width, int height, int *img){
	int  max, x, y;
	int red;
	char dummy[100];

	unsigned char *Image8bit = NULL;

	/*--< Read header information of 16bit ppm image from filept >--*/
	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);
	printf("%s %d %d %d\n", dummy, max, width, height);


	Image8bit = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image8bit, sizeof(unsigned char), width*height * 3, filept);

	int i = 0;
	/*--< Find maximum value in the image >--*/
	for (x = 0; x < width; x++){
		for (y = 0; y < height; y++){
			red = (int)Image8bit[(x + y*width) * 3];
			img[i] = red;
			i++;
		}
	}

	free(Image8bit);

}

// read a 16 bit color pgm image
inline void aux_read16pgm_1080p(FILE *filept, int *img)
{
	int  max, x, y;
	int red, green, blue, pixelmax;
	char dummy[100];

	bool divd = 0;

	unsigned short int *Image16bit = NULL;

	int width, height;

	//FILE* filept = fopen("007_007.ppm", "r");
	/*--< Read header information of 16bit ppm image from filept >--*/
	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);
	//printf("%s\n", dummy);

	/*--< Really 16bit ppm? Check it >--*/
	if (strncmp(dummy, "P5", 2) != 0) {
		fprintf(stderr, "Error: The input data is not PGM.\n");
		exit(1);
	}
	//if (max == 65535){
	//	divd = 1;
	//}

	Image16bit = (unsigned short int *)malloc(width*height* sizeof(unsigned short int));

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image16bit, sizeof(unsigned short int), width*height, filept);

	int i = 0;
	/*--< Find maximum value in the image >--*/
	pixelmax = 0;
	/* UNSW inverse depth already in 4K resolution, just crop to 1080p */
	for (x = 960; x < 960 + 1920; x++){
		for (y = 540; y < 540 + 1080; y++){
			//if (y>209 + 540 && x > 85 + 960 && y < 209 + 540 + 1080 && x < 85 + 960 + 1920)
				{
					red = Image16bit[(x + y*width)];
					//green = Image16bit[(x + y*width) * 3 + 1];
					//blue = Image16bit[(x + y*width) * 3 + 2];

					// Exhange upper 8bit and lower 8bit for Intel x86
					red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
					//green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
					//blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);

					if (pixelmax < red) pixelmax = red;
					//if (pixelmax < green) pixelmax = green;
					//if (pixelmax < blue) pixelmax = blue;


					//if (divd) // fix for 16bit to 10bit
					//{
					//	red = red >> 6;
					//	//green = green >> 6;
					//	//blue = blue >> 6;
					//}

					img[i] = red;
					//img[i + height*width] = green;
					//img[i + 2 * height*width] = blue;

					if (0){ //debug stuff
						fprintf(stderr, "sample %i\n", img[i]);
						/*fprintf(stderr, "sample %i\t%i\t%i\n", img[i] >> 6, img[i + height*width] >> 6, img[i + 2 * height*width] >> 6);*/
					}

					i++;

				}
		}
	}
	free(Image16bit);
}

inline void aux_read16ppm(FILE *filept, int width, int height, unsigned short *img)
{
	int  max, x, y;
	int red, green, blue, pixelmax;
	char dummy[100];

	unsigned short int *Image16bit = NULL;

	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);

	if (strncmp(dummy, "P6", 2) != 0) {
		fprintf(stderr, "Error: The input data is not binary PPM.\n");
		exit(1);
	}

	Image16bit = (unsigned short int *)malloc(width*height * 3 * sizeof(unsigned short int));

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image16bit, sizeof(unsigned short int), width*height * 3, filept);

	int i = 0;
	/*--< Find maximum value in the image >--*/
	pixelmax = 0;
	for (x = 0; x < width; x++){
		for (y = 0; y < height; y++){
			red = Image16bit[(x + y*width) * 3];
			green = Image16bit[(x + y*width) * 3 + 1];
			blue = Image16bit[(x + y*width) * 3 + 2];

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
			green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
			blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);

			if (pixelmax < red) pixelmax = red;
			if (pixelmax < green) pixelmax = green;
			if (pixelmax < blue) pixelmax = blue;


			img[i] = red;
			img[i + height*width] = green;
			img[i + 2 * height*width] = blue;

			if (0){ //debug stuff
				fprintf(stderr, "sample %i\t%i\t%i\n", img[i], img[i + height*width], img[i + 2 * height*width]);
				fprintf(stderr, "sample %i\t%i\t%i\n", img[i] >> 6, img[i + height*width] >> 6, img[i + 2 * height*width] >> 6);
			}

			i++;


		}
	}

	free(Image16bit);

}

// read a 16 bit color ppm image
inline void aux_read16ppm(FILE *filept, int width, int height, int *img)
{
	int  max, x, y;
	int red, green, blue, pixelmax;
	char dummy[100];

	bool divd = 0;

	unsigned short int *Image16bit = NULL;

	//FILE* filept = fopen("007_007.ppm", "r");
	/*--< Read header information of 16bit ppm image from filept >--*/
	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);
	//printf("%s\n", dummy);

	/*--< Really 16bit ppm? Check it >--*/
	if (strncmp(dummy, "P6", 2) != 0) {
		fprintf(stderr, "Error: The input data is not binary PPM.\n");
		exit(1);
	}
	if (max == 65535){
		//fprintf(stderr, "Warning: The input data is 16bit PPM.\n");
		divd = 1;
		//exit(1);
	}
	//if (max != 1023){
	//    fprintf(stderr, "Error: The input data is not 10bit PPM.\n");
	//    exit(1);
	//}

	Image16bit = (unsigned short int *)malloc(width*height * 3 * sizeof(unsigned short int));

	/*--< Read 16bit ppm image from filept >--*/
	fread(Image16bit, sizeof(unsigned short int), width*height * 3, filept);

	int i = 0;
	/*--< Find maximum value in the image >--*/
	pixelmax = 0;
	for (x = 0; x < width; x++){
		for (y = 0; y < height; y++){
			red = Image16bit[(x + y*width) * 3];
			green = Image16bit[(x + y*width) * 3 + 1];
			blue = Image16bit[(x + y*width) * 3 + 2];

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
			green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
			blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);

			if (pixelmax < red) pixelmax = red;
			if (pixelmax < green) pixelmax = green;
			if (pixelmax < blue) pixelmax = blue;


			if (divd) // fix for 16bit to 10bit
			{
				red = red >> 6;
				green = green >> 6;
				blue = blue >> 6;
			}

			img[i] = red;
			img[i + height*width] = green;
			img[i + 2 * height*width] = blue;

			if (0){ //debug stuff
				fprintf(stderr, "sample %i\t%i\t%i\n", img[i], img[i + height*width], img[i + 2 * height*width]);
				fprintf(stderr, "sample %i\t%i\t%i\n", img[i] >> 6, img[i + height*width] >> 6, img[i + 2 * height*width] >> 6);
			}

			i++;


		}
	}

	free(Image16bit);

}

// read a 16 bit color ppm image
//inline void aux_read16ppm(FILE *filept, int width, int height, int *img){
//    int  max, x, y;
//    int red, green, blue, pixelmax;
//    char dummy[100];
//    
//    unsigned short int *Image16bit=NULL;
//    
//    //FILE* filept = fopen("007_007.ppm", "r");
//    /*--< Read header information of 16bit ppm image from filept >--*/
//    fscanf(filept, "%s", dummy);
//    fscanf(filept, "%d %d\n", &width, &height);
//    fscanf(filept, "%d", &max);
//    fgetc(filept);
//    //printf("%s\n", dummy);
//    
//    /*--< Really 16bit ppm? Check it >--*/
//    if (strncmp(dummy, "P6", 2) != 0) {
//        fprintf(stderr, "Error: The input data is not binary PPM.\n");
//        exit(1);
//    }
//    if (max != 1023){
//        fprintf(stderr, "Error: The input data is not 10bit PPM.\n");
//        exit(1);
//    }
//    
//    Image16bit = (unsigned short int *)malloc(width*height*3*sizeof(unsigned short int));
//    
//    /*--< Read 16bit ppm image from filept >--*/
//    fread(Image16bit, sizeof(unsigned short int), width*height*3, filept);
//    
//    int i = 0;
//    /*--< Find maximum value in the image >--*/
//    pixelmax=0;
//    for(x=0;x<width;x++){
//        for(y=0;y<height;y++){
//            red  = Image16bit[(x+y*width)*3  ];
//            green= Image16bit[(x+y*width)*3+1];
//            blue = Image16bit[(x+y*width)*3+2];
//            
//            // Exhange upper 8bit and lower 8bit for Intel x86
//            red   = ((red   & 0x00ff)<<8)|((red   & 0xff00)>>8);
//            green = ((green & 0x00ff)<<8)|((green & 0xff00)>>8);
//            blue  = ((blue  & 0x00ff)<<8)|((blue  & 0xff00)>>8);
//            
//            if (pixelmax < red  ) pixelmax=red;
//            if (pixelmax < green) pixelmax=green;
//            if (pixelmax < blue ) pixelmax=blue;
//            
//            img[i] = red;
//            img[i+height*width] = green;
//            img[i+2*height*width] = blue;
//            i++;
//        }
//    }
//    
//    free(Image16bit);
//    
//}

template <class T>
void aux_write16pgm(FILE *fp, int width, int height, T *img)
{
	unsigned char *p;
	int i, tmp, j;

	unsigned short* img16bit = (unsigned short*)malloc(height*width*sizeof(unsigned short));
	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = (unsigned short)img[j + i*height];
			lin_ind = lin_ind + 1;
		}
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++){
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
	}

	fprintf(fp, "P5\n%d %d\n65535\n", width, height);
	fwrite(img16bit, sizeof(unsigned short int), width*height, fp);
	free(img16bit);
}

inline void aux_write16ppm_16(FILE *fp, int width, int height, unsigned short int *img){
	unsigned char *p;
	int i, tmp, j;

	//unsigned short int maxi = 0;

	unsigned short* img16bit = (unsigned short*)malloc(height*width * 3 * sizeof(unsigned short));
	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = img[j + i*height];
			img16bit[lin_ind + 1] = img[j + i*height + width*height];
			img16bit[lin_ind + 2] = img[j + i*height + 2 * width*height];
			lin_ind = lin_ind + 3;
		}
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++){
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
	}
	fprintf(fp, "P6\n%d %d\n65535\n", width, height);
	fwrite(img16bit, sizeof(unsigned short int), width*height * 3, fp);
	free(img16bit);
}

// write a 16 bit color ppm image
inline void aux_write16ppm(FILE *fp, int width, int height, unsigned short int *img){
	unsigned char *p;
	int i, tmp, j;

	unsigned short* img16bit = (unsigned short*)malloc(height*width * 3 * sizeof(unsigned short));
	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = img[j + i*height];
			img16bit[lin_ind + 1] = img[j + i*height + width*height];
			img16bit[lin_ind + 2] = img[j + i*height + 2 * width*height];
			lin_ind = lin_ind + 3;
		}
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++){
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
	}
	fprintf(fp, "P6\n%d %d\n1023\n", width, height);
	fwrite(img16bit, sizeof(unsigned short int), width*height * 3, fp);
	free(img16bit);
}

// This function transforms 3D lenslet image with size e.g. 13*434, 13*646, 3
// to an LFgamma structure with indices LFgamma[iview, ir, ic, icomp] where iview is the linear index of the view.
// The order is determined by the earlier introduced variables irMat01, icMat01 and Neigh9_165.
// LIremap can be zero.
inline void aux_Construct_LFgamma(int LI[], double LIremap[], int LFgamma[])
{
	// Local variables for CC_Construct_LFgamma
	int CMI_L1, CMI_L2, MIr, MIc, icglob, irglob, ig, in, in_lens_ok, inLgam, icomp, ival;
	int nr = 0;
	int nc = 0;
	int nvr = 0;
	int nvc = 0;
	aux_read_header(&nr, &nc, &nvr, &nvc);
	int nviews = nvr*nvc;
	int nrr = nr*nvr;
	int ncc = nc*nvc;

	int* GridCoordsYrows = (int*)malloc(nr*nc*sizeof(int));
	int* GridCoordsXcols = (int*)malloc(nr*nc*sizeof(int));
	for (int i = 0; i < nr; ++i) {
		for (int j = 0; j < nc; ++j) {
			GridCoordsXcols[i + j*nr] = 7 + j * 13;
			GridCoordsYrows[i + j*nr] = 7 + i * 13;
		}
	}

	for (MIr = 0; MIr < nr; MIr++)
		for (MIc = 0; MIc < nc; MIc++)
		{
		in = MIr + MIc*nr;
		CMI_L1 = GridCoordsXcols[in] - 1;
		CMI_L2 = GridCoordsYrows[in] - 1;
		// Check that centers are inside lenslett
		if ((CMI_L1 >= 0) && (CMI_L1 < ncc) && (CMI_L2 >= 0) && (CMI_L2 < nrr))
		{
			// all indices in the current MI
			for (ig = 0; ig < nviews; ig++)  //nviews
			{
				irglob = CMI_L2 + irMat01p_const[ig];
				icglob = CMI_L1 + icMat01p_const[ig];
				if ((icglob >= 0) && (icglob < ncc) && (irglob >= 0) && (irglob < nrr))
				{
					// in_lens_ok = sub2ind([nrr,ncc],irglob(inok),icglob(inok));
					for (icomp = 0; icomp < 3; icomp++)
					{
						// linear index LFgamma(ig,MIr,MIc,icomp) :  ((icomp*nc+MIc)*nr+MIr)*nviews+ig
						// linear index LI(irglob,icglob,icomp)  : (icomp*ncc+icglob)*nrr+irglob
						//LFgamma[ig+1, MIr+1, MIc+1, icomp+1] = LI[in_lens_ok+icomp];

						ival = LI[(icomp*ncc + icglob)*nrr + irglob];
						if (ival > 1023)
							ival = 1023;
						if (ival < 0)
							ival = 0;
						if (LIremap != 0) {
							int irM = irMat01p_const[ig] + 4;
							int icM = icMat01p_const[ig] + 4;
							//                                if( irM >=0 && icM >=0 && irM < 9 && icM < 9 ) {
							//                                    double gamma_ival = pow(((double)ival), 1.0/2.2);
							//                                    LIremap[(MIr*9+irM) + (MIc*9+icM)*9*nr + icomp*(9*nr*9*nc)] = gamma_ival;
							//                                }
						}

						LFgamma[((icomp*nc + MIc)*nr + MIr)*nviews + ig] = ival;
					}
				}
			}
		}
		}

	return;
}

// Function to go from LFgamma-structure to LI-structure.
// reconsLI is the output and LFgamma is the input.
// hatLI is used to initialize reconsLI. It can be zeros in the rectified case.
inline void aux_Inverse_LFgamma_TO_LI(int reconsLI[], int LFgamma[], int hatLI[])
{
	// Local variables same as in CC_Construct_LFgamma
	int CMI_L1, CMI_L2, MIr, MIc, icglob, irglob, ig, in, in_lens_ok, inLgam, icomp, i, ival;
	// read header
	int nr = 0;
	int nc = 0;
	int nvr = 0;
	int nvc = 0;
	aux_read_header(&nr, &nc, &nvr, &nvc);
	int nviews = nvr*nvc;
	int nrr = nr*nvr;
	int ncc = nc*nvc;
	int ncols = 3;

	int* GridCoordsYrows = (int*)malloc(nr*nc*sizeof(int));
	int* GridCoordsXcols = (int*)malloc(nr*nc*sizeof(int));
	for (int i = 0; i < nr; ++i) {
		for (int j = 0; j < nc; ++j) {
			GridCoordsXcols[i + j*nr] = 7 + j * 13;
			GridCoordsYrows[i + j*nr] = 7 + i * 13;
		}
	}


	// initialize reconsLI to the base JP2 hatLI
	for (i = 0; i < nrr*ncc * 3; i++)
		reconsLI[i] = hatLI[i];

	for (MIr = 0; MIr < nr; MIr++)
		for (MIc = 0; MIc < nc; MIc++)
		{
		in = MIr + MIc*nr;
		CMI_L1 = GridCoordsXcols[in] - 1;
		CMI_L2 = GridCoordsYrows[in] - 1;
		// Check that centers are inside lenslett
		if ((CMI_L1 >= 0) && (CMI_L1 < ncc) && (CMI_L2 >= 0) && (CMI_L2 < nrr))
		{
			// all indices in the current MI
			for (ig = 0; ig < nviews; ig++)  //ds.nviews
			{
				irglob = CMI_L2 + irMat01p_const[ig];
				icglob = CMI_L1 + icMat01p_const[ig];
				if ((icglob >= 0) && (icglob < ncc) && (irglob >= 0) && (irglob < nrr))
				{
					// in_lens_ok = sub2ind([nrr,ncc],irglob(inok),icglob(inok));
					for (icomp = 0; icomp < 3; icomp++)
					{
						// linear index LFgamma(ig,MIr,MIc,icomp) :  ((icomp*nc+MIc)*nr+MIr)*ds.nviews+ig
						// linear index LI(irglob,icglob,icomp)  : (icomp*ncc+icglob)*nrr+irglob
						//LFgamma[ig+1, MIr+1, MIc+1, icomp+1] = LI[in_lens_ok+icomp];
						ival = LFgamma[((icomp*nc + MIc)*nr + MIr)*nviews + ig];
						if (ival > 1023)
							ival = 1023;
						if (ival < 0)
							ival = 0;

						if (ival != 0) {
							reconsLI[(icomp*ncc + icglob)*nrr + irglob] = ival;
						}


						//							reconsLI[(icomp*ncc+icglob)*nrr+irglob] = LFgamma[((icomp*nc+MIc)*nr+MIr)*ds.nviews+ig];
					}
				}
			}
		}
		}
}

// Function to go from LFgamma-structure to LI-structure.
// reconsLI is the output and LFgamma is the input.
// hatLI is used to initialize reconsLI. It can be zeros in the rectified case.
inline void aux_Inverse_LFgamma_TO_LI(unsigned short int reconsLI[], int LFgamma[], int hatLI[])
{
	// Local variables same as in CC_Construct_LFgamma
	int CMI_L1, CMI_L2, MIr, MIc, icglob, irglob, ig, in, in_lens_ok, inLgam, icomp, i, ival;
	// read header
	int nr = 0;
	int nc = 0;
	int nvr = 0;
	int nvc = 0;
	aux_read_header(&nr, &nc, &nvr, &nvc);
	int nviews = nvr*nvc;
	int nrr = nr*nvr;
	int ncc = nc*nvc;
	int ncols = 3;

	int* GridCoordsYrows = (int*)malloc(nr*nc*sizeof(int));
	int* GridCoordsXcols = (int*)malloc(nr*nc*sizeof(int));
	for (int i = 0; i < nr; ++i) {
		for (int j = 0; j < nc; ++j) {
			GridCoordsXcols[i + j*nr] = 7 + j * 13;
			GridCoordsYrows[i + j*nr] = 7 + i * 13;
		}
	}


	// initialize reconsLI to the base JP2 hatLI
	for (i = 0; i < nrr*ncc * 3; i++)
		reconsLI[i] = hatLI[i];

	for (MIr = 0; MIr < nr; MIr++)
		for (MIc = 0; MIc < nc; MIc++)
		{
		in = MIr + MIc*nr;
		CMI_L1 = GridCoordsXcols[in] - 1;
		CMI_L2 = GridCoordsYrows[in] - 1;
		// Check that centers are inside lenslett
		if ((CMI_L1 >= 0) && (CMI_L1 < ncc) && (CMI_L2 >= 0) && (CMI_L2 < nrr))
		{
			// all indices in the current MI
			for (ig = 0; ig < nviews; ig++)  //ds.nviews
			{
				irglob = CMI_L2 + irMat01p_const[ig];
				icglob = CMI_L1 + icMat01p_const[ig];
				if ((icglob >= 0) && (icglob < ncc) && (irglob >= 0) && (irglob < nrr))
				{
					// in_lens_ok = sub2ind([nrr,ncc],irglob(inok),icglob(inok));
					for (icomp = 0; icomp < 3; icomp++)
					{
						// linear index LFgamma(ig,MIr,MIc,icomp) :  ((icomp*nc+MIc)*nr+MIr)*ds.nviews+ig
						// linear index LI(irglob,icglob,icomp)  : (icomp*ncc+icglob)*nrr+irglob
						//LFgamma[ig+1, MIr+1, MIc+1, icomp+1] = LI[in_lens_ok+icomp];
						ival = LFgamma[((icomp*nc + MIc)*nr + MIr)*nviews + ig];
						if (ival > 1023)
							ival = 1023;
						if (ival < 0)
							ival = 0;

						if (ival != 0) {
							reconsLI[(icomp*ncc + icglob)*nrr + irglob] = ival;
						}


						//							reconsLI[(icomp*ncc+icglob)*nrr+irglob] = LFgamma[((icomp*nc+MIc)*nr+MIr)*ds.nviews+ig];
					}
				}
			}
		}
		}
}

// allocates an int vector of size m
inline int* alocaVector(int m)
{
	int *vector, i;
	vector = (int*)malloc(m*sizeof(int));
	for (i = 0; i < m; i++)
	{
		vector[i] = 0;
	}
	return vector;
}

// allocates a double vector of size m
inline double* alocaDoubleVector(int m)
{
	double *vector;
	int i;
	vector = (double*)malloc(m*sizeof(double));
	for (i = 0; i < m; i++)
	{
		vector[i] = 0;
	}
	return vector;
}

inline long aux_GetFileSize(char* filename)
{
	struct stat stat_buf;
	int rc = stat(filename, &stat_buf);
	return rc == 0 ? stat_buf.st_size : -1;
}
#endif
