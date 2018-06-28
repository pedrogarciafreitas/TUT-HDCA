# TUT-HDCA

## Introduction

This is the TUT plenoptic image compression software.

## Installation instructions

This software has been developed and tested on Windows 7,10 and has mostly been developed using Visual Studio. However, a makefile for compiling on Linux is provided and it has been tested on Ubuntu 14.04 LTS.

Some of the encoding is currently done by JPEG2000 and for that we use the Kakadu Software

### Kakadu installation on Linux

Download Kakadu for Linux http://kakadusoftware.com/downloads/ and see Kakadu's README.txt for instructions regarding **LD_LIBRARY_PATH**. 

If you encounter the error *"kakadu/kdu_compress: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.21'* not found" do the following,

*sudo add-apt-repository ppa:ubuntu-toolchain-r/test*
 
*sudo apt-get update*

*sudo apt-get install* 

*sudo apt-get install libstdc++6*

## Running the encoder

## Running the decoder

## References

This work is based on academic research and any research based on this software should cite the following papers:

P. Astola, I. Tabus, *Light Field Compression of HDCA Images Combining Linear Prediction and JPEG 2000*, EUSIPCO 2018

I. Tabus, P. Helin, P. Astola, *Lossy Compression of Lenslet Images from Plenoptic Cameras Combining Sparse Predictive Coding and JPEG 2000*, ICIP 2017