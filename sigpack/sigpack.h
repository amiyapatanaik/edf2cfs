// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Version
//  1.01	Claes Rolén		2014-11-30	First version
//  1.02	Claes Rolén		2015-01-11	Added 'angle','specgram','fd_filter','gplot'
//                                      Changed file structure
//  1.03	Claes Rolén		2015-04-26	Added 'parser' class, 'err_handler','wrn_handler'
//                                      'freqz','phasez'
//  1.04	Claes Rolén		2015-08-01	Added FFTW class, 'unwrap', 'update_coeffs' + commenting for Doxygen
//  1.05	Claes Rolén		2015-10-11	Added plot to file in gplot 
//  1.06	Claes Rolén		2015-12-30	Added support for importing/exporting Wisdom in FFTW 


#ifndef ARMA_INCLUDES
#include <armadillo>
#endif

#ifndef SIGPACK_H
#define SIGPACK_H
#include "base/base.h"
#include "window/window.h"
#include "filter/filter.h"
#include "resampling/resampling.h"
#include "spectrum/spectrum.h"
#include "timing/timing.h"
//#include "gplot/gplot.h"
#include "parser/parser.h"
#include "fftw/fftw.h"
#endif

/// \mainpage notitle
///
/// \section intro_sec General
/// \tableofcontents
/// SigPack is a C++ signal processing library using the Armadillo library as a base. 
/// The API will be familiar for those who has used IT++ and Octave/Matlab. The intention 
/// is to keep it small and only implement the fundamental signal processing algorithms.
///
/// \section features_sec Features
/// \li Easy to use, based on Armadillo library
/// \li API similar to Matlab/Octave and IT++
/// \li FIR/IIR filter
/// \li Window functions - Hanning, Hamming, Bartlett, Kaiser ...
/// \li Spectrum and spectrogram
/// \li Timing/Delay
/// \li Gnuplot support
/// \li Up/Downsampling
/// \li Config file parser
/// \li FFTW support
///
/// \section install_sec Installation
/// Download Armadillo and SigPack and install/extract them to your install directory. 
/// Armadillo-4.320.2 version is used in the examples hereafter.
/// For Windows 64bit users: add the \<Armadillo install dir\>\\examples\\lib_win64 
/// to your path in your environment variables.

