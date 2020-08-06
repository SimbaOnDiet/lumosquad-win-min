#pragma once
#ifndef FFT_H
#define FFT_H

#include "stdafx.h"
#include "kiss_fft.h"
#include "kiss_fftnd.h"
#include "_kiss_fft_guts.h"
#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////
/// \brief Wrapper class for FFTW
////////////////////////////////////////////////////////////////////
class FFT
{
public:
	FFT();
	virtual ~FFT();
	/// \brief convolve image and filter using FFTW
	///
	/// \param source       source image
	/// \param kernel       convolution kernel
	/// \param xSource      width of source image
	/// \param ySource      height of source image
	/// \param xKernel      width of kernel
	/// \param yKernel      height yidth of kernel
	///
	/// \return Returns the convolved image in the 'image' array. If the convolve fails, returns false
	static bool convolve(float* source, float* kernel, int xSource, int ySource, int xKernel, int yKernel);
};

#endif
