#include "stdafx.h"
#include "FFT.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FFT::FFT()
{

}

FFT::~FFT()
{

}

bool FFT::convolve(float* source, float* kernel, int xSource, int ySource, int xKernel, int yKernel)
{
	int x, y, index;

	// get normalization params
	float maxCurrent = 0.0f;
	for (x = 0; x < xSource * ySource; x++)
		maxCurrent = (maxCurrent < source[x]) ? source[x] : maxCurrent;
	float maxKernel = 0.0f;
	for (x = 0; x < xKernel * yKernel; x++)
		maxKernel = (maxKernel < kernel[x]) ? kernel[x] : maxKernel;
	float maxProduct = maxCurrent * maxKernel;

	// retrieve dimensions
	int xHalf = xKernel / 2;
	int yHalf = yKernel / 2;
	int xResPadded = xSource + xKernel;
	int yResPadded = ySource + yKernel;

	if (xResPadded != yResPadded)
		(xResPadded > yResPadded) ? yResPadded = xResPadded : xResPadded = yResPadded;

	// create padded field
	const int ndims=2;
	int dims[ndims] = { xResPadded ,yResPadded };
	kiss_fftnd_cfg st;
	kiss_fft_cpx* padded = (kiss_fft_cpx*)malloc(xResPadded * yResPadded * sizeof(kiss_fft_cpx));
	st = kiss_fftnd_alloc(dims, ndims, 0, NULL, NULL);

	if (!padded)
	{
		cout << " IMAGE: Not enough memory! Try a smaller final image size." << endl;
		return false;
	}

	// init padded field
	for (index = 0; index < xResPadded * yResPadded; index++)
		padded[index].r = padded[index].i = 0.0f;
	index = 0;
	for (y = 0; y < ySource; y++)
		for (x = 0; x < xSource; x++, index++)
		{
			int paddedIndex = (x + xKernel / 2) + (y + yKernel / 2) * xResPadded;
			padded[paddedIndex].r = source[index];
		}

	// create padded filter
	kiss_fft_cpx* filter = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * xResPadded * yResPadded);
	if (!filter)
	{
		cout << " FILTER: Not enough memory! Try a smaller final image size." << endl;
		return false;
	}

	// init padded filter
	#pragma omp parallel for
	for (index = 0; index < xResPadded * yResPadded; index++)
		filter[index].r = filter[index].i = 0.0f;

	// quadrant IV
	#pragma omp parallel for
	for (y = 0; y < (yHalf + 1); y++)
		for (x = 0; x < (xHalf + 1); x++)
		{
			int filterIndex = x + xHalf + y * xKernel;
			int fieldIndex = x + (y + yResPadded - (yHalf + 1)) * xResPadded;
			filter[fieldIndex].r = kernel[filterIndex];
		}

	// quadrant I
	#pragma omp parallel for
	for (y = 0; y < yHalf; y++)
		for (x = 0; x < (xHalf + 1); x++)
		{
			int filterIndex = (x + xHalf) + (y + yHalf + 1) * xKernel;
			int fieldIndex = x + y * xResPadded;
			filter[fieldIndex].r = filter[fieldIndex].i = kernel[filterIndex];
		}

	// quadrant III
	#pragma omp parallel for
	for (y = 0; y < (yHalf + 1); y++)
		for (x = 0; x < xHalf; x++)
		{
			int filterIndex = x + y * xKernel;
			int fieldIndex = (x + xResPadded - xHalf) + (y + yResPadded - (yHalf + 1)) * xResPadded;
			filter[fieldIndex].r = filter[fieldIndex].i = kernel[filterIndex];
		}

	// quadrant II
	#pragma omp parallel for
	for (y = 0; y < yHalf; y++)
		for (x = 0; x < xHalf; x++)
		{
			int filterIndex = x + (y + yHalf + 1) * xKernel;
			int fieldIndex = (x + xResPadded - xHalf) + y * xResPadded;
			filter[fieldIndex].r= filter[fieldIndex].i = kernel[filterIndex];
		}

	// perform forward FFT on field
	kiss_fft_cpx* paddedTransformed = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * xResPadded * yResPadded);
	if (!paddedTransformed)
	{
		cout << " T-IMAGE: Not enough memory! Try a smaller final image size." << endl;
		return false;
	}
	kiss_fftnd(st, padded, paddedTransformed);

	// perform forward FFT on filter
	kiss_fft_cpx* filterTransformed = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * xResPadded * yResPadded);
	if (!filterTransformed)
	{
		cout << " T-FILTER: Not enough memory! Try a smaller final image size." << endl;
		return false;
	}
	kiss_fftnd(st, filter, filterTransformed);

	// apply frequency space filter
	#pragma omp parallel for
	for (index = 0; index < xResPadded * yResPadded; index++)
	{
		float newReal = paddedTransformed[index].r * filterTransformed[index].r -
			paddedTransformed[index].i * filterTransformed[index].i;
		float newIm = paddedTransformed[index].r * filterTransformed[index].i +
			paddedTransformed[index].i * filterTransformed[index].r;
		paddedTransformed[index].r = newReal;
		paddedTransformed[index].i = newIm;
	}

	// transform back
	kiss_fftnd_cfg st2 = kiss_fftnd_alloc(dims, ndims, 1, NULL, NULL);
	kiss_fftnd(st2, paddedTransformed, padded);

	// copy back into padded
	index = 0;
	#pragma omp parallel for
	for (y = 0; y < ySource; y++)
		for (x = 0; x < xSource; x++, index++)
		{
			int paddedIndex = (x + xKernel / 2) + (y + yKernel / 2) * xResPadded;
			source[index] = padded[paddedIndex].r;
		}

	// clean up
	free(padded);
	free(paddedTransformed);
	free(filter);
	free(filterTransformed);

	// if normalization is exceeded, renormalize
	float newMax = 0.0f;
	#pragma omp parallel for
	for (x = 0; x < xSource * ySource; x++)
		newMax = (newMax < source[x]) ? source[x] : newMax;
	if (newMax > maxProduct)
	{
		float scale = maxProduct / newMax;
		for (x = 0; x < xSource * ySource; x++)
			source[x] *= scale;
	}

	return true;
}
