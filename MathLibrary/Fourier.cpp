
#include <cmath>
#include "Fourier.h"

using namespace MATH;

void Fourier::fft(float* data, const int length, const DIRECTION direction){
	const float PI = 3.141592653589793238462643383279502884197f;
	int i, j;
	int n, mmax, m, istep;
	float tempreal, tempimg, wreal, wtemp, wpreal, wpimg, wimg;
	float theta, sign;

	sign = float(direction);
	int nn = length / 2;/// As mentioned above nn is the length of the data sample
	n = nn << 1;		/// Cute trick - too sexy for me

	/// This part is called the "bit-reversal" - it does the exchange of two complex numbers
	j = 1;
	for (i = 1; i < n; i += 2) {  
		if (j > i) {
			wtemp = data[j - 1], data[j - 1] = data[i - 1], data[i - 1] = wtemp;///swap (data[j - 1], data[i - 1]);
			wtemp = data[j], data[j] = data[i], data[i] = wtemp;				///swap (data[j], data[i]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	/// This is the Danielson-Lanczos part of the FFT process
	mmax = 2;
	while (n > mmax) {
		istep = mmax << 1;
		
		theta = sign*(2.0f*PI / mmax);
		wtemp = sinf(0.5f*theta);
		wpreal = -2.0f*wtemp*wtemp;
		wpimg = sinf(theta);
		wreal = 1.0f;
		wimg = 0.0f;
		for (m = 1; m<mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempreal = wreal*data[j - 1] - wimg*data[j];
				tempimg = wreal*data[j] + wimg*data[j - 1];
				data[j - 1] = data[i - 1] - tempreal;
				data[j] = data[i] - tempimg;
				data[i - 1] += tempreal;
				data[i] += tempimg;
			}
			wreal = (wtemp = wreal)*wpreal - wimg*wpimg + wreal;
			wimg = wimg*wpreal + wtemp*wpimg + wimg;
		}
		mmax = istep;
	}


} /// End of namespace MATH
