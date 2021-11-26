#ifndef FOURIER_H
#define FOURIER_H
namespace MATH {

/// This is my C++ interpretation of the original Fortran algorithm - SSF
/// The input is a 2*n (n being the number of sample points) array of floating
/// point values. The int "length" the actual length of the entire array (2*n).  
/// The real and imaginary components of the input are data[0] - real and data[1] -
/// imaginary; data[2] - real, data[3] - imaginary. 
/// In math terms it is: data[even] + i*data[odd]
/// If the data is a sampled set, then let data[even] be those values and let data[odd]
/// be zero (0.0f). 

class Fourier {
	
	public:
	enum DIRECTION : char {
		FORWARD = 1,
		REVERSE = -1
	};	
	void static fft(float* data, const int length, const DIRECTION direction);
};

} /// End of namespace MATH

#endif