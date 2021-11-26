
#include <time.h>
#include <stdlib.h>
#include <cmath>
#include "Randomizer.h"

using namespace MATH;

Randomizer::Randomizer(){
	seed = static_cast<long>(time(nullptr));
	ran1idnum = -seed;
	qd2idnum = -seed;
}

Randomizer::Randomizer(long seed_){
	seed = seed_;
	ran1idnum = -seed_;
	qd2idnum = -static_cast<int>(seed_);
}

Randomizer::~Randomizer(){
	seed = 0; 
	ran1idnum = 0;
	qd2idnum = 0;
}


double Randomizer::rand(const double low, const double high){
	return (static_cast<double>(ran1())*(high - low)) + low;
}

float Randomizer::rand(const float low, const float high) {
	return (ran1()*(high - low)) + low;
}

int Randomizer::rand(const int low, int high){
	return int(rand(double(low), double(high)));
}



/// Based on the work of Stephen K. Park and Keith W. Miller
float Randomizer::ran1(){
		#define IA 16807		/// primitive root modulo M31 7^5
		#define IM 2147483647	/// One of the Mersenne primes: 2^31-1
		#define AM (1.0f/static_cast<float>(IM))
		#define IQ 127773
		#define IR 2836
		#define NTAB 32
		#define NDIV (1+(IM-1)/NTAB)
		#define EPS 1.2e-7f
		#define RNMX (1.0f-EPS)
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	float temp;
	if (ran1idnum <= 0 || !iy) {
		if (-(ran1idnum) < 1) ran1idnum = 1;
		else ran1idnum = -(ran1idnum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (ran1idnum) / IQ;
			ran1idnum = IA*(ran1idnum - k*IQ) - IR*k;
			if (ran1idnum < 0) ran1idnum += IM;
			if (j < NTAB) iv[j] = ran1idnum;
		}
		iy = iv[0];
	}
	k = (ran1idnum) / IQ;
	ran1idnum = IA*(ran1idnum - k*IQ) - IR*k;
	if (ran1idnum < 0) ran1idnum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = ran1idnum;
	temp = AM * static_cast<float>(iy);
	if (temp > RNMX){
		return RNMX;
	} else {
		return temp;
	}
		#undef IA
		#undef IM
		#undef AM
		#undef IQ
		#undef IR
		#undef NTAB
		#undef NDIV
		#undef EPS
		#undef RNMX
}

/// G.E.P. Box and Mervin E. Muller, A Note on the Generation of Random Normal Deviates, 
/// "The Annals of Mathematical Statistics" 1958
double Randomizer::box_muller(const double mean, const double std_devation)	{
	int idnum = -1;
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last){ /// Use the value from the previous call - its a static 
		y1 = y2;
		use_last = 0;
	}
	else{
		do{
			x1 = 2.0 * ran1() - 1.0;
			x2 = 2.0 * ran1() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 *log(w)) / w);
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return(mean + y1 * std_devation);
}

