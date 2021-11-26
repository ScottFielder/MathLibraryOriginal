#ifndef RANDOMIZER_H
#define RANDOMIZER_H
namespace MATH {

class Randomizer{
	Randomizer(const Randomizer&) = delete;
	Randomizer(Randomizer&&) = delete;
	Randomizer& operator=(const Randomizer&) = delete;
	Randomizer& operator=(Randomizer&&) = delete;

private:
	long seed;
	long ran1idnum;
	int qd2idnum;
public:
	/// This will use the time in milliseconds as the seed value
	Randomizer();

	/// Using the same seed value will generate the same list of pseudo-random numbers
	/// each instance
	Randomizer(const long seed);

	~Randomizer();

	/// Generate numbers between the low (smallest) and high values
	double rand(const double low, const double high);

	/// Generate numbers between the low (smallest) and high values
	float rand(const float low, const float high);

	/// Generate numbers between the low (smallest) and high values
	int rand(const int low, const int high);

	/// This is the pseudo-random number generator 
	float ran1();

	/// Very fast but a bit sloppy (statistically) randomizer, I just included this just for fun
	inline double rand2(){
		unsigned long itemp;
		static unsigned long jflone = 0x3f800000;
		static unsigned long jflmsk = 0x007fffff;
		qd2idnum = 1664525L * qd2idnum + 1013904223L;
		itemp = jflone | (jflmsk & qd2idnum);
		return (double)(*(float *)&itemp) - 1.0;
	}

	/// This method developed by George Box and Mervin Muller in 1958 
	/// (the year I was born), is a mathematical way of computing a gaussian
	/// distribution  (zero expectation, unit variance) of random numbers, given 
	/// a source of standard linear random numbers. The distribution is around 
	/// the mean and spread by the standard deviation. 
	/// ref: G. E. P. Box and Mervin E. Muller, A Note on the Generation of Random Normal Deviates,
	/// The Annals of Mathematical Statistics (1958), Vol. 29, No. 2 pp. 610–611
	double box_muller(const double mean, const double std_devation);

	inline float box_muller(const float mean, const float std_devation) {
		return static_cast<float>(box_muller(static_cast<double>(mean), static_cast<double>(std_devation)));
	}
};

}/// End of namespace MATH

#endif