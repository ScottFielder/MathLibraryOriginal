#ifndef EMATH_H
#define EMATH_H
#include "Euler.h"
#include "Matrix.h"

namespace MATH {

	class EMath {
	public:
		static Matrix3 EulerToMatirix3(const Euler& e);
		static Euler Matrix3ToEuler(const Matrix3 &m);
	};

}
#endif