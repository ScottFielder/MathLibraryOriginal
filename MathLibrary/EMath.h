#ifndef EMATH_H
#define EMATH_H
#include "Euler.h"
#include "Matrix.h"
#include "Quaternion.h"

namespace MATH {

	class EMath {
	public:
		static Euler toEuler(const Matrix3 &m);
		static Euler toEular(const Quaternion& q);
	};
}
#endif