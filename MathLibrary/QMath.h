#ifndef QMATH_H
#define QMATH_H
#include <cmath>
#include "Quaternion.h"
#include "Matrix.h"
#include "Euler.h"

namespace MATH{
	
	class QMath {
	public:

		static Quaternion inverse(const Quaternion& q);
		static  Quaternion conjugate(const Quaternion& q);
		static float magnitude(const Quaternion& q);
		static Quaternion pow(const Quaternion& q, float exponent);

		static Quaternion angleAxisRotation(const float degrees, const Vec3& axis);
		static Vec3 rotate(const Vec3& v, const Quaternion& q);

		static  Quaternion fromEuler(const Euler& e);
		static  Euler fromQuaternion(const Quaternion& q);
		static  Matrix3 toMatrix3(const Quaternion& q);
		static  Matrix4 toMatrix4(const Quaternion& q);

		/// Not needed after C17
		//static float clamp(float x, float minVal, float maxVal) { return std::min(std::max(x, minVal), maxVal); }
	};
}
#endif