#ifndef QMATH_H
#define QMATH_H
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

		/// 2022-02-12 A  Dr. Umer Noor edit. Putting in a normalize method
		static Quaternion normalize(const Quaternion& q);

		static Quaternion lookAt(const Vec3& lookingAt, const Vec3& up = Vec3(0.0f,1.0f,0.0f));
		
		static Quaternion angleAxisRotation(const float degrees, const Vec3& axis);
		static Vec3 rotate(const Vec3& v, const Quaternion& q);

		

		/// 2022-04-04 A quaternion dot product
		inline static float dot(const Quaternion &a, const Quaternion &b){
			return((a.w * b.w) + (a.ijk.x * b.ijk.x) + (a.ijk.y * b.ijk.y) + (a.ijk.z * b.ijk.z));
		}


		static Quaternion slerp(const Quaternion& q1, const Quaternion& q2, float t);

		/// Not needed after C17 it is now std::clamp <algorithm>

		///static float clamp(float x, float minVal, float maxVal) { 
		///		return std::min(std::max(x, minVal), maxVal); 
		///}


		static  Quaternion toQuaternion(const Matrix3& m);
		static  Quaternion toQuaternion(const Matrix4& m);
		static  Quaternion toQuaternion(const Euler& e);

	};
}
#endif