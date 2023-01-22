#ifndef EMATH_H
#define EMATH_H
#include <algorithm> /// std::clamp
#include "Euler.h"
#include "Matrix.h"
#include "Quaternion.h"

namespace MATH {

	class EMath {
		public:
		static Euler toEuler(const Quaternion& q) {
			Euler result;
			result.roll = atan2(2.0f * (q.ijk.x * q.ijk.y + q.w * q.ijk.z), q.w * q.w + q.ijk.x * q.ijk.x - q.ijk.y * q.ijk.y - q.ijk.z * q.ijk.z);
			result.pitch = atan2(2.0f * (q.ijk.y * q.ijk.z + q.w * q.ijk.x), q.w * q.w - q.ijk.x * q.ijk.x - q.ijk.y * q.ijk.y + q.ijk.z * q.ijk.z);
			result.yaw = asin(std::clamp(-2.0f * (q.ijk.x * q.ijk.z - q.w * q.ijk.y), -1.0f, 1.0f));
			return result * RADIANS_TO_DEGREES;
		}


		static Euler toEuler(const Matrix3 &m) {
			Euler e;
			e.xAxis = atan2(m[5], m[8]);

			float cos1 = cos(e.xAxis);
			float cos2 = sqrt((m[0] * m[0]) + (m[1] * m[1]));
			float sin1 = sin(e.xAxis);

			e.yAxis = atan2(-m[2], cos2);
			e.zAxis = atan2(sin1*m[6] - cos1*m[3], cos1*m[4] - sin1*m[7]);

			e.xAxis *= RADIANS_TO_DEGREES;
			e.yAxis *= RADIANS_TO_DEGREES;
			e.zAxis *= RADIANS_TO_DEGREES;
			return e;
		}
	};
}
#endif