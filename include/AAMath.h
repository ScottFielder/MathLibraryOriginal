#ifndef AAMATH_H
#define AAMATH_H
#include "AxisAngle.h"
#include "Matrix.h"
#include "VMath.h"
namespace MATH {
	class AAMath {
	public:
		static AxisAngle toAxisAngle(const Matrix4& m){
			Vec3 axis( (m[5] - m[7]), (m[6] - m[2]), (m[1] - m[3]) );
			float angle = acos((m[0] + m[4] + m[8] - 1.0f) / 2.0f);
			return  AxisAngle(angle * RADIANS_TO_DEGREES, VMath::normalize(axis));
		}
		static AxisAngle toAxisAngle(const Matrix3& m){
			Vec3 axis((m[6] - m[9]), (m[8] - m[2]), (m[1] - m[4]));
			float angle = acos( (m[0] + m[5] + m[10] - 1.0f) / 2.0f);
			return  AxisAngle(angle * RADIANS_TO_DEGREES, VMath::normalize(axis));
		}
	};
}
#endif
