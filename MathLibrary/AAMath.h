#ifndef AAMATH_H
#define AAMATH_H
#include "AxisAngle.h"
#include "Matrix.h"
namespace MATH {
	class AAMath {
	public:
		static AxisAngle toAxisAngle(const Matrix4& pureRotationMatrix);
		static AxisAngle toAxisAngle(const Matrix3& pureRotationMatrix);
	};
}
#endif
