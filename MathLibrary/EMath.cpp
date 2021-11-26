#include "EMath.h"
#include "MMath.h"

using namespace MATH;

/// Convert Eular angles to a 3x3 rotation matrix
Matrix3 EMath::EulerToMatirix3(const Euler& e) {
	/// Note: If you want to multiply xaxis, yaxis,zaix in that order. I think
	/// it should be m = x * z * y <- reading right to left. .
	Matrix3 m = Matrix3(MMath::rotate(e.xAxis, Vec3(1.0f, 0.0f, 0.0f)) *
						MMath::rotate(e.zAxis, Vec3(0.0f, 0.0f, 1.0f)) *
						MMath::rotate(e.yAxis, Vec3(0.0f, 1.0f, 0.0f)));
	return m;
}



#define RADIANS2DEGREES 180.0f/M_PI
	Euler EMath::Matrix3ToEuler(const Matrix3 &m) {
		Euler e;
		e.xAxis = atan2(m[5], m[8]);

		float cos1 = cos(e.xAxis);
		float cos2 = sqrt((m[0] * m[0]) + (m[1] * m[1]));
		float sin1 = sin(e.xAxis);

		e.yAxis = atan2(-m[2], cos2);
		e.zAxis = atan2(sin1*m[6] - cos1*m[3], cos1*m[4] - sin1*m[7]);

		e.xAxis *= RADIANS2DEGREES;
		e.yAxis *= RADIANS2DEGREES;
		e.zAxis *= RADIANS2DEGREES;
		return e;
	}
#undef RADIANS2DEGREES




