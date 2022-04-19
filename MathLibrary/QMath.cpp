
#include <algorithm>    // std::min, std::clamp
#include "QMath.h"
#include "Matrix.h"
#include "Euler.h"
using namespace MATH;

float QMath::magnitude(const Quaternion& q) {
	return(sqrt((q.w * q.w) + (q.ijk.x * q.ijk.x) + (q.ijk.y * q.ijk.y) + (q.ijk.z * q.ijk.z)));
}

Quaternion QMath::conjugate(const Quaternion& q) {
	return Quaternion(q.w, -q.ijk);
}

Quaternion QMath::inverse(const Quaternion& q) {
	float mag = magnitude(q);
	Quaternion conj = conjugate(q);
	return  conj / (mag * mag);
}

Quaternion QMath::pow(const Quaternion& q, float exponent) {
	/// Watch out, if q.w is +/- 1, acos(q.w) will be 0 and cause divide by zero exception
	if (fabs(q.w) > VERY_CLOSE_TO_ONE) return q;

	float alpha = acos(q.w);
	float newAlpha = alpha * exponent;
	return Quaternion(cos(newAlpha), q.ijk * sin(newAlpha) / sin(alpha));
}


Quaternion QMath::normalize(const Quaternion& q){
return q / magnitude(q);
}

Euler QMath::fromQuaternion(const Quaternion& q) {
	Euler result;
	result.roll = atan2(2.0f * (q.ijk.x * q.ijk.y + q.w * q.ijk.z), q.w * q.w + q.ijk.x * q.ijk.x - q.ijk.y * q.ijk.y - q.ijk.z * q.ijk.z);
	result.pitch = atan2(2.0f * (q.ijk.y * q.ijk.z + q.w * q.ijk.x), q.w * q.w - q.ijk.x * q.ijk.x - q.ijk.y * q.ijk.y + q.ijk.z * q.ijk.z);
	result.yaw = asin(std::clamp(-2.0f * (q.ijk.x * q.ijk.z - q.w * q.ijk.y), -1.0f, 1.0f));
	return result * RADIANS_TO_DEGREES;
}

Quaternion QMath::fromEuler(const Euler& e) {
	float cosX = cos(0.5f * e.xAxis * DEGREES_TO_RADIANS);
	float cosY = cos(0.5f * e.yAxis * DEGREES_TO_RADIANS);
	float cosZ = cos(0.5f * e.zAxis * DEGREES_TO_RADIANS);
	float sinX = sin(0.5f * e.xAxis * DEGREES_TO_RADIANS);
	float sinY = sin(0.5f * e.yAxis * DEGREES_TO_RADIANS);
	float sinZ = sin(0.5f * e.zAxis * DEGREES_TO_RADIANS);

	/// For XYZ order
	return Quaternion((cosX * cosY * cosZ) + (sinX * sinY * sinZ),
		Vec3((sinX * cosY * cosZ) - (cosX * sinY * sinZ),
			(cosX * sinY * cosZ) + (sinX * cosY * sinZ),
			(cosX * cosY * sinZ) - (sinX * sinY * cosZ)));
}

Matrix3 QMath::toMatrix3(const Quaternion& q) {
	/// This is the fastest way I know...
	return Matrix3((1.0f - 2.0f * q.ijk.y * q.ijk.y - 2.0f * q.ijk.z * q.ijk.z), (2.0f * q.ijk.x * q.ijk.y + 2.0f * q.ijk.z * q.w), (2.0f * q.ijk.x * q.ijk.z - 2.0f * q.ijk.y * q.w),
		(2.0f * q.ijk.x * q.ijk.y - 2.0f * q.ijk.z * q.w), (1.0f - 2.0f * q.ijk.x * q.ijk.x - 2.0f * q.ijk.z * q.ijk.z), (2.0f * q.ijk.y * q.ijk.z + 2.0f * q.ijk.x * q.w),
		(2.0f * q.ijk.x * q.ijk.z + 2.0f * q.ijk.y * q.w), (2.0f * q.ijk.y * q.ijk.z - 2 * q.ijk.x * q.w), (1.0f - 2.0f * q.ijk.x * q.ijk.x - 2.0f * q.ijk.y * q.ijk.y));
}

Matrix4 QMath::toMatrix4(const Quaternion& q) {
	/// This is the fastest way I know...
	return Matrix4((1.0f - 2.0f * q.ijk.y * q.ijk.y - 2.0f * q.ijk.z * q.ijk.z), (2.0f * q.ijk.x * q.ijk.y + 2.0f * q.ijk.z * q.w), (2.0f * q.ijk.x * q.ijk.z - 2.0f * q.ijk.y * q.w), 0.0f,
		(2.0f * q.ijk.x * q.ijk.y - 2.0f * q.ijk.z * q.w), (1.0f - 2.0f * q.ijk.x * q.ijk.x - 2.0f * q.ijk.z * q.ijk.z), (2.0f * q.ijk.y * q.ijk.z + 2.0f * q.ijk.x * q.w), 0.0f,
		(2.0f * q.ijk.x * q.ijk.z + 2.0f * q.ijk.y * q.w), (2.0f * q.ijk.y * q.ijk.z - 2 * q.ijk.x * q.w), (1.0f - 2.0f * q.ijk.x * q.ijk.x - 2.0f * q.ijk.y * q.ijk.y), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);

	/// ... but this is the coolest. My way is just a bit faster on single processor machines,
	/// this method is faster on parallel multicore machines. Multicore can calc m1 and m2
	/// on separate threads. Just saying. 

	//Matrix4 m1( q.w,  q.ijk.z,  -q.ijk.y,  q.ijk.x,
	//		-q.ijk.z,   q.w,   q.ijk.x,  q.ijk.y,
	//		q.ijk.y,  -q.ijk.x,   q.w,  q.ijk.z,
	//		-q.ijk.x,  -q.ijk.y,  -q.ijk.z,  q.w);
	//
	//Matrix4 m2( q.w,   q.ijk.z,  -q.ijk.y,  -q.ijk.x,
	//			-q.ijk.z,   q.w,   q.ijk.x,  -q.ijk.y,
	//			q.ijk.y,  -q.ijk.x,   q.w,  -q.ijk.z,
	//			-q.ijk.x,   q.ijk.y,   q.ijk.z,   q.w);
	//return m1 * m2;

}
Quaternion QMath::angleAxisRotation(const float degrees, const Vec3& axis) {
	Vec3 rotationAxis = VMath::normalize(axis);
	float theta = degrees * DEGREES_TO_RADIANS;
	float cosVal = cos(theta / 2.0f);
	float sinVal = sin(theta / 2.0f);
	Quaternion result = Quaternion(cosVal, rotationAxis * sinVal);
	return result;
}

/// Given  Vector v and Quaternion q, return a rotated Vector by q.  
Vec3 QMath::rotate(const Vec3& v, const Quaternion& q) {
	/// This is the beauty of Quaternions 
	Quaternion p(0.0, v); /// convert the incoming vector to a Quaternion
	Quaternion qInv = QMath::inverse(q); /// Get the inverse - duh.
	Quaternion result = q * p * qInv;
	return Vec3(result.ijk);

	/*** I got the idea to do it this way from glm. They say it is faster. 
	Vec3 quatVector(q.v);
	Vec3 uv = VMath::cross(q.v, v);
	Vec3 uuv = VMath::cross(q.v, uv);
	return v + ((uv * q.w) + uuv) * 2.0f;***/


}


Quaternion QMath::slerp(const Quaternion& qa, const Quaternion& qb, float t) {	
	Quaternion q1 = qa;
	Quaternion q2 = qb;
	float cosTheta = dot(q1, q2); /// if cosTheta is nearly 1.0 will cause divide by zero error

	if (cosTheta < 0.0f) {		/// if cosTheta is negative, the angle is oblique. The shortest path 
		q2 = -q2;				/// would be the other representation of the same angle -q2 
		cosTheta = cosTheta;
	}
	float c1, c2;
	///If cosTheta is very close to 1.0 just lerp it to prevent divide by zero
	if (cosTheta > VERY_CLOSE_TO_ONE) {
		c1 = 1.0f - t;
		c2 = t;
	} else {
		float theta = acos(cosTheta);
		float sinTheta = sin(theta); 
		/// or sqrt(1.0 - cosTheta * cosTheta)
		//float sinTheta = sqrt(1.0f - (cosTheta * cosTheta));
		//float theta = atan2(sinTheta, cosTheta);

		c1 = sin((1.0f - t) * theta) / sinTheta;
		c2 = sin(t * theta) / sinTheta;
	}
	Vec3 ijk(c1 * q1.ijk.x + c2 * q2.ijk.x,
			c1 * q1.ijk.y + c2 * q2.ijk.y,
			c1 * q1.ijk.z + c2 * q2.ijk.z);
	return Quaternion(c1 * q1.w + c2 * q2.w, ijk);
		
}



