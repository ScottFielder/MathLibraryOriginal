#ifndef QMATH_H
#define QMATH_H
#include <algorithm>    // std::min, std::clamp
#include "Quaternion.h"
#include "Matrix.h"
#include "Euler.h"
#include "VMath.h"
namespace MATH{
	
	class QMath {
	public:

		static float magnitude(const Quaternion& q) {
			return(sqrt((q.w * q.w) + (q.ijk.x * q.ijk.x) + (q.ijk.y * q.ijk.y) + (q.ijk.z * q.ijk.z)));
		}

		static Quaternion conjugate(const Quaternion& q) {
			return Quaternion(q.w, -q.ijk);
		}

		static Quaternion inverse(const Quaternion& q) {
			float mag = magnitude(q);
			Quaternion conj = conjugate(q);
			return  conj / (mag * mag);
		}

		static Quaternion pow(const Quaternion& q, float exponent) {
			/// Watch out, if q.w is +/- 1, acos(q.w) will be 0 and cause divide by zero exception
			if (fabs(q.w) > VERY_CLOSE_TO_ONE) return q;

			float alpha = acos(q.w);
			float newAlpha = alpha * exponent;
			return Quaternion(cos(newAlpha), q.ijk * sin(newAlpha) / sin(alpha));
		}

		/// 2022-02-12 A  Dr. Umer Noor edit. Putting in a normalize method
		static Quaternion normalize(const Quaternion& q){
			return q / magnitude(q);
		}

		/// 2022-04-04 A quaternion dot product
		static float dot(const Quaternion &a, const Quaternion &b){
			return((a.w * b.w) + (a.ijk.x * b.ijk.x) + (a.ijk.y * b.ijk.y) + (a.ijk.z * b.ijk.z));
		}

		static Quaternion lookAt(const Vec3& direction, const Vec3& up) {
			Matrix3 result;
			Vec3 dir = VMath::normalize(direction);
			result.setColumn(Matrix3::Column::two, -dir);
			result.setColumn(Matrix3::Column::zero, VMath::normalize(VMath::cross(up, -dir)));
			result.setColumn(Matrix3::Column::one, VMath::cross(-dir, result.getColumn(Matrix3::Column::zero)));
			return toQuaternion(result);
		}

		static Quaternion toQuaternion(const Matrix3& m) {
			float trace;
			Quaternion qResult;
			if (m[8] < 0.0f) {
				if (m[0] > m[4]) {
					trace = 1.0f + m[0] - m[4] - m[8];
					Vec3 ijk(trace, m[1] + m[3], m[6] + m[2]);
					qResult.set(m[5] - m[7], ijk);
				} else {
					trace = 1.0f - m[0] + m[4] - m[8];
					Vec3 ijk(m[1] + m[3], trace, m[5] + m[7]);
					qResult.set( m[6] - m[2], ijk);
				}
			} else {
				if (m[0] < -m[4]) {
					trace = 1.0f - m[0] - m[4] + m[8];
					Vec3 ijk(m[6] + m[2], m[5] + m[7], trace);
					qResult.set(m[1] - m[3], ijk);
				} else {
					trace = 1.0f + m[0] + m[4] + m[8];
					Vec3 ijk(m[5] - m[7], m[6] - m[2], m[1] - m[3]);
					qResult.set(trace, ijk);
				}
			}
			qResult = qResult * (0.5f / sqrt(trace));
			return qResult;
	}



		static Quaternion toQuaternion(const Euler& e) {
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


		static Quaternion angleAxisRotation(const float degrees, const Vec3& axis) {
			Vec3 rotationAxis = VMath::normalize(axis);
			float theta = degrees * DEGREES_TO_RADIANS;
			float cosVal = cos(theta / 2.0f);
			float sinVal = sin(theta / 2.0f);
			Quaternion result = Quaternion(cosVal, rotationAxis * sinVal);
			return result;
		}

		/// Given  Vector v and Quaternion q, return a rotated Vector by q.  
		static Vec3 rotate(const Vec3& v, const Quaternion& q) {
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

		
		/// Not needed after C17 it is now std::clamp <algorithm>
		/***static float clamp(float x, float minVal, float maxVal) { 
				return std::min(std::max(x, minVal), maxVal); 
		}***/

		static Quaternion slerp(const Quaternion& qa, const Quaternion& qb, float t) {	
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
				////or float sinTheta = sqrt(1.0f - (cosTheta * cosTheta));
				///float theta = atan2(sinTheta, cosTheta);

				c1 = sin((1.0f - t) * theta) / sinTheta;
				c2 = sin(t * theta) / sinTheta;
			}
			Vec3 ijk(c1 * q1.ijk.x + c2 * q2.ijk.x,
					c1 * q1.ijk.y + c2 * q2.ijk.y,
					c1 * q1.ijk.z + c2 * q2.ijk.z);
			return Quaternion(c1 * q1.w + c2 * q2.w, ijk);
		
		}

	};
}
#endif