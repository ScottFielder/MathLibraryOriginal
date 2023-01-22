#ifndef MMATH_H
#define MMATH_H

#include "Matrix.h"
#include "Plane.h"
#include "AxisAngle.h"
#include "Euler.h"
#include "Quaternion.h"
namespace  MATH {

	class MMath {
	public:
		static Matrix4 rotate(float degrees_, float x_, float y_, float z_) {
			float cosang, sinang, cosm;
			Vec3 rotAxis(x_, y_, z_);
			rotAxis = VMath::normalize(rotAxis);

			cosang = cos(degrees_ * DEGREES_TO_RADIANS);
			sinang = sin(degrees_ * DEGREES_TO_RADIANS);
			cosm = (1.0f - cosang);

			Matrix4 m;

			m[0] = (rotAxis.x * rotAxis.x * cosm) + cosang;
			m[1] = (rotAxis.x * rotAxis.y * cosm) + (rotAxis.z * sinang);
			m[2] = (rotAxis.x * rotAxis.z * cosm) - (rotAxis.y * sinang);
			m[3] = 0.0;
			m[4] = (rotAxis.x * rotAxis.y * cosm) - (rotAxis.z * sinang);
			m[5] = (rotAxis.y * rotAxis.y * cosm) + cosang;
			m[6] = (rotAxis.y * rotAxis.z * cosm) + (rotAxis.x * sinang);
			m[7] = 0.0;
			m[8] = (rotAxis.x * rotAxis.z * cosm) + (rotAxis.y * sinang);
			m[9] = (rotAxis.y * rotAxis.z * cosm) - (rotAxis.x * sinang);
			m[10] = (rotAxis.z * rotAxis.z * cosm) + cosang;
			m[11] = 0.0;
			m[12] = 0.0;
			m[13] = 0.0;
			m[14] = 0.0;
			m[15] = 1.0;
			return m;
		}

		static Matrix4 rotate(const float degrees_, const Vec3& axis_) {
			return MMath::rotate(degrees_, axis_.x, axis_.y, axis_.z);
		}


		static Matrix4 perspective(const float fovy_, const float aspect_, const float zNear_, const float zFar_) {
			float cot = 1.0f / tan(fovy_ * 0.5f * DEGREES_TO_RADIANS);
			/// Don't forget, this looks row centric but it really is a column matrix - right-hand rule rules
			Matrix4 result(cot / aspect_, 0.0f, 0.0f, 0.0f,
				0.0f, cot, 0.0f, 0.0f,
				0.0f, 0.0f, (zNear_ + zFar_) / (zNear_ - zFar_), -1.0f,
				0.0f, 0.0f, (2.0f * zNear_ * zFar_) / (zNear_ - zFar_), 0.0f);
			return result;
		}



		/// This creates a transform from Normalized Device Coordinates (NDC) to 
		/// screen coordinates. OpenGL uses NDC as the base corrdinate system.			 
		///	              ------------------------------
		///	             /|                           /|
		///	            / |                          / |
		///	           /  |                         /  |
		///	          /   |                        /   |
		///	         /    |                       /    |
		///	        /     |                      /     |
		///	       /      |                     /      |
		///	      /       |                    /       |
		///	     /        |                   /        |
		///	    /         |                  /         |
		///	   /----------------------------/ (1.0,1.0)|
		///	   |          |                 |          |
		///	   |          |                 |          |      +Y
		///	   |          |                 |          | 
		///	   |          |-----------------|----------|      ^
		///	   |         /                  |         /       |
		///	   |        /                   |        /        |       -Z
		///	   |       /                    |       /         |
		///	   |      /                     |      /          |     /
		///	   |     /                      |     /           |    /
		///	   |    /                       |    /            |   /
		///	   |   /                        |   /             |  /
		///	   |  /                         |  /              | /
		///	   | /                          | /               |/
		///	   |/ (-1.0,-1.0)               |/                ----------------> +X
		///	   ------------------------------
		static Matrix4 viewportNDC(int width_, int height_) {
			float minZ = -1.0f;
			float maxZ = 1.0f;

			Matrix4 m;
			Matrix4 m1 = scale(1.0f, -1.0f, 1.0f);
			Matrix4 m2 = scale(float(width_) / 2.0f, float(height_) / 2.0f, maxZ - minZ);
			Matrix4 m3 = translate(float(width_) / 2.0f, float(height_) / 2.0f, minZ);
			m = m3 * m2 * m1;

			///This might be slightly faster way but who cares we do it rarely 
			/***
			m[0] = float(width_)/2.0f;
			m[5] = -float(height_)/2.0f;
			m[10] =  maxZ - minZ;
			m[12] = float(width_)/2.0f;
			m[13] = float(height_)/2.0f;
			m[14] = minZ;
			m[15] = 1.0f; ***/

			return m;
		}


		static Matrix4 orthographic(float xMin_, float xMax_, float yMin_, float yMax_, float zMin_, float zMax_) {
			Matrix4 m;
			Matrix4 m1 = MMath::scale(2.0f / (xMax_ - xMin_), 2.0f / (yMax_ - yMin_), -2.0f / (zMax_ - zMin_));
			Matrix4 m2 = MMath::translate(-(xMax_ + xMin_) / (xMax_ - xMin_), -(yMax_ + yMin_) / (yMax_ - yMin_), -(zMax_ + zMin_) / (zMax_ - zMin_));
			m = m2 * m1;
			/*** Here's another way to do it
			m[0] = 2.0f / (xMax - xMin);
			m[5] = 2.0f / (yMax - yMin);
			m[10] = -2.0f / (zMax - zMin);
			m[12] = -((xMax + xMin) / (xMax - xMin));
			m[13] = -((yMax + yMin) / (yMax - yMin));
			m[14] = -((zMax + zMin) / (zMax - zMin));
			m[15] = 1.0f;
			***/
			return m;
		}

		/// The orthographic projection matrix is linear and affine but nothing else so there is has no inverse
		/// Therefore, it is labeled singular or non-invertable.
		/// I would still like to back transform from screen space to world space though
		/// Here's my unOrtho - It undoes what the orthographic matrix creates
		/// Multiply screen coordinates by this matrix and you should get x and y world coordinates
		/// It's just for fun
		static Matrix4 unOrtho(const Matrix4& ortho) {
			Matrix4 m;
			m[0] = 1.0f / ortho[0];
			m[5] = 1.0f / ortho[5];
			m[10] = 1.0f / ortho[10];
			m[12] = -ortho[12] * m[0];
			m[13] = -ortho[13] * m[5];
			m[14] = -ortho[14] * m[10];
			m[15] = 1.0f;
			return m;
		}

		/// At first glance, it might look like this matrix 
		/// is written left-handed or transposed, it has not. 
		/// Remember how memory is layed out. It is still column based.  
		/// Tested Feb 1 2013 SSF  
		static Matrix4 translate(float x_, float y_, float z_) {
			return Matrix4(1.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 1.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 1.0f, 0.0f,
				x_, y_, z_, 1.0f);
		}

		static Matrix4 translate(const Vec3& translate_) {
			return MMath::translate(translate_.x, translate_.y, translate_.z);
		}

		static Matrix4 scale(float x_, float y_, float z_) {
			return Matrix4(x_, 0.0f, 0.0f, 0.0f,
				0.0f, y_, 0.0f, 0.0f,
				0.0f, 0.0f, z_, 0.0f,
				0.0f, 0.0f, 0.0f, 1.0f);
		}

		static Matrix4 scale(const Vec3& scale) {
			return MMath::scale(scale.x, scale.y, scale.z);
		}

		///Tested Feb 1 2013 SSF
		static Matrix4 lookAt(float eyeX, float eyeY, float eyeZ,
			float atX, float atY, float atZ,
			float upX, float upY, float upZ) {

			Vec3 at(atX, atY, atZ);
			Vec3 up(upX, upY, upZ);
			Vec3 eye(eyeX, eyeY, eyeZ);

			Matrix4 result;

			Vec3 forward = VMath::normalize(at - eye);
			up = VMath::normalize(up);
			Vec3 side = VMath::normalize(VMath::cross(forward, up));
			up = VMath::cross(side, forward);

			result[0] = side.x;
			result[1] = side.y;
			result[2] = side.z;
			result[3] = 0.0;

			result[4] = up.x;
			result[5] = up.y;
			result[6] = up.z;
			result[7] = 0.0;

			result[8] = -forward.x;
			result[9] = -forward.y;
			result[10] = -forward.z;
			result[11] = 0.0;

			result[12] = -VMath::dot(side, eye);
			result[13] = -VMath::dot(up, eye);
			result[14] = VMath::dot(forward, eye);
			result[15] = 1.0;

			return result;
		}

		static Matrix4 lookAt(const Vec3& eye, const Vec3& at, const Vec3& up) {
			return lookAt(eye.x, eye.y, eye.z, at.x, at.y, at.z, up.x, up.y, up.z);
		}

		/// Take the transpose of a matrix, swap row with columns 
		/// Tested 2016
		static Matrix4 transpose(const Matrix4& m) {
			return Matrix4(m[0], m[4], m[8], m[12],
				m[1], m[5], m[9], m[13],
				m[2], m[6], m[10], m[14],
				m[3], m[7], m[11], m[15]);

		}

		static Matrix3 transpose(const Matrix3& m) {
			return Matrix3(m[0], m[3], m[6],
				m[1], m[4], m[7],
				m[2], m[5], m[8]);
		}

		static float determinate(const Matrix3& m) {
			float a = m[0] * (m[4] * m[8] - m[7] * m[5]);
			float b = m[1] * (m[3] * m[8] - m[6] * m[5]);
			float c = m[2] * (m[3] * m[7] - m[6] * m[4]);
			return a - b + c;
		}

		static float determinate(const Matrix4& m) {
			float a = m[0] * (m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10]);
			float b = m[1] * (-m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10]);
			float c = m[2] * (m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9]);
			float d = m[3] * (-m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9]);
			return a + b + c + d;
		}

		static Matrix3 inverse(const Matrix3& m) {
			Matrix3 inverseM;
			Matrix3 transposeM;
			Matrix3 adjointM;
			float determinate;
			float a = m[0] * (m[4] * m[8] - m[7] * m[5]);
			float b = m[1] * (m[3] * m[8] - m[6] * m[5]);
			float c = m[2] * (m[3] * m[7] - m[6] * m[4]);
			determinate = a - b + c;
			transposeM = transpose(m);
			adjointM[0] = transposeM[4] * transposeM[8] - transposeM[7] * transposeM[5];
			adjointM[1] = -transposeM[3] * transposeM[8] + transposeM[6] * transposeM[5];
			adjointM[2] = transposeM[3] * transposeM[7] - transposeM[6] * transposeM[4];
			adjointM[3] = -transposeM[1] * transposeM[8] + transposeM[7] * transposeM[2];
			adjointM[4] = transposeM[0] * transposeM[8] - transposeM[6] * transposeM[2];
			adjointM[5] = -transposeM[0] * transposeM[7] + transposeM[6] * transposeM[1];
			adjointM[6] = transposeM[1] * transposeM[5] - transposeM[4] * transposeM[2];
			adjointM[7] = -transposeM[0] * transposeM[5] + transposeM[3] * transposeM[5];
			adjointM[8] = transposeM[0] * transposeM[4] - transposeM[3] * transposeM[1];

#ifdef _DEBUG  /// If in debug mode let's worry about divide by zero or nearly zero!!! 
			if (fabs(determinate) < VERY_SMALL) {
				std::string errorMsg = __FILE__ + __LINE__;
				throw errorMsg.append(": Divide by nearly zero! ");
			}
#endif
			for (int i = 0; i < 9; i++) {
				inverseM[i] = adjointM[i] / determinate;
			}
			return inverseM;
		}


		/// 4x4 no way, this is tough stuff
		/// Tested 2013
		static Matrix4 inverse(const Matrix4& m) {
			Matrix4 inverseM;
			Matrix4 adjointM;
			float determinate;

			adjointM[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
			adjointM[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
			adjointM[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
			adjointM[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
			adjointM[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
			adjointM[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
			adjointM[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
			adjointM[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
			adjointM[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
			adjointM[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
			adjointM[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
			adjointM[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
			adjointM[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
			adjointM[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
			adjointM[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
			adjointM[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

			determinate = m[0] * adjointM[0] + m[1] * adjointM[4] + m[2] * adjointM[8] + m[3] * adjointM[12];

#ifdef _DEBUG  /// If in debug mode let's worry about divide by zero or nearly zero!!! 
			if (fabs(determinate) < VERY_SMALL) {
				std::string errorMsg = __FILE__ + __LINE__;
				throw errorMsg.append(": Divide by nearly zero! ");
			}
#endif
			for (int i = 0; i < 16; i++) {
				inverseM[i] = adjointM[i] / determinate;
			}
			return inverseM;
		}

		/// Convert Eular angles to a 3x3 rotation matrix
		static Matrix3 toMatrix3(const Euler& e) {
			/// Note: If you want to multiply xaxis, yaxis,zaix in that order. I think
			/// it should be m = x * z * y <- reading right to left. .
			Matrix3 m = Matrix3(MMath::rotate(e.xAxis, Vec3(1.0f, 0.0f, 0.0f)) *
				MMath::rotate(e.zAxis, Vec3(0.0f, 0.0f, 1.0f)) *
				MMath::rotate(e.yAxis, Vec3(0.0f, 1.0f, 0.0f)));
			return m;
		}

		static Matrix4 toMatrix4(const AxisAngle& axisAngle_) {
			return MMath::rotate(axisAngle_.angle, axisAngle_.axis.x, axisAngle_.axis.y, axisAngle_.axis.z);
		}

		static Matrix3 toMatrix3(const Quaternion& q) {
			/// This is the fastest way I know...
			return Matrix3((1.0f - 2.0f * q.ijk.y * q.ijk.y - 2.0f * q.ijk.z * q.ijk.z), (2.0f * q.ijk.x * q.ijk.y + 2.0f * q.ijk.z * q.w), (2.0f * q.ijk.x * q.ijk.z - 2.0f * q.ijk.y * q.w),
				(2.0f * q.ijk.x * q.ijk.y - 2.0f * q.ijk.z * q.w), (1.0f - 2.0f * q.ijk.x * q.ijk.x - 2.0f * q.ijk.z * q.ijk.z), (2.0f * q.ijk.y * q.ijk.z + 2.0f * q.ijk.x * q.w),
				(2.0f * q.ijk.x * q.ijk.z + 2.0f * q.ijk.y * q.w), (2.0f * q.ijk.y * q.ijk.z - 2 * q.ijk.x * q.w), (1.0f - 2.0f * q.ijk.x * q.ijk.x - 2.0f * q.ijk.y * q.ijk.y));
		}


		static Matrix4 toMatrix4(const Quaternion& q) {
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

	};

}
#endif