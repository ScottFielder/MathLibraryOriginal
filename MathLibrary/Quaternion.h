#ifndef QUATERTNION_H
#define QUATERTNION_H
#include <iostream>
#include "VMath.h"
#include "Matrix.h"

/***
A quaternion is a mathematical object that can be used to represent rotations in 3D space.
It is a 4-dimensional extension of the complex numbers, and it consists of a scalar component (w)
and a 3-vector component (ijk).

Quaternions are often used in computer graphics, robotics, and aerospace engineering to represent
and manipulate rotations in a more efficient and stable way than traditional rotation matrices. 
They can be used to interpolate between two orientations in 3D space, and to avoid problems such
as gimbal lock and singularities that can occur when using Euler angles.

Quaternions can be represented in the form q = w + xi + yj + zk, where w, x, y, and z are real numbers,
and i, j, and k are imaginary units that satisfy the following properties: i^2 = j^2 = k^2 = ijk = -1.

Quaternions can be added, subtracted, and multiplied, and the multiplication is not commutative. 
The multiplication of two quaternions can be used to represent a rotation of one quaternion followed
by a rotation of the other quaternion.

The norm of a quaternion is defined as the square root of the sum of the squares of the components and
is also known as the magnitude of the quaternion.

A unit quaternion is a quaternion with a magnitude of 1. Unit quaternions can be used to represent
rotations without scaling.

Quaternions have many advantages over other representations of rotations in 3D space, such as Euler angles and rotation matrices, 
including being more efficient to compute, more robust to numerical errors and gimbal lock, and more convenient to interpolate and integrate.
***/

/// A quarternion can be written as a scalar plus a 3-vector (Vec3) component

namespace  MATH {
	struct Quaternion {
		float w;
		Vec3 ijk; /// These are the ijk components of the Quaternion 

		/// Just a little utility to populate a quaternion
		inline void set(float w_, float x_, float y_, float z_) {
			w = w_; ijk.x = x_; ijk.y = y_; ijk.z = z_;
		}

		/// Another little utility to populate a quaternion
		inline void set(float w_, Vec3 ijk_) {
			w = w_; ijk.x = ijk_.x; ijk.y = ijk_.y; ijk.z = ijk_.z;
		}

		/// Default constructor that sets the quaternion to the unit quaternion 
		inline Quaternion() {
			set(1.0f, 0.0f, 0.0f, 0.0f);
		}

		/// Another utility function to set the values of the quaternion
        /// w_ is the scalar component and ijk_ is the vector component
		inline Quaternion(float w_, const Vec3& ijk_) {
			set(w_, ijk_.x, ijk_.y, ijk_.z);
		}

		/// A copy constructor
		inline Quaternion(const Quaternion& q) {
			set(q.w, q.ijk.x, q.ijk.y, q.ijk.z);
		}

		/// An assignment operator   
		inline Quaternion& operator = (const Quaternion& q) {
			set(q.w, q.ijk.x, q.ijk.y, q.ijk.z);
			return *this;
		}

		/// Take the negative of a Quaternion
		inline const Quaternion operator - () const {
			return Quaternion(-w, Vec3(-ijk.x, -ijk.y, -ijk.z));
		}

			// 2022-02-12 Umer Noor edit. I think there is a bug here
			// I'll change the Vec3 on the stack to be ijk_result and
			// see if that helps. Compiler might be grabbing the ijk member
			// variable in the return line rather than the stack variable
		/// Multiply a two quaternions - using the right-hand rule 
		/// 2022-02-12 Scott edit. Worried that Umer uncovered a bug in my code,
		/// I derived the multiply over again (this time less sexy) and it seems to work 
		/// correctly, Thanks Umer
		inline const Quaternion operator * (const Quaternion& q) const {
			Quaternion result;
			result.w = (w * q.w) - (ijk.x * q.ijk.x) - (ijk.y * q.ijk.y) - (ijk.z * q.ijk.z);
			result.ijk.x = (w * q.ijk.x) + (ijk.x * q.w) + (ijk.z * q.ijk.y) - (ijk.y * q.ijk.z);
			result.ijk.y = (w * q.ijk.y) + (ijk.y * q.w) + (ijk.x * q.ijk.z) - (ijk.z * q.ijk.x);
			result.ijk.z = (w * q.ijk.z) + (ijk.z * q.w) + (ijk.y * q.ijk.x) - (ijk.x * q.ijk.y);
			return result;
			
		}

		inline const Quaternion& operator *= (const Quaternion& q) {
			*this = q * *this;
			return *this;
		}

		inline const Quaternion operator + (const Quaternion q) const {
			return Quaternion(w + q.w, Vec3(ijk.x + q.ijk.x, ijk.y + q.ijk.y, ijk.z + q.ijk.z));
		}

		inline const Quaternion operator - (const Quaternion q) const {
			return Quaternion(w - q.w, Vec3(ijk.x - q.ijk.x, ijk.y - q.ijk.y, ijk.z - q.ijk.z));
		}

		inline const Quaternion operator * (const float scalar) const {
			return Quaternion(w * scalar, Vec3(ijk.x * scalar, ijk.y * scalar, ijk.z * scalar));
		}

		inline const Quaternion operator / (const float scalar) const {
			return Quaternion(w / scalar, Vec3(ijk.x / scalar, ijk.y / scalar, ijk.z / scalar));
		}

		/// Now we can use the Quaternion like an array but we'll need two overloads
		inline const float operator [] (int index) const {  /// This one is for reading the Quaternion as if where an array
			return *((float*)this + index);
		}

		inline float& operator [] (int index) {	/// This one is for writing to the Quaternion as if where an array.  
			return *((float*)this + index);
		}

		inline void print(const char* comment = nullptr) const {
			if (comment) printf("%s\n", comment);
			printf("%1.4f %1.4f %1.4f %1.4f\n", w, ijk.x, ijk.y, ijk.z);
		}


		/////////////////////////////////////////////////////////////////////////
		/// This is just for teaching purposes - Caution, I'm getting out of control
		/////////////////////////////////////////////////////////////////////////
		/// Multiply a quaternion by a Vec3 (Quaternion * Vec3) 
		inline const Vec3 operator * (const Vec3& v_) const {
			/// Promote the Vec3 to a Quaternion and set w to be 0.0
			Quaternion p(0.0, v_);
			/// Now just call the Quaternion * Quaternion operator
			Quaternion result = *this * p;
			return result.ijk;
		}

		/// Multiply a Vec3 by a Quaternion (Vec3 * Quaternion) 
		friend Vec3 operator * (const Vec3 v, const Quaternion& q) {
			Quaternion qv(0.0f, v);
			Quaternion result = qv * q;
			return result.ijk;
		}

		/// Seriously, the tilde ~ is the complement operator not the 
		///  conjugate - but it was for fun. 
		inline Quaternion operator~() { return Quaternion(w, -ijk); }
		
	};
}
#endif
