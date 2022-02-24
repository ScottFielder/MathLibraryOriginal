#ifndef QUATERTNION_H
#define QUATERTNION_H
#include <iostream>
#include "VMath.h"
#include "Matrix.h"

///A quarternion can be written as a scalar plus a 3-vector (Vec3) component

namespace  MATH {
	struct Quaternion {
		float w;
		Vec3 ijk; /// These are the ijk components of the Quaternion 

		/// Just a little utility to populate a quaternion
		inline void set(float w_, float x_, float y_, float z_) {
			w = w_; ijk.x = x_; ijk.y = y_; ijk.z = z_;

			/// The memory layout is understood Dr. R!			
		  /**static_cast<float*>(&w+0) = w_;
			*static_cast<float*>(&w+1) = x_;
			*static_cast<float*>(&w+2) = y_;
			*static_cast<float*>(&w+3) = z_;*/
		}

		/// This is the unit quaterion by definition 
		inline Quaternion() {
			set(1.0f, 0.0f, 0.0f, 0.0f);
		}

		inline Quaternion(float w_, float x_, float y_, float z_) {
			set(w_, x_, y_, z_);
		}
		inline Quaternion(float w_, const Vec3& v) {
			set(w_, v.x, v.y, v.z);
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
			return Quaternion(-w, -ijk.x, -ijk.y, -ijk.z);
		}

		/// Multiply a two quaternions - using the right-hand rule 
		/// 2022-02-12 Scott edit. Worried that Umer uncovered a bug in my code,
		/// I derived the multiply over again (this time less sexy) and it seems to work 
		/// correctly.

		// 2022-02-12 Umer Noor edit. I think there is a bug here
		// I'll change the Vec3 on the stack to be ijk_result and
		// see if that helps. Compiler might be grabbing the ijk member
		// variable in the return line rather than the stack variable
		inline const Quaternion operator * (const Quaternion& q) const {
			Quaternion result;
			result.w = (w * q.w) - (ijk.x * q.ijk.x) - (ijk.y * q.ijk.y) - (ijk.z * q.ijk.z);
			result.ijk.x = (w * q.ijk.x) + (ijk.x * q.w) + (ijk.z * q.ijk.y) - (ijk.y * q.ijk.z);
			result.ijk.y = (w * q.ijk.y) + (ijk.y * q.w) + (ijk.x * q.ijk.z) - (ijk.z * q.ijk.x);
			result.ijk.z = (w * q.ijk.z) + (ijk.z * q.w) + (ijk.y * q.ijk.x) - (ijk.x * q.ijk.y);
			return result;
			/// This appears to be broken
			///Vec3 ijk(w * q.ijk + q.w * ijk + VMath::cross(ijk, q.ijk));
			/// Vec3 ijk_result((w * q.ijk) + (q.w * ijk) + VMath::cross(ijk, q.ijk));
			///return Quaternion(w * q.w - VMath::dot(ijk, q.ijk), ijk_result);
		}


		inline const Quaternion operator + (const Quaternion q) const {
			return Quaternion(w + q.w, ijk.x + q.ijk.x, ijk.y + q.ijk.y, ijk.z + q.ijk.z);
		}

		inline const Quaternion operator - (const Quaternion q) const {
			return Quaternion(w - q.w, ijk.x - q.ijk.x, ijk.y - q.ijk.y, ijk.z - q.ijk.z);
		}

		inline const Quaternion operator * (const float scalar) const {
			return Quaternion(w * scalar, ijk.x * scalar, ijk.y * scalar, ijk.z * scalar);
		}

		inline const Quaternion operator / (const float scalar) const {
			return Quaternion(w / scalar, ijk.x / scalar, ijk.y / scalar, ijk.z / scalar);
		}

		/// Now we can use the Quaternion like an array but we'll need two overloads
		inline const float operator [] (int index) const {  /// This one is for reading the Quaternion as if where an array
			return *((float*)this + index);
		}

		inline float& operator [] (int index) {	/// This one is for writing to the Quaternion as if where an array.  
			return *((float*)this + index);
		}

		inline void print(const char* comment = nullptr) {
			if (comment) printf("%s\n", comment);
			printf("%1.8f %1.8f %1.8f %1.8f\n", w, ijk.x, ijk.y, ijk.z);
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

		/// Seriously, the tilde ~ in the complement operator not the 
		///  conjugate - but it was for fun. 
		inline Quaternion operator~() { return Quaternion(w, -ijk); }
		/////////////////////////////////////////////////////////////////////////


	};
}
#endif
