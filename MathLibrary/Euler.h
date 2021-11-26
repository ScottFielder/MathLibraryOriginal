#ifndef EULER_H
#define EULER_H
#include <iostream> 
namespace MATH {
	union Euler {

		struct {
			float xAxis, yAxis, zAxis; 
		};
		struct {
			float pitch, yaw, roll;
		};
		struct {
			float pitch, heading, bank;
		};

	public:
		/// Just a little utility to populate a Euler
		inline void set(float xAxis_, float yAxis_, float zAxis_) {
			xAxis = xAxis_; yAxis = yAxis_; zAxis = zAxis_;
		}

		inline Euler(float xAxis, float yAxis, float zAxis) {
			set(xAxis, yAxis, zAxis);
		}

		/// Here's a set of constructors
		inline Euler() {
			set(0.0f, 0.0f, 0.0f);
		}

		/// A copy constructor
		inline Euler(const Euler& e) {
			set(e.xAxis, e.yAxis, e.zAxis);
		}

		/// An assignment operator   
		inline Euler& operator = (const Euler& e) {
			set(e.xAxis, e.yAxis, e.zAxis);
			return *this;
		}

		/// Multiply a Euler by a scalar
		inline const Euler  operator * (const float s) const {
			return Euler(s*xAxis, s*yAxis, s*zAxis);
		}

		/// Multiply a scaler by a Vec3   It's the scalar first then the Vec3
		inline friend Euler operator * (const float s, const Euler& e) {
			return e * s;
		}

		/// Multiply a Vec3 by a scalar and assign it to itself
		inline Euler& operator *= (const float s) {
			xAxis *= s;
			yAxis *= s;
			zAxis *= s;
			return *this;
		}

		/// An assignment operator   
		inline void print (const char* comment = nullptr) {
			if (comment) printf("%s\n", comment);
			printf("%1.8f %1.8f %1.8f\n", xAxis, yAxis, zAxis);
		}
	};
}
#endif