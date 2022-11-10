#ifndef PLANE_H
#define PLANE_H
#include <iostream>
#include "VMath.h"

///http://www.songho.ca/math/plane/plane.html
///http://mathworld.wolfram.com/HessianNormalForm.html
namespace  MATH {
	struct Plane : public Vec3 {
		float d;

		/// Just a little utility to populate a Plane
		inline void set(float x_, float y_, float z_, float d_) {
			x = x_; y = y_; z = z_, d = d_;
		}

	  
		/// Here's a set of constructors:
		///http://www.songho.ca/math/plane/plane.html
		///http://mathworld.wolfram.com/HessianNormalForm.html
		/// If the plane is defined by a normal (that should normalized), 
		/// the equation of a Plane is Ax + By + Cz - D = 0; where 
		/// A,B,C are the values of the normal n <A,B,C> and D is the 
		/// distance from the origin to the plane.

		/// If the the Vec3 is not normalized call PMath::normalize() to fix it  
		Plane(Vec3 n, float d) {
			set(n.x, n.y, n.z, -d);
		}

		Plane() {
			set(0.0f, 0.0f, 0.0f, 0.0f);
		}

		/// If the plane is defined by three points v0, v1, v2 
		/// then the equation of a Plane is Ax + By + Cz + D = 0.
		Plane(const Vec3& v0, const Vec3& v1, const Vec3& v2) {
			Vec3 a = v1 - v0;
			Vec3 b = v2 - v0;
			Vec3 n = VMath::normalize(VMath::cross(a, b));
			set(n.x, n.y, n.z, VMath::dot(n, v0));
		}

		/// These just set numbers - be careful. This may not be in the Hessian normal form.  
		/// Call PMath::normalize() to fix it
		Plane(float x, float y, float z, float d) {
			set(x, y, z, d);
		}
		/// A copy constructor
		Plane(const Plane& p) {
			set(p.x, p.y, p.z, p.d);
		}

		/// An assignment operator   
		inline Plane& operator = (const Plane& p) {
			set(p.x, p.y, p.z, p.d);
			return *this;
		}

		
		/// print the values of the plane and add a comment if you wish
		void print(const char* comment = nullptr) const {
			if (comment) printf("%s\n", comment);
			printf("%f %f %f %f\n", x, y, z, d);
		}
	
	};
}


#endif