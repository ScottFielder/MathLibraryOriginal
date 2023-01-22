#include "PMath.h"
#include "VMath.h"
using namespace MATH;

/// Creates the Hessian normal form of the plane 
Plane PMath::normalize(const Plane &p) {
	float mag = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
	return Plane(p.x / mag, p.y / mag, p.z / mag, p.d / mag);
}

/// Assuming the Hessian normal form of the plane 
float PMath::distance(const Vec3 &v, const Plane &p) {
	Vec3 n = p; /// Extract the normal from the plane
	return VMath::dot(n, v) + p.d;
}


/// Reflect an  incident vector off the normal of a plane 
Vec3 PMath::reflect(const Vec3 &v, const Plane &p) {
	Vec3 n = p; /// Extract the normal from the plane
	return v - (2.0f * VMath::dot(n, v)) * n;
}
