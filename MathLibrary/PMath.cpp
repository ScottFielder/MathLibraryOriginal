#include "PMath.h"
#include "VMath.h"
using namespace MATH;


Plane PMath::normalize(const Plane &p) {
	float mag = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
	return Plane(p.x / mag, p.y / mag, p.z / mag, 
		p.d / mag);
}


float PMath::distance(const Vec3 &v, const Plane &p) {
	Vec3 n = p;
	float mag = VMath::mag(n);
	return (p.x*v.x + p.y*v.y + p.z*v.z + p.d) / mag;
}


/// Reflect an  incident vector off the normal of a plane 
Vec3 PMath::reflect(const Vec3 &v, const Plane &p) {
	float lamda = 2.0f * VMath::dot(p, v);
	return v + lamda * p;
}
