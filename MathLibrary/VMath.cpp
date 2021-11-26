#include "VMath.h"
#include <assert.h>
using namespace MATH;
 
/// Return a normalized Vec3
Vec3 VMath::normalize(const Vec3 &a) {
		float magnitude = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
#ifdef _DEBUG  /// If in debug mode let's worry about divide by zero or nearly zero!!! 	
	if (magnitude < VERY_SMALL) {
		std::string errorMsg = __FILE__ + __LINE__;
		throw errorMsg.append(": Divide by nearly zero! ");
	}
#endif
	return Vec3(a.x / magnitude, a.y / magnitude, a.z / magnitude);
}

Vec3 VMath::reflect(const Vec3 &v, const Vec3 &n){
	float lamda = -2.0f * VMath::dot(n, v);
	return v + lamda * n;
}

float VMath::distance(const Vec3 &a, const Vec3 &b){
	Vec3 r  = a - b;
	return(mag(r));
}

Vec3 VMath::lerp(const Vec3 &v1, const Vec3 &v2, float t) {
	return (v1 + t * (v2 - v1));
}
