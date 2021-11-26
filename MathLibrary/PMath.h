#ifndef PMATH_H
#define PMATH_H
#include "Plane.h"
namespace MATH {
	class PMath {
	public:
		static Plane normalize(const Plane &p);

		/// Get the distance form a point (Vec3) to a plane
		static float distance(const Vec3 &v, const Plane &p);
		static Vec3 reflect(const Vec3 &v, const Plane &p);
	};
}
#endif
