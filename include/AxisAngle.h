#ifndef AXISANGLE_H
#define AXISANGLE_H
#include "Vector.h"
namespace MATH{
	struct AxisAngle {
		AxisAngle(float angle_, Vec3 rotAxis_): angle(angle_), axis(rotAxis_){}
		float angle;
		Vec3 axis;

		inline void print(const char* comment = nullptr) const {
			if (comment) printf("%s\n", comment);
			printf("%1.8f (%1.8f %1.8f %1.8f)\n",angle, axis.x, axis.y, axis.z);
		}
	};
}
#endif
