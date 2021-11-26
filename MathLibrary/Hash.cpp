#include "Hash.h"



/// Custom specialization of std::hash will be injected into namespace std
namespace std {
	using namespace MATH;
	size_t hash<Vec2>::operator() (Vec2 const& v) const  {
		size_t seed = 0;
		hash<float> hasher;
		combineHashes(seed, hasher(v.x));
		combineHashes(seed, hasher(v.y));
		return seed;
	}

	size_t hash<Vec3>::operator() (Vec3 const & v) const  {
		size_t seed = 0;
		hash<float> hasher;
		combineHashes(seed, hasher(v.x));
		combineHashes(seed, hasher(v.y));
		combineHashes(seed, hasher(v.z));
		return seed;
	}

	size_t hash<Vec4>::operator() (Vec4 const & v) const {
		size_t seed = 0;
		hash<float> hasher;
		combineHashes(seed, hasher(v.x));
		combineHashes(seed, hasher(v.y));
		combineHashes(seed, hasher(v.z));
		combineHashes(seed, hasher(v.w));
		return seed;
	}

}