#pragma once
#include <functional>
#include "Vector.h"

namespace MATH {

	/// https://en.cppreference.com/w/cpp/utility/hash
	/// This also defined in the Boost library as hash_combine()
	/// I renamed it to avoid collisions
	static void combineHashes(size_t& seed, size_t hash) {
		hash += 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hash;
	}

	/// I'm inserting these functions into namespace MATH here instead of 
	/// in Vector.h because the == operator has no place floating point math. 
	/// The == operator is true only if the binary patterns are exactly identical 
	/// this would never happen in any real math calculations.
	/// These overloads are used for hashing purposing only.
	/// 
	/// Tricky bit is, std::map and std::unordered_map need the == operator to be overloaded
	/// but notice it needs to be a two arg function; therefore, it cannot be part of the Vec2, 
	/// Vec3 or the Vec4 classes it needs to live outside the struct but within the namespace. 
	inline constexpr bool operator == (Vec2 const& v1, Vec2 const& v2) {
		return (v1.x == v2.x) && (v1.y == v2.y);
	}

	inline constexpr bool operator == (Vec3 const& v1, Vec3 const& v2) {
		return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
	}

	inline constexpr bool operator == (Vec4 const& v1, Vec4 const& v2) {
		return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z) && (v1.w == v2.w);
	}
}


/// These custom specializations of std::hash will be inserted into the namespace std
/// https://en.cppreference.com/w/cpp/utility/hash

namespace std {
	using namespace MATH;

	template<> struct hash<Vec2> {
		size_t operator() (Vec2 const& v) const {
			size_t seed = 0;
			hash<float> hasher;
			combineHashes(seed, hasher(v.x));
			combineHashes(seed, hasher(v.y));
			return seed;
		}
	};

	template<> struct hash<Vec3> {
		size_t operator() ( Vec3 const& v) const{
			size_t seed = 0;
			hash<float> hasher;
			combineHashes(seed, hasher(v.x));
			combineHashes(seed, hasher(v.y));
			combineHashes(seed, hasher(v.z));
			return seed;
		}
	};

	template<> struct hash<Vec4> {
		size_t operator() (Vec4 const& v) const{
			size_t seed = 0;
			hash<float> hasher;
			combineHashes(seed, hasher(v.x));
			combineHashes(seed, hasher(v.y));
			combineHashes(seed, hasher(v.z));
			combineHashes(seed, hasher(v.w));
			return seed;
		}
	};
}