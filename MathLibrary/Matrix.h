#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>	
#include "Vector.h"

namespace  MATH {
	/// These are the default vectors of the eye (camera) according to the right hand rule
	///								+Y   -Z
	///	                             |  /
	///   origin(0.0f,0.0f,0.0f);    | /
	///   up(0.0f,1.0f,0.0f);        |/_____+X
	///   forward(0.0f,0.0f,-1.0f);  /
	///                             /
	///                           +Z
	///
	class Matrix4 {

		/// Let's just make sure that all is clear about how this matrix is layed out. 

		/// How a matrix is really layed out is pretty much abitrary but we need to agree
		/// and the world has agreed (except for Microsoft) on the right-hand rule. 
		/// First, the 4x4 matrix is really just an array of 16 numbers.  
		/// We need to think of the array as a matrix in the following way
		/// 4x4 matrix - COLUMN MAJOR - The OpenGL, science, physics, mathematics and engineering way. 
		///	0  4  8  12        [0][0]  [1][0]  [2][0]  [3][0]   
		///	1  5  9  13  (or)  [0][1]  [1][1]  [2][1]  [3][1]
		///	2  6  10 14        [0][2]  [1][2]  [2][2]  [3][2]
		///	3  7  11 15        [0][3]  [1][3]  [2][3]  [3][3]
	private:
		float  m[16];
	public:

		inline Matrix4& operator = (const Matrix4& m_) {
			this->m[0] = m_[0]; this->m[1] = m_[1]; this->m[2] = m_[2]; this->m[3] = m_[3];
			this->m[4] = m_[4]; this->m[5] = m_[5]; this->m[6] = m_[6]; this->m[7] = m_[7];
			this->m[8] = m_[8]; this->m[9] = m_[9]; this->m[10] = m_[10]; this->m[11] = m_[11];
			this->m[12] = m_[12]; this->m[13] = m_[13]; this->m[14] = m_[14]; this->m[15] = m_[15];
			return *this;
		}


		inline Matrix4(float x0, float x1, float x2, float x3,
			float y0, float y1, float y2, float y3,
			float z0, float z1, float z2, float z3,
			float w0, float w1, float w2, float w3) {
			m[0] = x0;   m[4] = y0;   m[8] = z0; m[12] = w0;
			m[1] = x1;   m[5] = y1;   m[9] = z1; m[13] = w1;
			m[2] = x2;   m[6] = y2;   m[10] = z2; m[14] = w2;
			m[3] = x3;   m[7] = y3;   m[11] = z3; m[15] = w3;
		}

		/// Create the unit matrix probably the most common way of initializing a matrix
		/// (If the argument within the constructor is pre-defined, as in this case 1.0f, the 
		/// zero argument constructor will use 1.0f as the argument. Anything other than
		/// 1.0f will be assigned to the float d)
		inline Matrix4() {
			loadIdentity();
		}


		/// Creates the identity matrix
		inline void loadIdentity() {
			m[0] = 1.0f;   m[4] = 0.0f;   m[8] = 0.0f;   m[12] = 0.0f;
			m[1] = 0.0f;   m[5] = 1.0f;   m[9] = 0.0f;   m[13] = 0.0f;
			m[2] = 0.0f;   m[6] = 0.0f;   m[10] = 1.0f;   m[14] = 0.0f;
			m[3] = 0.0f;   m[7] = 0.0f;   m[11] = 0.0f;   m[15] = 1.0f;
		}

		/// Multiply two 4x4 matricies. 
		/// Grrr, I never liked mulipling maticies - but it needs to be done. 
		/// Tested 3/7/2017 SSF

		inline const Matrix4 operator*(const Matrix4& n) const {
			
			/// This approach is about 8 nanoseconds faster, not because I unrolled the loops but because of the constructor, ask me. 
			return Matrix4(
				(m[0 * 4 + 0] * n[0 * 4 + 0]) + (m[1 * 4 + 0] * n[0 * 4 + 1]) + (m[2 * 4 + 0] * n[0 * 4 + 2]) + (m[3 * 4 + 0] * n[0 * 4 + 3]),
				(m[0 * 4 + 1] * n[0 * 4 + 0]) + (m[1 * 4 + 1] * n[0 * 4 + 1]) + (m[2 * 4 + 1] * n[0 * 4 + 2]) + (m[3 * 4 + 1] * n[0 * 4 + 3]),
				(m[0 * 4 + 2] * n[0 * 4 + 0]) + (m[1 * 4 + 2] * n[0 * 4 + 1]) + (m[2 * 4 + 2] * n[0 * 4 + 2]) + (m[3 * 4 + 2] * n[0 * 4 + 3]),
				(m[0 * 4 + 3] * n[0 * 4 + 0]) + (m[1 * 4 + 3] * n[0 * 4 + 1]) + (m[2 * 4 + 3] * n[0 * 4 + 2]) + (m[3 * 4 + 3] * n[0 * 4 + 3]),
				(m[0 * 4 + 0] * n[1 * 4 + 0]) + (m[1 * 4 + 0] * n[1 * 4 + 1]) + (m[2 * 4 + 0] * n[1 * 4 + 2]) + (m[3 * 4 + 0] * n[1 * 4 + 3]),
				(m[0 * 4 + 1] * n[1 * 4 + 0]) + (m[1 * 4 + 1] * n[1 * 4 + 1]) + (m[2 * 4 + 1] * n[1 * 4 + 2]) + (m[3 * 4 + 1] * n[1 * 4 + 3]),
				(m[0 * 4 + 2] * n[1 * 4 + 0]) + (m[1 * 4 + 2] * n[1 * 4 + 1]) + (m[2 * 4 + 2] * n[1 * 4 + 2]) + (m[3 * 4 + 2] * n[1 * 4 + 3]),
				(m[0 * 4 + 3] * n[1 * 4 + 0]) + (m[1 * 4 + 3] * n[1 * 4 + 1]) + (m[2 * 4 + 3] * n[1 * 4 + 2]) + (m[3 * 4 + 3] * n[1 * 4 + 3]),
				(m[0 * 4 + 0] * n[2 * 4 + 0]) + (m[1 * 4 + 0] * n[2 * 4 + 1]) + (m[2 * 4 + 0] * n[2 * 4 + 2]) + (m[3 * 4 + 0] * n[2 * 4 + 3]),
				(m[0 * 4 + 1] * n[2 * 4 + 0]) + (m[1 * 4 + 1] * n[2 * 4 + 1]) + (m[2 * 4 + 1] * n[2 * 4 + 2]) + (m[3 * 4 + 1] * n[2 * 4 + 3]),
				(m[0 * 4 + 2] * n[2 * 4 + 0]) + (m[1 * 4 + 2] * n[2 * 4 + 1]) + (m[2 * 4 + 2] * n[2 * 4 + 2]) + (m[3 * 4 + 2] * n[2 * 4 + 3]),
				(m[0 * 4 + 3] * n[2 * 4 + 0]) + (m[1 * 4 + 3] * n[2 * 4 + 1]) + (m[2 * 4 + 3] * n[2 * 4 + 2]) + (m[3 * 4 + 3] * n[2 * 4 + 3]),
				(m[0 * 4 + 0] * n[3 * 4 + 0]) + (m[1 * 4 + 0] * n[3 * 4 + 1]) + (m[2 * 4 + 0] * n[3 * 4 + 2]) + (m[3 * 4 + 0] * n[3 * 4 + 3]),
				(m[0 * 4 + 1] * n[3 * 4 + 0]) + (m[1 * 4 + 1] * n[3 * 4 + 1]) + (m[2 * 4 + 1] * n[3 * 4 + 2]) + (m[3 * 4 + 1] * n[3 * 4 + 3]),
				(m[0 * 4 + 2] * n[3 * 4 + 0]) + (m[1 * 4 + 2] * n[3 * 4 + 1]) + (m[2 * 4 + 2] * n[3 * 4 + 2]) + (m[3 * 4 + 2] * n[3 * 4 + 3]),
				(m[0 * 4 + 3] * n[3 * 4 + 0]) + (m[1 * 4 + 3] * n[3 * 4 + 1]) + (m[2 * 4 + 3] * n[3 * 4 + 2]) + (m[3 * 4 + 3] * n[3 * 4 + 3]));
			/*** This is how you do it with loops
				Matrix4 result;
				for (int i = 0; i < 4; ++i) {
					for (int j = 0; j < 4; ++j) {
						result[i*4+j] = (m[0*4+j] * n[i*4 +0]) + (m[1*4+j] * n[i*4+1]) + (m[2*4+j] * n[i*4+2]) + (m[3*4+j] * n[i*4+3]);
					}
				}
				return result;
			***/
		
		
		}

		/// Multipling a matrix by itself is probably the most commom
		/// ("this" is the address of the matrix object itself. "*this" de-references that address
		inline Matrix4& operator*=(const Matrix4& n) {
			*this = *this * n;
			return *this;
		}

		/// 2022 June, Multiply a Vec4 by this matrix and return the resulting Vec4, 
		/// removed the divide by w in the result. 
		inline  Vec4 operator* (const Vec4& v) const {
			float x = v.x * m[0] + v.y * m[4] + v.z * m[8] + v.w * m[12];
			float y = v.x * m[1] + v.y * m[5] + v.z * m[9] + v.w * m[13];
			float z = v.x * m[2] + v.y * m[6] + v.z * m[10] + v.w * m[14];
			float w = v.x * m[3] + v.y * m[7] + v.z * m[11] + v.w * m[15];
			return Vec4(x, y, z, w);
		}

		/// Multiply a Vec3 by this matrix and return the resulting Vec3
		/// Mathematicians would say this is impossible but this is just 
		/// code.  I will assume the w-component of the Vec3 is 1.0.
		inline  Vec3 operator* (const Vec3& v) const {
			float x = v.x * m[0] + v.y * m[4] + v.z * m[8] + 1.0f * m[12];
			float y = v.x * m[1] + v.y * m[5] + v.z * m[9] + 1.0f * m[13];
			float z = v.x * m[2] + v.y * m[6] + v.z * m[10] + 1.0f * m[14];
			float w = v.x * m[3] + v.y * m[7] + v.z * m[11] + 1.0f * m[15];
			return Vec3(x, y, z);
		}

		/// Now I can use the structure itself as an array.
		/// When overloading the [] operator you need to declair one
		/// to read the array and one to write to the array. 
		///  Returns a const - the rvalue
		inline const float operator [] (int index) const {
			return *(m + index);
		}

		/// This one is for writing to the structure as if where an array 
		/// it returns a modifiable lvalue
		inline float& operator [] (int index) {
			return *(m + index);
		}

		/// These allows the convertion from type Matrix to const float * without issues
		inline operator float* () { return static_cast<float*>(&m[0]); }
		inline operator const float* () const { return static_cast<const float*>(&m[0]); }

		/// Print the Matrix4 in column form (right-hand rule), add a comment if you wish 
		inline void print(const char* comment = nullptr) const {
			if (comment) printf("%s\n", comment);
			printf("%1.4f %1.4f %1.4f %1.4f\n%1.4f %1.4f %1.4f %1.4f \n%1.4f %1.4f %1.4f %1.4f \n%1.4f %1.4f %1.4f %1.4f \n\n",
				m[0], m[4], m[8], m[12],
				m[1], m[5], m[9], m[13],
				m[2], m[6], m[10], m[14],
				m[3], m[7], m[11], m[15]);
		}

		/// These functions set and return the columns of a Matrix4
		/// Since the library is right-handed, the colunms are refered to 
		/// as basis vectors.  
		enum  Colunm { zero = 0, one, two, three };
		inline Vec4 getColumn(Colunm index) {
			return Vec4(m[4 * index + 0], m[4 * index + 1], m[4 * index + 2], m[4 * index + 3]);
		}

		inline void setColumn(Colunm index, const Vec4 v) {
			m[4 * index + 0] = v[0];
			m[4 * index + 1] = v[1];
			m[4 * index + 2] = v[2];
			m[4 * index + 3] = v[3];
		}

	};

	/********************************************************************************************/
	class Matrix3 {
		/// 3x3 matrix - COLUMN MAJOR 
		///	0	3	6			
		///	1	4	7	 
		///	2	5	8			

	private:
		float  m[9];
	public:
		/// Constuctors
		inline Matrix3(float xx, float yx, float zx,
			float xy, float yy, float zy,
			float xz, float yz, float zz) {
			m[0] = xx;   m[3] = xy;   m[6] = xz;
			m[1] = yx;   m[4] = yy;   m[7] = yz;
			m[2] = zx;   m[5] = zy;   m[8] = zz;

		}

		/// Create the unit matrix probably the most common way of initializing a matrix
		inline Matrix3() {
			loadIdentity();
		}


		/// Now I can use the class itself as an array.
		/// When overloading the [] operator you need to declair one
		/// to read the array and one to write to the array. 
		///  Returns a const rvalue
		inline const float operator [] (int index) const {
			return *(m + index);
		}

		/// This one is for writing to the class as if where an array 
		/// it returns a lvalue
		inline float& operator [] (int index) {
			return *(m + index);
		}

		inline Matrix3& operator = (const Matrix3& m_) {
			this->m[0] = m_[0]; this->m[1] = m_[1]; this->m[2] = m_[2];
			this->m[3] = m_[3]; this->m[4] = m_[4]; this->m[5] = m_[5];
			this->m[6] = m_[6]; this->m[7] = m_[7]; this->m[8] = m_[8];
			return *this;
		}


		/// Creates the identity matrix
		inline void loadIdentity() {
			m[0] = 1.0f;   m[3] = 0.0f;   m[6] = 0.0f;
			m[1] = 0.0f;   m[4] = 1.0f;   m[7] = 0.0f;
			m[2] = 0.0f;   m[5] = 0.0f;   m[8] = 1.0f;
		}

		/// Multiply a Vec3 by this matrix and return the resulting Vec3
		inline  Vec3 operator* (const Vec3& v) const {
			float x = v.x * m[0] + v.y * m[3] + v.z * m[6];
			float y = v.x * m[1] + v.y * m[4] + v.z * m[7];
			float z = v.x * m[2] + v.y * m[5] + v.z * m[8];
			return Vec3(x, y, z);
		}

		/// Multiply two 3x3 matricies. 
		/// In the Vec4 class I unroll the loop and seek efficiency everywhere,
		/// probably never matters in reality. 
		/// Here, I just do with just loops. 
		inline const Matrix3 operator*(const Matrix3& n) const {
			Matrix3 result;
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					result[i * 3 + j] = (m[0 * 3 + j] * n[i * 3 + 0]) + (m[1 * 3 + j] * n[i * 3 + 1]) + (m[2 * 3 + j] * n[i * 3 + 2]);
				}
			}
			return result;
		}

		inline Matrix3& operator*=(const Matrix3& n) {
			*this = *this * n;
			return *this;
		}

		/// Extracts the inner 3x3 from a 4x4 matrix
		/// using the 3x3 constructor
		inline Matrix3(const Matrix4& m_) {
			m[0] = m_[0]; m[1] = m_[1]; m[2] = m_[2];
			m[3] = m_[4]; m[4] = m_[5]; m[5] = m_[6];
			m[6] = m_[8]; m[7] = m_[9]; m[8] = m_[10];
		}

		/// Extracts the inner 3x3 from a 4x4 matrix
		/// The upper 3x3 of a 4x4 contains only scale and rotation components. 
		/// This will prove useful in many cases. 
		/// This uses the assignment operator, inspite of my tenet of "least astonishment" 
		inline Matrix3& operator = (const Matrix4& m_) {
			m[0] = m_[0]; m[1] = m_[1]; m[2] = m_[2];
			m[3] = m_[4]; m[4] = m_[5]; m[5] = m_[6];
			m[6] = m_[8]; m[7] = m_[9]; m[8] = m_[10];
			return *this;
		}


		/// These allow me convert from type Matrix to const float * without issues
		inline operator float* () { return static_cast<float*>(&m[0]); }
		inline operator const float* () const { return static_cast<const float*>(&m[0]); }


		/// These functions set and return the columns of a Matrix3
		/// Since the library is right-handed, the colunms are refered to 
		/// as basis vectors.  
		enum class Column : int { zero = 0, one, two };
		inline Vec3 getColumn(Column index) {
			return Vec3(m[3 * (int)index + 0], m[3 * (int)index + 1], m[3 * (int)index + 2]);
		}

		inline void setColumn(Column index, const Vec3 v) {
			m[3 * (int)index + 0] = v[0];
			m[3 * (int)index + 1] = v[1];
			m[3 * (int)index + 2] = v[2];
		}

		/// Print the Matrix3 in column form (right-hand rule), add a comment if you wish
		inline void print(const char* comment = nullptr) const {
			if (comment) printf("%s\n", comment);
			printf("%1.4f %1.4f %1.4f\n%1.4f %1.4f %1.4f\n%1.4f %1.4f %1.4f\n\n",
				m[0], m[3], m[6],
				m[1], m[4], m[7],
				m[2], m[5], m[8]);
		}

	};
}

#endif
