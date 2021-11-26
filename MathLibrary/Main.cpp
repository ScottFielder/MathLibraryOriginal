/// This is the test bed for my math library - SSF
/// This is currently version 1.20 - SSF June 2021

#include <iostream>
#include <fstream>
#include "MMath.h"
#include "QMath.h"
#include "EMath.h"
#include "PMath.h"
#include "AAMath.h"
#include "Sphere.h"
#include "Fourier.h"
#include "Randomizer.h"
#include "Hash.h"


#include <glm/vec3.hpp> /// glm::vec3
#include <glm/vec4.hpp> /// glm::vec4, glm::ivec4
#include <glm/mat4x4.hpp> /// glm::mat4
#include <glm/gtc/matrix_transform.hpp> /// glm::translate, glm::rotate, glm::scale, glm::perspective
#include <glm/gtc/type_ptr.hpp> /// glm::value_ptr
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/hash.hpp>


void FFT_Test();
void LookAtTest();
void inverseTest();
void UnOrthoTest();
void RotationTest();
void Vec3MultiplyMat4Test();
void multiplyMatrixTest();
void viewportNDCTest();
void moveCopyConstructors();
void rotationIsOrthogonal();
void planeTest();
void quaternionTest();
void hashTest();
void determinantTest();

/// Utility print() calls for glm to Scott's math library format 
void glmPrintM4(glm::mat4  mat, const char* comment = nullptr);
void glmPrintQ(glm::quat q, const char* comment = nullptr);
void glmPrintV3(glm::vec3 v, const char* comment = nullptr);
void glmPrintV4(glm::vec4 v, const char* comment = nullptr);


using namespace MATH;
using namespace glm;
using namespace std;


int main(int argc, char*argv[]) {


	
}

void determinantTest(){
	/// These vectors should return a value of 30 - it does.
	/// Swap any two and the sign should change - it does
	Matrix4 m4;
	m4.setColunm(Matrix4::Colunm::zero,Vec4(1,0,2,-1));
	m4.setColunm(Matrix4::Colunm::one, Vec4(3, 0, 0, 5));
	m4.setColunm(Matrix4::Colunm::two, Vec4(2, 1, 4, -3));
	m4.setColunm(Matrix4::Colunm::three,Vec4(1, 0, 5, 0));
	
	printf("%f\n",MMath::determinate(m4));

	/// deternimant of the identity matrix = 1.0 - it is 
	Matrix3 m3;
	printf("%f\n", MMath::determinate(m3));
}

void hashTest(){
	Vec3 v1(1.1f, 1.0f, 1.0f);
	Vec3 v2(1.1f, 0.0f, 1.0f);
	Vec3 v3(1.1f, 0.0f, 1.0f);
	Vec3 v4(1.1f, 1.0f, 1.0f);

	vec3 glmV1(1.1f, 1.0f, 1.0f);
	vec3 glmV2(1.1f, 1.0f, 1.0f);
	bool t = (v1 == v2);
	bool glmt = (glmV1 == glmV2);
	hash<Vec3> hasher;
	size_t myHash = hasher(v1);
	
	hash<vec3> glmHasher;
	size_t glmHash = glmHasher(glmV1);
	printf("%x vs. %x\n", myHash, glmHash);

	std::unordered_map<Vec3, uint32_t> uniqueVerts;
	static uint32_t count = 0;
	if (uniqueVerts.count(v1) == 0) {
		uniqueVerts[v1] = count;
		++count;
	}
	printf("%d\n", uniqueVerts.size());

	if (uniqueVerts.count(v2) == 0) {
		uniqueVerts[v2] = count;
		++count;
	}
	printf("%d\n", uniqueVerts.size());

	if (uniqueVerts.count(v3) == 0) {
		uniqueVerts[v3] = count;
		++count;
	}
	printf("%d\n", uniqueVerts.size());


	if (uniqueVerts.count(v4) == 0) {
		uniqueVerts[v4] = count;
		++count;
	}
	printf("%d\n\n", uniqueVerts.size());

	/////////// glm ///////////////
	std::unordered_map<vec3, uint32_t> glmUniqueVerts;
	static uint32_t glmCount = 0;
	if(glmUniqueVerts.count(glmV1) == 0){
		glmUniqueVerts[glmV1] = glmCount;
		++glmCount;
	}
	printf("%d\n", glmUniqueVerts.size());
	if (glmUniqueVerts.count(glmV2) == 0) {
		glmUniqueVerts[glmV2] = glmCount;
		++glmCount;
	}
	printf("%d\n", glmUniqueVerts.size());


	

}

void quaternionTest() {
	/// glm version 
	vec3 eulerAngles(radians(90.0), radians(45.0), radians(0.0));
	quat myQuaternion = quat(eulerAngles);
	glmPrintQ(myQuaternion, "glm Quat");

	/// my version 
	
	Euler e9(90.0, 45.0, 0.0);
	Quaternion q3 = QMath::fromEuler(e9);
	q3.print("my Quat");

	eulerAngles = glm::eulerAngles(myQuaternion);
	glmPrintV3(eulerAngles, "glm from quat");

	e9 = QMath::fromQuaternion(q3);
	e9.print("from Quat");



	mat4 glmRot = glm::toMat4(myQuaternion);
	glmPrintM4(glmRot,"glm rotation Matrix");
	

	
	Matrix4 meMat = QMath::toMatrix4(q3);
	meMat.print("My matrix");


	/// Lets say I have a unit vector along the x-axis 1,0,0. 
	/// Rotate 45.0 degree around the z-axis
	/// The resulting vector should be 0.70711, 0.70711, 0.0 
	/// Let's test this in every way I can think of
	Vec3 v(1.0, 0.0, 0.0);
	Quaternion q = QMath::angleAxisRotation(90.0,Vec3(0.0,1.0,0.0));
	Euler e2 = QMath::fromQuaternion(q);
	e2.print("from Q");
	
	q.print("The rotation Quaternion");
	Euler e(0.0, 0.0, 45.0);

	Quaternion qe = QMath::fromEuler(e);
	
	
	

	qe.print("from Euler");
	Quaternion p(0.0, v);
	Vec3 v2 = q * v * ~q;
	v2.print("The slow way");



	Vec3 v3 = QMath::rotate(v, q);
	v3.print("faster way");

	Matrix3 m3 = QMath::toMatrix3(q);
	Vec3 v4 = m3 * v;
	v4.print("Mat3");

	Matrix4 m4 = QMath::toMatrix4(q);
	Vec3 v5 = m4 * v;
	v5.print("Mat4");

	Quaternion q2 = QMath::angleAxisRotation(90.0, Vec3(0.0, 0.0, 1.0));
	q2 = QMath::pow(q2, 0.5);
	Vec3 v6 = QMath::rotate(v, q2);
	v6.print("Using the pow function");

	printf("Magnitude of q \n%f\n", QMath::magnitude(q));
	Quaternion conj_q = QMath::conjugate(q);
	conj_q.print("conjugate of q");

	Quaternion inv_q = QMath::inverse(q);
	inv_q.print("inv of q");

	Quaternion q4 = q * inv_q;
	q4.print("q * q-1 is the identity");
}


void inverseTest(){
	Matrix4 rot = MMath::rotate(90.0f, Vec3(0.0f,1.0f,0.0f));
	Matrix4 invRot = MMath::inverse(rot);
	Matrix4 product = rot * invRot;
	product.print();

	Matrix3 rot3 = MMath::rotate(45.0f, Vec3(0.0f, 1.0f, 0.0f));
	Matrix3 invRot3 = MMath::inverse(rot);
	Matrix3 product3 = rot * invRot;
	product3.print();
}


void planeTest() {
	Plane p1(2.0f, -2.0f, 5.0f, 8.0f);
	p1.print();
	Vec3 v = Vec3(4.0f, -4.0f, 3.0f);
	v.print();
	float distance = PMath::distance(v, p1);
	printf("%f vs. 6.79\n", distance);
	Plane p2 = PMath::normalize(p1);
	p2.print();
	float distance2 = PMath::distance(v, p2);
	printf("%f vs. 6.79\n", distance2);

	Plane p3(1, 2, 2, -6);
	Vec3 v3(-1, -2, -3);
	float distance3 = PMath::distance(v3, p3);
	printf("%f vs. %f\n", distance3, -17.0/3.0);
	Plane p4 = PMath::normalize(p3);
	p4.print();
	float distance4 = PMath::distance(v3, p4);
	printf("%f vs. %f\n", distance4, -17.0 / 3.0);

	
	Plane p5(Vec3(1, 0, 0), 0);
	Vec3 v5(-5, 0, 0);
	float distance5 = PMath::distance(v5, p5);
	printf("%f\n", distance5);
	///Vec3 v2 =VMath::reflect(v, n);
}


void rotationIsOrthogonal() {
	Matrix4 M = MMath::rotate(180.0f, Vec3(0, 1, 0));
	Vec4 v0 = M.getColumn(Matrix4::Colunm::zero);
	Vec4 v1 = M.getColumn(Matrix4::Colunm::one);
	Vec4 v2 = M.getColumn(Matrix4::Colunm::two);
	Vec4 v3 = M.getColumn(Matrix4::Colunm::three);
	printf("%f\n", VMath::dot(v3, v0));
	printf("%f\n", VMath::dot(v3, v1));
	printf("%f\n", VMath::dot(v3, v2));
	printf("%f\n", VMath::dot(v3, v3));
	printf("If all the values are zero, the matrix is orthoganal\n");
}
void randomizerTest() {
	Randomizer r;
	ofstream myfile;
	myfile.open("data.csv");
	for (double i = 0.0; i < 512.0; i++) {
		double val = r.box_muller(0.0, 2.0);
		myfile << i << "," << val << "\n";
	}
	myfile.close();
}


void moveCopyConstructors() {
	
}

void viewportNDCTest() {
	Matrix4 m = MMath::viewportNDC(1024, 1024);
	m.print();
	Vec3 pos0(0, 0, 0);
	Vec3 result0 = m * pos0;
	result0.print();

	Vec3 pos1(-1, 1, 1);
	Vec3 result1 = m * pos1;
	result1.print();
}

void multiplyMatrixTest() {
	Matrix4 tmSSF = MMath::translate(10.0f, 10.0f, 10.0f);
	Matrix4 rmSSF = MMath::rotate(90.0f, 0.0f, 1.0f, 0.0f);
	Matrix4 smSSF = MMath::scale(0.75f, 0.75f, 0.75f);
	Matrix4 resultSSF = tmSSF * rmSSF * smSSF;
	resultSSF.print();

	glm::mat4 mt = glm::translate(glm::mat4(1.0f), glm::vec3(10.0f, 10.0f, 10.0f));
	glm::mat4 mr = glm::rotate(mat4(), glm::radians(90.0f), vec3(0, 1.0f, 0));
	glm::mat4 ms = glm::scale(mat4(), glm::vec3(0.75f, 0.75f, 0.75f));
	glm::mat4 result = mt * mr * ms;
	
	glmPrintM4(result);
}
void Vec3MultiplyMat4Test() {
	Matrix4 translate = MMath::rotate(90.0, 0.0, 1.0,0.0);
	Vec3 pos(5.0, 0.0, 0.0);
	Vec4 xxx = pos;
	xxx.print();
	Vec3 result = translate * pos;
	result.print();
}
void RotationTest(){
	mat4 rot2 = rotate(mat4(), 3.141592654f/2.0f, vec3(1.0f,0.0f,0.0f));
	float  m[16] = {0.0};

	const float *pSource = (const float*)glm::value_ptr(rot2);
	for (int i = 0; i < 16; ++i)
		m[i] = pSource[i];
	printf("%1.8f %1.8f %1.8f %1.8f\n%1.8f %1.8f %1.8f %1.8f\n%1.8f %1.8f %1.8f %1.8f\n%1.8f %1.8f %1.8f %1.8f\n\n",
				m[0], m[4], m[8],  m[12],
				m[1], m[5], m[9],  m[13],
				m[2], m[6], m[10], m[14],
				m[3], m[7], m[11], m[15]);

}
void UnOrthoTest() {
	/// Just seeing if I can deconstruct the the ortho matrix
	int w = 800, h = 600;
	Matrix4 ndc = MMath::viewportNDC(w,h);
	
	float xMax = 10.0, xMin = -10.0, yMax = 10.0, yMin = -10.0, zMax = 1.0, zMin = -10.0;
	Matrix4 ortho = MMath::orthographic(xMin, xMax, 
										yMin, yMax, 
										zMin, zMax);

	Matrix4 projection = ortho * ndc;
	projection.print();
	
	Matrix4 m;
	/// This is the ortho * ndc matrix broken down into its parts 
	Matrix4 m1 = MMath::scale(2.0f / (xMax - xMin), 2.0f / (yMax - yMin),-2.0f / (zMax - zMin));
	Matrix4 m2 = MMath::translate( -(xMax + xMin) / (xMax - xMin), -(yMax + yMin) / (yMax - yMin), -(zMax + zMin) / (zMax - zMin)); 
	Matrix4 m3 = MMath::scale(1.0f, -1.0f, 1.0f);
	Matrix4 m4 = MMath::scale(float(w)/2.0f, float(h)/2.0f, 1 - 0);
	Matrix4 m5 = MMath::translate(float(w)/2.0f,float(h)/2.0f, 0);

	/// Here they are in their inverse 
	Matrix4 m6 = MMath::inverse(m1);
	Matrix4 m7 = MMath::translate( (xMax + xMin) / (xMax - xMin), (yMax + yMin) / (yMax - yMin), (zMax + zMin) / (zMax - zMin)); 
	Matrix4 m8 = MMath::scale(1.0f, -1.0f, 1.0f);
	Matrix4 m9 = MMath::inverse(MMath::scale(float(w)/2.0f, float(h)/2.0f, 1 - 0));
	Matrix4 m10 = MMath::translate(-float(w)/2.0f,-float(h)/2.0f, 0);

	m = m1*m2*m3*m4*m5;  /// creates the ortho * ndc
	m *= m10 *m9 *m8 *m7 *m6; /// Now back out 
	m.print(); /// Should be an identity matrix
	/// It is!!!
		
	Matrix4 invOrtho = MMath::unOrtho(projection );
	invOrtho.print();
	Vec3 v1(400.0,300.0,0.0);
	Vec3 v2(10.0,0.0,10.0);
	Vec3 result1 = invOrtho  * v1;
	result1.print();	
}


void LookAtTest(){
	glm::mat4 mt = glm::lookAt(glm::vec3(0.0f, 0.0f, -10.0f),
								glm::vec3(0.0f, 0.0f, 0.0f), 
								glm::vec3(1.0f, 0.0f, 0.0f));
	
	Vec3 v = Vec3(0, 0, 0);
	Matrix4 lookat = MMath::lookAt(Vec3(0.0,0.0,-10.0), Vec3(0,0,0), Vec3(1,0,0));
	
	lookat.print();
	glmPrintM4(mt);
	Vec3 v1 = lookat * v;
	v1.print();
}


void FFT_Test(){
#define SAMPLE_SIZE 512
	FILE *fp;

	float data[2 * SAMPLE_SIZE];
	float orig_data[2 * SAMPLE_SIZE];
	float transformed[2 * SAMPLE_SIZE];

	/// Change this as you will, keep it under the Nyquist frequency (1/2*step)
	float freq = 2.0f;
	float theta = 0.0f;
	float step = 2.0f * M_PI / SAMPLE_SIZE;

	Randomizer r; /// I'll use this to create some noise

	//////////////////////////////////////////////////////////////////
	/// Create a data sample SAMPLE_SIZE long times 2 (real and imaginary components)
	for (int i = 0; i < 2 * SAMPLE_SIZE; i += 2){
		data[i] = cos(theta * freq) + 0.7f*cos(theta * freq * 3) + (float) r.box_muller(0.0, 0.5); /// real
		data[i + 1] = 0.0f;									  ///img
		theta += step;
	}
	//////////////////////////////////////////////////////////////////

	/// Just make a copy of the original data
	memcpy(orig_data, data, 2 * SAMPLE_SIZE * sizeof(float));


	/// Now do the FFT on the noisy data
	Fourier::fft(data, 2 * SAMPLE_SIZE, Fourier::DIRECTION::FORWARD);

	/// Keep a copy of the tranformed data
	memcpy(transformed, data, 2 * SAMPLE_SIZE * sizeof(float));

	//////////////////////////////////////////////////////////////////
	/// A cheezy version of a filter
	//for (int i = 0; i < 2 * SAMPLE_SIZE; i++){
	//if (abs(data[i] < 100.0f)) data[i] = 0.0f;
	//}

	//////////////////////////////////////////////////////////////////
	/// Now do the reverse transform then renormalize
	Fourier::fft(data, 2 * SAMPLE_SIZE, Fourier::DIRECTION::REVERSE);

	/// Re-normalize the data
	for (int i = 0; i < 2 * SAMPLE_SIZE; i++){
		data[i] *= 1.0f / float(SAMPLE_SIZE);
	}

	//////////////////////////////////////////////////////////////////
	/// Write it all out in files
	//////////////////////////////////////////////////////////////////
	if (fopen_s(&fp, "0.orig_data.csv", "w") != 0){
		printf("Can't open file\n");
		return;
	}
	for (int i = 0; i < 2 * SAMPLE_SIZE; i += 2){
		fprintf(fp, "%f, %f\n", orig_data[i], orig_data[i + 1]);
	}
	fclose(fp);



	if (fopen_s(&fp, "1.transformed.csv", "w") != 0){
		printf("Can't open file\n");
		return;
	}
	for (int i = 0; i < 2 * SAMPLE_SIZE; i += 2){
		fprintf(fp, "%f, %f\n", transformed[i], transformed[i + 1]);
	}
	fclose(fp);


	if (fopen_s(&fp, "2.orig&reverse.csv", "w") != 0){
		printf("Can't open file\n");
		return;
	}
	for (int i = 0; i < 2 * SAMPLE_SIZE; i += 2){
		fprintf(fp, "%f, %f, %f, %f\n", orig_data[i], orig_data[i + 1], data[i], data[i + 1]);
	}
	fclose(fp);


}
///////////////////////////////////////////////////////////////////////////////////////////////
/// These are print statements for glm - they don't have them  
///////////////////////////////////////////////////////////////////////////////////////////////
void glmPrintM4(glm::mat4  mat, const char* comment){
	
	int i, j;
	for (j = 0; j<4; j++) {
		for (i = 0; i<4; i++) {
			printf("%1.8f ", mat[i][j]);
		}
		printf("\n");
	}
}

void glmPrintQ(glm::quat q, const char* comment) {
	if (comment) printf("%s\n", comment);
	///                                    w     i     j     k
	printf("%1.8f %1.8f %1.8f %1.8f \n", q[3], q[0], q[1], q[2]);
}

void glmPrintV3(glm::vec3 v, const char* comment) {
	if (comment) printf("%s\n", comment);
	printf("%1.8f %1.8f %1.8f\n", v[0], v[1], v[2]);
}

void glmPrintV4(glm::vec4 v, const char* comment) {
	if (comment) printf("%s\n", comment);
	printf("%1.8f %1.8f %1.8f %1.8f\n", v[0], v[1], v[2], v[3]);
}