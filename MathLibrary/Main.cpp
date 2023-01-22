/// This is the test bed for my math library - SSF
/// This is currently version 1.30 - SSF Jan 2023

#include <iostream>
#include <fstream>
#include <algorithm> 

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
void slerpTest();

/// Utility print() calls for glm to Scott's math library format 
void glmPrintM4(glm::mat4  mat, const char* comment = nullptr);
void glmPrintM3(glm::mat3  mat, const char* comment = nullptr);

void glmPrintQ(glm::quat q, const char* comment = nullptr);
void glmPrintV3(glm::vec3 v, const char* comment = nullptr);
void glmPrintV4(glm::vec4 v, const char* comment = nullptr);


using namespace MATH;
using namespace glm;
using namespace std;




int main(int argc, char* argv[]) {
	planeTest();


}


void slerpTest() {
	Euler e1(90.0f, 0.0f, 0.0f);
	Quaternion q1 = QMath::toQuaternion(e1);
	q1.print("Start");

	Euler e2(0.0f, 90.0f, 0.0f);
	Quaternion q2 = QMath::toQuaternion(e2);
	q2.print("End");

	for (float t = 0.0f; t < 1.1f; t+=0.1f) {
		Quaternion q = QMath::slerp(q1, q2, t);
		q.print("slerping");
	}

}

void determinantTest(){
	/// These vectors should return a value of 30 - it does.
	/// Swap any two and the sign should change - it does
	Matrix4 m4;
	m4.setColumn(Matrix4::Colunm::zero,Vec4(1, 0, 2, -1));
	m4.setColumn(Matrix4::Colunm::one, Vec4(3, 0, 0, 5));
	m4.setColumn(Matrix4::Colunm::two, Vec4(2, 1, 4, -3));
	m4.setColumn(Matrix4::Colunm::three,Vec4(1, 0, 5, 0));
	
	printf("%f\n",MMath::determinate(m4));

	/// deternimant of the identity matrix = 1.0 - it is 
	Matrix3 m3;
	printf("%f\n", MMath::determinate(m3));


	Matrix4 perspectiveM = MMath::perspective(45.0f, (16.0f / 9.0f), 0.5f, 100.0f);
	printf("determinant of the perspective Matrix %f\n", MMath::determinate(perspectiveM));
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
	Quaternion qLookat = QMath::lookAt(Vec3(1.0f, 0.0f, 1.0f),Vec3(0.0f, 1.0f, 0.0));
	qLookat.print();
	
	glm::quat glmlookat = glm::quatLookAt(normalize(vec3(1.0f, 0.0f, 1.0)), vec3(0.0f, 1.0f, 0.0));
	glmPrintQ(glmlookat,"glm lookat");


	Matrix3 rm = Matrix3 (MMath::rotate(-270.0f, Vec3(1.0f, 0.0f, 0.0f)));
	rm.print("My rotation matrix");
	Quaternion qm = QMath::toQuaternion(rm);
	qm.print("Quaternion rotate");
	Matrix4 meMat = MMath::toMatrix4(qm);
	meMat.print("My matrix from Q");

	mat3 glmrot = mat3(rotate(glm::radians(-270.0f), vec3(1.0f,0.0f,0.0f)));
	glmPrintM4(glmrot, "GLM rot");
	glm::quat glmqm= glm::quat_cast(glmrot);
	glmPrintQ(glmqm,"glm matrix from Q");

	glm::mat4 glmRot = glm::toMat4(glmqm);
	glmPrintM4(glmRot,"glm rotation Matrix from Q");
	

	Quaternion q3 = QMath::angleAxisRotation(-90.0,Vec3(1.0f,0.0f,0.0f));
	q3.print("my Quat");
	glm::quat myQuaternion = glm::angleAxis(glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	glmPrintQ(myQuaternion, "glm Quat");


	
	/// Lets say I have a unit vector along the x-axis 1,0,0. 
	 /*Rotate 45.0 degree around the z-axis
	 The resulting vector should be 0.70711, 0.70711, 0.0 
	 Let's test this in every way I can think of*/
	Vec3 v(1.0, 0.0, 0.0);
	Quaternion q = QMath::angleAxisRotation(90.0,Vec3(0.0,1.0,0.0));
	Euler e2 = EMath::toEuler(q);
	e2.print("from Q");
	
	q.print("The rotation Quaternion");
	Euler e(0.0, 0.0, 45.0);

	Quaternion qe = QMath::toQuaternion(e);
	
	qe.print("from Euler");
	Vec3 v2 = qe * v * ~qe;
	v2.print("The slow way");



	Vec3 v3 = QMath::rotate(v, qe);
	v3.print("faster way");

	Matrix3 m3 = MMath::toMatrix3(qe);
	Vec3 v4 = m3 * v;
	v4.print("Mat3");

	Matrix4 m4 = MMath::toMatrix4(qe);
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
	/*Plane p1(2.0f, -2.0f, 5.0f, 8.0f);
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
	printf("%f vs. %f\n", distance4, -17.0 / 3.0);*/

	Matrix4 proj = MMath::perspective(52.8f, 1.0f, 1.0f, 10.0f) * MMath::lookAt(Vec3(0.0f, 0.0f, 8.0f), Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f));
	
	proj.print();

#define INDEX(column,row) (proj[(row-1) * 4 + (column-1)])

	
	Plane near(INDEX(4,1) + INDEX(3,1),
			INDEX(4,2) + INDEX(3,2),
			INDEX(4,3) + INDEX(3,3),
			INDEX(4,4) + INDEX(3,4));
	near = PMath::normalize(near);
	near.print();
	Plane far(INDEX(4,1) - INDEX(3,1),
			INDEX(4,2) - INDEX(3,2),
			INDEX(4,3) - INDEX(3,3),
			INDEX(4,4) - INDEX(3,4));
	far = PMath::normalize(far);
	far.print();
	
	Vec3 v5(0.0, 0.0, 0.0);
	float distance5 = PMath::distance(v5, near);
	printf("near %f\n", distance5);
	distance5 = PMath::distance(v5, far);
	printf("far %f\n", distance5);
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
	printf("If all the values are zero, the matrix is orthogonal\n");
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
	if (comment) printf("%s\n", comment);
	for (j = 0; j<4; j++) {
		for (i = 0; i<4; i++) {
			printf("%1.4f ", mat[i][j]);
		}
		printf("\n");
	}
}
void glmPrintM3(glm::mat3  mat, const char* comment){
	int i, j;
	if (comment) printf("%s\n", comment);
	for (j = 0; j<3; j++) {
		for (i = 0; i<3; i++) {
			printf("%1.4f ", mat[i][j]);
		}
		printf("\n");
	}
}

void glmPrintQ(glm::quat q, const char* comment) {
	if (comment) printf("%s\n", comment);
	///                                    w     i     j     k
	printf("%1.4f %1.4f %1.4f %1.4f \n", q[3], q[0], q[1], q[2]);
}

void glmPrintV3(glm::vec3 v, const char* comment) {
	if (comment) printf("%s\n", comment);
	printf("%1.4f %1.4f %1.4f\n", v[0], v[1], v[2]);
}

void glmPrintV4(glm::vec4 v, const char* comment) {
	if (comment) printf("%s\n", comment);
	printf("%1.4f %1.4f %1.4f %1.4f\n", v[0], v[1], v[2], v[3]);
}