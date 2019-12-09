/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################

	This class provides basic operations with 4x4 matrices
	Matrix is stored in row-major order:
		/ 0 1 2 \
		| 3 4 5 |
		\ 6 7 8 /
*/
#include "matrix_4x4.h"
#include <cstring>
/*
	C'tor of the 4x4 matrix
	=> x00 First element of the matrix
	=> x02 Second element (first row, second column)
	=> x33 Last element
	=> x03-x32 Rest of the elements
*/
matrix_4x4::matrix_4x4(
	double x00,double x01,double x02,double x03,
	double x10,double x11,double x12,double x13,
	double x20,double x21,double x22,double x23,
	double x30,double x31,double x32,double x33
){
	/* Simple assign to the internal matrix */
	m_data[0][0]=x00;m_data[0][1]=x01;m_data[0][2]=x02;m_data[0][3]=x03;
	m_data[1][0]=x10;m_data[1][1]=x11;m_data[1][2]=x12;m_data[1][3]=x13;
	m_data[2][0]=x20;m_data[2][1]=x21;m_data[2][2]=x22;m_data[2][3]=x23;
	m_data[3][0]=x30;m_data[3][1]=x31;m_data[3][2]=x32;m_data[3][3]=x33;
}
/*
	C'tor of identity matrix
*/
matrix_4x4::matrix_4x4(){}
/*
	Makes translate matrix
	=> x X value of the translation
	=> y Y value of the translation
	=> z Z value of the translation
*/
void matrix_4x4::translate(double x,double y,double z){
	/* Set X component of last column (weight) to the translation value */
	m_data[0][3] += x;
	/* same with both y */
	m_data[1][3] += y;
	/* and z */
	m_data[2][3] += z;
}
/*
	Transforms input vector (3 element) by this 4x4 matrix
	Result just equals to "(M*[v;0])(1:3)" in octave/matlab
	<=> input/output 3 element vector
*/
void matrix_4x4::transform_vector(double* v){
	/* Rename data - make code shorter and easier to read */
	auto& d=m_data;
	/* Calculate X component of the transformed vector */
	double x=v[0]*d[0][0]+v[1]*d[0][1]+v[2]*d[0][2]+d[0][3];
	/* Y the same way, this is just dot product of 2. matrix row
		and input vector */
	double y=v[0]*d[1][0]+v[1]*d[1][1]+v[2]*d[1][2]+d[1][3];
	/* and Z is dot product with 3. row */
	double z=v[0]*d[2][0]+v[1]*d[2][1]+v[2]*d[2][2]+d[2][3];
	/* and we save the new X value to the I/O vector */
	v[0]=x;
	/* same with y */
	v[1]=y;
	/* and z */
	v[2]=z;
}
/*
	Converts this matrix to its inverse 
	(WARNING! does not care if singular)
*/
void matrix_4x4::inverse(){
	/* This inversion is done by adjungate matrix and
		subdeterminant method - it seems to be fast,
		even when we know its dimension apriori.
		Hopefully no one will have to debug the following... */
	/* At first, we rename the data to make code shorer and clearer */
	auto& d=m_data;
	/* The following line is subdeterminant of following elements:
		            / . . . . \
		            | . . . . |
		A2323 = det | . . # # |
		            \ . . # # /
	*/
	double A2323=d[2][2]*d[3][3]-d[2][3]*d[3][2];
	/* And the following line is subdeterminant of following elements:
		            / . . . . \
		            | . . . . |
		A1323 = det | . # . # |
		            \ . # . # /
	*/
	double A1323=d[2][1]*d[3][3]-d[2][3]*d[3][1];
	/* This way we get all 2x2 subdeterminants we need ... */
	double A1223=d[2][1]*d[3][2]-d[2][2]*d[3][1];
	double A0323=d[2][0]*d[3][3]-d[2][3]*d[3][0];
	double A0223=d[2][0]*d[3][2]-d[2][2]*d[3][0];
	double A0123=d[2][0]*d[3][1]-d[2][1]*d[3][0];
	double A2313=d[1][2]*d[3][3]-d[1][3]*d[3][2];
	double A1313=d[1][1]*d[3][3]-d[1][3]*d[3][1];
	double A1213=d[1][1]*d[3][2]-d[1][2]*d[3][1];
	double A2312=d[1][2]*d[2][3]-d[1][3]*d[2][2];
	double A1312=d[1][1]*d[2][3]-d[1][3]*d[2][1];
	double A1212=d[1][1]*d[2][2]-d[1][2]*d[2][1];
	double A0313=d[1][0]*d[3][3]-d[1][3]*d[3][0];
	double A0213=d[1][0]*d[3][2]-d[1][2]*d[3][0];
	double A0312=d[1][0]*d[2][3]-d[1][3]*d[2][0];
	double A0212=d[1][0]*d[2][2]-d[1][2]*d[2][0];
	double A0113=d[1][0]*d[3][1]-d[1][1]*d[3][0];
	double A0112=d[1][0]*d[2][1]-d[1][1]*d[2][0];
	/* Then we compute the determinant (or its inverse),
		we use subdeterminants and Laplace expansion to make it "faster" */
	double inv_det=1/(
		+d[0][0]*(d[1][1]*A2323-d[1][2]*A1323+d[1][3]*A1223) 
		-d[0][1]*(d[1][0]*A2323-d[1][2]*A0323+d[1][3]*A0223) 
		+d[0][2]*(d[1][0]*A1323-d[1][1]*A0323+d[1][3]*A0123) 
		-d[0][3]*(d[1][0]*A1223-d[1][1]*A0223+d[1][2]*A0123)
	);
	double nd[4][4];
	/* And then we calculate each element of the inverted matrix,
		we use here rule of block matrices */
	nd[0][0]=inv_det*+(d[1][1]*A2323-d[1][2]*A1323+d[1][3]*A1223);
	nd[0][1]=inv_det*-(d[0][1]*A2323-d[0][2]*A1323+d[0][3]*A1223);
	nd[0][2]=inv_det*+(d[0][1]*A2313-d[0][2]*A1313+d[0][3]*A1213);
	nd[0][3]=inv_det*-(d[0][1]*A2312-d[0][2]*A1312+d[0][3]*A1212);
	nd[1][0]=inv_det*-(d[1][0]*A2323-d[1][2]*A0323+d[1][3]*A0223);
	nd[1][1]=inv_det*+(d[0][0]*A2323-d[0][2]*A0323+d[0][3]*A0223);
	nd[1][2]=inv_det*-(d[0][0]*A2313-d[0][2]*A0313+d[0][3]*A0213);
	nd[1][3]=inv_det*+(d[0][0]*A2312-d[0][2]*A0312+d[0][3]*A0212);
	nd[2][0]=inv_det*+(d[1][0]*A1323-d[1][1]*A0323+d[1][3]*A0123);
	nd[2][1]=inv_det*-(d[0][0]*A1323-d[0][1]*A0323+d[0][3]*A0123);
	nd[2][2]=inv_det*+(d[0][0]*A1313-d[0][1]*A0313+d[0][3]*A0113);
	nd[2][3]=inv_det*-(d[0][0]*A1312-d[0][1]*A0312+d[0][3]*A0112);
	nd[3][0]=inv_det*-(d[1][0]*A1223-d[1][1]*A0223+d[1][2]*A0123);
	nd[3][1]=inv_det*+(d[0][0]*A1223-d[0][1]*A0223+d[0][2]*A0123);
	nd[3][2]=inv_det*-(d[0][0]*A1213-d[0][1]*A0213+d[0][2]*A0113);
	nd[3][3]=inv_det*+(d[0][0]*A1212-d[0][1]*A0212+d[0][2]*A0112);
	/* The last step is to save the matrix */
	memcpy(m_data,nd,sizeof(double)*16);
}
/*
	Method which multiplies this matrix with another one (from right)
	=> mat Multiplicator
	<= Result of "M*mat"
*/
matrix_4x4 matrix_4x4::operator*(const matrix_4x4& mat){
	/* Creating new matrix (on stack - don't worry, this is overloaded
		operator) */
	matrix_4x4 res;
	/* For each row in this matrix */
	for(size_t a=0; a<4; a++){
		/* For each column in second matrix */
		for(size_t b=0; b<4; b++){
			/* we will calculate sum of the corresponding elements (we do dot
				product of this matrix row and second matrix column) */
			double sum = 0;
			/* For each one of the corresponding elements */
			for(size_t c=0; c<4; c++){
				/* We add the result of the multiplication */
				sum+=m_data[a][c]*mat.m_data[c][b];
			}
			/* Then we can set the value of the current element */
			res.m_data[a][b]=sum;
		}
	}
	/* Finally, we return the result */
	return res;
}
