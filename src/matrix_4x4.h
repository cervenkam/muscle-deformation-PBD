/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################

	This class provides basic operations with 4x4 matrices
*/
#pragma once
#include "point_3D.h"
class matrix_4x4{
	public:
		/*
			C'tor of the 4x4 matrix
			=> x00 First element of the matrix
			=> x02 Second element (first row, second column)
			=> x33 Last element
			=> x03-x32 Rest of the elements
		*/
		matrix_4x4(
			double x00,double x01,double x02,double x03,
			double x10,double x11,double x12,double x13,
			double x20,double x21,double x22,double x23,
			double x30,double x31,double x32,double x33
		);
		/*
			C'tor of identity matrix
		*/
		matrix_4x4();
		/*
			Makes translate matrix
			=> x X value of the translation
			=> y Y value of the translation
			=> z Z value of the translation
		*/
		void translate(double x,double y,double z);
		/*
			Transforms input vector (3 element) by this 4x4 matrix
			Result just equals to "(M*[v;0])(1:3)" in octave/matlab
			<=> input/output 3 element vector
		*/
		void transform_vector(double* v);
		/*
			Converts this matrix to its inverse 
			(WARNING! does not care if singular)
		*/
		void inverse();
		/*
			Method which multiplies this matrix with another one (from right)
			=> mat Multiplicator
			<= Result of "M*mat"
		*/
		matrix_4x4 operator*(const matrix_4x4& mat);
	private:
		double m_data[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
};
