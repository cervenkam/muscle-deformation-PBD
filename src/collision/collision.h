/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################
*/
#pragma once
#include <vector>
#include "../matrix_4x4.h"
#include "../graph.h"
/*
	This class provides interface for various types
	of collision detection algorithms
*/
class collision_detection{
	public:
		/*
			Method finds out if collision occurs or not
			=> x Source position
			=> p Target position
			<=> res "safe" position (only if collision occurs)
			<= Does line segment collide with current model
		*/
		bool find(const double* x,const double* p,double* res);
		/*
			Method finds out if collision occurs or not
			=> x Source position in collision detection space
			=> p Target position in collision detection space
			<=> res "safe" position in collision detection space
				(only if collision occurs)
			<= Does line segment collide with current model
			(in collision detection space)
		*/
		virtual bool find_internal(const double*,const double*,double*) = 0;
		/*
			Method generates collision detection according to the graph
			=> part Model ID
			=> graph Graph which is used to generate collision detection
		*/
		virtual void generate(size_t,graph*) = 0;
		/*
			D'tor
		*/
		virtual ~collision_detection();
		/*
			Sets the transform matrix for collision detection (we won't have
			to recalculate it in each iteration this way)
			=> mat Transformation matrix
		*/
		void set_transform(matrix_4x4 transform);
	protected:
		/* last used affine transform */
		matrix_4x4 m_transform;
		/* inverse of last used affine transform */
		matrix_4x4 m_inverse_transform;
		/* total transform from "real" space into collision detection space */
		matrix_4x4 m_total_transform;
		/* total inverse transform from "real" space into collision
			detection space */
		matrix_4x4 m_total_inverse_transform;
};
