/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################

	This class provides interface for various types
	of collision detection algorithms
*/
#include "collision.h"
/*
	Method finds out if collision occurs or not
	=> x Source position
	=> p Target position
	<=> res "safe" position (only if collision occurs)
	<= Does line segment collide with current model
*/
bool collision_detection::find(
	const double* x,const double* p,double* res
){
	/* at first, we prepare auxiary vectors */
	double x_trans[3],p_trans[3]; 
	/* then, we copy original position vector (to preserve the original one) */
	CPY_VEC(x_trans,x);
	/* same with working position vector */
	CPY_VEC(p_trans,p);
	/* then, we transform position vector into collision detection space
	(which does not change) */
	m_total_inverse_transform.transform_vector(x_trans);
	/* same with working position */
	m_total_inverse_transform.transform_vector(p_trans);
	/* then we can find the collision - implementation dependent,
		this is just interface */
	bool result = find_internal(x_trans,p_trans,res);
	/* finally, we have to transform result back from collision
		detection space */
	m_total_transform.transform_vector(res);
	/* and we return the result, if it collides */
	return result;
}
/*
	Sets the transform matrix for collision detection (we won't have
	to recalculate it in each iteration this way)
	=> mat Transformation matrix
*/
void collision_detection::set_transform(matrix_4x4 mat){
	/* at first, we store total transformation
		(from "real" space to "collison detection" space) */
	m_total_transform = mat*m_total_transform;
	/* then, we store last used transformation */
	m_transform = mat;
	/* we store the same matrix into total inverse transform */
	m_total_inverse_transform = m_total_transform;
	/* but we have to (of course) inverse this matrix (hopefully not singular,
		but we assume it is rigid transform so determinant should be 1 */
	m_total_inverse_transform.inverse();
	/* we do the same with last used transform, we store it */
	m_inverse_transform = m_transform;
	/* and inverse */
	m_inverse_transform.inverse();
}
/*
	D'tor
*/
collision_detection::~collision_detection(){}
