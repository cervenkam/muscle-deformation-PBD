/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################

	This class holds additional information about a single point of PBD
	algorithm (force, velocity working position, etc.)
	Position is stored in the graph directly, which significantly speeds up the
	whole calculation
*/
#include "point_3D.h"
#include <cstring>
#include <iostream>

/*
	Default c'tor
*/
point_3D::point_3D(){
	/* ID has to be noticeable, if we won't select one */
	init(666);
}
/*
	C'tor with point ID
	=> Point id (for debugging)
*/
point_3D::point_3D(size_t id){
	init(id);
}
/*
	Inits this point with given ID
	=> id New point ID
*/
void point_3D::init(size_t id){
	m_id = id;
}
/*
	Returns this point ID
	<= Point ID
*/
size_t& point_3D::id(){
	return m_id;
}
/*
	Returns if this point is fixed (cannot be moved by PBD algorithm,
	for example it is a muscle point attached to a bone)
	<= Can this point move or not?
*/
bool point_3D::is_fixed(){
	return this->m_fixed < std::numeric_limits<size_t>::max();
}
/*
	Set this point fixed, which means PBD cannot move with this point
	=> fixed Can be this point moved or not?
*/
void point_3D::set_fixed(size_t fixed){
	this->m_fixed = fixed;
}
/*
	Returns ID of the model to which this point is connected
	<= Attachment model ID
*/
size_t point_3D::get_fixed(){
	return this->m_fixed;
}
/*
	Returns current velocity of this point
	<= Velocity (corrdinates/seconds)
*/
double* point_3D::get_velocity(){
	return m_velo;
}
/*
	Returns current working positon
	In original PBD, this point is called "p", this is a candidate
	point, where should be this point moved in the next iteration,
	but has to be adjusted by collision detection etc.
	<= Working position
*/
double* point_3D::get_working_position(){
	return m_work;
}
/*
	Sets working position, for "working position", see "get_working_position"
	=> ptr Working position
*/
void point_3D::set_working_position(const double* ptr){
	memcpy(m_work,ptr,sizeof(double)*POINT_DIM);
}
/*
	Sets force applied on this point
	=> ptr Force to be applied on this point
*/
void point_3D::set_force(const double* ptr){
	memcpy(m_force,ptr,sizeof(double)*POINT_DIM);
}
/*
	Returns amount of force acting on the point
	<= Amount of force
*/
double* point_3D::get_force(){
	return m_force;
}
/*
	Returns inverse of mass of this point
	<= Inverse of mass (1/m [1/kg])
*/
double& point_3D::inverse_mass(){
	return m_inv_mass;
}
/*
	Resets accumulator (for accumulator, see "accumulate")
*/
void point_3D::reset_accumulator(){
	memset(m_accum,0,sizeof(double)*POINT_DIM);
}
/*
	Returns accumulator (for accumulator, see "acumulate")
	<= Accumulator
*/
double* point_3D::accumulator(){
	return m_accum;
}
/*
	Accumulate a vector to this point, each point has its accumulator,
	to which we can add more and more vectors and can be reset by
	"reset_accumulator". It is used for volume calculation, when
	each point has to remember a vector of surrounding points information
	=> ptr Vector added to the accumulator (3 elements)
*/
void point_3D::accumulate(double* ptr){
	ADD_VEC(m_accum,m_accum,ptr);
}
