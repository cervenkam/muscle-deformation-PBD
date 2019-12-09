/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################
*/
#pragma once
#include <cstddef>
#include <limits>
#define POINT_DIM 3
/* Macros which performs vector calculations (macros for speedup, sometimes
	you shouldn't trust compiler /especially MSVC :D/), notice that
	OUTPUT argument IS the FIRST one! */
/* Prints 3 element vector to standard output (debugging purpose) */
#define PRINT_VEC(X)\
	std::cout << (X)[0] << " ; " << (X)[1] << " ; " << (X)[2] << std::endl
/* Calculates dot product of two 3 element vectors */
#define DOT_VEC(v1,v2)\
	((v1)[0]*(v2)[0]+(v1)[1]*(v2)[1]+(v1)[2]*(v2)[2])
/* Casts each element of the vector to given type */
#define CAST_VEC(dest,src,type)\
	(dest)[0]=static_cast<type>(std::floor((src)[0]));\
	(dest)[1]=static_cast<type>(std::floor((src)[1]));\
	(dest)[2]=static_cast<type>(std::floor((src)[2]))
/* Calculates cross product of two 3 element vectors */
#define CROSS_VEC(dest,v1,v2)\
	(dest)[0]=(v1)[1]*(v2)[2]-(v1)[2]*(v2)[1];\
	(dest)[1]=(v1)[2]*(v2)[0]-(v1)[0]*(v2)[2];\
	(dest)[2]=(v1)[0]*(v2)[1]-(v1)[1]*(v2)[0]
/* Adds two vectors together */
#define ADD_VEC(dest,v1,v2) \
	(dest)[0]=(v1)[0]+(v2)[0];\
	(dest)[1]=(v1)[1]+(v2)[1];\
	(dest)[2]=(v1)[2]+(v2)[2]
/* Adds a scalar value to a vector */
#define ADD_VEC_SCALAR(dest,v1,val) \
	(dest)[0]=(v1)[0]+(val);\
	(dest)[1]=(v1)[1]+(val);\
	(dest)[2]=(v1)[2]+(val)
/* Subtract one vector from another */
#define SUB_VEC(dest,v1,v2) \
	(dest)[0]=(v1)[0]-(v2)[0];\
	(dest)[1]=(v1)[1]-(v2)[1];\
	(dest)[2]=(v1)[2]-(v2)[2]
/* Calculate normal to pair of input vectors */
#define NORM_VEC(dest,v0,v1)\
	CROSS_VEC((dest),(v0),(v1));\
	UNIF_VEC((dest),(dest))
/* Sets given vector it given values */
#define SET_VEC(dest,x,y,z)\
	(dest)[0] = (x);\
	(dest)[1] = (y);\
	(dest)[2] = (z)
/* Converts vector, so it will retain direction, but have length of 1 */
#define UNIF_VEC(dest,src) \
	{\
		double len = std::sqrt(DOT_VEC((src),(src)));\
		if(std::abs(len)>1e-10){\
			DIV_VEC((dest),(src),len);\
		}else{\
			SET_VEC(dest,0.,0.,0.);\
		}\
	}
/* Calculates normal vector of a triangle */
#define NORM_TRI(dest,tmp1,tmp2,v0,v1,v2)\
	SUB_VEC((tmp1),(v1),(v0));\
	SUB_VEC((tmp2),(v2),(v0));\
	NORM_VEC((dest),(tmp1),(tmp2))
/* Copies vector to another one */
#define CPY_VEC(dest,src)\
	(dest)[0] = (src)[0];\
	(dest)[1] = (src)[1];\
	(dest)[2] = (src)[2]
/* Multiplies vector by a scalar */
#define MULT_VEC(dest,src,var)\
	(dest)[0] = (src)[0]*(var);\
	(dest)[1] = (src)[1]*(var);\
	(dest)[2] = (src)[2]*(var)
/* Divides vector by a scalar value */
#define DIV_VEC(dest,src,var)\
	{ \
		double val = 1./(var);\
		MULT_VEC(dest,src,val);\
	}
/* Gets smaller of the two values */
#define MIN2(A,B) (((A)>(B))?(B):(A))
/* Gets bigger of the two values */
#define MAX2(A,B) (((A)<(B))?(B):(A))
/* Get the biggest value of the three */
#define MIN3(A,B,C) (((A)<(B))?(((A)<(C))?(A):(C)):(((B)<(C)?(B):(C))))
/* Get the smallest value of the three */
#define MAX3(A,B,C) (((A)>(B))?(((A)>(C))?(A):(C)):(((B)>(C)?(B):(C))))
/*
	This class holds additional information about a single point of PBD
	algorithm (force, velocity working position, etc.)
	Position is stored in the graph directly, which significantly speeds up the
	whole calculation
*/
class point_3D{
	public:
		/*
			Default c'tor
		*/
		point_3D();
		/*
			C'tor with point ID
			=> Point id (for debugging)
		*/
		point_3D(size_t id);
		/*
			Returns current velocity of this point
			<= Velocity (corrdinates/seconds)
		*/
		double* get_velocity();
		/*
			Returns current working positon
			In original PBD, this point is called "p", this is a candidate
			point, where should be this point moved in the next iteration,
			but has to be adjusted by collision detection etc.
			<= Working position
		*/
		double* get_working_position();
		/*
			Returns amount of force acting on the point
			<= Amount of force
		*/
		double* get_force();
		/*
			Sets working position (for "working position", see
			"get_working_position")
			=> ptr Working position
		*/
		void set_working_position(const double*);
		/*
			Sets force applied on this point
			=> ptr Force to be applied on this point
		*/
		void set_force(const double*);
		/*
			Resets accumulator (for accumulator, see "accumulate")
		*/
		void reset_accumulator();
		/*
			Returns accumulator (for accumulator, see "acumulate")
			<= Accumulator
		*/
		double* accumulator();
		/*
			Accumulate a vector to this point, each point has its accumulator,
			to which we can add more and more vectors and can be reset by
			"reset_accumulator". It is used for volume calculation, when each
			point has to remember a vector of surrounding points information
			=> ptr Vector added to the accumulator (3 elements)
		*/
		void accumulate(double* pnt);
		/*
			Returns inverse of mass of this point
			<= Inverse of mass (1/m [1/kg])
		*/
		double& inverse_mass();
		/*
			Returns this point ID
			<= Point ID
		*/
		size_t& id();
		/*
			Returns if this point is fixed (cannot be moved by PBD algorithm,
			for example it is a muscle point attached to a bone)
			<= Can this point move or not?
		*/
		bool is_fixed();
		/*
			Set this point fixed, which means PBD cannot move with this point
			=> fixed Can be this point moved or not?
		*/
		void set_fixed(size_t);
		/*
			Returns ID of the model to which this point is connected
			<= Attachment model ID
		*/
		size_t get_fixed();
	private:
		/*
			Inits this point with given ID
			=> id New point ID
		*/
		void init(size_t id);
		/* ID of the current point */
		size_t m_id = 666;
		/* velocity of this node (coordinates/seconds) */
		double m_velo[POINT_DIM] = {0};
		/* current force applied */
		double m_force[POINT_DIM] = {0};
		/* working position ("p" in PBD) */
		double m_work[POINT_DIM] = {0};
		/* accumulator (see "accumulate" method) */
		double m_accum[POINT_DIM] = {0};
		/* inverse mass (1/m [1/kg]) */
		double m_inv_mass = 1;
		/* ID of the model to which this point is attached */
		size_t m_fixed = std::numeric_limits<size_t>::max();
};
