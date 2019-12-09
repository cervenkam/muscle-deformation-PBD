/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################
*/
#pragma once
#include "collision.h"
#ifndef PARAMETER
#define PARAMETER 64
#endif
/* Content of each voxel */
struct voxel_element{
	/* is this voxel inside the model */
	uint8_t inside:1;
	/* is this voxel on the border of the model */
	uint8_t border:1;
};
/*
	This class implements collision detection with
	space division algoithm known as "voxelization". Space
	is divided uniformly into multiple blocks 
*/
class voxel: public collision_detection{
	public:
		/*
			Default c'tor
		*/
		voxel();
		/*
			Method generates collision detection
			=> part Model ID
			=> graph Data about given model
		*/
		void generate(size_t,graph* graph);
		/*
			Method finds out if collision occurs or not
			=> x Source position in collision detection space
			=> p Target position in collision detection space
			<=> ret "safe" position in collision detection space
				(only if collision occurs)
			<= Does line segment collide with current model (in collision
				detection space)
		*/
		bool find_internal(const double*,const double*,double*);
	private:
		/*
			Computes AABB (bounding box) of a triangle
			=> p1 First point of a triangle
			=> p2 Second point of a triangle
			=> p3 Third point of a triangle
			<=> aabb Bounding box ([x_min y_min z_min x_max y_max z_max])
		*/
		void computeAABB(const double*,const double*,const double*,double*);
		/*
			Calculates, if given triangle intersects selected voxel
			=> p1 First point of the triangle
			=> p2 Second point of the triangle
			=> p3 Third point of the triangle
			=> n Triangle normal (can be computed internally, but 
				it is faster to just pass it)
			=> x Voxel coordinates
			<= Does triangle collide with given voxel
		*/
		bool box_collide(
			const double* p1,const double* p2,const double* p3,
			const double* n,const size_t* x
		);
		/*
			Method returns if point collides with this model
			=> vec Point to check
			<= Does point collide?
		*/
		bool is_inside(const double* vec);
		/*
			Returns voxel indicies from original coordinates
			This is the inverse of "get_original" method
			=> pnt Original coordinates
			<=> vc Voxel indicies
		*/
		void get_subblock(const double* pnt,double* vc);
		/*
			Returns orignal coordinates from voxel indicies
			This is inverse of "get_subblock" method
			=> pnt Voxel indicies
			<=> vc Original coordinates
		*/
		void get_original(const double* pnt,double* vc);
		/* all voxels data */
		struct voxel_element m_data[PARAMETER][PARAMETER][PARAMETER];		
		/* AABB bounding box of this model (with block size as well) */
		double m_bounds[9];
		/* size of the largest diagonal of the block */
		double m_diagonal;
};
