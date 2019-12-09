/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################

	This class implements collision detection with
	space division algoithm known as "voxelization". Space
	is divided uniformly into multiple blocks 
*/
#include "voxel.h"
#include <fstream>
#include <cstring>
#include <cmath>
/*
	Default c'tor
*/
voxel::voxel(){}
/*
	Returns orignal coordinates from voxel indicies
	This is inverse of "get_subblock" method
	=> pnt Voxel indicies
	<=> vc Original coordinates
*/
void voxel::get_original(const double* pnt,double* vc){
	/* for each coordinate, we just multiply */
	vc[0]=pnt[0]*m_bounds[6]+m_bounds[0];
	/* the coordinate with given scale and */
	vc[1]=pnt[1]*m_bounds[7]+m_bounds[1];
	/* add minimum coordinate value */
	vc[2]=pnt[2]*m_bounds[8]+m_bounds[2];
}
/*
	Returns voxel indicies from original coordinates
	This is the inverse of "get_original" method
	=> pnt Original coordinates
	<=> vc Voxel indicies
*/
void voxel::get_subblock(const double* pnt,double* vc){
	/* here we just subtract smallest value */
	vc[0]=(pnt[0]-m_bounds[0])/m_bounds[6];
	/* and divide the result by voxel size */
	vc[1]=(pnt[1]-m_bounds[1])/m_bounds[7];
	/* for each coordinate (X,Y,Z) */
	vc[2]=(pnt[2]-m_bounds[2])/m_bounds[8];
}
/*
	Method generates collision detection
	=> part Model ID
	=> graph Data about given model
*/
void voxel::generate(size_t part,graph* graph){
	/* first of all, we reset all voxels */
	memset(m_data,0,sizeof(m_data));
	/* then, we calculate AABB of the model */
	graph->boundary(part,
		m_bounds[0],m_bounds[3],
		m_bounds[1],m_bounds[4],
		m_bounds[2],m_bounds[5]
	);
	/* next we need to calculate one voxel size, so we subtract minimum
		coordinates from maximum ones */
	SUB_VEC(m_bounds+6,m_bounds+3,m_bounds);
	/* and divide it by number of voxel in each axis */
	DIV_VEC(m_bounds+6,m_bounds+6,PARAMETER);
	/* diagonal size is sqrt(3) (diagonal in a cube of edge length 1) */
	m_diagonal = std::sqrt(3)+1e-3; 
	/* then, we get number of polygons in the model */
	size_t len = graph->polys(part);
	/* and prepare AABB for the model */
	double aabb[8];
	/* starting voxel position */
	double start[3];
	/* and ending voxel position */
	double end[3];
	/* for each polygon/triangle in the model */
	for(size_t a=0; a<len; a++){
		/* we get all three indices of the triangle */
		std::array<size_t,3>& pnts = graph->poly(part,a);
		/* so we can get coordinates of the first */
		const double* p1 = graph->get_position(pnts[0]);
		/* second */
		const double* p2 = graph->get_position(pnts[1]);
		/* and third point forming current triangle */
		const double* p3 = graph->get_position(pnts[2]);
		/* then we calculate AABB of the triangle */
		computeAABB(p1,p2,p3,aabb);
		/* we transform the coordinates into voxel indicies */
		get_subblock(aabb,start);
		/* the whole AABB will be transformed here */
		get_subblock(aabb+3,end);
		/* we add one to the end (to make sure we iterate over every single
			voxel */
		ADD_VEC_SCALAR(end,end,1);
		/* here we acquire normal vector of the triangle (in "real" space) */
		const double* n = graph->get_normal(part,a);
		/* next, we need temporary vector, which will be used for "walking"
			through all voxels */
		size_t p[3];
		/* so we "walk" through all X coordinates, which surrounds triangle */
		for(p[0]=start[0]; p[0]<end[0]; p[0]++){
			/* same with Y */
			for(p[1]=start[1]; p[1]<end[1]; p[1]++){
				/* and Z coordinates */
				for(p[2]=start[2]; p[2]<end[2]; p[2]++){
					/* if given triangle collides with current voxel */
					if(box_collide(p1,p2,p3,n,p)){
						/* and the voxel coordinates are inside the voxelized
							area */
						if(p[0]<PARAMETER && p[1]<PARAMETER && p[2]<PARAMETER){
							/* we can tell there is a border */
							m_data[p[0]][p[1]][p[2]].border=1;
						}
					}
				}
			}
		}
	}
	/* There the whole space will be declared as "inside",
		will be fixed later, so we "walk" through all X coordinates */
	for(size_t x=0; x<PARAMETER; x++){
		/* through Y coordinates */
		for(size_t y=0; y<PARAMETER; y++){
			/* and even Z ones */
			for(size_t z=0; z<PARAMETER; z++){
				/* the voxel is inside by default */
				m_data[x][y][z].inside=1;
			}
		}
	}
	/* CRAZY SEEDFILL follows (running out of time, don't judge ;)) */
	/* we remember all voxels, which was found */
	std::vector<std::array<int,3> > found;
	/* into this vector we put the first one (which can't be inside...) */
	found.push_back({0,0,0});
	/* and until the vector is empty */
	while(!found.empty()){
		/* we get found voxel */
		std::array<int,3> p = found.back();
		/* remove it from the array */
		found.pop_back();
		/* and get its value */
		struct voxel_element& el = m_data[p[0]][p[1]][p[2]];	
		/* if this voxel is not interesting (we already have been there) */
		if(!el.inside || el.border){
			/* we continue with next one */
			continue;
		}
		/* If there is no border, we can tell this part is outside */
		m_data[p[0]][p[1]][p[2]].inside=0;
		/* and now we start searching around this voxel to each side by one */
		for(int a=-1; a<2; a+=2){
			/* we get X coordinate of an adjacent voxel */
			int c = p[0]+a;
			/* and if it is in voxelized space */
			if(c<PARAMETER && c>0){
				/* we can try to set it as outside */
				found.push_back({c,p[1],p[2]});
			}
			/* we get Y coordinate of an adjacent voxel */
			c = p[1]+a;
			/* and if it is in voxelized space */
			if(c<PARAMETER && c>0){
				/* we can try to set it as outside */
				found.push_back({p[0],c,p[2]});
			}
			/* we get Z coordinate of an adjacent voxel */
			c = p[2]+a;
			/* and if it is in voxelized space */
			if(c<PARAMETER && c>0){
				/* we can try to set it as outside */
				found.push_back({p[0],p[1],c});
			}
		}
	}
}
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
bool voxel::box_collide(
	const double* p1,const double* p2,const double* p3,
	const double* n,const size_t* x
){
	/* at first, we create temporary vectors */
	double p1s[3],p2s[3],p3s[3];
	/* we transform coordinates into "voxel" coordinates, first point */
	get_subblock(p1,p1s);
	/* second point */
	get_subblock(p2,p2s);
	/* and third point of the triangle too */
	get_subblock(p3,p3s);
	/* then we need more temporary vectors */
	double tmp1[3],tmp2[3],w[3];
	/* with which we calculate normal vector of the triangle (in voxel space) */
	NORM_TRI(w,tmp1,tmp2,p1s,p2s,p3s);
	/* then we need to get distance of the given voxel to a triangle plane */
	double diff = DOT_VEC(w,p1s)-DOT_VEC(w,x);
	/* if the distance is too large (more than half of a diagonal) */
	if(std::abs(diff)>m_diagonal/2){
		/* then triangle is certainly not intersecting given voxel */
		return false;
	}
	/* exact test follows, we need one more temporary vector */
	double sub[3];
	/* at first, we get directional vector between first and second point */
	SUB_VEC(sub,p2s,p1s);
	/* then we get normal of a plane perpedicular to triangle plane and
		containing one of the three edges */
	NORM_VEC(w,sub,n);
	/* then we calculate distance to this plane */
	diff = DOT_VEC(w,x)-DOT_VEC(w,p1s);
	/* if the distance is too high */
	if(diff>m_diagonal){
		/* it cannot collide */
		return false;
	}
	/* similar procedure follows for edge formed by third and second point
		at first, we calculate directional vector */
	SUB_VEC(sub,p3s,p2s);
	/* then we need perpendicular plane to triangle plane, which also
		contains edge formed by third and second point (like before) */	
	NORM_VEC(w,sub,n);
	/* so we calculate distance of given voxel to this plane */
	diff = DOT_VEC(w,x)-DOT_VEC(w,p2s); 
	/* and if it is too high */
	if(diff>m_diagonal){ 
		/* it cannot collide */
		return false;
	}
	/* finally, we test the last edge, so we subtratct third point from the
		first one */
	SUB_VEC(sub,p1s,p3s);
	/* then we get the last perpendicular plane containing one of the edges
		in it */
	NORM_VEC(w,sub,n);
	/* also we need distance to this plane */
	diff = DOT_VEC(w,x)-DOT_VEC(w,p3s);
	/* it given voxel is too far */
	if(diff>m_diagonal){
		/* it cannot collide */
		return false;
	}
	/* we tested (and excluded) everything, thus it collides */
	return true;
}
/*
	Computes AABB (bounding box) of a triangle
	=> p1 First point of a triangle
	=> p2 Second point of a triangle
	=> p3 Third point of a triangle
	<=> aabb Bounding box ([x_min y_min z_min x_max y_max z_max])
*/
void voxel::computeAABB(
	const double* p1,const double* p2,const double* p3,double* aabb
){
	for(size_t a=0; a<3; a++){
		aabb[a] = MIN3(p1[a],p2[a],p3[a]);
		aabb[a+3] = MAX3(p1[a],p2[a],p3[a]);
	}
}
/*
	Method returns if point collides with this model
	=> vec Point to check
	<= Does point collide?
*/
bool voxel::is_inside(const double* vec){
	/* we prepare temporary vector of different data type */
	size_t idx[3];
	/* then we cast the original one (we perform rounding down) */
	CAST_VEC(idx,vec,size_t);
	/* and finally we return if we are both in bounds and inside object */
	return idx[0]<PARAMETER && idx[1]<PARAMETER && idx[2]<PARAMETER &&
		m_data[idx[0]][idx[1]][idx[2]].inside;
}
/*
	Method finds out if collision occurs or not
	=> x Source position in collision detection space
	=> p Target position in collision detection space
	<=> ret "safe" position in collision detection space
		(only if collision occurs)
	<= Does line segment collide with current model (in collision
		detection space)
*/
bool voxel::find_internal(
	const double* x,const double* p,double* ret
){
	/* first, we prepare temporary vectors */
	double X[3],P[3],diff[3];
	/* then, we transform origin coordinates into voxel indicies */
	get_subblock(x,X);
	/* same with target point */
	get_subblock(p,P);
	/* then, we calculate difference (in voxels) */
	SUB_VEC(diff,X,P);
	/* its length */
	double len = std::sqrt(DOT_VEC(diff,diff));
	/* and lets make it uniform */
	DIV_VEC(diff,diff,len);
	/* presumption of innocence, it does not collide... */
	bool retval = false;
	/* if it collides at its origin, it means bone has moved */
	if(is_inside(X)){
		/* so we move original point to the return value */
		CPY_VEC(ret,x);
		/* do the same transformation as bone has */
		m_transform.transform_vector(ret);
		/* and we tell it collides */
		retval = true;
	}else{
		/* if it is OK in original position, we go in reverse order,
			from target point to the origin..., at first, we need to know
			how many steps we made */
		size_t a;
		/* then, we walk from target point to the origin */
		for(a=0; a<(len+1); a++){
			/* if we collide */
			if(is_inside(P)){
				/* its bad, we have to set the flag */
				retval=true;
			/* it it is OK now */
			}else{
				/* we can break the "walking" loop */
				break;
			}
			/* there we add one step to the target point towards origin */
			ADD_VEC(P,P,diff);	
		}
		/* if we end up going back to the origin */
		if(a>len+1){
			/* we do not need to calculate anything, just return origin */
			CPY_VEC(ret,x);
		/* in case we stopped somewhere on the way */
		}else{
			/* we calculate the position from voxel indices */
			get_original(P,ret);
		}
	}
	/* and finally we return if line segment collides */
	return retval;
}
