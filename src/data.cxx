/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################

	This class is used as an interface between PBD algorithm
	and VTK. It holds data from VTK (passed from "main" module),
	but it is independent on it.
*/
#include "point_3D.h"
#include "graph.h"
#include "data.h"
#include <cmath>
/*
	C'tor for data interface
	=> points Point array (x1, y1, z1, x2, y2, z2, x3, ...)
	=> num_points Number of points in "point" array (points.len/3)
	=> polys Polygons/Triangles array (3 i j k 3 l m n 3 ...)
	=> num_polys Number of triangles in "polys" array (polys.len/4)
	=> scalars Scalar data (e.g. colors) assigned to each point
		(structure is the same as "points" vector)
	=> num_scalars Number of scalar data entries
	=> points_lines Coordinates of the point forming muscle fibres
		(structure is the same as is in "points" array)
	=> num_points_lines Number of points in "points_lines" array
	=> lines Lines array (for example 2 i j 2 k l 2 m n 2 o p)
	=> num_lines Number of lines in "lines" array
*/
data::data(
	double* points,
	long long int num_points,
	long long int* polys,
	long long int num_polys,
	unsigned char* scalars,
	long long int num_scalars,
	double* points_lines,
	long long int num_points_lines,
	long long int* lines,
	long long int num_lines
){
	/* Just passing the parameters into the object */
	m_points = points;
	m_num_points = num_points;
	m_lines = lines;
	m_num_lines = num_lines;
	m_polys = polys;
	m_num_polys = num_polys;
	m_scalars = scalars;
	m_num_scalars = num_scalars;
	m_points_lines = points_lines;
	m_num_points_lines = num_points_lines;
	m_lines = lines;
	m_num_lines = num_lines;
	/* and allocation space for closest points (each point determining line
		references a signle (closest) point of the model) */
	m_closest_points = new long long int[num_points_lines];
}
/*
	D'tor
*/
data::~data(){
	/* Frees the closes point reference array */
	delete[] m_closest_points;
}
/*
	Converts this data to the graph structure
	=> points Number of point in each model
	=> classes Number of models in each class
	=> limit Upper bound for full matrix creation
*/
graph* data::create_graph(
	std::vector<size_t> points,
	std::vector<size_t> classes,
	size_t limit
){
	/* At first, we crate an empty graph */
	graph* ret = new graph(points,classes,limit);
	/* then, we put points into graph */
	ret->set_positions(m_points);
	/* and information about colors */
	ret->set_colors(m_scalars);
	/* Next, we have to iterate over all triangles in the model */
	FOREACH_DATA(this,polys){
		/* at first, we prepare vector, which holds all point references 
			in the given triangles (probably should be std::array, but
			this is more general) */
		std::vector<size_t> points;
		/* then, for each point in actual triangle */
		FOREACH_CELL {
			/* we get next index of the triangle (modulo) - triangle 1 2 3:
				cell_value: 1 => next_index: 2
				cell_value: 2 => next_index: 3
				cell_value: 3 => next_index: 1
				this allows to create edges
			*/
			long long int next_index =
				data_array[(cell_array_index<array_index+cell_size-1)
					?(cell_array_index+1)
					:(array_index+1)];
			/* Then, we prepare temporary array for position difference */
			double tmp[POINT_DIM];
			/* we acquire position of the current point */
			const double* p1 = ret->get_position(cell_value);
			/* and the position of the following point on triangle */
			const double* p2 = ret->get_position(next_index);
			/* we calculate its difference */
			SUB_VEC(tmp,p1,p2);
			/* and length of the edge (difference) */
			double dist= std::sqrt(DOT_VEC(tmp,tmp));
			/* in the next step, we need to get edge stiffness (it varies
				because of anisotropy) */
			double stiffness = find_stiffness(p1,tmp);
			/* next, we create an actual edge containing previous information */
			ret->edge(cell_value,stiffness,next_index) = dist;
			/* and push point reference into vector */
			points.push_back(static_cast<size_t>(cell_value));
	
		} END_FOREACH
		/* Then, we try to add the whole polygon into the graph,
			this is not trivial, because we have to distinguish
			different models (bones/muscles) */
		/* At first, we get number of points in the current polygon
			(should be 3 for triangles, but this is more general) */
		size_t len = points.size();
		/* We prepare variable, which tells us which model we are in */
		size_t part = static_cast<size_t>(-1);
		/* then, we iterate over all points in the polygon/triangle */
		for(size_t a=0; a<len; a++){
			/* we get current model */
			size_t curr_part = ret->get_part_from_node(points[a]);
			/* if this is not the first point (then previous model ID does
				not mind) and the model IDs are not the same) */
			if(a && curr_part!=part){
				/* then this is not a valid triangle, so we get rid of it */
				goto do_not_add;
			/* if this is a valid triangle so far */
			}else{
				/* we adjust the current model ID for the next iteration */
				part=curr_part;
			}
		}
		/* if the polygon is valid, we can add it into the graph */
		ret->add_poly(part,{points[0],points[1],points[2]});
		/* or we jump here otherwise */
		do_not_add:;
	} END_FOREACH
	/*for(int a=0; a<6931; a+=10){
		for(int b=0; b<6931; b+=10){
			if(a==b){
				continue;
			}
			double tmp[POINT_DIM];
			const double* p1 = ret->get_position(a);
			const double* p2 = ret->get_position(b);
			SUB_VEC(tmp,p1,p2);
			double dist= std::sqrt(DOT_VEC(tmp,tmp));
			double stiffness = find_stiffness(p1,tmp);
			ret->edge(a,stiffness,b) = dist;
		}
	}*/
	/* TODO */
	/*		double tmp[POINT_DIM];
			const double* p1 = ret->get_position(795);
			const double* p2 = ret->get_position(556);
			SUB_VEC(tmp,p1,p2);
			double dist= std::sqrt(DOT_VEC(tmp,tmp));
			double stiffness = find_stiffness(p1,tmp);
			ret->edge(795,stiffness,556) = dist;
	*/
	/* finally, we can return this graph (PBD d'tor does free this instance) */
	return ret;
}
/*
	Method returns point coordinates which forms fibres
	<= Fibre points coordinates
*/
double* data::get_points_lines(){
	return m_points_lines;
}
/*
	Method returns number of points which forms fibres
	<= Number of points in fibres
*/
long long int data::get_number_of_points_lines(){
	return m_num_points_lines;
}
/*
	Method returns point coordinates of all models (bones/muscles)
	<= All points coordinates
*/
double* data::get_points(){
	return m_points;
}
/*
	Method returns total number of points
	<= Total point count
*/
long long int data::get_number_of_points(){
	return m_num_points;
}
/*
	Method returns array of point reference which together forms fibres
	<= Fibre point reference array
*/
long long int* data::get_lines(){
	return m_lines;
}
/* 
	Method returns number of line segments
	<= Number of line segments
*/
long long int data::get_number_of_lines(){
	return m_num_lines;
}
/*
	Method returns model points reference closest to each of
	the fibre node point
	<= Closest fibre points to the model points
*/
long long int* data::get_closest_points(){
	return m_closest_points;
}
/*
	Method returns polygons/triangles reference array
	<= Polygons/triangles reference array (e.g. [3] 0 1 2 [3] 3 4 5)
*/
long long int* data::get_polys(){
	return m_polys;
}
/*
	Method returns number of polygons/triangles forming this graph
	<= Number of triangles/polygons
*/
long long int data::get_number_of_polys(){
	return m_num_polys;
}
/*
	Method returns scalar data (color data) for each point
	<= Color array for each point
*/
unsigned char* data::get_scalars(){
	return m_scalars;
}
/*
	Method returns number of scalar data
	<= Number of scalar/color data
*/
long long int data::get_number_of_scalars(){
	return m_num_scalars;
}
/*
	Calculates stiffness of a single edge respecting anisotropy
	=> p1 Current model point
	=> direction Directional vector of the model edge
	<= Edge stiffness
*/
double data::find_stiffness(const double* p1,const double* direction){
	/* at the beggining, the stiffness will be
		maxed out (it can be betwenn 0 and 1) */
	double stiffness=1;
	/* we have to know, which point of the model is the closest one */
	size_t min_index = MAX;
	/* so we need to get store its min. distance to each fibre point */
	double min_dist = std::numeric_limits<double>::max();
	/* this is a temporary vector, which we use to calculate distance */
	double tmp[3];
	/* and then we store pointer to the closest point coordinates */
	double* best_pnt = nullptr;
	/* So we iterate over all fibres point */
	FOREACH_POINT(this,points_lines){
		/* for each fibre point, we calculate direction vector from given
			model point */
		SUB_VEC(tmp,pt,p1);	
		/* then, we calculate its distance (L2^2) */
		double dst = DOT_VEC(tmp,tmp);
		/* it this point is the closest one */
		if(dst<min_dist){
			/* we declare so */
			min_dist = dst;
			/* we store index to this point */
			min_index = current_point/3;
			/* and the coordinates, too */
			best_pnt = pt;
		}
	} END_FOREACH
	/* If we found closest point, then we can continue */
	if(min_index<MAX){
		/* we prepare auxiary vectors */
		double unif_dir[3],diff_fib[3];
		/* and best stiffnes factor (should be in <-1;1> with max abs. value */
		double best_stiff = 0;
		/* then, we need to make direction uniform, so we can calculate proper
			dot product in <-1;1> */
		UNIF_VEC(unif_dir,direction);
		/* For each line segment */
		FOREACH_DATA(this,lines){
			/* and for each point in given segment */
			FOREACH_CELL {
				/* if selected point is in this segment */
				if(static_cast<size_t>(cell_value) == min_index){
					/* then we iterate once more over given segment */
					FOREACH_CELL {
						/* there we should skip the zero vector, but we don't
							care, because dot product will be zero, so we skip
							the condition */
						double* pnt = m_points_lines+3*cell_array_index;
						/* at this point, we calculate difference
							(can be zero vector) */
						SUB_VEC(diff_fib,pnt,best_pnt);
						/* we make this vector of length 1
							(zero vector stays zero) */
						UNIF_VEC(diff_fib,diff_fib);
						/* then, we calculate stiffnes with dot product */
						double stiff = DOT_VEC(unif_dir,diff_fib);
						/* and if this edge is better then the previous ones */
						if(std::abs(stiff)>best_stiff){
							/* we store its value */
							best_stiff = stiff;
						}
					} END_FOREACH
				}
			} END_FOREACH
		} END_FOREACH
		/* Finally, we have to negate the value */
		stiffness = 1-std::abs(best_stiff);
	}
	/* and we can return the stiffnes value */
	return stiffness;
}
