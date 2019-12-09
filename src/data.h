/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################
*/
#pragma once
#include "graph.h"
#include <vector>
#define __FOREACH_DATA(ARR,LEN) \
long long int cell_size = 0; \
long long int* data_array = (ARR); \
long long int vbgrbw = (LEN); \
for(long long int array_index=0,cell_index=0; cell_index<vbgrbw; \
		array_index+=cell_size,cell_index++){ \
	cell_size = data_array[array_index]+1;
/* This macro iterates over all values in single cell, 2 variables can be used:
	cell_array_index: Index to the cell array, where the element is stored,
	cell_value: Value of the current element
	!!! Has to be closed by END_FOREACH !!! */
#define FOREACH_CELL \
for(long long int cell_array_index=array_index+1; \
		cell_array_index<array_index+cell_size; cell_array_index++) { \
	[[maybe_unused]] long long int cell_value = data_array[cell_array_index];
#define END_FOREACH }
/* This macro iterates over all point coordinates (of point/lines/polygons etc.)
	2 variables can be used:
	pt: coordinates of each point (3 element vector)
	current_point: index of the current point
	!!! Has to be closed by END_FOREACH !!! */
#define FOREACH_DATA(DATA,Y) \
	__FOREACH_DATA((DATA)->get_##Y(),(DATA)->get_number_of_##Y())

#define __FOREACH_POINT(ARR,LEN) \
long long int vqasept = 3*(LEN); \
double* point_array = (ARR);\
for(long long int current_point = 0; current_point<vqasept; current_point+=3){\
	double* pt = point_array+current_point;
#define FOREACH_POINT(DATA,Y) \
	__FOREACH_POINT((DATA)->get_##Y(),(DATA)->get_number_of_##Y())
constexpr size_t MAX = std::numeric_limits<size_t>::max();
class graph;
/*
	This class is used as an interface between PBD algorithm
	and VTK. It holds data from VTK (passed from "main" module),
	but it is independent on it.
*/
class data{
	public:
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
		data(
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
		);
		/*
			D'tor
		*/
		~data();
		/*
			Converts this data to the graph structure
			=> points Number of point in each model
			=> classes Number of models in each class
			=> limit Upper bound for full matrix creation
		*/
		graph* create_graph(std::vector<size_t>,std::vector<size_t>,size_t);
		/*
			Method returns point coordinates of all models (bones/muscles)
			<= All points coordinates
		*/
		double* get_points();
		/*
			Method returns point coordinates which forms fibres
			<= Fibre points coordinates
		*/
		double* get_points_lines();
		/*
			Method returns total number of points
			<= Total point count
		*/
		long long int get_number_of_points();
		/*
			Method returns number of points which forms fibres
			<= Number of points in fibres
		*/
		long long int get_number_of_points_lines();
		/*
			Method returns array of point reference which together forms fibres
			<= Fibre point reference array
		*/
		long long int* get_lines();
		/* 
			Method returns number of line segments
			<= Number of line segments
		*/
		long long int get_number_of_lines();
		/*
			Method returns polygons/triangles reference array
			<= Polygons/triangles reference array (e.g. [3] 0 1 2 [3] 3 4 5)
		*/
		long long int* get_polys();
		/*
			Method returns number of polygons/triangles forming this graph
			<= Number of triangles/polygons
		*/
		long long int get_number_of_polys();
		/*
			Method returns scalar data (color data) for each point
			<= Color array for each point
		*/
		unsigned char* get_scalars();
		/*
			Method returns number of scalar data
			<= Number of scalar/color data
		*/
		long long int get_number_of_scalars();
		/*
			Method returns model points reference closest to each of
			the fibre node point
			<= Closest fibre points to the model points
		*/
		long long int* get_closest_points();
	private:
		/*
			Calculates stiffness of a single edge respecting anisotropy
			=> p1 Current model point
			=> direction Directional vector of the model edge
			<= Edge stiffness
		*/
		double find_stiffness(const double*,const double*);
		/* stores model point coordinates */
		double* m_points = nullptr;
		/* stores number of points in the stored model */
		long long int m_num_points = 0;
		/* stores fibre segment array (for example: [2] 0 1 [2] 2 3) */
		long long int* m_lines = nullptr;
		/* stores number of fibre segments */
		long long int m_num_lines = 0;
		/* stores polygons/triangles (for example [3] 0 1 2 [3] 3 4 5) */
		long long int* m_polys = nullptr;
		/* stores number of polygons/triangles */
		long long int m_num_polys = 0;
		/* stores scalar values for each model point (colors) */
		unsigned char* m_scalars = nullptr;
		/* stores number of scalar values */
		long long int m_num_scalars = 0;
		/* stores fibre coordinates */
		double* m_points_lines = nullptr;
		/* stores number of fibre points */
		long long int m_num_points_lines = 0;
		/* stores reference to the closest point from each fibre point
			to one of the model point (which one is stored here) */
		long long int* m_closest_points = nullptr;
};
