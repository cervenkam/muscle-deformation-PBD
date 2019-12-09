/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################
*/
#pragma once
#include "point_3D.h"
#include "data.h"
#include "collision/collision.h"
#include <array>
#include <iostream>
class data;
class collision_detection;
/*
	This class implements operations on graph structure, it
	works with nodes, edges, and even faces.	
*/
class graph{
	public:
		/* struct which holds information about single edge */
		struct edge{
			/* edge length */
			double value;
			/* edge stiffnes (used in simulation) */
			double stiffness;
			/* edge stiffness (without power) */
			double abs_stiffness;
			/* target point */
			size_t end;
		};
		/*
			Graph c'tor
			=> nodes Number of nodes in each model
			=> classes Number of models in each class (bone/muscle)
			=> limit Maximum number of nodes for which dense matrix
			will be formed
		*/
		graph(std::vector<size_t> nodes,std::vector<size_t> classes,size_t);
		/*
			Graph d'tor
		*/
		~graph();
		/*
			Method creates edge and returns its value
			=> start Starting node ID
			=> stiffnes Stiffness of the new edge
			=> end Ending node ID
			<=> Edge length reference
		*/
		double& edge(size_t start,double stiffness,size_t end);
		/*
			Generates collision detection on a single model
			=> part Model ID
		*/
		void generate_collision_detection(size_t);
		/*
			Returns node by its ID
			=> node Node ID
			<= Node information (velocity, force, etc.)
		*/
		point_3D& node(size_t node);
		/*
			Method updates stiffness of each node to specified power
			=> power Power of the stiffness
		*/
		void update_stiffness(double power);
		/*
			Returns total number of nodes in the whole graph
			<= Total number of nodes
		*/
		size_t nodes();
		/*
			Affine transforms model according to specified transformation matrix
			=> part Transformed model ID
			=> trans Affine transform matrix
		*/
		void transform(size_t,matrix_4x4);
		/*
			Returns volume of the selected model
			=> part Model ID
			<= Volume of the model
		*/
		double& volume(size_t);
		/*
			Precalculates model volume and saves it for later
			=> part Model ID
		*/
		void calculate_volume(size_t);
		/*
			Returns starting position and ending position of
			all nodes in given class
			=> cls Class ID
			<= start Node start ID
			<= end Node end ID
		*/
		void nodes_from_class(size_t cls,size_t& start,size_t& end);
		/*
			Returns stating and ending index of nodes in specified model
			=> part Model ID
			<= start Node starting ID
			<= end Node ending ID
		*/
		void nodes_from_part(size_t cls,size_t& start,size_t& end);
		/*
			Returns stating and ending index of models in specified class
			(for example muscles, bones, ...)
			=> cls Class ID
			<= start Model starting ID
			<= end Model ending ID
		*/
		void parts_from_class(size_t cls,size_t& start,size_t& end);
		/*
			Calculates AABB of single given model in this graph
			=> part ID of the model
			=> min_x Minimum X coordinate
			=> max_x Maximum X coordinate
			=> min_y Minimum Y coordinate
			=> max_y Maximum Y coordinate
			=> min_z Minimum Z coordinate
			=> max_z Maximum Z coordinate
		*/
		void boundary(size_t,double&,double&,double&,double&,double&,double&);
		/*
			Returns edge from sparse structure (vector), this is linear,
			we should address dence matrix "m_arr", if we have enough memory
			=> start Start node ID
			=> end End node ID
			<=> Edge information
		*/
		struct edge& slow_at(size_t start,size_t end);
		/*
			Returns edge lenght between two specified nodes
			=> start Starting node ID
			=> end Ending node ID
			<=> Edge length reference
		*/
		double& at(size_t start,size_t end);
		/*
			Returns edges outgoing from specified node
			=> start Outgoing node ID
			<=> Edges outgoing from given node
		*/
		std::vector<struct edge>* get_edges(size_t start);
		/*
			Method adds polygon/triangle into this graph
			=> part Model ID, where triangle belongs
			=> polys Triangle point references
		*/
		void add_poly(size_t,std::array<size_t,3> polys);
		/*
			Method returns given triangle
			=> part Model ID
			=> index Triangle ID
		*/
		std::array<size_t,3>& poly(size_t,size_t index);
		/*
			Method returns number of triangles in the given model
			=> part Model ID
			<= Number of triangles
		*/
		size_t polys(size_t);
		/*
			Method returns collision detection algorithm of single model
			=> part Model ID
			<= Collision detection for given model
		*/
		collision_detection* get_collision_detection(size_t); 
		/*
			Method returns number of models in this graph
			<= Number of models
		*/
		size_t get_number_of_parts();
		/*
			Returns model ID from node ID
			=> node Node ID
			<= Model ID
		*/
		size_t get_part_from_node(size_t node);
		/*
			Returns normal vector of specified face/triangle
			=> part Model ID
			=> index Triangle/face ID
			<= Normal vector (length 1) of the triangle
		*/
		double* get_normal(size_t,size_t index);
		/*
			Returns angle vector between face pairs
			=> part Model ID
			=> index Triangle id
			<= Angles to each neighbouring triangle (3 elements)
		*/
		double* get_phi(size_t,size_t index);
		/*
			Returns indices to adjacent triangles
			=> part Model ID
			=> index Triangle ID to which we want adjacent triangle indices
			<= Indices to adjacent triangles
		*/
		std::array<size_t,3>& poly_neighbours(size_t,size_t index);
		/*
			This method calculates angle between two triangles
			=> part Model ID
			=> tri1 Triangle 1 ID
			=> tri2 Triangle 2 ID
			<= Angle between these two triangles
		*/
		double calculate_angle(size_t,size_t tri1, size_t tri2);
		/*
			This method calculates angle between adjoining triangles
			=> part Model ID
		*/
		void calculate_neighbours(size_t);
		/*
			This method looks for adjacent triangle to the given one
			=> part Model ID
			=> v1 Point ID, which is shared between both triangles
			=> v2 Second point ID, which is shared between both triangles
			=> exclude_poly Triangle ID, which we dont want
			(we want the second one)
			<= ID of the third point of the second triangle
		*/
		size_t poly_search(size_t,size_t v1,size_t v2,size_t exclude_poly);
		/*
			Method updates single normal for given triangle in certain model
			=> part Model ID
			=> poly_idx Triangle ID
		*/
		void update_normal(size_t part,size_t poly_idx);
		/*
			Method updates all normals of all triangles in given model
			=> part Model ID
		*/
		void update_normals(size_t part);
		/*
			Sets working positions for all nodes
			=> pos Working positions of all nodes
		*/
		void set_positions(double*);
		/*
			Sets colors of each node in this graph
			=> col Color values array (R1 G1 B1 R2 G2 B2 R3 ...)
		*/
		void set_colors(unsigned char*);
		/*
			Returns current position of a single node (constant one)
			=> idx Node ID
			<= Node current position
		*/
		const double* get_position(size_t);
		/*
			Returns current position of a single node
			=> idx Node ID
			<= Node current position
		*/
		double* change_position(size_t idx);
		/*
			Returns color of a single node
			=> idx Node ID
			<= Color of the given node (R G B)
		*/
		unsigned char* get_color(size_t);
	private:
		/*
			Method calculates volume of a single model specified by its ID
			=> part Model ID
			<= Volume of a model
			(assuming it is manifold triangle surface model)
		*/
		double calc_volume(size_t);
		/* edge dense matrix for each model */
		double** m_arr = nullptr;
		/* node array */
		point_3D* m_data;
		/* node coordinate vector (x1 y1 z1 x2 y2 z2 x3 ...) */
		double* m_data_position;
		/* color of the nodes (r1 g1 b1 r2 g2 b2 r3 ...) */
		unsigned char* m_data_color;
		/* edges outgoing from each node */
		std::vector<struct edge>* m_edges;
		/* triangle indices in each model */
		std::vector<std::array<size_t,3> >* m_polys;
		/* adjacent triangles of each triangle in each model */
		std::vector<std::array<size_t,3> >* m_poly_neighbours;
		/* normal vectors for each triangle in each model (+ angles) */
		std::vector<std::array<double,POINT_DIM+3> >* m_normals;
		/* total number of nodes */
		size_t m_nodes;
		/* number of nodes in each model */
		std::vector<size_t> m_parts;
		/* number of model in each class */
		std::vector<size_t> m_classes;
		/* vector of collision detection (for each model) */
		collision_detection** m_collision = nullptr;
		/* volume of each model */
		double* m_volume = nullptr;
		/* Variable, if something went wrong */
		static double FAILURE;
		/* Zero edge */
		static struct edge ZERO;
};
