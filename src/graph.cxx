/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################

	This class implements operations on graph structure, it
	works with nodes, edges, and even faces.	
*/
#include "point_3D.h"
#include "graph.h"
#include "collision/collision.h"
#include "collision/voxel.h"
#include <iostream>
#include <limits>
#include <cmath>
/* Variable, if something went wrong */
double graph::FAILURE = std::numeric_limits<double>::quiet_NaN();
/* Zero edge */
struct graph::edge graph::ZERO;
/*
	Graph c'tor
	=> nodes Number of nodes in each model
	=> classes Number of models in each class (bone/muscle)
	=> limit Maximum number of nodes for which dense matrix will be formed
*/
graph::graph(
	std::vector<size_t> nodes,
	std::vector<size_t> classes,
	size_t limit
){
	/* At first, we save number of nodes in each model */
	m_parts = nodes;
	/* then we save number of models in each class */
	m_classes = classes;
	/* next we acquire number of models */
	size_t len = nodes.size();
	/* we also prepare for each model a dense matrix of edge values */
	m_arr = new double*[len];
	/* and each model has its own collision detection */
	m_collision = new collision_detection*[len];
	/* Each model has its own volume, too */
	m_volume = new double[len];
	/* Next, each model has its own list of triangles */
	m_polys =           new std::vector<std::array<size_t,3> >[len];
	/* and each triangle has a list of neighbouring triangles */
	m_poly_neighbours = new std::vector<std::array<size_t,3> >[len];
	/* and precomputed normals */
	m_normals =         new std::vector<std::array<double,POINT_DIM+3> >[len];
	/* Next, we need to count all nodes in all models */
	size_t sum = 0;
	/* so we iterate over all models */
	for(size_t a=0; a<len; a++){
		/* we get number of nodes in the single model */
		size_t nds = nodes[a];
		/* and we add it to the total number of nodes */
		sum+=nds;
		/* If the current model is sufficiently small */
		if(nds<=limit){
			/* Then we create dense matrix (for edges) */
			m_arr[a] = new double[nds*nds]();
		/* else if the model is huge */
		}else{
			/* we won't create dense matrix (edges will be saved sparse) */
			m_arr[a] = nullptr;
		}
		/* and we initiate collision detection (to null, because it is no
			point to create it at this point) */
		m_collision[a] = nullptr;
	}
	/* Then, we create edge vector (sparse representation) */
	m_edges = new std::vector<struct edge>[sum]();
	/* prepare points for all nodes (velocities, forces, etc.) */
	m_data = new point_3D[sum]();
	/* and finally, we save total number of nodes */
	m_nodes = sum;
}
/*
	Graph d'tor
*/
graph::~graph(){
	/* at first, we need to know number of models stored in this graph */
	size_t len = m_parts.size();
	/* we iterate over them */
	for(size_t a=0; a<len; a++){
		/* if the current model has dense matrix allocated */
		if(m_arr[a]){
			/* we free it */
			delete[] m_arr[a];
		}
		/* if it has collision detection */
		if(m_collision[a]){
			/* we free it as well */
			delete[] m_collision[a];
		}
	}
	/* We remove array of dense matrices */
	delete[] m_arr;
	/* array of collision deteciton algorithm */
	delete[] m_collision;
	/* array of edge vectors */
	delete[] m_edges;
	/* Then, we free all point information (velocities, forces, etc.) */
	delete[] m_data;
	/* all faces */
	delete[] m_polys;
	/* precomputed normals */
	delete[] m_normals;
	/* model volumes */
	delete[] m_volume;
	/* and neighbouring faces */
	delete[] m_poly_neighbours;
}
/*
	Method creates edge and returns its value
	=> start Starting node ID
	=> stiffnes Stiffness of the new edge
	=> end Ending node ID
	<=> Edge length reference
*/
double& graph::edge(size_t start,double stiffness,size_t end){
	/* at first, we put the new edge int sparse structure */
	m_edges[start].push_back({0.,stiffness,stiffness,end});
	/* and next, we try to acquire the length from dense matrix if possible */
	return at(start,end);
}
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
void graph::boundary(size_t part,
	double& min_x,double& max_x,
	double& min_y,double& max_y,
	double& min_z,double& max_z
){
	/* at first, we initialize minimum values to double's max */
	min_x=min_y=min_z=std::numeric_limits<double>::max();
	/* and for maximum values other way around */
	max_x=max_y=max_z=std::numeric_limits<double>::lowest();
	/* next, we prepare start and end variable */
	size_t start,end;
	/* and fill it with start and end node ID of given model ID */
	nodes_from_part(part,start,end);
	/* we precalculate end index */
	end*=POINT_DIM;
	/* and then we iterate over all points in the given model */
	for(size_t a=start*POINT_DIM; a<end; a+=POINT_DIM){
		/* we get coordinate of this node */
		double* p = m_data_position+a;
		/* if X coordinate is lower than minimum, we rewrite it */
		if(p[0]<min_x) min_x=p[0];	
		/* if Y coordinate is lower than minimum, we rewrite it */
		if(p[1]<min_y) min_y=p[1];	
		/* if Z coordinate is lower than minimum, we rewrite it */
		if(p[2]<min_z) min_z=p[2];	
		/* if X coordinate is higher than maximum, we rewrite it */
		if(p[0]>max_x) max_x=p[0];	
		/* if Y coordinate is higher than maximum, we rewrite it */
		if(p[1]>max_y) max_y=p[1];	
		/* if Z coordinate is higher than maximum, we rewrite it */
		if(p[2]>max_z) max_z=p[2];	
	}
}
/*
	Method returns collision detection algorithm of single model
	=> part Model ID
	<= Collision detection for given model
*/
collision_detection* graph::get_collision_detection(size_t part){
	return m_collision[part];
}
/*
	Method calculates volume of a single model specified by its ID
	=> part Model ID
	<= Volume of a model (assuming it is manifold triangle surface model)
*/
double graph::calc_volume(size_t part){
	/* Temporary vector is created */
	double res[POINT_DIM];
	/* and volume is zero */
	double vol=0;
	/* then, we iterate over all triangles of the model */
	for(std::array<size_t,3>& poly: m_polys[part]){
		/* and get coordinates of the first point on the triangle */
		const double* p1=get_position(poly[0]);
		/* second point */
		const double* p2=get_position(poly[1]);
		/* and third point */
		const double* p3=get_position(poly[2]);
		/* then, we calculate cross product (we get volume of a parallelogram
			formed by the origin and the triagle) */
		CROSS_VEC(res,p1,p2);
		/* and we add the volume to the total */
		vol += DOT_VEC(res,p3);	
	}
	/* To be accurate, we divide the total volume by six to get actual volume,
		but teoretically speaking if we are just comparing volumes, there is
		no need to do so */
	return vol/6;
}
/*
	Precalculates model volume and saves it for later
	=> part Model ID
*/
void graph::calculate_volume(size_t part){
	m_volume[part] = calc_volume(part);
}
/*
	Returns volume of the selected model
	=> part Model ID
	<= Volume of the model
*/
double& graph::volume(size_t part){
	return m_volume[part];
}
/*
	Sets working positions for all nodes
	=> pos Working positions of all nodes
*/
void graph::set_positions(double* pos){
	/* at first, we store the pointer */
	m_data_position = pos;
	/* then, we prepare number of nodes in total */
	size_t total = nodes();
	/* and we iterate over all nodes */
	for(size_t a=0; a<total; a++){
		/* and we set its working position */
		m_data[a].set_working_position(pos+a*POINT_DIM);
	}
}
/*
	Sets colors of each node in this graph
	=> col Color values array (R1 G1 B1 R2 G2 B2 R3 ...)
*/
void graph::set_colors(unsigned char* col){
	m_data_color = col;
}
/*
	Returns current position of a single node
	=> idx Node ID
	<= Node current position
*/
double* graph::change_position(size_t idx){
	return m_data_position+POINT_DIM*idx;
}
/*
	Returns current position of a single node (constant one)
	=> idx Node ID
	<= Node current position
*/
const double* graph::get_position(size_t idx){
	return m_data_position+POINT_DIM*idx;
}
/*
	Returns color of a single node
	=> idx Node ID
	<= Color of the given node (R G B)
*/
unsigned char* graph::get_color(size_t idx){
	return m_data_color+POINT_DIM*idx;
}
/*
	Generates collision detection on a single model
	=> part Model ID
*/
void graph::generate_collision_detection(size_t part){
	/* it the collision detection object does not exist */
	if(!m_collision[part]){
		/* we create one */
		m_collision[part] = new voxel();
	}
	/* and then we just regenerate it */
	m_collision[part]->generate(part,this);
}
/*
	Returns node by its ID
	=> node Node ID
	<= Node information (velocity, force, etc.)
*/
point_3D& graph::node(size_t node){
	return m_data[node];
}
/*
	Returns stating and ending index of nodes in specified model
	=> part Model ID
	<= start Node starting ID
	<= end Node ending ID
*/
void graph::nodes_from_part(size_t part,size_t& start,size_t& end){
	/* at first, we initiate starting position to zero */
	start=0;
	/* loop control variable has to be declared here (see below) */
	size_t a=0;
	/* and we will iterate over all models before the given one */
	for(; a<part; a++){
		/* and we add number of nodes to the start */
		start+=m_parts[a];
	}
	/* if the selected model is the last one */
	if(part==m_parts.size()-1){
		/* then end is the same as total number of nodes */
		end = m_nodes;
	/* but if it is not the last one */
	}else{
		/* we can get the ending position by adding number of nodes
			in selected part to the starting position */
		end = start+m_parts[a];
	}
	
}
/*
	Returns stating and ending index of models in specified class
	(for example muscles, bones, ...)
	=> cls Class ID
	<= start Model starting ID
	<= end Model ending ID
*/
void graph::parts_from_class(size_t cls,size_t& start,size_t& end){
	/* at first, we initiate starting position to zero */
	start = 0;
	/* loop control variable has to be declared here (see below) */
	size_t a;
	/* and we will iterate over all classes before the given one */
	for(a=0; a<cls; a++){
		/* and we add number of models to the start */
		start+=m_classes[a];
	}
	/* if the selected class is the last one */
	if(a==m_classes.size()){
		/* then end is the same as total number of models */
		end = m_parts.size();
	/* but if it is not the last one */
	}else{
		/* we can get the ending position by adding number of models
			in selected class to the starting position */
		end = start+m_classes[a];
	}
}
/*
	Returns starting position and ending position of all nodes in given class
	=> cls Class ID
	<= start Node start ID
	<= end Node end ID
*/
void graph::nodes_from_class(size_t cls,size_t& start,size_t& end){
	/* at first, we prepare count of models from the start */
	size_t sum = 0;
	/* we will iterate over all classes before the given one */
	for(size_t a=0; a<cls; a++){
		/* and we add number of models to the start */
		sum+=m_classes[a];
	}
	/* then, we set starting position to zero */
	start = 0;
	/* and prepare loop control variable to zero */
	size_t a = 0;
	/* and for each model which doesn't belongs to given class */
	for(; a<sum; a++){
		/* we add number of point to the start */
		start+=m_parts[a];
	}
	/* end will be same as the start */
	end = start;
	/* but we have to add all nodes in current class */
	for(size_t b=0; b<m_classes[cls]; b++){
		/* so we do so */
		end+=m_parts[a++];
	}
}
/*
	Returns total number of nodes in the whole graph
	<= Total number of nodes
*/
size_t graph::nodes(){
	return m_nodes;
}
/*
	Returns normal vector of specified face/triangle
	=> part Model ID
	=> index Triangle/face ID
	<= Normal vector (length 1) of the triangle
*/
double* graph::get_normal(size_t part,size_t index){
	return m_normals[part][index].data();
}
/*
	Returns angle vector between face pairs
	=> part Model ID
	=> index Triangle id
	<= Angles to each neighbouring triangle (3 elements)
*/
double* graph::get_phi(size_t part,size_t index){
	return m_normals[part][index].data()+POINT_DIM;
}
/*
	Returns edge lenght between two specified nodes
	=> start Starting node ID
	=> end Ending node ID
	<=> Edge length reference
*/
double& graph::at(size_t start,size_t end){
	/* at first, we have to move to correct model, so we get number of
	nodes in each model */
	size_t len = m_parts.size();
	/* then, we have to know how many node we have already skipped */
	size_t sum = 0;
	/* for each model */
	for(size_t a=0; a<len; a++){
		/* we get number of nodes */
		size_t curr = m_parts[a];
		/* and if both start and end is in this model */
		if(start>=sum && start<sum+curr && end>=sum && end<sum+curr){
			/* we can return the edge length, if we have dense matrix */
			if(m_arr[a]){
				/* we return just single matrix element */
				return m_arr[a][start-sum+(end-sum)*curr];
			/* if we do not have dense matrix */
			}else{
				/* we have to find in linearly in the sparse structure */
				return slow_at(start,end).value;
			}
		}
		/* this was not right model, we have to count the point of this model */
		sum+=curr;
	}
	/* start and end does not address same model, which is bad */
	return FAILURE;
}
/*
	Returns edge from sparse structure (vector), this is linear, we should
	address dence matrix "m_arr", if we have enough memory
	=> start Start node ID
	=> end End node ID
	<=> Edge information
*/
struct graph::edge& graph::slow_at(size_t start,size_t end){
	/* at first, we get all edges which are associated with starting node */
	std::vector<struct graph::edge>& edges = m_edges[start];
	/* then we get number of outgoing edges */
	size_t len = edges.size();
	/* and we iterate over each one */
	for(size_t a=0; a<len; a++){
		/* if the edge is the one we search */
		if(edges[a].end==end){
			/* we return it */
			return edges[a];
		}
	}	
	/* it there is no edge, we return default one */
	return ZERO;
}
/*
	Returns model ID from node ID
	=> node Node ID
	<= Model ID
*/
size_t graph::get_part_from_node(size_t node){
	/* we get number of models */
	size_t len = m_parts.size();
	/* and we prepare variable for number of discovered points so far */
	size_t sum = 0;
	/* then, we iterate over all models */
	for(size_t a=0; a<len; a++){
		/* and we add number of point in current node to the total sum */
		sum+=m_parts[a];
		/* if we have passed the point we wanted to find */
		if(sum>node){
			/* then we return ID of the model we are in right now */
			return a;
		}
	}
	/* if we weren't successfull, we return something reasonable */
	return len-1;
}
/*
	Returns edges outgoing from specified node
	=> start Outgoing node ID
	<=> Edges outgoing from given node
*/
std::vector<struct graph::edge>* graph::get_edges(size_t start){
	return m_edges+start;
}
/*
	Method adds polygon/triangle into this graph
	=> part Model ID, where triangle belongs
	=> polys Triangle point references
*/
void graph::add_poly(size_t part,std::array<size_t,3> polys){
	/* we push the references back */
	m_polys[part].push_back(polys);
	/* and allocate some space for corresponding normal vector
		and angle between triangles */
	m_normals[part].push_back({0,0,0,0,0,0});
	/* finally we calculate the normal properly */
	update_normal(part,m_polys[part].size()-1);
}
/*
	Method updates all normals of all triangles in given model
	=> part Model ID
*/
void graph::update_normals(size_t part){
	/* at first, we get number of triangles forming the model */
	size_t count = m_polys[part].size();
	/* and for each one */
	for(size_t a=0; a<count; a++){
		/* we calculate normal vector */
		update_normal(part,a);
	}
}
/*
	Method updates single normal for given triangle in certain model
	=> part Model ID
	=> poly_idx Triangle ID
*/
void graph::update_normal(size_t part,size_t poly_idx){
	/* first of all, we get references to the 3 points forming triangle */
	std::array<size_t,3> ply = poly(part,poly_idx);
	/* then, we prepare two auxiary vectors */
	double tmp1[POINT_DIM],tmp2[POINT_DIM];
	/* next, we need coordinates of the first point */
	double* p0 = m_data[ply[0]].get_working_position();
	/* second point */
	double* p1 = m_data[ply[1]].get_working_position();
	/* and third point forming the triangle */
	double* p2 = m_data[ply[2]].get_working_position();
	/* then, we get pointer, where normal vector will be stored */
	double* nrm = get_normal(part,poly_idx);
	/* and the normal is calculated */
	NORM_TRI(nrm,tmp1,tmp2,p0,p1,p2);
}
/*
	This method calculates angle between adjoining triangles
	=> part Model ID
*/
void graph::calculate_neighbours(size_t part){
	/* at first, we acquire number of triangles in given model */
	size_t poly_count = polys(part);
	/* for each triangle */
	for(size_t a=0; a<poly_count; a++){
		/* we get 3 point indicies forming a triangle */
		std::array<size_t,3>& ply = poly(part,a);
		/* next, we find point which does not belong to the current
			triangle, but forms adjoining one */
		size_t x = poly_search(part,ply[1],ply[2],a);
		/* these points */
		size_t y = poly_search(part,ply[2],ply[0],a);
		/* are three on a manifold triangle surface */
		size_t z = poly_search(part,ply[0],ply[1],a);
		/* so now we can store point references */
		m_poly_neighbours[part].push_back({x,y,z});
		/* finally, we have to calculate angles into following memory */
		double* phi = get_phi(part,a);
		/* so we use dedicated method */
		phi[0] = calculate_angle(part,x,a);
		/* to do so */
		phi[1] = calculate_angle(part,y,a);
		/* which calculates all three angles */
		phi[2] = calculate_angle(part,z,a);
	}
}
/*
	This method calculates angle between two triangles
	=> part Model ID
	=> tri1 Triangle 1 ID
	=> tri2 Triangle 2 ID
	<= Angle between these two triangles
*/
double graph::calculate_angle(size_t part,size_t tri1, size_t tri2){
	/* this is pretty straightforward, we acquire first normal */
	double* n1 = get_normal(part,tri1);
	/* and second normal */
	double* n2 = get_normal(part,tri2);
	/* the angle is then only acos of a dot product of these two normals */
	return std::acos(DOT_VEC(n1,n2));
}
/*
	This method looks for adjacent triangle to the given one
	=> part Model ID
	=> v1 Point ID, which is shared between both triangles
	=> v2 Second point ID, which is shared between both triangles
	=> exclude_poly Triangle ID, which we dont want (we want the second one)
	<= ID of the third point of the second triangle
*/
size_t graph::poly_search(size_t part,size_t v1,size_t v2,size_t exclude_poly){
	/* so, we get number of polygons */
	size_t poly_count = polys(part);
	/* and for each one */
	for(size_t a=0; a<poly_count; a++){
		/* we check if we are not in the same triangle as the input is */
		if(a==exclude_poly){
			/* if these triangles are the same (we don't want this
			situation), we continue with next triangle */
			continue;
		}
		/* we get triangle indices */
		std::array<size_t,3>& ply = poly(part,a);
		/* if two points are shared, we found the correct triangle */
		if(
			(ply[0] == v1 && ply[1] == v2) || (ply[0] == v2 && ply[1] == v1) ||
			(ply[1] == v1 && ply[2] == v2) || (ply[1] == v2 && ply[2] == v1) ||
			(ply[2] == v1 && ply[0] == v2) || (ply[2] == v2 && ply[0] == v1)
		){
			/* so we can return this triangle ID */
			return a;
		}
	}
	/* if it fails, we return some nonsence ;) */
	return std::numeric_limits<size_t>::max();
}
/*
	Method returns given triangle
	=> part Model ID
	=> index Triangle ID
*/
std::array<size_t,3>& graph::poly(size_t part,size_t index){
	return m_polys[part].at(index);
}
/*
	Method returns number of triangles in the given model
	=> part Model ID
	<= Number of triangles
*/
size_t graph::polys(size_t part){
	return m_polys[part].size();
}
/*
	Returns indices to adjacent triangles
	=> part Model ID
	=> index Triangle ID to which we want adjacent triangle indices
	<= Indices to adjacent triangles
*/
std::array<size_t,3>& graph::poly_neighbours(size_t part,size_t index){
	return m_poly_neighbours[part].at(index);
}
/*
	Method returns number of models in this graph
	<= Number of models
*/
size_t graph::get_number_of_parts(){
	return m_parts.size();
}
/*
	Affine transforms model according to specified transformation matrix
	=> part Transformed model ID
	=> trans Affine transform matrix
*/
void graph::transform(size_t part,matrix_4x4 trans){
	/* at first, we need to know where given model starts and ends */
	size_t start,end;
	/* this we get with "nodes_from_part" method */
	nodes_from_part(part,start,end);
	/* We will traverse array of triplets */
	end*=POINT_DIM;
	/* so we iterate over all points */
	for(size_t a=start*POINT_DIM; a<end; a+=POINT_DIM){
		/* we get its position */
		double* p = m_data_position+a;
		/* and this position will be trasnformed */
		trans.transform_vector(p);
	}
	/* then we need to find nodes, which are attached to given model,
		these has to be trasformed too. So we will iterate over all points */
	for(size_t a=0; a<m_nodes; a++){
		/* if current point is attached to this model */
		if(m_data[a].get_fixed()==part){
			/* we get position of this point */
			double* x = change_position(a);
			/* and we transform it */
			trans.transform_vector(x);
		}
	}
	/* finally, we tell collision detection it should transform itself too */
	if(m_collision[part]){
		m_collision[part]->set_transform(trans);
	}
}
/*
	Method updates stiffness of each node to specified power
	=> power Power of the stiffness
*/
void graph::update_stiffness(double power){
	/* for each node in the graph */
	for(size_t a=0; a<m_nodes; a++){
		/* we get edge outgoing from each node */
		std::vector<struct edge>* edges = get_edges(a);
		/* and we get number of this edges */
		size_t len = edges->size();
		/* for every edge (in the graph) */
		for(size_t b=0; b<len; b++){
			/* we get edge structure */
			struct edge& e = edges->at(b);
			/* and set stiffness according to the power */
			e.stiffness = std::pow(e.abs_stiffness,power);
		}
	}
}
