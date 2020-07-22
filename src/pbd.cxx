/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################

	This class calculates PBD deformation. It uses data
	provided by "graph" class.
	This file is highly influened by following paper:
		M. Muller, B. Heidelberger, M. Hennix, J. Ratcliff.
		Position based dynamics. Journal of Visual Communication and Image
		Representation. vol. 18. p. 109-118.
		https://doi.org/10.1016/j.jvcir.2007.01.005 (2007)
*/
#include <cmath>
#include <chrono>
#include <iomanip>
#include <chrono>
#include <thread>
#include <cstring>
#include <algorithm>
#include "pbd.h"
#define MUSCLE_CLASS 0
#define BONE_CLASS 1
/* Loops over all verticies in graph, don't forget to close it with END macro */
#define FOREACH_VERTEX_IN_GRAPH_AS(PNT,IDX) \
	for(size_t (IDX)=m_start; (IDX)<m_end; (IDX)++) {\
		point_3D& (PNT) = m_graph->node(IDX);
/* Loops over all verticies in a signle model, end it with END macro */
#define FOREACH_VERTEX_IN_PART_AS(PART,PNT,IDX) \
	size_t strt,nd; \
	m_graph->nodes_from_part(PART,strt,nd);\
	for(size_t (IDX)=strt; (IDX)<nd; (IDX)++) {\
		[[maybe_unused]] point_3D& (PNT) = m_graph->node(IDX);
/* Loops over all adjacent triangles, end it with END_NEIGHBOURS */
#define FOREACH_POLY_NEIGHBOURS(PART,CURR,IDX1,IDX2) \
	size_t poly_count = m_graph->polys((PART));\
	for(size_t (IDX1)=0; (IDX1)<poly_count; (IDX1)++){\
		std::array<size_t,3>& idxs = m_graph->poly_neighbours((PART),(IDX1));\
		for(size_t (CURR)=0; (CURR)<3; (CURR)++){ \
				size_t (IDX2) = idxs[(CURR)]; \
				std::array<size_t,3>& idxs_neigh = \
					m_graph->poly_neighbours((PART),(IDX2));\
				double phi = m_graph->get_phi((PART),(IDX1))[(CURR)];
/* End of FOREACH_VERTEX_IN_GRAPH_AS and FOREACH_VERTEX_IN_PART_AS macros */
#define END }
/* End of FOREACH_POLY_NEIGHBOURS macro */
#define END_NEIGHBOURS }}
/* PI constant (this should be enought precision) */
constexpr double PI = 3.141592653589793238463;
/*
	Method inits PBD algorithm
	=> grp Graph containing structured data for PBD algorith,
	=> data Raw data
	=> scenario ID of the dynamic scenario (it has to be implemented this way,
		we can force user to write his own scenario and compile it, but who
		would do that. We also can create a scripting language, but who would
		like to learn in, so I dont care too
*/
void pbd::init(graph* grp,data* dta,size_t scenario){
	/* at first, we retype function pointer to something more readable */
	typedef void (pbd::*scenarios)(size_t);
	/* and now, we create array of function pointers, which point to all
		implemented scenarios */
	scenarios scrs[] = {
		&pbd::scenario_0,
		&pbd::scenario_1,
		&pbd::scenario_2,
		&pbd::scenario_3,
		&pbd::scenario_4,
		&pbd::scenario_5,
		&pbd::scenario_6
	};
	/* at this point, we select scenario and save its function pointer */
	m_scenario = scrs[scenario];
	/* then, we save structured data */
	m_graph = grp;
	/* and raw data, too */
	m_data = dta;
}
/*
	Returns structured data (graph) associated with this PBD algorithm
	=> Graph/structured data
*/
graph* pbd::get_graph(){
	return m_graph;
}
/*
	This method does one "outer" iteration of the PBD algorithm
	=> iter Iteration number
	<= Return state (0 - OK, iteartion has been done, else FAIL)
*/
int pbd::execute(size_t iter){
	/* if the anisotropy has been changed between iterations */
	if(std::abs(m_anisotropy-m_last_anisotropy)>1e-9){
		/* we save the last one */
		m_last_anisotropy = m_anisotropy;
		/* and recalculate stiffnes according to anisotropy level */
		m_graph->update_stiffness(m_anisotropy);
	}
/* In case we do not have HUD (Head-Up Display,
	in this case sliders and labels) */
#ifdef NO_HUD
	/* and if we are in first iteration */
	if(!iter){
		/* we have to start it immediately, there is no way user can start
			the simulation if no slider is visible */
		m_pause=0;
	}
#endif
	/* If pause is set */
	if(m_pause>0.5){
		/* we sleep so FPS will be slightly less than 30 FPS */
		std::this_thread::sleep_for(std::chrono::milliseconds(33));
		/* and return that we does not anything useful :) */
		return 1;
	}
	/* at first, we get time at the start of this iteration */
	CAPTURE_TIME(t1);
	/* then, we check if we are allowed to do this iteration */
	if(iter>=m_iter){
		/* if we are not, we end with error code */
		return 1;
	}
	/* at this point, we need node IDs for muscles */
	m_graph->nodes_from_class(MUSCLE_CLASS,m_start,m_end);
	/* so model IDs for muscles */
	m_graph->parts_from_class(MUSCLE_CLASS,m_muscle_start,m_muscle_end);
	/* and model IDs for bones */
	m_graph->parts_from_class(BONE_CLASS,m_bone_start,m_bone_end);
	/* if we are in first iteration */
	if(!iter){
		/* we initialize the PBD (precompute everything possible, that's why
			first iteration take so long) */
		initialize_first_iteration();
	}
	/* then, we specify gravity force, which is exponential function of slider
		value */
	double grav = std::pow(10,m_gravity)-1;
	/* for each muscle node */
	for(size_t a=m_start; a<m_end; a++){
		/* we get force applied to this node */
		double* f = m_graph->node(a).get_force();
		/* and set this force to the gravity */
		SET_VEC(f,0,0,-grav);
	}
	/* follows scenario specific commands */
	(this->*m_scenario)(iter);
	/* now, we can do inner interation (gauss-seidel sovler) */
	for(subiter=0; subiter<m_constr; subiter++){
		/* we get time before we apply forces */
		CAPTURE_TIME(t2);
		/* then, we move all nodes by given external forces */
		move_points_by_force();
		/* and get spend time */
		CAPTURE_TIME(t3);
		/* if we can apply vertex (distance) constraints */
		if(!m_no_edge){
			/* we apply them */
			project_vertex_constraints();
		}
		/* and calculate spend time of this operation */
		CAPTURE_TIME(t41);
		/* if we are allowed to preserve model volume */
		if(!m_no_volume){
			/* we prepare array for each model (PROBABLY PREALLOCATE !!!)*/
			double* numerators = new double[m_muscle_end-m_muscle_start];
			/* this array is passed, so values of this array are being set */
			calculate_volume_constraints(numerators);
			/* then, we project this constraint (change node positions) */
			project_volume_constraints(numerators);
			/* and we have to clean memory */
			delete[] numerators;
		}
		/* we once more calculate, how long does it take to calculate
			volume constraint */
		CAPTURE_TIME(t42);
		/* if we can preserve model local shape */
		if(!m_no_bending){
			/* we do so */
			project_bending_constraints();
		}
		/* and save time, which did it take */
		CAPTURE_TIME(t43);
		/* Now we prepare structure, which contains various information
			about collisions occured */
		collisions coll;
		/* if we can avoid collisions */
		if(!m_no_collision){
			/* we find out where collisions occur */
			find_collisions(coll);
		}
		/* and we save current time, so we know how long did it take to
			find collisions */
		CAPTURE_TIME(t44);
		/* if we can avoid collisions */
		if(!m_no_collision){
			/* we try to fix all collisions */
			project_collision_constraints(coll);
		}
		/* and get time of collision fixing */
		CAPTURE_TIME(t5);
		/* then, we update all nodes velocities */
		update_velocities();
/* if we are interested in time measurement */
#ifdef TIME_MEASUREMENT
		/* we get final time */
		CAPTURE_TIME(t6);
		/* we print out iteration ID */
		std::cout << "ITERATION " << iter << " (" << subiter << "/" <<
			m_constr << ")" << std::endl;
		/* and also print various time differences */
		PRINT_DIFF("Initialization:           ",tmp1 ,t2 ,t1 );
		PRINT_DIFF("New position calculation: ",tmp2 ,t3 ,t2 );
		PRINT_DIFF("Applying constraints:     ",tmp4 ,t5 ,t3 );
		PRINT_DIFF(" |- Edge constraints:     ",tmp41,t41,t3 );
		PRINT_DIFF(" |- Volume constraint:    ",tmp42,t42,t41);
		PRINT_DIFF(" |- Bending constraint:   ",tmp43,t43,t42);
		PRINT_DIFF(" |- Collision constraint: ",tmp44,t5 ,t43);
		PRINT_DIFF("     | - Find:            ",tmp45,t44,t43);
		PRINT_DIFF("     | - Correct:         ",tmp5 ,t5, t44);
		PRINT_DIFF("Update velocities:        ",tmp6 ,t6 ,t5 );
#endif
	}
	/* one iteration has been done successfully */
	return 0;
}
/*
	D'tor
*/
pbd::~pbd(){
	/* We have to delete graph, we adpoted it :) */
	delete m_graph;
}
/*
	This method is called during first iteration
	It is supposed to precalculate as much as possible
	This is the reason first iteration take so much time
*/
void pbd::initialize_first_iteration(){
	/* at first, we set correct colors */
	reset_colors();
	/* then, we iterate over all bone models */
	for(size_t a=m_bone_start; a<m_bone_end; a++){
		/* and generate collision detection */
		m_graph->generate_collision_detection(a);
	}
	/* next we iterate over all muscle models */
	for(size_t a=m_muscle_start; a<m_muscle_end; a++){
		/* we calculate its initial volume */
		m_graph->calculate_volume(a);
		/* and adjacent triangles */
		m_graph->calculate_neighbours(a);
	}
	/* next, we get pointer to closest point (for each fibre node
		to a signle node on model surface) */
	long long int* closest = m_data->get_closest_points();
	/* for each fibre node */
	FOREACH_POINT(m_data,points_lines){
		/* we prepare minimum index */
		size_t min_index = MAX;
		/* and minimum distance */
		double min_dist = std::numeric_limits<double>::max();
		/* for each "muscle" node */
		for(size_t a=m_start; a<m_end; a++) {
			/* we get position of this node */
			const double* ptl = m_graph->get_position(a);
			/* prepare temporary vector */
			double diff[3];
			/* which we fill with difference vector */
			SUB_VEC(diff,pt,ptl);
			/* we then calculate distance between muscle and fibre node */
			double dist = DOT_VEC(diff,diff);	
			/* it it is the smallest so far */
			if(dist<min_dist){
				/* we have to rewrite smallest distance */
				min_dist=dist;
				/* and we save index of this node */
				min_index=a;
			}
		} 
		/* finally, we can save index of closest point into array */
		closest[current_point/3]=min_index;
	} END_FOREACH
}
/*
	This method moves each point according to external force
*/
void pbd::move_points_by_force(){
/* if we want parallel computing (for huge muscles) */
#ifdef PARALLEL
/* we calculate it in parallel */
#pragma omp parallel for
#endif
	/* for each vertex in graph */
	FOREACH_VERTEX_IN_GRAPH_AS(vertex,vertex_index){
		/* we change velocity with respect to external forces */
		update_velocity_with_external_forces(vertex);
		/* then we apply force loss */
		damp_velocity(vertex);
		/* next, we get actual position */
		const double* x = m_graph->get_position(vertex_index);
		/* and compute destination position */
		calculate_new_position(x,vertex);
	}END
}
/*
	This method finds all colisions between bones and muscles
	<=> Information about collisions
*/
void pbd::find_collisions(collisions& coll){
	/* for each muscle model */
	for(size_t part=m_muscle_start; part<m_muscle_end; part++){
		/* we find collision with one muscle only */
		find_collision_constraints(coll,part,m_bone_start,m_bone_end);
	}
}
/*
	Method projects (corrects) edge length constraints (preserves distance)
*/
void pbd::project_vertex_constraints(){
	std::vector<double> retval;
/* if we want parallel computing (for huge muscles) */
#ifdef PARALLEL
/* we calculate it in parallel */
#pragma omp parallel for 
#endif
	/* for each vertex in graph */
	FOREACH_VERTEX_IN_GRAPH_AS(vertex,vertex_index){
		/* we fix distance between all points */
		project_vertex_constraints(vertex,vertex_index,retval);
	} END
	if(m_measurement_vertex){
		std::sort(retval.begin(),retval.end());
		size_t l = retval.size();
		if(subiter==0)std::cout << retval[l/10] << "," << retval[l/4] << "," <<
			retval[l/2] << "," << retval[3*l/4] << "," << retval[9*l/10] << std::endl;
	}
}
/*
	Method projects (corrects) collisions
	=> coll Information about collisions
*/
void pbd::project_collision_constraints(collisions& coll){
	/* at first, we need to know number of collisions occurs */
	size_t len = coll.size();
/* if we want parallel computing (for huge muscles) */
#ifdef PARALLEL
/* we calculate it in parallel */
#pragma omp parallel for
#endif
	/* for each collision */
	for(size_t a=0; a<len; a++){
		/* we get target point, where node wanted to go */
		double* p = m_graph->node(coll[a].vec_dynamic).get_working_position();
		/* and we pass it into correction method */
		pbd::constraints::project_collision_constraint(p,coll[a]);
	}
}
/*
	Method calculates volume constraints (calculates current volume)
		=> numerators Memory, where volume numerators can be stored
*/
void pbd::calculate_volume_constraints(double* numerators){
	/* for each muscle model */
	for(size_t part=m_muscle_start; part<m_muscle_end; part++){
		/* we get single space for one numerator */
		double& num = numerators[part-m_muscle_start];
		/* we init it to zero */
		num = 0;
		/* and we get number of polygons/triangles in current model */
		size_t poly_count = m_graph->polys(part);
/* if we want parallel computing (for huge muscles) */
#ifdef PARALLEL
/* we calculate it in parallel */
#pragma omp parallel for reduction(+:num) firstprivate(part)
#endif
		/* for each triangle */
		for(size_t idx=0; idx<poly_count; idx++){
			/* we calculate volume of the parallelogram */
			calculate_volume_constraint(idx,part,num);
		}
		/* volume correction from parallellogram to tetrahedron */
		num/=6;
	}
}
/*
	This method project volume constraints to each muscle node
	=> numerator Volume numerator for each model
*/
void pbd::project_volume_constraints(double* numerators){
	/* for each muscle model */
	for(size_t part=m_muscle_start; part<m_muscle_end; part++){
		/* we prepare volume denominator */
		double denominator=0;
		/* for each node in this model */
		FOREACH_VERTEX_IN_PART_AS(part,vertex,vertex_index){
			/* we add accumulated value */
			double* accu = vertex.accumulator();
			/* to the denominator */
			denominator+=DOT_VEC(accu,accu);
		} END
		/* numerator is just the computed one compared with 
			pressure and original volume */
		double num = (numerators[part-m_muscle_start]
			-m_pressure*m_graph->volume(part))/denominator;
		if(m_measurement_volume){
			if(subiter==0)std::cout << numerators[part-m_muscle_start]/m_graph->volume(part) << std::endl;
		}
		{
			/* at this point, we can project volume constraint of
				each node individually */
			FOREACH_VERTEX_IN_PART_AS(part,vertex,vertex_index){
				pbd::constraints::project_volume_constraint(
					vertex,
					num
				);
			} END
		}
	}
}
/*
	Projects bending contraint (preserving local shape) to each node
	of each muscle model
*/
void pbd::project_bending_constraints(){
	std::vector<double> retval;
	/* for each muscle model */
	for(size_t part=m_muscle_start; part<m_muscle_end; part++){
		/* we update normals of each triangle in given muscle model */
		m_graph->update_normals(part);
		/* for each adjacent triangle */
		FOREACH_POLY_NEIGHBOURS(part,curr,idx1,idx2){
			/* we prepare index, which tells us which node is in the other
				side */
			size_t second_side = 0;
			/* if the 2. node (with index 1) is the same as triangle index */
			if(idxs_neigh[1] == idx1){
				/* we know that the 2. node is the same (index 1) */
				second_side = 1;
			/* or if the 3. node (with index 2) is the same as triangle index */
			}else if(idxs_neigh[2] == idx1){
				/* finally we know that the 3. node is the same (index 2) */
				second_side = 2;
			}
			/* if the triangle index is lower than the other triangle (this
				forces only one check for each pair) and if both indices are
				lower than number of triangles (overflow check) */
			if(idx1>idx2 && idx1<poly_count && idx2<poly_count){
					/* then we can project this constraint on each triangle */
					pbd::constraints::project_bending_constraint(
						m_graph,
						phi,
						curr,
						second_side,
						part,
						idx1,
						idx2,
						m_bending,
						retval
					);
			}
		} END_NEIGHBOURS
	}
	if(m_measurement_bending){
		std::sort(retval.begin(),retval.end());
		size_t l = retval.size();
		if(subiter==0)std::cout << retval[l/10] << "," << retval[l/4] << "," <<
			retval[l/2] << "," << retval[3*l/4] << "," << retval[9*l/10] << std::endl;
	}
}
/*
	This method update velocities of each node in each muscle model
*/
void pbd::update_velocities(){
/* if we want parallel computing (for huge muscles) */
#ifdef PARALLEL
/* we calculate it in parallel */
#pragma omp parallel for
#endif
	/* for each vertex in graph */
	FOREACH_VERTEX_IN_GRAPH_AS(vertex,vertex_index){
		/* we at first check, if we  can move with this point */
		if(!m_graph->node(vertex_index).is_fixed()){
			/* we need source position */
			double* pos = m_graph->change_position(vertex_index);
			/* so we can update velocity */
			update_velocity(pos,vertex);
			/* and we can also update position */
			update_position(pos,vertex);
		}
	}END
}
/*
	Calculate volume constraint of single triangle
	=> idx Triangle ID
	=> part Model ID
	<=> numerator Volume numerator
*/
void pbd::calculate_volume_constraint(
	size_t idx,
	size_t part,
	double& numerator
){
	/* at first, we get triangle indices */
	std::array<size_t,3>& poly = m_graph->poly(part,idx);
	/* then, we get coordinates of the first */
	const double* p0 = m_graph->node(poly[0]).get_working_position();
	/* second */
	const double* p1 = m_graph->node(poly[1]).get_working_position();
	/* and third triangle forming node */
	const double* p2 = m_graph->node(poly[2]).get_working_position();
	/* then we prepare temporary vector */
	double tmp[POINT_DIM];
	/* which we fill with cross product of one of the point pair */
	CROSS_VEC(tmp,p0,p1);
	/* this value is cumulated into node, which has not been involved */
	m_graph->node(poly[2]).accumulate(tmp);
	/* and this value is added to numerator
		(assuming point with index 0 is the current one) */
	numerator+=DOT_VEC(tmp,p2);
	/* then, we calculate cross product vector for edge 0-2 */
	CROSS_VEC(tmp,p2,p0);
	/* same as above, acumulate to node 1, because we've calculated edge 0-2 */
	m_graph->node(poly[1]).accumulate(tmp);
	/* finally, we calculate cross product for edge 1-2 */
	CROSS_VEC(tmp,p1,p2);
	/* once more, this value is cumulated into node, which has not
		been involved */
	m_graph->node(poly[0]).accumulate(tmp);
}
/*
	Method updates velocity with external force applied to a single point
	=> vertex Vertex on which the force is applied
*/
void pbd::update_velocity_with_external_forces(point_3D& vertex){
	/* at first, we calculate weight (it depends on time difference too) */
	double w = vertex.inverse_mass()*m_time_diff;
	/* then, we can get current velocity */
	double* v = vertex.get_velocity();
	/* and force */
	double* f = vertex.get_force();
	/* We prepare temporary vector */
	double tmp[POINT_DIM];
	/* and this vector is filled with weighted force - velocity */
	MULT_VEC(tmp,f,w);
	/* this velocity is added to the node */
	ADD_VEC(v,v,tmp);
}
/*
	Performs velocity loss
	=> vertex Node, which loses a bit of velocity
*/
void pbd::damp_velocity(point_3D& vertex){
	/* so we get velocity of the given point */
	double* v = vertex.get_velocity();
	/* and multiply it by a constant, that's it ;) */
	MULT_VEC(v,v,m_damp);
}
/*
	Method calculates new destination position according to its velocity
	=> x Current position
	=> vertex Moving node with current and destination position
*/
void pbd::calculate_new_position(const double* x,point_3D& vertex){
	/* at first, we get node velocity */
	double* v = vertex.get_velocity();
	/* then we need one auxiary vector */
	double tmp[POINT_DIM];
	/* which will be filled with direction of the node */
	MULT_VEC(tmp,v,m_time_diff);
	/* this direction is added to the current position, so we get destination */
	ADD_VEC(tmp,tmp,x);
	/* and destination is simply copied into node */
	vertex.set_working_position(tmp);
}
/*
	Finds collision for each node in a signle muscle model
	<=> coll Information about collision
	=> muscle_part Muscle model ID
	=> start First bone model ID
	=> end Last bone model ID
*/
void pbd::find_collision_constraints(
	collisions& coll,
	size_t muscle_part,
	size_t start,
	size_t end
){
	/* for each muscle model */
	for(size_t a=start; a<end; a++){
		/* we calculate collision with single bone */
		find_collision_constraints_one_bone(coll,muscle_part,a);
	}
}
/*
	Calculates collision between single muscle and single bone
	<=> coll Information about collisions occured
	=> muscle_part Muscle model ID
	=> bone_part Bone model ID
*/
void pbd::find_collision_constraints_one_bone(
	collisions& coll,
	size_t muscle_part,
	size_t bone_part
){
	/* first of all, we get collision detection algorithm associated
		with current bone model */
	collision_detection* coldet = m_graph->get_collision_detection(bone_part);
	/* then, we will need node range */
	size_t start,end;
	/* of the current muscle */
	m_graph->nodes_from_part(muscle_part,start,end);
	/* we won't do anything fancy here, we just go through all nodes
		(no spartial division etc.) */
	for(size_t a=start; a<end; a++){
		/* if this node can be moved */
		if(!m_graph->node(a).is_fixed()){
			/* we prepare vector for new position */
			double new_pos[3];
			/* if the point moved through bone */
			if(coldet->find(
				m_graph->get_position(a),
				m_graph->node(a).get_working_position(),
				new_pos
			)){
/* if we want parallel computing (for huge muscles) */
#ifdef PARALLEL
/* we have to be sure we won't add to vector at same time
	from multiple threads */
#pragma omp critical
#endif
				/* we find the best position, so it won't collide */
				coll.push_back({muscle_part,bone_part,a,
					new_pos[0],new_pos[1],new_pos[2]});
				/*unsigned char* c = m_graph->get_color(a);
				c[0]=0;
				c[1]=0;
				c[2]=255;*/
			}
		}
	}
}
/*
	This method resets color of all models
	Muscles will be red, bones light gray (nearly white)
*/
void pbd::reset_colors(){
	/* for each muscle model */
	for(size_t part=m_muscle_start; part<m_muscle_end; part++){
		/* for each node in muscle model */
		FOREACH_VERTEX_IN_PART_AS(part,pnt,a){ 
			/* we get color information */
			unsigned char* clr = m_graph->get_color(a);
			/* and we change it */
			clr[0] = 255;
			/* to color (255,0,0) */
			clr[1] = 0;
			/* which is RED */
			clr[2] = 0;
		} END
	}
	/* for each bone model */
	for(size_t part=m_bone_start; part<m_bone_end; part++){
		/* for each node in bone model */
		FOREACH_VERTEX_IN_PART_AS(part,pnt,a){ 
			/* we get color information */
			unsigned char* clr = m_graph->get_color(a);
			/* and we change it */
			clr[0] = 200;
			/* to color (200,200,200) */
			clr[1] = 200;
			/* which is LIGHT GRAY */
			clr[2] = 200;
		} END
	}
}
/*
	Method project vertex constraints (edge lenght preservation) to
	each edge outgoing from specified node
	=> vertex Point, which has outgoing edges
	=> index Muscle model ID
*/
void pbd::project_vertex_constraints(point_3D& vertex,size_t index,std::vector<double>& retval){
	/* at first, we get all edges we are interested in */
	std::vector<struct graph::edge>* edges = m_graph->get_edges(index);
	/* we get its count */
	size_t count = edges->size();
	/* and for each one */
	for(size_t a=0; a<count; a++){
		/* we get target point index */
		size_t target_index = edges->at(a).end;
		/* then, we get stiffness, if anisotropy is concerned, it has
			been already calculated, else it is just 1 */
		double stiffness = m_no_anisotropy?1:edges->at(a).stiffness;
		/* we get target node */
		point_3D& target = m_graph->node(target_index);
		/* and we get initial length of the edge */
		double len = m_graph->at(index,target_index);
		/* now we can project distance to the single edge */
		pbd::constraints::project_distance_constraint(
			vertex,
			target,
			len,
			stiffness,
			retval
		);
	}
}
/*
	Method updates velocity according to origin and target position
	=> x Origin position
	=> vertex Moving node
*/
void pbd::update_velocity(const double* x,point_3D& vertex){
	/* at first, we get target position */
	double* px = vertex.get_working_position();
	/* and current velocity */
	double* v = vertex.get_velocity();
	/* we can subtract positions to get "velocity" */
	SUB_VEC(v,px,x);
	/* but for velocity we have to divide it with time difference */
	DIV_VEC(v,v,m_time_diff);
}
/*
	Method updates target position of the node
	=> pos New target position
	=> vertex Node to change its target position
*/
void pbd::update_position(double* pos,point_3D& vertex){
	/* so we get memory, where old target position is stored */
	double* work = vertex.get_working_position();
	/* and we changed it according to the input */
	memcpy(pos,work,sizeof(double)*POINT_DIM);	
}

void pbd::set_parameters(std::vector<double>& params){
	double* arr[] = {
		&m_time_diff,&m_damp,&m_bending,&m_pressure,
		&m_constr,&m_gravity,&m_pause,&m_anisotropy
	};
	size_t len = sizeof(arr)/sizeof(arr[0]);
	if(params.size()<len){
		len=params.size();
	}
	for(size_t a=0; a<len; a++){
		*(arr[a]) = params[a];	
	}
}
/* CONSTRAINTS SECTION BEGINS - THERE ARE METHODS WHICH
	SOLVES CONSTRAINTS OF SINGLE ELEMENTS */

/*
	This method compute volume constraint for single point
	=> pnt Point, which target position should be fixed
	=> mult Volume of the tetrahedron
*/
void pbd::constraints::project_volume_constraint(
	point_3D& pnt,
	double mult
){
	/* at first, we get memory where target position is stored */
	double* p = pnt.get_working_position();
	/* we get accumulated position */
	double* accu = pnt.accumulator();
	/* we multiply accumulator with volume */
	MULT_VEC(accu,accu,mult);
	/* and subtract the target position by this vector */
	SUB_VEC(p,p,accu);
	/* finally, we do not need accumulated value */
	pnt.reset_accumulator();
}
/*
	This method projects bending constraint to a single triangle
	=> graph Structured data about all models
	=> phi_0 Original dihedral angle between two faces
	=> n_th_neightbour ID of the node, which is different on first triangle
		considering second triangle
	=> second_side_neighbourID of the node, which is different on second
		triangle considering first triangle
	=> part Muscle model ID
	=> first First polygon ID
	=> second Second polygon ID
	=> bending Bending factor
*/
void pbd::constraints::project_bending_constraint(
	graph* graph,
	double& phi_0,
	size_t& n_th_neighbour,
	size_t& second_side_neighbour,
	size_t& part,
	size_t& first,
	size_t& second,
	double& bending,
	std::vector<double>& retval
){
	/* if 0th neighbour, then share points 1,2 */
	/* if 1st neighbour, then share points 0,2 */
	/* if 2nd neighbour, then share points 0,1 */
	/* at first we prepare bunch of temporary vectors */
	double q1[POINT_DIM],q2[POINT_DIM],q3[POINT_DIM],q4[POINT_DIM];
	/* we get first triangle indices */
	std::array<size_t,3>& vecs1 = graph->poly(part,first);
	/* and second too */
	std::array<size_t,3>& vecs2 = graph->poly(part,second);
	/* there we acquire first point (which have both triangles in common) */
	point_3D& p1_p = graph->node(vecs1[(n_th_neighbour+1)%3]);
	/* second point (which have both triangles in common too) */
	point_3D& p2_p = graph->node(vecs1[(n_th_neighbour+2)%3]);
	/* third point (this belong only to second triangle) */
	point_3D& p3_p = graph->node(vecs2[second_side_neighbour]);
	/* and fourth point (this belongs only to first triangle) */
	point_3D& p4_p = graph->node(vecs1[n_th_neighbour]);
	/* now, we need th get position of the first */
	double* p1_abs = p1_p.get_working_position();
	/* second */
	double* p2_abs = p2_p.get_working_position();
	/* third */
	double* p3_abs = p3_p.get_working_position();
	/* and fourth point */
	double* p4_abs = p4_p.get_working_position();
	/* We prepare vectors for moved points by -p1 */
	double p2[POINT_DIM],p3[POINT_DIM],p4[POINT_DIM];
	/* And we calculate relative coordinates of the second point to
		the first one  */
	SUB_VEC(p2,p2_abs,p1_abs);
	/* the same with third */
	SUB_VEC(p3,p3_abs,p1_abs);
	/* and fourth point */
	SUB_VEC(p4,p4_abs,p1_abs);
	/* at this point, we need to calculate normals of relative coordinates */
	double n1[POINT_DIM],n2[POINT_DIM];
	/* so this can be done using only two points (first is zero) */
	CROSS_VEC(n1,p2,p3);
	/* we need to get lenght of this vector */
	double cross23 = std::sqrt(DOT_VEC(n1,n1));
	/* and make this normal of lenght one */
	DIV_VEC(n1,n1,cross23);
	/* the same with second normal, we use only two points */
	CROSS_VEC(n2,p2,p4);
	/* calculate its length */
	double cross24 = std::sqrt(DOT_VEC(n2,n2));
	/* and divide the normal by its length */
	DIV_VEC(n2,n2,cross24);
	/* now, we calculate "angle" between both triangles */
	double dot = DOT_VEC(n1,n2);
	/* if this is nonsence, or nearly coplanar */
	if(dot>1-1e-10 || dot < -1+1e-10){
		/* then this task has bad conditionality, so we rather do anything */
		return;
	}
	/* follows computation from the PBD paper
		(referenced on top of this source file) */
	/* we begin with q3 */
	CROSS_VEC(q3,p2,n2);
	CROSS_VEC(q1,n1,p2);
	MULT_VEC(q1,q1,dot);
	ADD_VEC(q3,q3,q1);
	DIV_VEC(q3,q3,cross23);
	/* then q4 follows */
	CROSS_VEC(q4,p2,n1);
	CROSS_VEC(q1,n2,p2);
	MULT_VEC(q1,q1,dot);
	ADD_VEC(q4,q4,q1);
	CROSS_VEC(q1,p2,p4);
	DIV_VEC(q4,q4,cross24);
	/* next is q2 */
	CROSS_VEC(q2,p3,n2);
	CROSS_VEC(q1,n1,p3);
	MULT_VEC(q1,q1,dot);
	ADD_VEC(q2,q2,q1);
	CROSS_VEC(q1,p2,p3);
	DIV_VEC(q2,q2,std::sqrt(DOT_VEC(q1,q1)));
	double tmp[POINT_DIM];
	CROSS_VEC(tmp,p4,n1);
	CROSS_VEC(q1,n2,p4);
	MULT_VEC(q1,q1,dot);
	ADD_VEC(tmp,tmp,q1);
	CROSS_VEC(q1,p2,p4);
	DIV_VEC(tmp,tmp,std::sqrt(DOT_VEC(q1,q1)));
	ADD_VEC(q2,q2,tmp);
	/* the following line (only) helps compute q1 */
	SUB_VEC(q1,q2,q3);
	MULT_VEC(q2,q2,-1);
	/* and finally q1 */
	SUB_VEC(q1,q1,q4);
	retval.push_back(std::abs(std::acos(dot)-(PI-phi_0)));
	/* we compute numerator (see PBD paper included) */
	double num = std::sqrt(1-dot*dot)*(std::acos(dot)-(PI-phi_0));
	/* and denominator (we make sure it is not zero) */
	double denom = MAX2(1e-4,
		DOT_VEC(q1,q1)+DOT_VEC(q2,q2)+DOT_VEC(q3,q3)+DOT_VEC(q4,q4)
	);
	/* then we need weighted sum, so we add up all masses together */
	double sum_mass = p1_p.inverse_mass()+p2_p.inverse_mass()+
		p3_p.inverse_mass()+p4_p.inverse_mass();
	/* and calculate constant ratio, according to the weighed sum */
	double ratio = -4*num*bending/(denom*sum_mass);
	/* if this changes behaviour drastically */
	if(ratio<1e-1 || ratio >1e1){
		/* we'd rather let it be, do not touch it it it works ;) */
		return;
	}
	/* now, we finally finish weighted sum for first point */
	MULT_VEC(q1,q1,ratio*p1_p.inverse_mass());
	/* second */
	MULT_VEC(q2,q2,ratio*p2_p.inverse_mass());
	/* third */
	MULT_VEC(q3,q3,ratio*p3_p.inverse_mass());
	/* and fourth too */
	MULT_VEC(q4,q4,ratio*p4_p.inverse_mass());
	/* and at the end of the day we can change node coordinates */
	ADD_VEC(p1_abs,p1_abs,q1);
	/* second */
	ADD_VEC(p2_abs,p2_abs,q2);
	/* third */
	ADD_VEC(p3_abs,p3_abs,q3);
	/* and fourth, too */
	ADD_VEC(p4_abs,p4_abs,q4);
}
/*
	Project edge length preserving constraint on a pair of nodes
	=> source Source node ID (from edge point of wiew)
	=> target Destination node ID (from edge point of view)
	=> len Original length of this edge
	=> stiffness Stiffness of the edge
*/
void pbd::constraints::project_distance_constraint(
	point_3D& source,
	point_3D& target,
	double len,
	double stiffness,
	std::vector<double>& retval
){
	/* at first, we get source node working position */
	double* sx = source.get_working_position();
	/* and target node working position */
	double* tx = target.get_working_position();
	/* then we need to know first node mass */
	double weight_source = source.inverse_mass();
	/* second one too */
	double weight_target = target.inverse_mass();
	/* we will do weighted sum, so we need to know total weight */
	double weight_sum = weight_source+weight_target;
	/* it weight is too small */
	if(weight_sum<1e-5){
		/* we will divide by "zero", so we won't rather do anything */
		return;
	}
	/* we need two auxiary vectors */
	double dx[POINT_DIM],ratio_dx[POINT_DIM];
	/* we fill one with difference between two point coordinates */
	SUB_VEC(dx,sx,tx);
	/* and then we calculate current length of the edge
		(squared - dot product) */
	double curr_len_squared = DOT_VEC(dx,dx);
	/* now, we calculate constant multiplicator (this is independent to
		node masses */
	retval.push_back(std::sqrt(curr_len_squared)/len);
	double ratio = (1-len/std::sqrt(curr_len_squared))*stiffness;
	/* now we multiply edge length by this factor (it gets shorter/longer) */
	MULT_VEC(ratio_dx,dx,ratio);
	/* now we calculate weight of first node */
	double normalized_source_weight = weight_source/weight_sum;
	/* and the second node */
	double normalized_target_weight = weight_target/weight_sum;
	/* Now, we perform weighted sum for first node */
	MULT_VEC(dx,ratio_dx,normalized_source_weight);
	/* and this we subtract from first node position */
	SUB_VEC(sx,sx,dx);
	/* same for second node */
	MULT_VEC(dx,ratio_dx,normalized_target_weight);
	/* but now we have to add it (we have only one vector, but we need opposite
		direction) */
	ADD_VEC(tx,tx,dx);
}
/*
	Method projects collision constraint on single node
	=> position Target position of the node
	=> col Information about collision which occurs
*/
void pbd::constraints::project_collision_constraint(
	double*	position,
	struct collision& col
){
	/* only what we will do is set the position to the best
		one collision detection found */
	SET_VEC(position,col.x,col.y,col.z);
}
/* SCENARIOS SECTION, THERE ARE HARD CODED DYNAMIC SCENARIOS, IF YOU 
	WANT TO ADD ONE, FEEL FREE TO, BUT DON'T FORGET TO ADD FUNCTION
	POINTER INTO "init" METHOD */

/*
	Default, empty scenario, does nothing usefull
	=> iter Iteration value
*/
void pbd::scenario_0([[maybe_unused]] size_t iter){}
/*
	Scenario for real data
	Pelvis and femur, femur is moved
	=> iter Iteration value
*/
void pbd::scenario_1(size_t iter){
	/* if we are in first iteration */
	if(!iter){
		/* we perform out own collision detection */
		collisions c;
		/* which tells us which points collides at the start */
		find_collisions(c);
		/* we fix these points, so they wont be able to move by PBD */
		size_t len = c.size();
		/* so for each colliding point */
		for(size_t a=0; a<len; a++){
			/* we fix its position */
			m_graph->node(c[a].vec_dynamic).set_fixed(c[a].part_static);
		}
	}
	/* if we are under 200th iteration */
	if(iter<200){
		/* we set rotation angle of the femur bone */
		double angle = (iter<100)?-0.01:0.01;
		/* we need to know cosine of this angle */
		double sin = std::sin(angle);
		/* and sine too */
		double cos = std::cos(angle);
		/* we transform coordinate system to the center of femur-pelvis joint */
		matrix_4x4 trans(1,0,0,146,0,1,0,238,0,0,1,-431,0,0,0,1);
		/* we rotate femur with center of rotaion set in joint */
		matrix_4x4 rot(1,0,0,0,0,cos,-sin,0,0,sin,cos,0,0,0,0,1);
		/* and we tranform back */
		matrix_4x4 inv_trans(1,0,0,-146,0,1,0,-238,0,0,1,431,0,0,0,1);
		/* finally, we make one single matrix from those three */
		matrix_4x4 res = trans*rot*inv_trans;
		/* and we apply this tranform to the third model (with index 2) */
		m_graph->transform(2,res);
	}
}
/*
	Second scenario, this is used for "squasing" muscle between two plates
	=> iter Interation number
*/
void pbd::scenario_2(size_t iter){
	/* if we are under 200th iteration */
	if(iter<200){
		/* we declare move direction and value (if <100, we sqash the muscle,
			else we release it) */
		double move = (iter<100)?-0.9:0.9;
		/* we create transformation matrix accoring to the move */
		matrix_4x4 trans(1,0,0,0,0,1,0,move,0,0,1,0,0,0,0,1);
		/* and we apply transformation to second model (with index 1) */
		m_graph->transform(1,trans);
	}
}
/*
	Third scenario - tetrahedron preserving distance between nodes
	=> iter Iteration number from start of simulation
*/
void pbd::scenario_3([[maybe_unused]] size_t iter){
	/* we prepare force, which will be applied to the "muscle" tetrahedron */
	double force[3] = {-1e4,0,0};
	/* we apply this force to first node */
	m_graph->node(0).set_force(force);
	/* second */
	m_graph->node(1).set_fixed(true);
	/* and third node are fixed, which creates sort of pendulum */
	m_graph->node(2).set_fixed(true);
	/* in this test we are not intersed in collisions */
	m_no_collision=true;
	/* volume preservation */
	m_no_volume=true;
	/* neither holding local shape */
	m_no_bending=true;
	/* neither anisotropy */
	m_no_anisotropy=true;
}
/*
	Scenario 4 - tetrahedron which preserves local shape
	=> iter Interation number
*/
void pbd::scenario_4([[maybe_unused]] size_t iter){
	/* at first, we make sure bending is maxed */
	m_bending=1;
	/* then we prepare force */
	double force[3] = {0,0,1e4};
	/* which we apply to the fourth node of the tetrahedron */
	m_graph->node(3).set_force(force);
	/* and we fix first node, this should cause tetrahedron to expand */
	m_graph->node(0).set_fixed(true);
	/* we do not allow collisions */
	m_no_collision=true;
	/* and volume */
	m_no_volume=true;
	/* and edge length presevration */
	m_no_edge=true;
	/* neither anisotropy */
	m_no_anisotropy=true;
}
/*
	Scenario 5 - tetrahedron which preserves its volume
	=> iter Iteration number
*/
void pbd::scenario_5([[maybe_unused]] size_t iter){
	/* first of all, we prepare force */
	double force[3] = {0,0,1e4};
	/* and we apply it to the forth node */
	m_graph->node(3).set_force(force);
	/* now, we fix first */
	m_graph->node(0).set_fixed(true);
	/* and second node, so tetrahedron should converge to plane */
	m_graph->node(1).set_fixed(true);
	/* we do not allow collisions */
	m_no_collision=true;
	/* nor bending */
	m_no_bending=true;
	/* nor edge length preservation */
	m_no_edge=true;
	/* neither anisotropy */
	m_no_anisotropy=true;
}
void pbd::scenario_6(size_t iter){
	/* if we are in first iteration */
	if(!iter){
		matrix_4x4 trans(1,0,0,0,0,1,0,-5,0,0,1,0,0,0,0,1);
		m_graph->transform(0,trans);
		size_t start,end;
		m_graph->nodes_from_part(0,start,end);
		for(size_t a=start; a<end; a++){
			double* realpos = m_graph->change_position(a);
			m_graph->node(a).set_working_position(realpos);
		}
		/* we perform out own collision detection */
		collisions c;
		/* which tells us which points collides at the start */
		find_collisions(c);
		/* we fix these points, so they wont be able to move by PBD */
		size_t len = c.size();
		/* so for each colliding point */
		for(size_t a=0; a<len; a++){
			/* we fix its position */
			m_graph->node(c[a].vec_dynamic).set_fixed(c[a].part_static);
		}
	}
	/*std::vector<int> fx_1({3634,6959,3627,3654,3676,7269,3691});
	std::vector<int> fx_2({8367,8214,8213,8297,8303,8438,8159});
	for(int& i:fx_1){
		m_graph->node(i).set_fixed(1);
	}
	for(int& i:fx_2){
		m_graph->node(i).set_fixed(2);
	}*/
	if(iter<200){
		double angle = (iter<100)?-0.01:0.01;
		double sin = std::sin(angle);
		double cos = std::cos(angle);
		matrix_4x4 trans(1,0,0,146,0,1,0,238,0,0,1,-431,0,0,0,1);
		matrix_4x4 rot(1,0,0,0,0,cos,-sin,0,0,sin,cos,0,0,0,0,1);
		matrix_4x4 inv_trans(1,0,0,-146,0,1,0,-238,0,0,1,431,0,0,0,1);
		matrix_4x4 res = trans*rot*inv_trans;
		m_graph->transform(2,res);
	}
}
