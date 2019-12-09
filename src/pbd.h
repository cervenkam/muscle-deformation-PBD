/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################
*/
#pragma once
#include "point_3D.h"
#include "graph.h"
/* Following only if we are interested in time measurement */
#ifdef TIME_MEASUREMENT
#include <iostream>
#include <iomanip>
/* This macro capture acutal time */
#define CAPTURE_TIME(X) auto (X) = std::chrono::high_resolution_clock::now()
/* This macro prints time difference between two timestamps */
#define PRINT_DIFF(X,TMP,D1,D2) \
	auto (TMP) = std::chrono::duration_cast<std::chrono::nanoseconds>( \
		(D1)-(D2)).count(); \
	std::cout << (X) << std::right << std::setw(10) << (TMP) << \
	" (limit:  " << std::right << std::setw(10) << static_cast<size_t>(1e9/ \
	static_cast<double>(TMP)) << " FPS )" << std::endl
#else
/* This macro does not anything (we do not want time measurement) */
#define CAPTURE_TIME(X)
/* This macro does not anything (we do not want time measurement) */
#define PRINT_DIFF(X,TMP,D1,D2)
#endif
/*
	Struct which holds information about a single collision
*/
struct collision{
	/* Muscle model ID */
	size_t part_dynamic;
	/* Bone model ID */
	size_t part_static;
	/* Muscle node ID */
	size_t vec_dynamic;
	/* Muscle node best X coordinate, which avoids collision */
	double x;
	/* Muscle node best Y coordinate, which avoids collision */
	double y;
	/* Muscle node best Z coordinate, which avoids collision */
	double z;
};
/* We make notation for collision shorter */
typedef std::vector<struct collision> collisions;
/*
	This class calculates PBD deformation. It uses data
	provided by "graph" class.
*/
class pbd{
	public:
		/*
			This method does one "outer" iteration of the PBD algorithm
			=> iter Iteration number
			<= Return state (0 - OK, iteartion has been done, else FAIL)
		*/
		int execute(size_t);
		/*
			This class provides elementary constraint solving
		*/
		class constraints{
			public:
				/*
					Project edge length preserving constraint on a pair of nodes
					=> source Source node ID (from edge point of wiew)
					=> target Destination node ID (from edge point of view)
					=> len Original length of this edge
					=> stiffness Stiffness of the edge
				*/
				static void project_distance_constraint(
					point_3D&,
					point_3D&,
					double,
					double,
					std::vector<double>&
				);
				/*
					This method compute volume constraint for single point
					=> pnt Point, which target position should be fixed
					=> mult Volume of the tetrahedron
				*/
				static void project_volume_constraint(
					point_3D& pnt,
					double volume
				);
				/*
					This method projects bending constraint to a single node
					=> graph Structured data about all models
					=> phi_0 Original dihedral angle between two faces
					=> n_th_neightbour ID of the node, which is different on
						first triangle considering second triangle
					=> second_side_neighbourID of the node, which is different
						on second triangle considering first triangle
					=> part Muscle model ID
					=> first First polygon ID
					=> second Second polygon ID
					=> bending Bending factor
				*/
				static void project_bending_constraint(
					graph*,
					double&,
					size_t&,
					size_t&,
					size_t&,
					size_t&,
					size_t&,
					double&,
					std::vector<double>&
				);
				/*
					Method projects collision constraint on single node
					=> position Target position of the node
					=> col Information about collision which occurs
				*/
				static void project_collision_constraint(
					double*,
					struct collision&
				);
		};
		/*
			D'tor
		*/
		~pbd();
		/*
			Method inits PBD algorithm
			=> grp Graph containing structured data for PBD algorith,
			=> data Raw data
			=> scenario ID of the dynamic scenario (it has to be implemented
				this way, we can force user to write his own scenario and
				compile it, but who would do that. We also can create a
				scripting language, but who would like to learn in.
		*/
		void init(graph* grp,data* dta,size_t);
		/*
			Returns structured data (graph) associated with this PBD algorithm
			=> Graph/structured data
		*/
		graph* get_graph();
	private:
		/*
			This method resets color of all models
			Muscles will be red, bones light gray (nearly white)
		*/
		void reset_colors();
		/*
			Finds collision for each node in a signle muscle model
			<=> coll Information about collision
			=> muscle_part Muscle model ID
			=> start First bone model ID
			=> end Last bone model ID
		*/
		void find_collision_constraints(collisions&,size_t,size_t,size_t);
		/*
			Calculates collision between single muscle and single bone
			<=> coll Information about collisions occured
			=> muscle_part Muscle model ID
			=> bone_part Bone model ID
		*/
		void find_collision_constraints_one_bone(collisions&,size_t,size_t);
		/*
			Method updates velocity with external force applied
			to a single point
			=> vertex Vertex on which the force is applied
		*/
		void update_velocity_with_external_forces(point_3D& vertex);
		/*
			Performs velocity loss
			=> vertex Node, which loses a bit of velocity
		*/
		void damp_velocity(point_3D& vertex);
		/*
			Method calculates new destination position according to its velocity
			=> x Current position
			=> vertex Moving node with current and destination position
		*/
		void calculate_new_position(const double*,point_3D& vertex);
		/*
			Method project vertex constraints (edge lenght preservation) to
			each edge outgoing from specified node
			=> vertex Point, which has outgoing edges
			=> index Muscle model ID
		*/
		void project_vertex_constraints(point_3D& vertex,size_t,std::vector<double>&);
		/*
			Method updates velocity according to origin and target position
			=> x Origin position
			=> vertex Moving node
		*/
		void update_velocity(const double*,point_3D& vertex);
		/*
			Method updates target position of the node
			=> pos New target position
			=> vertex Node to change its target position
		*/
		void update_position(double*,point_3D& vertex);
		/*
			Calculate volume constraint of single triangle
			=> idx Triangle ID
			=> part Model ID
			<=> numerator Volume numerator
		*/
		void calculate_volume_constraint(size_t,size_t,double&);
		/*
			This method is called during first iteration
			It is supposed to precalculate as much as possible
			This is the reason first iteration take so much time
		*/
		void initialize_first_iteration();
		/*
			This method moves each point according to external force
		*/
		void move_points_by_force();
		/*
			This method finds all colisions between bones and muscles
			<=> Information about collisions
		*/
		void find_collisions(collisions& coll);
		/*
			Method projects (corrects) edge length constraints
			(preserves distance)
		*/
		void project_vertex_constraints();
		/*
			Method projects (corrects) collisions
			=> coll Information about collisions
		*/
		void project_collision_constraints(collisions& coll);
		/*
			Method calculates volume constraints (calculates current volume)
				=> numerators Memory, where volume numerators can be stored
		*/
		void calculate_volume_constraints(double* numerator);
		/*
			This method project volume constraints to each muscle node
			=> numerator Volume numerator for each model
		*/
		void project_volume_constraints(double* numerator);
		/*
			Projects bending contraint (preserving local shape) to each node
			of each muscle model
		*/
		void project_bending_constraints();
		/*
			This method update velocities of each node in each muscle model
		*/
		void update_velocities();
		/* SCENARIOS */
		/*
			Default, empty scenario, does nothing usefull
			=> iter Iteration value
		*/
		void scenario_0(size_t);
		/*
			Scenario for real data
			Pelvis and femur, femur is moved
			=> iter Iteration value
		*/
		void scenario_1(size_t);
		/*
			Second scenario, this is used for "squasing" muscle between
			two plates
			=> iter Interation number
		*/
		void scenario_2(size_t);
		/*
			Third scenario - tetrahedron preserving distance between nodes
			=> iter Iteration number from start of simulation
		*/
		void scenario_3(size_t);
		/*
			Scenario 4 - tetrahedron which preserves local shape
			=> iter Interation number
		*/
		void scenario_4(size_t);
		/*
			Scenario 5 - tetrahedron which preserves its volume
			=> iter Iteration number
		*/
		void scenario_5(size_t);
		void scenario_6(size_t);
	public:
		/* time difference between two "outer" interations */
		double m_time_diff = 1e-3;
		/* velocity damp <0;1> */
		double m_damp = 0.99;
		/* bending factor <0;1> */
		double m_bending = 0.9;
		/* pressure factor (0;+\infty), rather (0.5;2) */
		double m_pressure = 1;
		/* number of "internal" interation */
		double m_constr = 3;
		/* maximum velocity allowed (protection agains "explosion") */
		double m_max_velocity = 1;
		/* gravity factor */
		double m_gravity = 0;
		/* pause <0;1> (binary, threshold 0.5) */
		double m_pause = 1;
		/* anisotropy factor <0;+\infty), rather <0;10> */
		double m_anisotropy = 0.5;
		/* anisotropy used in previous iteration */
		double m_last_anisotropy = m_anisotropy;
		/* is edge length preservation disabled */
		bool m_no_edge=false;
		/* should we not care about local shape */
		bool m_no_bending=false;
		/* and what about anisotropy */
		bool m_no_anisotropy=false;
		/* and volume preservation */
		bool m_no_volume=false;
		/* and finally collisions */
		bool m_no_collision=false;
		/* dynamic scenario used (zero "empty" by default) */
		void (pbd::*m_scenario)(size_t) = &pbd::scenario_0;
	private:
		const bool m_measurement_bending = false;
		const bool m_measurement_volume = false;
		const bool m_measurement_vertex = false;
		/* number of iteration, after which simulation stops */
		double m_iter = std::numeric_limits<double>::max();
		/* Model ID of the first muscle */
		size_t m_muscle_start;
		/* Model ID of the last muscle +1 */
		size_t m_muscle_end;
		/* Model ID of the first bone */
		size_t m_bone_start;
		/* Model ID of the last bone +1 */
		size_t m_bone_end;
		/* First muscle node ID */
		size_t m_start;
		/* Last muscle node ID +1 */
		size_t m_end;
		/* Structured data about models */
		graph* m_graph;
		/* Raw loaded data */
		data* m_data;
		size_t subiter;
};
