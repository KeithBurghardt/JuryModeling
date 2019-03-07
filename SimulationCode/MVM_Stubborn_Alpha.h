
/*
******************************************************
Competing Infection Header for CompetingInfections_Poisson_Temporal.c 

Created By: Keith Burghardt
******************************************************
*/
#ifndef _CompetingInfections_Poisson_h
#define _CompetingInfections_Poisson_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort_int.h>

	
////////////////////////////////////////////////////////////////////////////////

typedef struct
{
    // Node properties
    int                 	ID;
    int                 	degree;
    int                 	current_number_of_neighbors;
    int				tri_vertex_edges;

    // Dynamics on nodes
    double			time_infected;//time with infection
    double			time_last_checked;//time with infection
    double			current_stubbornness_p;//time with infection
    int 			num_active_edges;
    int 			num_contrary_edges;
    double			frequency_exposed;//number of times exposed to infection
    double			s; //re-scaling of infection rate
    int				initial_preference;//prefer up state/downstate

    // Node arrays
    int * __restrict__		neighbors;
    int * __restrict__ 		wait_time;//wait times
    int * __restrict__ 		wait_tau;
 
    // Dynamics on nodes arrays
    unsigned int * __restrict__	infection; //infection value: 1, 2, 3, ...
    int * __restrict__		active_neighbors;
    int * __restrict__		contrary_neighbors;
    double * __restrict__	Bj;
} node;

typedef struct
{
    // Graph Properties
    int                 	NumNodes;
    int                 	NumEdges;
    double              	average_degree;
    double              	assort_p;
    double              	alpha;

    //Dynamic properties on the graph
    int				Q;
    int				Ni;
    int				NumActiveEdges;
    int				SumBij;
    int				NumGuilty;

    //Graph arrays
    node * __restrict__ 	nodes;
    int * __restrict__  	degree_sum;
    int * __restrict__		Regions;
    double * __restrict__	Rho_i;
    double * __restrict__	v0;

    //Infection algorithm arrays
    int * __restrict__  	InfectiousNodes;
    int * __restrict__  	CumSumActiveEdges;
    double * __restrict__  	CumSumBi;

    // Properties of infection
    long double 		m;
    long double 		m2;
    long double 		m4;
    long double * __restrict__  MeanFractionInfected;
    double 			ConsensusTime;
    int * __restrict__   	FreqExposureUntilInfectedHist;
    int * __restrict__   	TimeUntilInfectedHist;

} graph;

typedef struct
{
    // Thread properties
    int                 	thread;
    gsl_rng *           	r;
    graph               	G;

    // Record of dynamics within thread
    int				Q;

    double * __restrict__   	FractionInfectedVersusTime;
    double * __restrict__	UVersusTime;
    double * __restrict__  	Mean;
    double * __restrict__  	m;
    double * __restrict__  	m2;
    double * __restrict__  	m4;
    double * __restrict__	Correlation;    
    double * __restrict__	Rho_i;
    double * __restrict__	ConsensusTimes;
    double * __restrict__	FreezeTime;
    double * __restrict__	DynamicsP;
    double * __restrict__ 	v0;
    
    int * __restrict__   	T_distribution;
    int * __restrict__   	TOTAL_TimeUntilInfectedHist;
    int * __restrict__   	TOTAL_FreqExposureUntilInfectedHist;
    int * __restrict__  	whos_infected;
    int * __restrict__  	NumsGuilty;

    /*// Files to make graphs
    FILE * 			s,
    FILE * 			f*/
} thread_args;

void free_thread_args (	thread_args * 		args);

void free_graph (       graph * __restrict__ 	G);

void free_node (        node * __restrict__   	N);

static __inline__ void make_node (
		        node * __restrict__   	n,
                        int         	  	ID,
                        int         	  	deg,
			double		  	fract_tri_vertices,
			int			Q);

void make_graph (	graph * __restrict__   	G,
                        int                	N,
                        double            	alpha,
                        double            	k_av,
                        double             	assort_p,
                        double             	fract_tri_vertices,
			int			Q,
                 	FILE **                 s_deg,
                 	FILE **                 s_graph,
                        gsl_rng *          	r);

static __inline__ void make_file_nodes (  
			graph * __restrict__    G,
                        FILE **                 s_deg,
                        FILE **                 s_graph);

static __inline__ void make_grid(         
			graph * __restrict__    G,
                        int                     N,
			double			assort_p);

static __inline__ void make_SFdegree_nodes ( 
			graph * __restrict__ 	G,
                        int                  	N,
                        double               	alpha,
                        int                  	k_max,
                        double               	k_av,
                        double               	assort_p,
                        double               	fract_tri_vertices,
                        const gsl_rng *		r);

static __inline__ void make_Kregular_nodes (
			  graph * __restrict__ 	G,
                          int 			N,
                          double 		av_deg,
                          double 		fract_tri_vertices);

static __inline__ void make_complete_nodes (
			  graph * __restrict__ 	G,
                          int 			N);


static __inline__ void make_poisson_degree_nodes (
			graph * __restrict__ 	G,
			int	   		N,
			double	   		av_deg,
			double	   		fract_tri_vertices,
			double	   		assort_p,
			gsl_rng *  		r);

static __inline__ void connect_nodes (    
			graph * __restrict__ 	G,
                        const gsl_rng * 	r);

static __inline__ void spatial_connect_nodes (    
			graph * __restrict__ 	G,
                        const gsl_rng * 	r);

static __inline__ void connect_node (     
			graph * __restrict__    G,
	                int                     x,
                        int                     y,
			const gsl_rng * 	r);

static __inline__ double dist (
			int 			i, 
			int 			j,
			int			L);

static __inline__ void find_random_edges (
			graph * __restrict__ 	G,
                        int 			random_nodes [],
                        const gsl_rng * 	r);

static __inline__ void pick_edge (	graph * __restrict__    G,
                	int                     edge_array [],
	        	int                     index,
	        	const gsl_rng *        	r);
/*
void pick_infected_edge (
             		graph * __restrict__    G,
                	int * __restrict__      edge_array,
			double			tau,
			double			beta,
			double			t,
                	const gsl_rng *         r);
*/
static __inline__ int node_search_int(
			int * __restrict__ 	degree_sum,
			int 			random_edge,	
			int 			N);

static __inline__ void sort (		
			graph * __restrict__ 	G, 
			int * 			node_array);

static __inline__ void rewire (   	graph * __restrict__ 	G,
                	int 			random_nodes[],
                	int		 	ordered_nodes[],
                	double 			p,
                	const gsl_rng * 	r);

static __inline__ void rewire_links (	
			graph * __restrict__ 	G, 
			int * __restrict__	new_neighbors, 
			int * __restrict__	old_neighbors);

static __inline__ int already_neighbors (	graph * __restrict__ 	G, 
			int * __restrict__	new_neighbors);

static __inline__ void replace_links(	graph * __restrict__ 	G, 
			int * __restrict__	old_neighbors, 
			int * __restrict__	new_neighbors);

static __inline__ int     find_neighbor_from_list (
			int 			n1, 
			int 			neighbor_list []);

static __inline__ int     find_neighbor_pos (
			int 			neighbor, 
			node * __restrict__ 	n);

static __inline__ int     find_neighbor(  graph * __restrict__ 	G,
                        int             	node,
                        int             	j,
                        int * __restrict__ 	multi_edge,
                        const gsl_rng * 	r);

void find_tri_neighbors(     
			graph * __restrict__ 	G,
                        int             	node,
                        int             	j,
                        int * __restrict__  	multi_edge1,
                        int * __restrict__  	multi_edge2,
                        int * __restrict__  	neighbor1,
                        int * __restrict__  	neighbor2,
                        const gsl_rng * 	r);

static __inline__ void find_distribution( 
			double             	alpha,
                        int                	N,
                        double             	k_av,
                        int * __restrict__ 	k_min,
                        int               	k_max,
                        int * __restrict__ 	n);

static __inline__ void find_distribution_fast( 
			double             alpha,
                        double             	k_av,
                        int * __restrict__ 	k_min,
			double * __restrict__  	P_k_min_minus_1);

static __inline__ double av_deg (
		        double             	alpha,
                        int                	k_min,
                        int                	k_max);

static __inline__ double av_deg_fast ( 	
			double  	   	alpha,
 	       	        double     	   	k_min);

static __inline__ double finite_size_poly_sum(
			double         		alpha,
                        int            		k_min,
                        int            		k_max);

static __inline__ int ran_scale_free (
			const gsl_rng *         r,
                        int                     x_min,
                        double                  alpha,
                        int			N);

static __inline__ double ran_ref_1D_walk (
                        double                  T,
                        const gsl_rng *         r);

void make_thread_args ( thread_args *	   	args,
                        int                	thread,
                        int                	N,
			FILE **			s_deg,
			FILE **			s_graph,
                        double             	alpha,
                        double             	k_av,
                        double             	assort_p,
                        double             	fract_tri_vertices,
			int			Q);



void link_dynamics(     graph * __restrict__    G,
                        int                     node,
                        double                  average_wait_tau,
                        double                  wait_tau_alpha,
                        gsl_rng *         	r);

static __inline__ void resample_wait_tau (
			graph * __restrict__    G,
                        int			node,
                        int			neighbor,
                        double			tau_av,
                        double			alpha,
                        gsl_rng *		r);

int rand_SFtau (	double 		   	tau_av, 
			double 		   	alpha, 
			int			N, 
			gsl_rng * 	   	r);

static __inline__ void synchronize_link ( 
			graph * __restrict__    G,
                        int                     node,
                        int                     neighbor);

//infection function definitions

static __inline__ void count_infections ( 
			graph * __restrict__    G,
                        int * __restrict__      NUMBER_INFECTED,
                        double               	average_wait_tau,
                        double                  wait_tau_alpha,
			gsl_rng * 		r);
/*
static __inline__ void find_2nd_mag_moment(
                        graph * __restrict__    G,
                        double * __restrict__   X);

static __inline__ void find_4th_mag_moment(
                        graph * __restrict__            G,
                 	long double * __restrict__	U);
*/

static __inline__ void copy_infection(    
			graph * __restrict__    G,
                        int                     neighbor,
                        int                     n,
			int			is_infected,
			int			Q,
			double			t);

static __inline__ void infect(
	  		node * __restrict__   	n, 
			int 			strain,
			double			t);

static __inline__ int is_susceptible ( 	
			node * __restrict__    	n);

static __inline__ int distinct_infection (
			node * __restrict__    	n, 
			node * __restrict__    	neighbor);

static __inline__ int find_strain (	
			node * __restrict__ 	n,
			int			Q);

static __inline__ unsigned int BIT_POS (
			int 			var, 
			int 			pos);

static __inline__ void randomize (	
			int * __restrict__ 	list, 
			int 			N, 
			const gsl_rng * 	r);

static __inline__ void infect_and_recover (
                graph * __restrict__    G,
                int * __restrict__      NUMBER_INFECTED,
                double                  beta,
                double                  gamma,
                double                  temp,
                double                  mvm_p,
                double                  tau,
                double                  tau_hung,
                double                  stop_rate,
                double                  mu,
                double                  average_wait_tau,
                double                  wait_tau_alpha,
                double                  hung_p_innocent,
                double                  hung_p_guilty,
                double                  S,
                gsl_rng *               r);


/*static __inline__ void add_time_infected (
			graph * __restrict__    G,
                        double                  dt);*/

static __inline__ void infect_and_recover_neutral
                    (   graph * __restrict__	G,
			int * __restrict__	NUMBER_INFECTED,
                        double                  beta,
                        double                  gamma,
                        double                  temp,
			double			tau,
     	                gsl_rng *               r);

static __inline__ void infect_and_recover_inout (
                graph * __restrict__    G,
                int * __restrict__      NUMBER_INFECTED,
                double                  beta,
                double                  gamma,
                double                  temp,
                double                  mvm_p,
                double                  tau,
                double                  tau_hung,
                double                  stop_rate,
                double                  mu,
                double                  average_wait_tau,
                double                  wait_tau_alpha,
                double                  hung_p_innocent,
                double                  hung_p_guilty,
                double                  S,
                gsl_rng *               r);

static __inline__ void infect_outward (   
			graph * __restrict__ 	G,
			int * __restrict__	NUMBER_INFECTED,
                        double                  gamma,
                        double                  beta,
                        double                  tau,
                        double                  mu,
                        double                  t,
			gsl_rng * 		r);

static __inline__ void infect_inward (    graph *                 G,
                        int * __restrict__      NUMBER_INFECTED,
                        double                  beta,
                        double                  mvm_p,
                        double                  tau,
                        double                  tau_hung,
                        double                  stop_rate,
                        double                  mu,
                        double                  t,
                        int	                total_vote,
                        gsl_rng *               r);


static __inline__ void infect_susceptible (
                        graph * __restrict__    G,
			int * __restrict__	NUMBER_INFECTED,
                        int                     n,
                        int              	neighbor,
                        double                  beta,
                        double                  t,
                        gsl_rng *               r);

static __inline__ void become_turncoat(   
			graph * __restrict__    G,
			int * __restrict__	NUMBER_INFECTED,
                        int                     n,
                        int                     neighbor,
                        double                  beta,
                        double                  tau,
                        double                  t,
                        gsl_rng *               r);

static __inline__ void find_turncoat_p (  
			graph * __restrict__    G,
                        double * __restrict__   turncoat_p,
                        int                     n,
                        int                     neighbor,
                        double                  beta,
                        double                  tau,
                        double                  delta_t);

static __inline__ void recover (
		        graph * __restrict__    G,
			int * __restrict__	NUMBER_INFECTED,
                        double                  tau,
                        double                  t,
                        gsl_rng *               r);
/*void recover (
			graph * __restrict__	G,
			gsl_rng *		r);*/

static __inline__ void pick_infectious_node (
                        graph * __restrict__    G,
                        int * __restrict__	node,
                        gsl_rng *               r);

static __inline__ void infect_active_edge (
                        graph * __restrict__    G,
			int * __restrict__	NUMBER_INFECTED,
                        double                  tau,
			double			t,
                        gsl_rng *               r);

static __inline__ void 	pick_active_edge (
                    	graph * __restrict__    G,
                        int * __restrict__	infectious_node,
                        int * __restrict__	susceptible_node,
                        gsl_rng *               r);

static __inline__ void infect_contrary_edge (
                        graph * __restrict__    G,
			int * __restrict__	NUMBER_INFECTED,
                        double                  tau,
			double			t,
                        gsl_rng *               r);

static __inline__ void 	pick_contrary_edge (
                        graph * __restrict__    G,
                        int * __restrict__	contrary_node1,
                        int * __restrict__	contrary_node2,
                        gsl_rng *               r);

static __inline__ void calc_mean_fraction_infected (
                        graph * __restrict__    G,
                        double * __restrict__   FRACTION_INFECTED,
                        double                  dt);

static __inline__ void normalize_mean_fraction_infected (
                        graph * __restrict__    G,
                        double                  delta_t);

static __inline__ int one_infection_survives (
                        int * __restrict__	NUMBER_INFECTED,
                        int                     Q,
                        int                     N,
                        int * __restrict__	surviving_infection);

static __inline__ void update_infection_states_recover (
                        graph * __restrict__    G,
			int * __restrict__	NUMBER_INFECTED,
                        int                     node,
                        double                  tau,
                        double                  t);

static __inline__ void update_infection_states_infect (
                        graph * __restrict__    G,
			int * __restrict__	NUMBER_INFECTED,
                        int                     node,
                        int                     strain,
                        double                  tau,
                        double			t);

static __inline__ void update_infection_states_contrary (
                        graph * __restrict__    G,
			int * __restrict__	NUMBER_INFECTED,
                        int                     node,
                        int                     strain,
                        double                  tau,
                        double                  t);

static __inline__ void add_links (        graph * __restrict__    G,
                        int                     node,
                        int                     neighbor,
                        int                     strain);

static __inline__ void add_active_links ( 
			graph * __restrict__    G,
                        int                     node,
                        int                     neighbor);

static __inline__ void add_active_links_to_neighbors (
                        graph * __restrict__    G,
                        int                     node);

static __inline__ void add_contrary_links (
                        graph * __restrict__    G,
                        int                     node,
                        int                     neighbor,
                        int                     strain);

static __inline__ void remove_contrary_links (
                        graph * __restrict__    G,
                        int                     node);

static __inline__ void remove_active_links (
                        graph * __restrict__    G,
                        int                     node);

static __inline__ void remove_list_int (  
			int * __restrict__	*list,
             	        int * __restrict__	list_size);

static __inline__ void remove_list_double (       
			double * __restrict__	*list,
                        int * __restrict__	list_size);

static __inline__ void add_element_int (  
			int * __restrict__	* list,
                        int * __restrict__	list_size,
                 	int                     val);

static __inline__ void add_element_double (
                        double * __restrict__	*list,
                        int * __restrict__	list_size,
                        double                  val);

static __inline__ void remove_element_int (
                	int * __restrict__	*list,
                	int * __restrict__	list_size,
                	int                     value,
                	int *      		pos);

static __inline__ void remove_element_double (
                	double * __restrict__	*list,
                	int * __restrict__	list_size,
                	int			pos);

static __inline__ int find_pos_int (
                	int * __restrict__	list,
                	int                     list_size,
                	int                     value);

static __inline__ int find_pos_double (
                	double * __restrict__	list,
                	int	                list_size,
                	double                  value);

static __inline__ void CalcCumSumBi(
			graph * __restrict__ 	G,
			double			tau,
			double			t);

static __inline__ void CalcCumSumActiveEdges(graph * __restrict__ G);

static __inline__ void mutate (
		   	node * __restrict__	node,
                	double                  mu,
			int			Q,
			double			t,
                	gsl_rng *               r);

static __inline__ int emp_rnd (           
			double * __restrict__	cdf,
                        double                  InfNum,
	                gsl_rng *               r);

static __inline__ void hist2cdf (         
			int * __restrict__	hist,
                        int                     length,
                        double * __restrict__   cdf);

static __inline__ void frac_of_total(     
			double * __restrict__   FRACTION_INFECTED,
                        int * __restrict__      NUMBER_INFECTED,
                        int                     N,
			int			Q);

//static __inline__ 
void setup_experiment(  
			graph * __restrict__    G,
                      	double                  average_wait_tau,
                      	double                  wait_tau_alpha,
			double			tau,
			double			beta,
            	 	double             	S,
            	 	double             	delta_inf,
            	 	double             	binom,
                      	int * __restrict__      NUMBER_INFECTED,
                      	double * __restrict__   FRACTION_INFECTED,
                      	int * __restrict__      T_distribution,
                      	gsl_rng *               r);

static __inline__ void seed_infection (   
			graph * __restrict__    G,
                        int * __restrict__	NUMBER_INFECTED,
                        double                  tau,
                        double			beta,
            	 	double             	S,
            	 	double             	delta_inf,
            	 	double             	binom,
                        gsl_rng *         	r);

static __inline__ void initializing_state (    
				graph * __restrict__    G,
                                int * __restrict__	NUMBER_INFECTED);

static __inline__ void init_time_infected_dist (
			int * __restrict__ 	T_distribution,
			int			Q);

static __inline__ void init_link_wait_time (	
			graph * __restrict__ 	G, 
			double 			average_wait_tau, 
			double			wait_tau_alpha, 
			gsl_rng * 		r);

static __inline__ void make_susceptible (
			node * __restrict__ 	n,
			int			Q);

static __inline__ void add_stubbornness(  
			graph * __restrict__    G,
			double			S,
			double			beta,
                        gsl_rng *               r);

static __inline__ int node_already_picked (
                        int                     i_ran,
                        int * __restrict__	nodes,
                        int                     i);

static __inline__ void flip_spin (        
			graph * __restrict__    G,
			int * __restrict__	NUMBER_INFECTED,
			double			tau,
			double			t,
                        gsl_rng *               r);	

static __inline__ void make_flip  (
			node * __restrict__     n,
			int * __restrict__	NUMBER_INFECTED,
                        int                     Q,
			double			t,
                        gsl_rng *               r);

static __inline__ void swap(
			int *	__restrict__ 	list, 
			int			length,
			int	 		pos, 
			int 			new_val);

static __inline__ int find(
			int * __restrict__	list,
			int			length,
			int			value);

static __inline__ void create_random_strain_list( 
			int * __restrict__ 	random_strain_list, 
			int			InfNum,
			double			delta_inf,
			double			binom,
			int			Q,
			gsl_rng * 		r);

static __inline__ void store_number_of_infected_nodes (   
			double * __restrict__   FractionInfectedVersusTime,
			double * __restrict__   UVersusTime,
                        double * __restrict__   FRACTION_INFECTED,
                        graph *  __restrict__ 	G,
                        int *    __restrict__	T_distribution,
                        double                  tau,
                        double                  t);

static __inline__ double sum_double (
			double * __restrict__	y,
			int			n);

static __inline__ double sum_int (
			int * __restrict__	y,
			double			n);

void data_analysis ( 	thread_args * 		t, 
			int * 			count);

void determine_regional_rho ( 
			graph * __restrict__ 	G);

static __inline__ double determine_regional_distance (
			int			i,
			int			j);

static __inline__ void cross_correlation_r (      
			graph * __restrict__    G,
			int                     r,
			double * __restrict__   cross_correlation,
			int * __restrict__	num_regions_r);

static __inline__ double spatial_correlation (
                        graph * __restrict__    G,
                        int                     r,
                        double                  cross_correlation,
                        int                     num_regions_r);

void run_experiment (
                     graph * __restrict__	G,
                     double * __restrict__	FractionInfectedVersusTime,
                     double * __restrict__	UVersusTime,
                     long long                  start,
                     int * __restrict__         T_distribution,
                     double                     beta,
                     double                     S,
                     double                     delta_inf,
                     double                     gamma,
                     double                     temp,
                     double                     mu,
                     double                     mvm_p,
                     double                     tau,
                     double                     tau_hung,
                     double                     stop_rate,
                     double                     stop_rate_hung,
                     double                     average_wait_tau,
                     double                     wait_tau_alpha,
                     double                     hung_p_innocent,
                     double                     hung_p_guilty,
                     gsl_rng *                  r);
/*
static __inline__ void run_experiment (
                     graph * __restrict__       G,
                     double * __restrict__      FractionInfectedVersusTime,
                     double * __restrict__      UVersusTime,
                     long long                  start,
                     int * __restrict__         T_distribution,
                     double                     beta,
                     double                     S,
                     double                     delta_inf,
                     double                     gamma,
                     double                     temp,
                     double                     mu,
                     double                     mvm_p,
                     double                     tau,
                     double                     tau_hung,
                     double                     stop,
                     double                     stop_hung,
                     double                     average_wait_tau,
                     double                     wait_tau_alpha,
                     double                     hung_p_innocent,
                     double                     hung_p_guilty,
                     gsl_rng *                  r);
*/

void* perform_trials (	void *      	    	trial);

void write_graph_to_file (
			thread_args * 		trial, 
			FILE * 			p);

void write_simulation_data_to_file(
			thread_args * 		trial, 
			int 			i, 
			int 			infection, 
			double			tau, 
			FILE * 			p);

void write_data_to_file(pthread_t * 		threads, 
			thread_args * 		trial);

void create_threads (	pthread_t * 		threads, 
			thread_args * 		trial);

void behavior_space ();

void time_and_run_behavior_space ();

void seed_rand(         int         	    thread_n,
                        gsl_rng *   	    r);

#endif
