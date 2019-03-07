//Written By: Keith Burghardt
/*
 This is code models juries on an a large variety of networks. 
 - parameter COMPLETE forces network to be all-to-all and ignores 
   other network parameters like degree.
 - Mu in our paper is 1/Tau in the simulation. Tau represents 
   the timescale of stubbornness.

 NOTES:
 - compile with "$ gcc -Wall -O3 -pthread
   -Wl,-rpath=/path/to/gsl.1.16/lib
   -L/path/to/gsl.1.16/lib -lgsl -lgslcblas -lm
   -I/path/to/gsl.1.16/include -funroll-loops
   [maybe -fvariable-expansion-in-unroller] -o Simulation MVM_Stubborn_Alpha.c

 - allocate memory in parallel with Hoard, or similar
 - run ./Simulation
 - Number of Threads = NUMBER_OF_TRIALS * (number of degree distributions) * (number of clustering coefficients)

  output:
 "name.dat": " AVERAGE_DEGREE_MIN = , AVERAGE_DEGREE_MAX, AVERAGE_DEGREE_STEP_SIZE= ,
 gamma= ... , beta= ... ,..., initial number infected = ... # Trials =__  \n"
 
 then array: "Average Fraction Infected \t Variance \t Coefficient for 1/f^alpha noise \n"
 for all trials, then all beta's, then all gamma's, and then all degrees.
 - then we output time and final vote

 - Behavior Space determined by definitions below.
 
 

 */
////////////////////////////////////////////////////////////////////////////////

#include "MVM_Stubborn_Alpha.h"

/********************************************************/
////////////////    Code Parameters     //////////////////
/********************************************************/

// explicitly state all steps
#define VERBOSE			(0)

// record the number of times exposed before node changes strain
#define RECORD_EXPOSURE		(0)

//record the length of time infected OVER TIME
#define RECORD_TIME_INFECTED	(0)

//record the time to reach equilibrium
//IF THIS IS 1 "RECORD_LONG_CONSENSUS_TIME" MUST BE 1 AS WELL
#define RECORD_EQUILIBRIUM_TIME	(0)

//record the time to reach consensus for the minority
#define RECORD_CONSENSUS_TIME	(0)

//Same a above, but we increase speed and decrease memory 
//but the only data that can be recorded is consensus times and final vote
#define RECORD_LONG_CONSENSUS_TIME (1)

// Record the number in majority when simulation ends
//USE THIS WHEN RECORDING VOTE AND TIME TOGETHER
#define RECORD_DIST 		   (1)

// USE THIS WHEN RECORDING VOTES OVER TIME FOR ALL SIMULATIONS
#define RECORD_VOTE_TIME 	   (0)

// USE THIS TO REMOVE HUNG CONDITIONS
#define REMOVE_HUNG_CONDITIONS	   (0)

// USE THIS TO REMOVE STOP HUNG CONDITION
#define REMOVE_STOP_HUNG_CONDITION	   (0)

// USE THIS TO REMOVE STUBBORNNESS HUNG CONDITION
#define REMOVE_STUBBORNNESS_HUNG_CONDITION (1)

// USE THIS TO REMOVE Q TIME DEPENDENCE (DEPENDENCE ON CURRENT VOTE)
#define REMOVE_Q_TIME_DEPENDENCE (0)

// USE THIS TO REMOVE STUBBORN TIME DEPENDENCE (DEPENDENCE ON TIME WITH LATEST OPINION)
#define REMOVE_STUBBORN_TIME_DEPENDENCE (0)

// USE THIS WHEN RECORDING THE TIME FROZEN AND FINAL P (versus V_g)
#define RECORD_FREEZE_P         (0)

// Threshold between hanging and not:
//CA, OR (civil): 0.75
//OR (criminal): 0.833333 (~0.82)
//WA/NE: 1.0 (~0.99)
#define VOTE_THRESHOLD			(0.82)

//find v/v_0 = v Q / N given Q and N from real data
#define REAL_DATA_DISTRIBUTION	(0)
#define INFECT_CANDIDATES	(0)
#define NUMBER_OF_ELECTIONS	(1)


//record what nodes are infected on a graph
#define RECORD_INFECTIONS_ON_GRAPH (0)

//record the graph itself
#define RECORD_GRAPH 		(0)

//record the difference in the infections when all nodes are first infected
#define FIND_INTIAL_INFECTION_DIFFERENCE	(0)

//sampling small portion of network
#define SAMPLE			(0)
#define SAMPLE_FRACTION		(1.0)

// a small value
#define eps             	(0.00000001)

//way in which infection takes place
#define INFECT_INWARD		(1)
#define INFECT_NEUTRAL		(0)
#define INFECT_OUTWARD		(0)

// definitions
#define INFECTED                (1)
#define SUSCEPTIBLE             (0)

/////////////behavior space parameters////////////////////


/********************************************************/
/////////////     Network Parameters     /////////////////
/********************************************************/


// Order of the network (N)
  
#define MIN_NUMBER_OF_NODES     (12)
#define MAX_NUMBER_OF_NODES     (12)
#define STEP_SIZE_NUMBER_OF_NODES (4.0)// Must be greater than 1


// Average degree (<k>) for a random network

#define AVERAGE_DEGREE_MIN	 (10.0)
#define AVERAGE_DEGREE_MAX	 (10.0)
#define AVERAGE_DEGREE_STEP_SIZE (2.0)


///////////////////////////////////////////////////
///////////     DEGREE PARAMETERS     /////////////
///////////////////////////////////////////////////

//defines whether graph is K-regular
#define K_REGULAR               (0)

//K-regular annealed graph (no constant links)
#define K_REGULAR_ANNEALED	(0)

//alpha coefficient for SF networks

#define ALPHA_MIN               (0.0)
#define ALPHA_MAX               (0.0)
#define ALPHA_STEP_SIZE         (0.4)

//fraction of nodes maximum degree can cover
//alpha = 2.5
// for Finland,1995, MAKE THIS 0.8
//alpha = 2.01:
// for Italy1972, MAKE THIS 0.5
#define K_CUT_OFF               (0.7)

// probability of rewiring multiple edges or self loops
#define RETRY_P                 (0.1)
///////////////////////////////////////////////////


///////////////////////////////////////////////////
/////////// CLUSTERING / ASSORTATIVITY ////////////
///////////////////////////////////////////////////

//clustering

#define FRACT_TRI_VERTEX_MIN    (0.0)
#define FRACT_TRI_VERTEX_MAX    (0.0)
#define FRACT_TRI_VERTEX_STEP_SIZE (0.2)

//assortativity

#define ASSORT_P_MIN    	(0.0)
#define ASSORT_P_MAX    	(0.0)
#define ASSORT_P_STEP_SIZE    	(0.2)
///////////////////////////////////////////////////



///////////////////////////////////////////////////
////////////     GRID PARAMETERS     //////////////
///////////////////////////////////////////////////
//2D Grid, 4 neighbors, periodic BC
#define GRID_GRAPH 		(0)

//determine regional correlation
#define CORRELATION 		(0)

//election district regions
#define REGION_DIMENSION	(5) // < sqrt(N)

///////////////////////////////////////////////////

///////////////////////////////////////////////////
////////////     SPATIAL PARAMETERS     //////////////
///////////////////////////////////////////////////
// length scale of spatial connections
// set to zero for no spatial dependence
#define R_C			(0.0)

#define SPATIALLY_CLOSE		(0)

///////////////////////////////////////////////////


///////////////////////////////////////////////////
////////////     COMPLETE GRAPH?     //////////////
///////////////////////////////////////////////////

// Defines whether graph is all-to-all (i.e. complete)
// Ignore all previous parameters if so
#define COMPLETE		(1)

///////////////////////////////////////////////////


/********************************************************/
//////////    TEMPORAL LINKS ON THE NETWORK     //////////
/********************************************************/

#define AVERAGE_WAIT_TAU_MIN		(0.0)
#define AVERAGE_WAIT_TAU_MAX		(0.0)
#define AVERAGE_WAIT_TAU_STEP_SIZE 	(4.0)
 
#define WAIT_TAU_ALPHA_MIN		(0.0)
#define WAIT_TAU_ALPHA_MAX		(0.0)
#define WAIT_TAU_ALPHA_STEP_SIZE	(0.2)

/********************************************************/
//////////////     Infection Parameters    ///////////////
/********************************************************/

// number of infection strains
#define NUMBER_OF_INFECTIONS		(2)

//number of strains in infection [i]
// assume 1 strain per node: 2^32 strains possible
// in the future we may allow for multiple strains per node
#define STRAINS_PER_SITE		(4294967296)//(8*sizeof(unsigned int))

//Voter Model / Invasion Process ?
#define VOTER_MODEL			(0)


///////////////////////////////////////////////////
////////////    VIRALITY PARAMETER   //////////////
///////////////////////////////////////////////////

//infection probabilities
#define BETA_MIN                (1.0)
#define BETA_MAX                (1.0)
#define BETA_STEP_SIZE          (0.5)

//difference in beta between infection 1 and 2
#define DELTA_BETA              (0.0)

#define SCALE_FREE_S            (0)
#define TOP_HAT_S               (0)

///////////////////////////////////////////////////
////////////      MVM PARAMETER      //////////////
///////////////////////////////////////////////////

// Majority voter model parameter
// agent chooses majority opinion with probability MVM_P
// minority opinion with probability 1-MVM_P
#define MVM_P_MIN		(0.89)
#define MVM_P_MAX		(0.93)
#define MVM_P_STEP_SIZE		(0.005)


///////////////////////////////////////////////////
////////////  RANDOM FLIP PARAMETER  //////////////
///////////////////////////////////////////////////

//add temperature to model
#define T_MIN                   (0.0)
#define T_MAX                   (0.0)
#define T_STEP_SIZE             (0.4)


///////////////////////////////////////////////////
////////////    RECOVERY PARAMETER   //////////////
///////////////////////////////////////////////////

//recovery probability

#define GAMMA_MIN               (0.0)
#define GAMMA_MAX               (0.0)
#define GAMMA_STEP_SIZE         (0.1)

#define STUBBORN_GAMMA          (0)

///////////////////////////////////////////////////
/////////    STRAIN MUTATION PARAMETER   //////////
///////////////////////////////////////////////////

//mutation probability
// This is not the same as mu as the paper
// The paper has a variable "mu" which, in this code is 1/tau

#define MU_MIN          	(0.0)
#define MU_MAX           	(0.0)
#define MU_STEP_SIZE     	(0.4)


///////////////////////////////////////////////////
////////    TIME STUBBORNNESS PARAMETERS  /////////
///////////////////////////////////////////////////

// record delta t distribution for tau
#define RECORD_TAU              (0)

// Make deccelerating dynamics heaviside step function
#define HEAVISIDE_FREEZE        (0)

//Time-scale of "Stubbornness"/"Stickyness"
// in our paper Mu = 1/Tau. This is not the same as mu shown later in the code

#define TAU_MIN                 (1.0)//set to 0 if recording consensus
#define TAU_MAX                 (10)//set to 1 if recording consensus
#define TAU_STEP_SIZE           (1.0)

///////////////////////////////////////////////////////////////////////
////////////   Stop:Stop_hung ratio or Tau:Tau_hung ratio   ///////////
///////////////////////////////////////////////////////////////////////

// Hung jury: fraction of non-hung jury tau
// EX: tau_hung = 0.1, tau = 100 
// 	-> hung jury stubbornness rate moves from 0.01 -> 0.1
#define HUNG_RATIO_MIN            (0.0)//keep at 0 to remove hung juries
#define HUNG_RATIO_MAX            (0.2)
#define HUNG_RATIO_STEP_SIZE      (0.05)

#define Q_0              	  (0.3)

///////////////////////////////////////////////////
////////////    STOPPING PARAMETER    /////////////
///////////////////////////////////////////////////
#define STOP_MIN            	(0.01)
#define STOP_MAX            	(0.02)
#define STOP_STEP_SIZE      	(0.002)

#define RANDOM_WALK_ALPHA	(0)

///////////////////////////////////////////////////
////////////    NODE BIAS PARAMETER   /////////////
///////////////////////////////////////////////////

//SET STUBBORNS TO 0, IF NOT IN USE
//Variable stubbornness
#define NUM_STUBBORN_NODES	(0)//equally split stubborn nodes between strains

//prefer initial condition
#define PREFER_IC		(0)

//RANDOM WALK AND STUBBORN NODES ONLY
//Beta parameterizes "Scale Free S"
//average stubbornness value
//s= 0 -> not stubborn
#define S_MIN			(0.0)
#define S_MAX			(0.0)
#define S_STEP_SIZE		(2000000.0)

//Random Walk bias
//	-Time is in hours
// 	-If time is 0.0, we ignore
#define JURY_TRIAL_TIME		(0)

//STOP if values above 75% or below 25%
// stop with probability = beta/sqrt(S) * hung_p OR beta * hung_p (if NOT JURY_TRIAL_TIME)
// where 1/sqrt(S) is the timescale 
#define HUNG_P_GUILTY_MIN	(0.0)
#define HUNG_P_GUILTY_MAX	(0.0)
#define HUNG_P_GUILTY_STEP_SIZE	(10.0)

#define HUNG_P_INNOCENT_MIN	(0.0)
#define HUNG_P_INNOCENT_MAX	(0.0)
#define HUNG_P_INNOCENT_STEP_SIZE (10.0)

//OLD: s2-s1
#define DELTA_S			(0.0)	// 12 man jury: 0.025584/2					// 6  man jury: 0.09316/2

//Variable stubbornness
#define RANDOM_S		(0)

///////////////////////////////////////////////////
////////////    INITIAL CONDITIONS   //////////////
///////////////////////////////////////////////////

//initial fraction of infected nodes

#define INIT_FRACTION_INFECTED	(1.0)

// Initial condition delta  rho, between 2 infections
// DELTA_INFECTED = 0.1 => 55% infection 1, 45% infection 2
#define DELTA_INFECTED_MIN	(0.0)
#define DELTA_INFECTED_MAX	(0.0)
#define DELTA_INFECTED_STEP_SIZE (0.1)

// BINOMIAL probability of picking guilty vs innocent
#define BINOM_INIT_CONDIT_MIN	(0.5904)
#define BINOM_INIT_CONDIT_MAX	(0.5904)
#define BINOM_INIT_CONDIT_STEP_SIZE (0.01)

// probability of picking guilty vs innocent 
#define REAL_DATA_INIT_CONDIT	(0)

//for each trial, keep the infection on a random vertex with the degree AVERAGE_DEGREE
#define KEEP_ON_AVERAGE_DEGREE 	(0)

//Seed the same nodes for every Beta
#define SEED_SAME_NODES		(0)

/********************************************************/
/////////////////   Time Parameters    ///////////////////
/********************************************************/


//maximum number of ticks

#define RUN_TIME                (20000)

#define IGNORE_TIME             (19999)


/********************************************************/
/////////////////   Number Of Trials    //////////////////
/********************************************************/

#define NUM_RUNS		(5000)// number of runs per network (doesn't effect # threads)
#define NUMBER_OF_TRIALS 	(32)// number of networks (proportional to # threads)



/********************************************************
 ********************************************************
 
 PROGRAM STARTS HERE
 
 ********************************************************
 ********************************************************/


int main(int argc, char** argv) 
{

    if(ALPHA_MIN>eps)
    {
	if(VERBOSE)
   	    printf("\nWe are implimenting networks with a SCALE FREE distribution\n");
    }
    else
    {
        if(VERBOSE)
   	    printf("\nWe are implimenting networks with a POISSON distribution\n");
        if(ALPHA_MAX > 0.0)
	{
	    if(VERBOSE)
   	        printf("\nAnd we are implimenting networks with a SCALE FREE distribution\n");
	}
    }
    if(KEEP_ON_AVERAGE_DEGREE)
    {
	printf("Seeding nodes with degree = <K>\n\n");
	if((int)(INIT_FRACTION_INFECTED*MIN_NUMBER_OF_NODES)/(double)MIN_NUMBER_OF_NODES > 0.05)
	    printf("NOTE: High likelyhood we are seeding the same nodes twice. MAY EFFECT MODEL RESULTS!!\n\n");
	if(INIT_FRACTION_INFECTED > 1.0)
	{
	    printf("ERROR: Infecting more nodes than there are in the network.");
	    exit(1);
	}
    }

    time_and_run_behavior_space();
    
    return 0;
}

void time_and_run_behavior_space()

{
    //Time code, in clock cycles and time of dat
    
    //start time, in clock cycles
    clock_t start = clock();//timer
    clock_t end;
    clock_t time;

    struct timeval actual_start;
    struct timeval actual_end;
    int actual_time;
    //start time in UTC
    gettimeofday(&actual_start,NULL);

    //THIS IS THE SIMULTATION
    behavior_space();

    // end time, in clock cycles
    end=clock();
    time = end-start;

    //end time, in UTC
    gettimeofday(&actual_end,NULL);
    actual_time = actual_end.tv_sec-actual_start.tv_sec;

    char * file = "/Users/keithb/Desktop/test_log.txt";//MultiMu0.1Trial_log.dat";
    FILE * p = fopen(file,"w");
    if(p==NULL)
    {
     	printf("\nError opening file: %s\n",file);
        exit(1);
    }
    fprintf(p,"Total Time: %li s (clock cycles) or %i sec (actual)\n",time/(CLOCKS_PER_SEC/1000000),actual_time);
    fclose(p);

    if(VERBOSE)
        printf("Total Time: %li s (clock cycles) or %i sec (actual)\n",time/(CLOCKS_PER_SEC/1000000),actual_time);
    
}

void behavior_space () 
{
    //We determine the total number of threads (i.e. total number of graphs)
    int Alpha_steps = (ALPHA_MAX-ALPHA_MIN)/ALPHA_STEP_SIZE +1+eps;
    int K_steps = (AVERAGE_DEGREE_MAX-AVERAGE_DEGREE_MIN)/AVERAGE_DEGREE_STEP_SIZE +1;
    int Tri_steps = (FRACT_TRI_VERTEX_MAX-FRACT_TRI_VERTEX_MIN)/FRACT_TRI_VERTEX_STEP_SIZE +1+eps;
    int Assort_steps = (ASSORT_P_MAX-ASSORT_P_MIN)/ASSORT_P_STEP_SIZE +1+eps;
    int Node_steps = log((double)MAX_NUMBER_OF_NODES/MIN_NUMBER_OF_NODES)/log((double)STEP_SIZE_NUMBER_OF_NODES)+1+eps;
    if(MAX_NUMBER_OF_NODES == MIN_NUMBER_OF_NODES) Node_steps = 1;

    int NumThreads = NUMBER_OF_TRIALS*Alpha_steps*K_steps*Tri_steps*Assort_steps*Node_steps;
    if(REAL_DATA_DISTRIBUTION)
	NumThreads = NUMBER_OF_ELECTIONS * NUMBER_OF_TRIALS;
    if(VERBOSE)
        printf("NumThreads: %i\n",NumThreads);

    pthread_t threads [NumThreads];
    thread_args trial [NumThreads];

    //run threads
    create_threads(threads,trial);
    //writes data to file (i.e. the "export" folder)
    write_data_to_file(threads,trial);

}

void create_threads (	pthread_t * 	threads, 
			thread_args * 	trial) 
{
    int rc,thread, count=0,num_nodes;
    double alpha,K_av,fract_tri_vertices, assort_p;
    FILE * s_deg;
    FILE * s_graph; 

 
    //OTHERWISE: run simulations on with pre-defined Q and N
        //scale free coefficient
        for(alpha = ALPHA_MIN;
            alpha <= ALPHA_MAX+eps;
            alpha += ALPHA_STEP_SIZE)
        {
 	    //number of nodes
            for(num_nodes = MIN_NUMBER_OF_NODES;
                num_nodes <= MAX_NUMBER_OF_NODES;
                num_nodes *= STEP_SIZE_NUMBER_OF_NODES)
	    {
	        //average degree
     	        for(K_av = AVERAGE_DEGREE_MIN;
                    K_av <= AVERAGE_DEGREE_MAX+eps;
                    K_av += AVERAGE_DEGREE_STEP_SIZE)
                {
		    //fraction of triangles
                    for(fract_tri_vertices=FRACT_TRI_VERTEX_MIN;
                        fract_tri_vertices<=FRACT_TRI_VERTEX_MAX + eps;
                        fract_tri_vertices +=FRACT_TRI_VERTEX_STEP_SIZE)
                    {
		        //assortativity/dis-assortativity parameter
             	        for(assort_p = ASSORT_P_MIN;
                            assort_p <= ASSORT_P_MAX + eps;
                            assort_p += ASSORT_P_STEP_SIZE)
                        {
	 	 	   //number of trials, holding graph parameters fixed
                            for (thread=count;
                                 thread<count + NUMBER_OF_TRIALS;
                                 ++thread)
                            {
				
			        //make the arguments for each thread
                     	        make_thread_args(&trial[thread],
                                                 thread,
                                                 num_nodes,
				    		 &s_deg,
				    		 &s_graph,
                                                 alpha,
                                                 K_av,
                                                 fract_tri_vertices,
                                                 assort_p,
			 	 	         NUMBER_OF_INFECTIONS);
			        //create threads
				
                                rc = pthread_create(&threads[thread],
                                                    NULL,
                                                    perform_trials,
                                                    (void *) &trial[thread]);
			        //assert no errors
                                assert(0 == rc);
							        
                                if(VERBOSE)
                                    printf("Starting Thread #: %i\n",thread);
			
                            }
                            count += NUMBER_OF_TRIALS;
		        }
		    }
	        }
            }
        }
    
}

void write_data_to_file (pthread_t * threads, thread_args * trial) 
{
    
    FILE * p = NULL;
    //file name/directory

    const char * file = "/Volumes/2TB Drive/2StrainJuryDelib-1min;Beta=1;N=12;Binomial_p=0.5904_1_0.5904;MVM_p=0.89_0.005_0.93;Alpha=0.01_0.002_0.02;Alpha_hung=0_1_0;Tau=1_1_10;HungRatio=0.0_0.05_0.2;NumRuns=5000;NumTrials=32;Q0=0.3;CompleteGraph;AllInfected;OR12_Criminal_NoStubbornnessHung.dat";


    int rc,thread, count=0, i, infection, num_nodes;
    //variables used to define total write space
    int gamma_steps = (GAMMA_MAX-GAMMA_MIN)/GAMMA_STEP_SIZE + 1 + eps;
    int beta_steps = (BETA_MAX - BETA_MIN)/BETA_STEP_SIZE + 1 + eps;
    int s_steps = (S_MAX - S_MIN)/S_STEP_SIZE + 1 + eps;
    int delta_infected_steps = (DELTA_INFECTED_MAX - DELTA_INFECTED_MIN)/DELTA_INFECTED_STEP_SIZE + 1 + eps;
    int hung_p_guilty_steps = (HUNG_P_GUILTY_MAX - HUNG_P_GUILTY_MIN)/HUNG_P_GUILTY_STEP_SIZE + 1 + eps;
    int hung_p_innocent_steps = (HUNG_P_INNOCENT_MAX - HUNG_P_INNOCENT_MIN)/HUNG_P_INNOCENT_STEP_SIZE + 1 + eps;
    int stop_steps = (STOP_MAX - STOP_MIN)/STOP_STEP_SIZE + 1 + eps;
    int binom_steps = (BINOM_INIT_CONDIT_MAX - BINOM_INIT_CONDIT_MIN)/BINOM_INIT_CONDIT_STEP_SIZE + 1 + eps;
    int mvm_p_steps = (MVM_P_MAX - MVM_P_MIN)/MVM_P_STEP_SIZE + 1 + eps;
    int tau_steps = (TAU_MAX - TAU_MIN)/TAU_STEP_SIZE + 1 + eps;
    int hung_steps = (HUNG_RATIO_MAX - HUNG_RATIO_MIN)/HUNG_RATIO_STEP_SIZE + 1 + eps;
    int mu_steps = (MU_MAX - MU_MIN)/MU_STEP_SIZE + 1 + eps;
    int temp_steps = fabs(T_MAX - T_MIN)/T_STEP_SIZE + 1 + eps;
    int average_wait_tau_steps = (AVERAGE_WAIT_TAU_MAX - AVERAGE_WAIT_TAU_MIN)/AVERAGE_WAIT_TAU_STEP_SIZE + 1 + eps;
    int wait_tau_alpha_steps = (WAIT_TAU_ALPHA_MAX - WAIT_TAU_ALPHA_MIN)/WAIT_TAU_ALPHA_STEP_SIZE + 1 + eps;

    int write_space = NUM_RUNS * gamma_steps * temp_steps * beta_steps * s_steps * delta_infected_steps * hung_p_guilty_steps * hung_p_innocent_steps * stop_steps * binom_steps * mvm_p_steps * hung_steps * tau_steps * mu_steps * average_wait_tau_steps * wait_tau_alpha_steps * NUMBER_OF_INFECTIONS;
    if (VERBOSE)
        printf("WRITE SPACE: %i\n",write_space);
    double K_av,alpha,beta,gamma,fract_tri_vertices,assort_p,wait_tau_alpha,wait_tau,mvm_p,mu,temp,S,hung_p_guilty, hung_p_innocent,delta_inf,hung,T,stop_rate,binom;
    int run;

    if(VERBOSE)
        printf("\nOpening file %s \n\n",file);
    //open file we use to write data to
    p = fopen(file, "w");//in order to write everything on ONE text file, we open the file here
    if(p==NULL)
    {
     	printf("\nError opening file: %s\n",file);
        exit(1);
    }
    else if(VERBOSE)
	printf("file opened\n");
    
    //Statement of Behavior Space Selected
    
    if(VERBOSE)
    {
    printf( "Behavior space recorded in the following order: \n\nReal Data? %i\n Number of elections: %i\n VM? %i\n Stubborn(0)/Initial Preference(1): %i\nNumber of stubborn nodes: %i\nDependent on jury trial time? %i\tP(stop if NOT hung) 1-SUPER_MAJORITY (%f) threshold %f:%f:%f, SUPER_MAJORITY (%f) threshold %f:%f:%f\n Infect in(-1)/neutral(0)/out(1)?: %i\n Number of nodes=%i:%f:%i\n sample fraction=%f\n Grid? %i\tCorrelation? %i\tR_c = %f\talpha=%f:%f:%f\n <K>=%f:%f:%f\n fraction of links that are edges of triangles: %f:%f:%f\n assort_p=%f:%f:%f\n alpha coefficient for wait time = %f:%f:%f\n wait time = %f:%f:%f\n mutation probability = %f:%f:%f\n stubborn recovery? %i \trecovery probability=%f:%f:%f\n flip rate = %f:%f:%f \n infection probability=%f:%f:%f\n Stubbornness (S) = %f:%f:%f \tscale free s? %i\n Heaviside freeze (0 = No/1 = Yes): %i \t MVM_p = %f:%f:%f \t stubborn time-scale=%f:%f:%i \t hung ratio =%f:%f:%f \tstopping rate if not hung = |1/2 - v_g/N|*(%f:%f:%f), binomial_p: %f:%f:%f \n infection #:%i:%i:%i\n initial fraction infected = %f\tinfect candidates? %i\n delta infected = %f:%f:%f\n Number of Runs=%i \n Trials=%i\n transient time = %i\n run time = %i\n\n\n",
               REAL_DATA_DISTRIBUTION, NUMBER_OF_ELECTIONS,
               VOTER_MODEL,
               PREFER_IC,
	       NUM_STUBBORN_NODES,
	       JURY_TRIAL_TIME,
	       1-VOTE_THRESHOLD,
	       HUNG_P_INNOCENT_MIN,HUNG_P_INNOCENT_STEP_SIZE,HUNG_P_INNOCENT_MAX,
	       VOTE_THRESHOLD,
	       HUNG_P_GUILTY_MIN, HUNG_P_GUILTY_STEP_SIZE,HUNG_P_GUILTY_MAX,
               -1 * INFECT_INWARD + 0 * INFECT_NEUTRAL + INFECT_OUTWARD,
               MIN_NUMBER_OF_NODES, STEP_SIZE_NUMBER_OF_NODES, MAX_NUMBER_OF_NODES,
               SAMPLE_FRACTION,
	       GRID_GRAPH,
	       CORRELATION,
	       R_C,
               ALPHA_MIN, ALPHA_STEP_SIZE, ALPHA_MAX,
               AVERAGE_DEGREE_MIN, AVERAGE_DEGREE_STEP_SIZE, AVERAGE_DEGREE_MAX,
               FRACT_TRI_VERTEX_MIN, FRACT_TRI_VERTEX_STEP_SIZE, FRACT_TRI_VERTEX_MAX,
               ASSORT_P_MIN, ASSORT_P_STEP_SIZE, ASSORT_P_MAX,
               WAIT_TAU_ALPHA_MIN, WAIT_TAU_ALPHA_STEP_SIZE, WAIT_TAU_ALPHA_MAX,
               AVERAGE_WAIT_TAU_MIN, AVERAGE_WAIT_TAU_STEP_SIZE, AVERAGE_WAIT_TAU_MAX,
               MU_MIN, MU_STEP_SIZE, MU_MAX,
               STUBBORN_GAMMA,GAMMA_MIN, GAMMA_STEP_SIZE, GAMMA_MAX,
               T_MIN, T_STEP_SIZE, T_MAX,
               BETA_MIN, BETA_STEP_SIZE, BETA_MAX,
               S_MIN, S_STEP_SIZE, S_MAX,
	       SCALE_FREE_S,
	       HEAVISIDE_FREEZE,
               MVM_P_MIN, MVM_P_STEP_SIZE, MVM_P_MAX,
               TAU_MIN, TAU_STEP_SIZE, TAU_MAX,
               HUNG_RATIO_MIN, HUNG_RATIO_STEP_SIZE, HUNG_RATIO_MAX,
               STOP_MIN, STOP_STEP_SIZE, STOP_MAX,
	       BINOM_INIT_CONDIT_MIN, BINOM_INIT_CONDIT_STEP_SIZE, BINOM_INIT_CONDIT_MAX,
               1, 1, NUMBER_OF_INFECTIONS,
               INIT_FRACTION_INFECTED,
	       INFECT_CANDIDATES,
               DELTA_INFECTED_MIN,DELTA_INFECTED_STEP_SIZE, DELTA_INFECTED_MAX,
               NUM_RUNS,
               NUMBER_OF_TRIALS,
               IGNORE_TIME,
               RUN_TIME);
    }

    //parameter summary, written into file
    
    fprintf(p,
               "Behavior space recorded in the following order: \n\nReal Data? %i\n Number of elections: %i\n VM? %i\n Stubborn(0)/Initial Preference(1): %i\tNumber of stubborn nodes: %i\tDependent on jury trial time? %i\tP(stop if NOT hung) 1-SUPER_MAJORITY (%f) threshold %f:%f:%f, SUPER_MAJORITY (%f) threshold %f:%f:%f\n Infect in(-1)/neutral(0)/out(1)?: %i\n Number of nodes = %i : %f : %i\n sample fraction = %f \n Grid? %i\tCorrelation? %i\tR_c = %f\talpha = %f : %f : %f\n <K> = %f : %f : %f\n fraction of links that are edges of triangles = %f : %f : %f\n assort_p = %f : %f : %f\n alpha coefficient for wait time = %f : %f : %f\n wait time = %f : %f : %f\n mutation probability = %f : %f : %f\n stubborn recovery? %i \trecovery probability = %f : %f : %f\n flip rate = %f : %f : %f \n infection probability = %f : %f : %f \n stubbornness = %f : %f : %f\tscale free s? %i \n Heaviside freeze (0 = No/1 = Yes): %i \t MVM_p = %f:%f:%f \t stubborn time-scale=%f:%f:%i \t Remove ALL hung conditions %i \t Remove only alpha hung conditions %i \t remove only stubbornness hung conditions %i \t remove q time dependence %i \t remove s time dependence %i \t hung ratio =%f:%f:%f  \tstopping rate if not hung = |1/2 - v_g/N|*(%f:%f:%f), binomial_p: %f:%f:%f \n infection #: %i:%i:%i \n initial fraction infected = %f\tinfect candidates? %i \n delta infected = %f:%f:%f \n Number of Runs= %i \n Trials= %i \n transient time = %i \n run time = %i\n\n\n",
	       REAL_DATA_DISTRIBUTION, NUMBER_OF_ELECTIONS,
               VOTER_MODEL,
               PREFER_IC,
	       NUM_STUBBORN_NODES,
	       JURY_TRIAL_TIME,
	       1-VOTE_THRESHOLD,
	       HUNG_P_INNOCENT_MIN,HUNG_P_INNOCENT_STEP_SIZE,HUNG_P_INNOCENT_MAX,
	       VOTE_THRESHOLD,
	       HUNG_P_GUILTY_MIN,HUNG_P_GUILTY_STEP_SIZE,HUNG_P_GUILTY_MAX,
	       -1 * INFECT_INWARD + 0 * INFECT_NEUTRAL + INFECT_OUTWARD,
	       MIN_NUMBER_OF_NODES, STEP_SIZE_NUMBER_OF_NODES,	MAX_NUMBER_OF_NODES,
	       SAMPLE_FRACTION,
	       GRID_GRAPH,
	       CORRELATION,
	       R_C,
               ALPHA_MIN, ALPHA_STEP_SIZE, ALPHA_MAX,
               AVERAGE_DEGREE_MIN, AVERAGE_DEGREE_STEP_SIZE, AVERAGE_DEGREE_MAX,
               FRACT_TRI_VERTEX_MIN, FRACT_TRI_VERTEX_STEP_SIZE, FRACT_TRI_VERTEX_MAX,
               ASSORT_P_MIN, ASSORT_P_STEP_SIZE, ASSORT_P_MAX,
               WAIT_TAU_ALPHA_MIN, WAIT_TAU_ALPHA_STEP_SIZE, WAIT_TAU_ALPHA_MAX,
               AVERAGE_WAIT_TAU_MIN, AVERAGE_WAIT_TAU_STEP_SIZE, AVERAGE_WAIT_TAU_MAX,
       	       MU_MIN, MU_STEP_SIZE, MU_MAX,
               STUBBORN_GAMMA,GAMMA_MIN, GAMMA_STEP_SIZE, GAMMA_MAX,
               T_MIN, T_STEP_SIZE, T_MAX,
               BETA_MIN, BETA_STEP_SIZE, BETA_MAX,
               S_MIN, S_STEP_SIZE, S_MAX,
	       SCALE_FREE_S,
	       HEAVISIDE_FREEZE,
               MVM_P_MIN, MVM_P_STEP_SIZE, MVM_P_MAX,
               TAU_MIN, TAU_STEP_SIZE, TAU_MAX,
	       REMOVE_HUNG_CONDITIONS,
	       REMOVE_STOP_HUNG_CONDITION,
	       REMOVE_STUBBORNNESS_HUNG_CONDITION,
	       REMOVE_Q_TIME_DEPENDENCE,
	       REMOVE_STUBBORN_TIME_DEPENDENCE,
               HUNG_RATIO_MIN, HUNG_RATIO_STEP_SIZE, HUNG_RATIO_MAX,
               STOP_MIN, STOP_STEP_SIZE, STOP_MAX,
	       BINOM_INIT_CONDIT_MIN, BINOM_INIT_CONDIT_STEP_SIZE, BINOM_INIT_CONDIT_MAX,
               1, 1, NUMBER_OF_INFECTIONS,
               INIT_FRACTION_INFECTED,
	       INFECT_CANDIDATES,
               DELTA_INFECTED_MIN,DELTA_INFECTED_STEP_SIZE, DELTA_INFECTED_MAX,
	       NUM_RUNS,
               NUMBER_OF_TRIALS,
  	       IGNORE_TIME,
  	       RUN_TIME);


    /***********************************************************************/
    ///////////// We Record Simulation Data Here /////////////
    /***********************************************************************/
 
    
    //IF NOT ELECTION DATA: this code writes data to file
	//scale-free degree distribution coefficient
        for (alpha = ALPHA_MIN;
            alpha <= ALPHA_MAX+eps;
            alpha += ALPHA_STEP_SIZE)
        {
	    //number of nodes
       	    for(num_nodes = MIN_NUMBER_OF_NODES;
                num_nodes <= MAX_NUMBER_OF_NODES;
                num_nodes *= STEP_SIZE_NUMBER_OF_NODES)
            {
		//average degree
                for (K_av = AVERAGE_DEGREE_MIN;
                     K_av <= AVERAGE_DEGREE_MAX+eps;
                     K_av += AVERAGE_DEGREE_STEP_SIZE)
            	{
		    //fraction of links that are trangles (clustering)
                    for(fract_tri_vertices=FRACT_TRI_VERTEX_MIN;
                        fract_tri_vertices<=FRACT_TRI_VERTEX_MAX + eps;
                        fract_tri_vertices +=FRACT_TRI_VERTEX_STEP_SIZE)
                    {
			//assortativity/dis-assortativity parameter
             	        for(assort_p = ASSORT_P_MIN;
                            assort_p <= ASSORT_P_MAX +eps;
                            assort_p += ASSORT_P_STEP_SIZE)
                        {
			    //scale-free coefficient for wait-times
                            i = 0;
                            //fprintf(p,"\n");
	                    for(wait_tau_alpha = WAIT_TAU_ALPHA_MIN; wait_tau_alpha <= WAIT_TAU_ALPHA_MAX + eps; wait_tau_alpha += WAIT_TAU_ALPHA_STEP_SIZE)
		            {
				//average waiting time before links are available
                                //fprintf(p,"\n");
                 	        for(wait_tau = AVERAGE_WAIT_TAU_MIN; wait_tau <=AVERAGE_WAIT_TAU_MAX +eps; wait_tau += AVERAGE_WAIT_TAU_STEP_SIZE)
			        {
				    //probability of random mutations to a new opinion
                    	            //fprintf(p,"\n");
	       	       	            for(mu = MU_MIN; mu <= MU_MAX+eps; mu += MU_STEP_SIZE)
			            {
					//recovery probability
                    	 	        //fprintf(p,"\n");
                    		        for (gamma=GAMMA_MIN; gamma<=GAMMA_MAX+eps;gamma+=GAMMA_STEP_SIZE)
                    		        {
					    //spin-flip rate ("temperature")
                     		            //fprintf(p,"\n");
                    		            for (temp=T_MIN; temp<=T_MAX+eps;temp+=T_STEP_SIZE)
                    		            {
						//P simulation stops when below the threshold
                     		                //fprintf(p,"\n");
                            		        for (delta_inf=DELTA_INFECTED_MIN; delta_inf<=DELTA_INFECTED_MAX+eps; delta_inf+=DELTA_INFECTED_STEP_SIZE)
                            			{
					 	    //P simulation stops when below the threshold
                     		                    //fprintf(p,"\n");
                            		            for (hung_p_innocent=HUNG_P_INNOCENT_MIN; hung_p_innocent<=HUNG_P_INNOCENT_MAX+eps; hung_p_innocent+=HUNG_P_INNOCENT_STEP_SIZE)
                            			    {

						        //P simulation stops when above the threshold
                     		                        //fprintf(p,"\n");
                            			        for (hung_p_guilty=HUNG_P_GUILTY_MIN; hung_p_guilty<=HUNG_P_GUILTY_MAX+eps; hung_p_guilty+=HUNG_P_GUILTY_STEP_SIZE)
                            			        {
							    //infection probability
                     		                	    //fprintf(p,"\n");
                        	                	    for (beta=BETA_MIN; beta<=BETA_MAX+eps; beta+=BETA_STEP_SIZE)
                        	                	    {

 	 					                //node-based bias/stubbornness
                     		                                //fprintf(p,"\n");
                            			                for (S=S_MIN; S<=S_MAX+eps; S+=S_STEP_SIZE)
                            			                {
						                    //MVM_p: probability one adopts majority vs minority opinion
                            		                            //fprintf(p,"\n");
                            		                            for(mvm_p = MVM_P_MIN; mvm_p<=MVM_P_MAX+eps; mvm_p += MVM_P_STEP_SIZE)
                            		                            {

 						                        //stubbornness time-scale
                             		                                //fprintf(p,"\n");
                            		                                for(T = TAU_MIN; T<=TAU_MAX+eps; T += TAU_STEP_SIZE)
                            		                                {
						                            //hung stubbornness time-scale
                            		                                    //fprintf(p,"\n");
                            		                                    for(hung = HUNG_RATIO_MIN; hung<=HUNG_RATIO_MAX+eps; hung += HUNG_RATIO_STEP_SIZE)
                            		                                    {
 						                                // stop rate: stop_rate * |1/2 - v_g/N|, v/N >= VOTE_THRESHOLD, or v/N <=  1-VOTE_THRESHOLD
                             		                                        //fprintf(p,"\n");
                            		                                        for(stop_rate = STOP_MIN; stop_rate<=STOP_MAX+eps; stop_rate += STOP_STEP_SIZE)
                            		                                        {	
						                                    //hung stubbornness time-scale
                            		                                            //fprintf(p,"\n");
										    for(binom = BINOM_INIT_CONDIT_MIN; binom<=BINOM_INIT_CONDIT_MAX+eps; binom += BINOM_INIT_CONDIT_STEP_SIZE)
                            		                                            {
						              	                        //run number, holding the network fixed
                            		                                                //fprintf(p,"\n");
							   	                        for(run = 0; run < NUM_RUNS; run++)
							   	                        {
							                                    //infection #
                             		                                                    for(infection = 0; infection<NUMBER_OF_INFECTIONS; infection ++)
                                	                                                    {
							                                        // trial number, holding parameters constant
                                    		                                                for (thread = count;
                                        	                                                thread < count + NUMBER_OF_TRIALS;
                                        	                                                ++thread)
                                    		                                                {
  				       		                                                    //We wait for the thread to finish for each trial
				        	                                                    if(i==0)
			  	        	                                                    {
								                                        //join thread
		                            		                                                rc = pthread_join(threads[thread],(void *)&trial[thread]);
                		            		                                                assert(0 == rc);
				  	    	 	                                                if(VERBOSE)
                		                	                                                    printf("\n\nThread #%i finished:\n\tAlpha: %f, K: %f, Assort_p:%0.5f, Fract_tri:%0.5f, Trial: %i\n\n",thread,alpha,K_av,assort_p,fract_tri_vertices,thread);
 				         	                                                    }
							   	                                    //WRITE DATA TO FILE HERE
				        	                                                    write_simulation_data_to_file(&trial[thread],i,infection,T,p);
                                    		                                                 }
                                    		                                                 i++;
  					                                                         if(i==(write_space))
					                                                         {
				       		                                                     for (thread = count;
                                        		                                             thread < count + NUMBER_OF_TRIALS;
                                        		                                             ++thread)
                                        	                                                     {
								                                         //WRITE GRAPH TO FILE once all data recorded
							                                                 if(RECORD_GRAPH)
				            		                                                     write_graph_to_file(&trial[thread],p);
					   		                                                 fprintf(p,"\n");
		                            		                                                 free_thread_args(&trial[thread]);
										                     }
						                                                }
									                    }
											}
										    }
										}
									    }
									}
							            }
							        }
							    }
						        }
						    }
					        }
					    }
				        }
			            }
                                }
                            }
                            count+=NUMBER_OF_TRIALS;
		        }
		    }
                }
            }
        }
    
       
    if(VERBOSE)
        printf("\nWritten Successfully to: %s\n",file);

    fclose(p);
}

void write_simulation_data_to_file (	thread_args * 	trial,
					int 		i,
					int 		infection,
					double 		tau,
					FILE * 		p) 
{

    int T,K,hist_val,j,t_infected,N = trial->G.NumNodes,now,end,Q=trial->Q;
    double mean,m,m2,m4,rho_t,time;

    if(!RECORD_LONG_CONSENSUS_TIME)
    {
        //mean density of infection_i
        mean = trial->Mean[i];
        m = trial->m[i/Q];
        m2 = trial->m2[i/Q];
        m4 = trial->m4[i/Q];
        fprintf(p,"%.15lf\t%.15lf\t%.15lf\t%.15lf\n\n",mean,m,m2,m4);
	if(VERBOSE)
	    printf("%.15lf\t%.15lf\t%.15lf\t%.15lf\n\n",mean,m,m2,m4);
    }

    if(CORRELATION)
    {
	int r,k;
	double rho_i;
        int l = REGION_DIMENSION;
        int L = sqrt((double)(N)) + eps;
        int num_regions = (L/l) * (L/l);

        //bin width
        int r_bin_width = 1;

        //maximum r
        int r_max = (L/l)/sqrt(2) + 1;

	double c_r;

	fprintf(p,"\n");
        //starting from nearest neighbor, going outward
	for(r = 1; r < r_max; r+= r_bin_width)
	{
	    c_r = trial -> Correlation[i * r_max + r];
	    fprintf(p,"%lf\n",c_r);
	}
	
	fprintf(p,"\n");
	for(k = 0; k < num_regions; k++)
	{
	    //NOTE this only captures the results of 1 TRIAL
	    //This is a test of the process NOT the final code
	    rho_i = trial -> Rho_i[i * num_regions + k];
	    fprintf(p,"%lf\n",rho_i);
	}
    }

    if(RECORD_EXPOSURE)
    {
	fprintf(p,"\n");
        for(T = 0; T < TAU_MAX * 10; ++T)
        {
	    hist_val = trial -> TOTAL_TimeUntilInfectedHist[i * (TAU_MAX * 10 + 1) + T];
            fprintf(p,"\t %i",hist_val);
        }
	fprintf(p,"\n");
	int av_deg = AVERAGE_DEGREE_MAX;
        for(K = 0; K < (int)(av_deg * TAU_MAX * 10); ++K)
        {
	    hist_val = trial -> TOTAL_FreqExposureUntilInfectedHist[i * 10 * TAU_MAX * av_deg + K];
            fprintf(p,"\t %i",hist_val);
        }
	fprintf(p,"\n");

    }

    if(RECORD_LONG_CONSENSUS_TIME && i % NUMBER_OF_INFECTIONS == 0)
    {
        //record initial vote here: t0 = 1, 2, 3, 4, 5

        if(RECORD_VOTE_TIME)
        {
	    double v,t0 = 2.5,epsilon = 1/(2*((double)N));
            for (t0=1.0;t0<6.0;t0++)
            {
                v = trial->v0[i/Q + (int)t0-1];
                fprintf(p,"%f\n",v);
            }
        }


        // record final time here
	time = trial->ConsensusTimes[i/NUMBER_OF_INFECTIONS];
	//if(RECORD_TIME_VOTE)
	fprintf(p,"\n %f",time);
	if(RECORD_DIST)
	{
            int votes = trial->NumsGuilty[i/NUMBER_OF_INFECTIONS];
	    //if(RECORD_TIME_VOTE)
	    fprintf(p,"\t%i",votes);
	}
	if(RECORD_FREEZE_P)
        {
            double Freeze = trial->FreezeTime[i/NUMBER_OF_INFECTIONS];
            double P = trial->DynamicsP[i/NUMBER_OF_INFECTIONS];
            fprintf(p,"\t%f\t%f",Freeze,P);
        }
    }

    //record rho(t,t') before equilibrium
    if(TAU_MAX>0)
    { //histogram of time_infected
	if(RECORD_CONSENSUS_TIME)
	{
	    if(mean < 0.5)
	    {
		       
		end = (RUN_TIME - IGNORE_TIME) * (TAU_MAX + 1) * i + (RUN_TIME - IGNORE_TIME - 1) * (TAU_MAX + 1) + 1;

                for(time = RUN_TIME - IGNORE_TIME - 1; time >= 0; -- time)
                {
		    now = (RUN_TIME - IGNORE_TIME) * (TAU_MAX + 1) * i + (time) * (TAU_MAX + 1) + 1;
		    //even if there is no tau, we set tau to 1 in order to determine the consensus time
		    if(tau < 1.0) tau = 1.0;
		    if(sum_int(&trial->T_distribution[now],tau) > sum_int(&trial->T_distribution[end],tau))
		    {
			if(time == RUN_TIME - IGNORE_TIME - 2)
			    break;

  		        fprintf(p,"\n %f",time + IGNORE_TIME - 1);
		        break;
		    }
		}
	    }
	}

	if(RECORD_TIME_INFECTED)
	{
	    //we make sure to record the number of infected individuals even when there is no tau

            for(time = 0; time< RUN_TIME - IGNORE_TIME; ++ time)
            {
		now = (RUN_TIME - IGNORE_TIME) * i + time;

                rho_t = trial->FractionInfectedVersusTime[now];
                fprintf(p,"%0.10f\t",rho_t);		

		if(RECORD_TAU)
		{
		    now = (RUN_TIME - IGNORE_TIME) * (TAU_MAX + 1) * i + (time) * (TAU_MAX + 1) + 1;
             	    fprintf(p,"\n");
	            if(tau < 1.0) tau = 1;
                    for(t_infected = 0; t_infected < tau; ++t_infected)
                    {
                        fprintf(p,"%i\t",trial->T_distribution[now + t_infected]);
		    }
	        }
            }
        }
    }

    if(RECORD_INFECTIONS_ON_GRAPH && infection==0)
    {
	fprintf(p,"\n\n");
	for(j=0; j<N; ++j)
	{
	    fprintf(p,"%i\t",trial->whos_infected[i/NUMBER_OF_INFECTIONS * N + j]);
	}
    }
    fprintf(p,"\n\n\n");
}

void write_graph_to_file (	thread_args * 	trial, 
				FILE * 		p) 
{
    int i,j;
    fprintf(p,"\n\n");

    //a different graph for each trial
    for (i=0; i<trial->G.NumNodes; ++i)
    {
        for(j=0; j<trial->G.nodes[i].degree;++j)
        {
       	    fprintf(p,"%i ",trial->G.nodes[i].neighbors[j]);
        }
        fprintf(p,"\n");
    }
}

/****************************************************************
following function was partly barrowed from stackoverflow forum
****************************************************************/

void seed_rand(int thread_n, gsl_rng* r) 
{
    struct timeval tv;
    
    gettimeofday(&tv, NULL);
    gsl_rng_set(r,tv.tv_sec * thread_n + tv.tv_usec);
}
////////////////////////////////////////////////////////////

void* perform_trials(void * trial) 
{

    thread_args* t = (thread_args *) trial;//equate pointers
    
    double gamma=0.0,beta=0.0,temp=0.0,wait_tau_alpha=0.0,average_wait_tau=0.0,tau=0.0,mu=0.0,S=0.0,hung_p_guilty = 0.0, hung_p_innocent = 0.0, delta_inf = 0.0;
    double mvm_p = 0.0, hung = 0.0, stop_rate = 0.0, binom = 0.0;
    int num_runs, count = 0;
    long long Beginning_Of_Array=0;
    //For each graph/thread run simulations of the following:

    for (wait_tau_alpha = WAIT_TAU_ALPHA_MIN; wait_tau_alpha <= WAIT_TAU_ALPHA_MAX+eps; wait_tau_alpha +=WAIT_TAU_ALPHA_STEP_SIZE)
    {
	for(average_wait_tau = AVERAGE_WAIT_TAU_MIN; average_wait_tau <= AVERAGE_WAIT_TAU_MAX+eps; average_wait_tau +=AVERAGE_WAIT_TAU_STEP_SIZE)
	{
	    for(mu = MU_MIN; mu <= MU_MAX + eps; mu += MU_STEP_SIZE)
	    {
                for (gamma = GAMMA_MIN; gamma<=GAMMA_MAX+eps;gamma+=GAMMA_STEP_SIZE)
                {
                    for (temp=T_MIN; temp<=T_MAX+eps; temp+=T_STEP_SIZE)
		    {
                        for (delta_inf=DELTA_INFECTED_MIN; delta_inf<=DELTA_INFECTED_MAX+eps; delta_inf+=DELTA_INFECTED_STEP_SIZE)
                        {
                            for (hung_p_innocent=HUNG_P_INNOCENT_MIN; hung_p_innocent<=HUNG_P_INNOCENT_MAX+eps; hung_p_innocent+=HUNG_P_INNOCENT_STEP_SIZE)
                            {
                                for (hung_p_guilty=HUNG_P_GUILTY_MIN; hung_p_guilty<=HUNG_P_GUILTY_MAX+eps; hung_p_guilty+=HUNG_P_GUILTY_STEP_SIZE)
                                {
                        	    for (beta=BETA_MIN; beta<=BETA_MAX+eps; beta+=BETA_STEP_SIZE)
				    {
                                        for (S=S_MIN; S<=S_MAX+eps; S+=S_STEP_SIZE)
  	 		                {
                                            for (mvm_p = MVM_P_MIN; mvm_p <= MVM_P_MAX+eps; mvm_p += MVM_P_STEP_SIZE)
                                            {
                                                for (tau = TAU_MIN; tau <= TAU_MAX+eps; tau += TAU_STEP_SIZE)
                                                {
 	                                            for (hung = HUNG_RATIO_MIN; hung <= HUNG_RATIO_MAX+eps; hung += HUNG_RATIO_STEP_SIZE)
        	                                    {
							for(binom = BINOM_INIT_CONDIT_MIN; binom<=BINOM_INIT_CONDIT_MAX+eps; binom += BINOM_INIT_CONDIT_STEP_SIZE)
                                            		{
                                            		    for (stop_rate = STOP_MIN; stop_rate <= STOP_MAX+eps; stop_rate += STOP_STEP_SIZE)
                                            		    {

		                                                Beginning_Of_Array = count * (RUN_TIME - IGNORE_TIME) * (TAU_MAX + 1);

 	 		      			                for(num_runs = 0; num_runs < NUM_RUNS; num_runs ++)
	  					                {
                                                		    run_experiment (
		                                                        &t->G,
					       	                        t->FractionInfectedVersusTime,
			       			           	        t->UVersusTime,
                                            			        Beginning_Of_Array,
						            	        t->T_distribution,
     		                	                       	        beta,
                                        		    	        S,
                                            	        		delta_inf,
                         		                   	        gamma,
							    	        temp,
				            	        		mu,
									mvm_p,
        		                                    	        tau,
        		                                    	        hung,
									stop_rate,
									binom,
					 	            	        average_wait_tau,
			 			            	        wait_tau_alpha,
					      			        hung_p_innocent,
					      			        hung_p_guilty,
                                                		        t->r);
		 		    	
					
		 		                                    data_analysis(t,&count);
							        }   
							    }   
							}   
						    }   
						}   
					    }   
					}   
				    }   
			        }   
			    }   
			}   
                    }
		}
            }
        }
    }
}

void data_analysis (thread_args * t, int * count) 
{
    ////////////////////VARIABLES FOR DATA ANALYSIS/////////////////////
    int i, N = t->G.NumNodes, Q = t->Q;
    double mean;
    ////////////////////////////////////////////////////////////////////
    /////////////////////// DATA ANALYSIS /////////////////////

    // on a 2D grid, we can define "regions",
    // whose candidate preference we average over
    if(CORRELATION)
    {
	int r,n;
	int l = REGION_DIMENSION;
        int L = sqrt((double)(N)) + eps;

        //bin width
        int r_bin_width = 1;

        //maximum r
        int r_max = (L/l)/sqrt(2) + 1;

	//number of regions
	int num_regions = (L/l) * (L/l);

	double cross_correlation;
	int num_regions_r;

	//determine mean candidate preference for region
	determine_regional_rho(&(t->G));
	for(n= 0; n<num_regions;++n)
	{
	    t->Rho_i[(*count) * num_regions + n] = t->G.Rho_i[n];
	}

	// starting from nearest neighbor, going outward
	// we find the correlation, C(r)
	for(r = 1; r < r_max; r+= r_bin_width)
	{
	    cross_correlation_r (&(t->G), r, &cross_correlation, &num_regions_r);
	    t -> Correlation[ (*count) * r_max + r] = spatial_correlation (&(t->G), r, cross_correlation, num_regions_r);
	}
    }

    //record infection for each node at a time "t"
    if(RECORD_INFECTIONS_ON_GRAPH)
    {
        int infection_count = (*count)/Q * N, j,strain;
        for(j=0; j<N; ++j)
        {
	    strain = t->G.nodes[j].infection[0];
	    t->whos_infected[infection_count + j] = strain;
        }  
    }

    if(RECORD_VOTE_TIME)
    {
        double t0 = 2.5,epsilon = 1/(2*((double)N));
        for (t0=1.0;t0<6.0;t0++)
        {
            t->v0[(*count)/Q + (int)t0-1] = t->G.v0[(int)t0-1];
        }
    }


    //we record consensus time AND NOTHING ELSE
    if(RECORD_LONG_CONSENSUS_TIME && (*count) % Q == 0)
    {
        //printf("%i\n",(*count)/Q);
	t->ConsensusTimes [(*count)/Q] = t->G.ConsensusTime;
	if(RECORD_DIST)
 	    t->NumsGuilty [(*count)/Q] = t->G.NumGuilty;
	if(RECORD_FREEZE_P)
        {
            t->FreezeTime[(*count)/Q] = t->G.nodes[0].time_last_checked;
            t->DynamicsP[(*count)/Q] = t->G.nodes[0].current_stubbornness_p;
        }
    }
    if(!RECORD_LONG_CONSENSUS_TIME)
    {
        //recording <m^4>
    	if((*count) % Q == 0)
    	{
            //Umean = gsl_stats_mean (t->UVersusTime, STRIDE, n);
	
	    t->m[(*count)/Q] = t->G.m;
	    t->m2[(*count)/Q] = t->G.m2;
	    t->m4[(*count)/Q] = t->G.m4;
    	}
    }

    for(i=0; i<Q; ++i)
    {
    	if(!RECORD_LONG_CONSENSUS_TIME)
	{
	    mean = 	t->G.MeanFractionInfected[i];

	    t->Mean[*count] = mean;
	    if(REAL_DATA_DISTRIBUTION)
	        t->Mean[(*count)] = mean * (t->G.Q);
            //if(VERBOSE)
	    //    printf("Infection: %i, Mean: %.8lf\n",i+1,mean);
	}
        (*count)++;
    }   
}

void determine_regional_rho ( graph * __restrict__ G) 
{
    // determine the average preference to a candidate
    // for each (l x l) region (alike to "counties", etc)

    int N = G->NumNodes;
    int l = REGION_DIMENSION;
    int L = sqrt((double)(N)) + eps;
    int num_regions = (L/l) * (L/l);
    int i,node,region;
    double inf_1;

    //start at G->Rho_i = 0
    memset(G->Rho_i, 0, num_regions * sizeof(G->Rho_i[0]));

    for( region = 0; region< num_regions; region ++)
    {
	for (i = 0; i< l*l; ++i)
	{
            node = G->Regions[region * l*l + i];
	
            inf_1 = G->nodes[node].infection[0];
	    if(inf_1 > 0)
                G->Rho_i[region]+=(inf_1-1);
	}
        G->Rho_i[region] = G->Rho_i[region]/(l*l);
    }
}

double determine_regional_distance (int i, int j) 
{
    //////////////////////////////////////////////
    ///    Finding Distance Between Regions    ///
    //////////////////////////////////////////////
    /*
        We find the spatial distance between region centers
    	(Euclidian distance of centers of region i -> j)
    */

    int N = MIN_NUMBER_OF_NODES;
    int l = REGION_DIMENSION;
    int L = sqrt((double)(N)) + eps;

    double center_i_x, center_i_y,center_j_x, center_j_y,d_sq,delta_x,delta_y;

    //center coordinates of region i
    center_i_x = i % (L/l);
    center_i_y = (int)(i/(L/l));
	
    //center coordinates of region j
    center_j_x = j % (L/l);
    center_j_y = (int)(j/(L/l));
	
    delta_x = fabs(center_i_x - center_j_x);
    if(delta_x > 0.5 * (L/l))
	delta_x = (L/l) - delta_x;

    delta_y = fabs(center_i_y - center_j_y);
    if(delta_y > 0.5 * (L/l))
	delta_y = (L/l) - delta_y;

    //square distance
    d_sq = delta_x * delta_x + delta_y * delta_y;

    //absolute distance
    return sqrt(d_sq);
}

void cross_correlation_r ( 	graph * __restrict__ 	G,
				int 			r,
				double * __restrict__ 	cross_correlation,
				int * __restrict__	num_regions_r) 
{
    ////////////////////////////////////////////////////////////////
    /// Find Cross-correlation for All Regions Within r +/- dr ///
    ////////////////////////////////////////////////////////////////
    /*
        r = Euclidian distance of centers of region i -> j
       Bin = r +/- dr, dr = l/2
    */

    int N = G->NumNodes;
    int l = REGION_DIMENSION;
    int L = sqrt((double)(N)) + eps;
    int num_regions = (L/l) * (L/l);
    int i,j;
    // delta r
    double dr = 0.5, dist;
	
    //cross correlation within r
    (*cross_correlation) = 0;

    //count number of regions within r +/- dr
    (*num_regions_r) = 0;

    for (i = 0; i< num_regions; ++ i)
    {
        for(j = 0; j < num_regions; ++j)
        {
            // if d_ij within r +/- dr
	    dist = determine_regional_distance(i,j);

            if(dist < ((double)r + dr) && dist > ((double)r - dr))
            {
               //find cross correlation between all i and j s.t. d_ij = r +/- dr
               (*cross_correlation) += G->Rho_i[i] * G->Rho_i[j];
               (*num_regions_r)++;
            }
        }
    }
}

double spatial_correlation (
			graph * __restrict__	G,
			int 			r,
			double 			cross_correlation,
			int 			num_regions_r) 
{
    /////////////////////////////////////////////
    /// Finding The Spatial Correlation ///
    /////////////////////////////////////////////
    /*
        Definition of correlation:
    	    -> C (r) = (<rho_i rho_j | d_ij = r> - <rho>^2)/(<rho^2> - <rho>^2);
	    -> r = Euclidian distance of centers of region i -> j
	    -> Bin = r +/- dr, dr = l/2
    */
    int N = G->NumNodes;
    int l = REGION_DIMENSION;
    int L = sqrt((double)(N)) + eps;
    int num_regions = (L/l) * (L/l);
    int i;

    double rho_av_sq = 0, rho_sq_av = 0, sigma_rho;

    cross_correlation = cross_correlation/num_regions_r;

    for(i = 0; i< num_regions; i++)
    {
	rho_av_sq += G->Rho_i[i];
	rho_sq_av += G->Rho_i[i] * G->Rho_i[i];
    }
    rho_av_sq = (rho_av_sq * rho_av_sq)/(num_regions * num_regions);
    
    rho_sq_av = rho_sq_av / num_regions;
    sigma_rho = rho_sq_av - rho_av_sq;
    //printf("\nrho = %f\nsigma_rho = %f\n",sqrt(rho_av_sq),sigma_rho);

    double C_r = (cross_correlation - rho_av_sq)/sigma_rho;
    //printf("C(%i) = %f\n\n",r,C_r);
    return C_r;
}

void run_experiment (
                     graph * __restrict__ 	G,
		     double * __restrict__ 	FractionInfectedVersusTime,
		     double * __restrict__ 	UVersusTime,
                     long long 			start,
		     int * __restrict__ 	T_distribution,
                     double 			beta,
                     double 			S,
                     double 			delta_inf,
                     double 			gamma,
                     double 			temp,
		     double 			mu,
                     double 			mvm_p,
                     double 			tau,
                     double 			hung,
                     double 			stop_rate,
                     double 			binom,
		     double 			average_wait_tau,
		     double	 		wait_tau_alpha,
		     double			hung_p_innocent,
		     double			hung_p_guilty,
                     gsl_rng * 			r) 
{
    int Q = G->Q;
    int NUMBER_INFECTED [Q];
    double FRACTION_INFECTED [Q];

    int * __restrict__ T_dist;

    //if(RECORD_TAU)
        T_dist = &T_distribution[start];

    setup_experiment(	G,
			average_wait_tau,
			wait_tau_alpha,
			tau,
			beta,
			S,
			delta_inf,
			binom,
                        NUMBER_INFECTED,
                        FRACTION_INFECTED,
			T_dist,
                        r);


    infect_and_recover(G, NUMBER_INFECTED,beta, gamma,temp,mvm_p,tau,hung,stop_rate,mu,average_wait_tau,wait_tau_alpha,hung_p_innocent,hung_p_guilty,S,r);
}

double sum_double (	double * __restrict__ 	y, 
			int 			n) 
{
    int i;
    double s=0;

    for (i=0; i<n; ++i) s+=y[i];

    return s;
}

double sum_int (int * __restrict__ 	y, 
		double 			n) 
{
    int i;
    double s=0;

    for (i=0; i<n; ++i) s+=y[i];

    return s;
}


void store_number_of_infected_nodes ( 	double * __restrict__ 	FractionInfectedVersusTime,
					double * __restrict__ 	UVersusTime,
					double * __restrict__ 	FRACTION_INFECTED,
					graph * __restrict__ 	G,
					int * __restrict__	T_distribution,
					double 			tau,
					double 			t) 
{
    int i,delta_t,n,N=(G->NumNodes),Q = (G->Q),strain;
    if( t >= IGNORE_TIME && t < RUN_TIME && !RECORD_LONG_CONSENSUS_TIME)
    {

        UVersusTime[(int)(t) - IGNORE_TIME] = (1-2*FRACTION_INFECTED[0]) * (1-2*FRACTION_INFECTED[0])*(1-2*FRACTION_INFECTED[0]) * (1-2*FRACTION_INFECTED[0]);
        for(i = 0; i < Q; ++i)
	{
            FractionInfectedVersusTime[i * (RUN_TIME - IGNORE_TIME) + (int)(t) - IGNORE_TIME] = FRACTION_INFECTED[i];
	}
	if(RECORD_TAU)
	{
	    //for all nodes...
	    for(n = 0; n < N; ++n)
 	    {
	    	//if node is infected/opinionated...
	    	if(!is_susceptible(&(G->nodes[n])))
	    	{
	   	    delta_t = t - (G->nodes[n].time_infected);

	            //record which infection strain we have
		    strain = G->nodes[n].infection[0];

		    // if we are not in a frozen state, we record the time a node has been infected
		    if (delta_t < tau)
		    {
                        T_distribution[(strain - 1) * (RUN_TIME - IGNORE_TIME) * (TAU_MAX+1) + ((int)(t)-IGNORE_TIME)*(TAU_MAX+1) + (int)(delta_t)]++;
		    }
		
		    //if in the frozen state, we add it to the "t' >= tau" bin
            	    else if (delta_t > tau && tau >= 1.0)
		    {
                        T_distribution[(strain - 1) * (RUN_TIME - IGNORE_TIME) * (TAU_MAX+1) + ((int)(t)-IGNORE_TIME)*(TAU_MAX+1) + (int)(tau)]++;
		    }
		    // if tau < 1, then there is only 1 bin
		    else if (delta_t > tau && tau < 1.0)
		    {
                        T_distribution[(strain - 1) * (RUN_TIME - IGNORE_TIME) * (TAU_MAX+1) + ((int)(t)-IGNORE_TIME)*(TAU_MAX+1) + 1]++;
		    }
		}
	    }
	}
    }
}




void setup_experiment(graph * __restrict__ 	G,
		      double 			average_wait_tau,
		      double 			wait_tau_alpha,
		      double 			tau,
		      double 			beta,
		      double 			S,
		      double 			delta_inf,
		      double 			binom,
                      int * __restrict__	NUMBER_INFECTED,
                      double * __restrict__	FRACTION_INFECTED,
                      int * __restrict__	T_distribution,
                      gsl_rng * 		r) 
{
    int Q = G->Q, N = G->NumNodes;

    /////////////////////////// INITIALIZING STATE OF NODES ///////////////////////////
    initializing_state(G,NUMBER_INFECTED);

    //////////////////////////// INFECTING INITIAL NODES /////////////////////////////
    seed_infection(G,NUMBER_INFECTED,tau,beta,S,delta_inf,binom,r);

    //////////////////////////// FRACTION OF NODES INFECTED ///////////////////////////
    frac_of_total(FRACTION_INFECTED, NUMBER_INFECTED,N,Q);

    ///////////////////////// INITIALIZING TIME_INFECTED DISTRIBUTION /////////////////
    init_time_infected_dist(T_distribution,Q);

    //////////////////////// INITIALIZING WAITING_TIME FOR EACH LINK ////////////////
    init_link_wait_time(G,average_wait_tau,wait_tau_alpha,r);

}

void seed_infection ( 	graph * __restrict__ 	G,
			int * __restrict__ 	NUMBER_INFECTED,
			double 			tau,
			double 			beta,
			double 			S,
			double 			delta_inf,
			double 			binom,
			gsl_rng * 		r) 
{
    int i, i_ran, infection_strain, N=G->NumNodes,Q = G->Q, InfNum;

    ////////////////////////// DETERMINING NUMBER OF INFECTED NODES ///////////////////////////
    if(INFECT_CANDIDATES)
	InfNum = Q;
    else
    {
        if (Q < (int)(N*INIT_FRACTION_INFECTED))
	    InfNum = (int)(N*INIT_FRACTION_INFECTED);
	else
	{
	    InfNum = Q;
	    printf("NOTE: Q > N * Fraction Seeded = %i\nSimulation will run as desired but with a larger initial fraction infected.",(int)(N*INIT_FRACTION_INFECTED));
	}
    } 
	
    /////////////////////// CREATING RANDOMIZED LIST OF STRAINS ////////////////////////
    int random_strain_list [(long int)InfNum];
    create_random_strain_list(random_strain_list,InfNum,delta_inf,binom,Q,r);
       
    /////////////////////////////// INFECTING NODES ///////////////////////////////
    for(i = 0; i < InfNum; ++i)
    {
	//strain to infect node with
        infection_strain = random_strain_list[i];
	
	//finding susceptible node
	//printf("gsl_rng_uniform(): seed_infection()\n");
	do{
	    i_ran = gsl_rng_uniform_int(r,N);
        }while(!is_susceptible(&(G->nodes[i_ran])));
	
	//link-dynamics (neutral) case
	if(INFECT_NEUTRAL)
	{
	    int t0 = 0.0;
	    update_infection_states_infect(G,NUMBER_INFECTED,i_ran,infection_strain,tau,t0);
	}
	else
	{
	    //infect node
            infect(&(G->nodes[i_ran]),infection_strain,0.0);

	    //number infected ++ for strain
            NUMBER_INFECTED[infection_strain - 1] ++;
	    //added to list of infectious nodes
	    add_element_int(&(G->InfectiousNodes),&(G->Ni),i_ran);
	}

    } 

    ////////////////////// ADDING STUBBORNNESS ///////////////////////
    add_stubbornness(G,S,beta,r);

}

void initializing_state (	graph * __restrict__ 	G,
				int * __restrict__	NUMBER_INFECTED) 
{
    int Q = G->Q, N = G->NumNodes,i;

    ///////////////////////// INITIALIZING NUMBER OF INFECTIONS ////////////////////////
    G->Ni = 0;
    G->InfectiousNodes = malloc(0);
    memset(NUMBER_INFECTED, 0, Q * sizeof(NUMBER_INFECTED[0]));
    //printf("N1(0) = %i\n\n",NUMBER_INFECTED[0]);
    ///////////////////////      INITIALIZING CONSENSUS TIME     ///////////////////////
    G->ConsensusTime = 0.0;	   

    ///////////////////////      INITIALIZING NUMBER GUILTY      ///////////////////////
    G->NumGuilty = 0;	   

    ////////////////////// MAKING ALL NODES INITIAL SUSCEPTIBLE  ///////////////////////
    for(i = 0; i < N; i++)
    {
	make_susceptible(&(G->nodes[i]),Q);
    }

    /////////////////////// INITIALIZING CumSumBi, CumSumActiveEdges, NumActiveEdges, AND SumBij ////////////////////////
    if(INFECT_NEUTRAL)
    {
        G->SumBij = 0.0;
	// Set Cumulative Sums to 0's
	memset(G->CumSumBi, 0, N * sizeof(G->CumSumBi[0]));
        G->NumActiveEdges = 0;
	memset(G->CumSumActiveEdges, 0, N * sizeof(G->CumSumActiveEdges[0]));
    }

    /////////////////////// SETTING MEAN FRACTION INFECTED, MAGNETISM MOMENTS TO 0 ////////////////////////
    memset(G->MeanFractionInfected, 0, Q * sizeof(G->MeanFractionInfected[0]));
    G->m = 0.0;
    G->m2 = 0.0;
    G->m4 = 0.0;

    /////////////////////////////////////// RECORDING INITIAL VOTES ////////////////////////////////////////
    if (RECORD_VOTE_TIME)
    {
        int num_times = 5;//timesteps: 1.0,2.0,3.0,4.0,5.0
        G->v0 = malloc(sizeof(double) * num_times);
    }
}

void init_time_infected_dist (	int * __restrict__ 	T_distribution, 
				int 			Q) 
{
   // Set myArray to all 0's
    if(RECORD_TAU)
        memset(T_distribution, 0, Q * (RUN_TIME - IGNORE_TIME) * (TAU_MAX + 1) * sizeof(int));
}

void init_link_wait_time ( 	graph * __restrict__ 	G,
				double 			average_wait_tau,
				double 			wait_tau_alpha,
				gsl_rng * 		r) 
{
    if(AVERAGE_WAIT_TAU_MAX > eps)
    {
        int i,j,neighbor,wait_tau, N = G->NumNodes;
        for(i=0; i<N; ++i)
        {
            for(j=0; j<G->nodes[i].degree; ++j)
            {
                neighbor = G->nodes[i].neighbors[j];
                if(neighbor >= 0)
                {
		    if(AVERAGE_WAIT_TAU_MIN > eps)
		    {
	 	        resample_wait_tau(G,i,j,average_wait_tau,wait_tau_alpha,r);
                        wait_tau = G->nodes[i].wait_tau[j];

                        if(wait_tau > eps)
                        {
                            G->nodes[i].wait_time[j] = gsl_rng_uniform_int(r,wait_tau);

                            //make sure wait_time is between 1 and tau
                            G->nodes[i].wait_time[j] ++;
		        }
                    }
     	        }
	    }
        }
    }
}

void make_susceptible ( node * __restrict__ 	n,
			int 			Q) 
{
    // Set infection to SUSCEPTIBLE
    //memset(n->infection, SUSCEPTIBLE, (Q / STRAINS_PER_SITE + 1) * sizeof(n->infection[0]));
    n->infection[0] = 0;
    n->time_infected = 0.0;
    n->frequency_exposed = 0;
    n->time_last_checked = 0.0;
    n->current_stubbornness_p = 1.0;
}


void add_stubbornness ( graph * __restrict__ 	G,
			double 			S,
			double 			beta,
                        gsl_rng * 		r) 
{
    if(S > 0.0)
    {
        int i,
	    i_ran,
	    N = G->NumNodes,
	    // force the first strain picked to be "2"
	    // then the second strain picked will be "1"
	    strain = 1;
	int n=0,s_min=0;

    	double alpha = 2.5,
		s_av = 1/beta,
		//arbitrarily large value. no significant effect on speed
		s_max = 1000000000000000;

        double s = 0.0;
	//printf("\ns_av = %f \n",s_av);
	//printf("s_av = %f \n s_max = %i\n n = %i",s_av,s_max,n);
	if(SCALE_FREE_S)	
    	    find_distribution(alpha,N,s_av,&s_min,s_max,&n);

        for(i = 0; i < N; ++i)
        {
	    if(JURY_TRIAL_TIME)
	    {
                s = ran_ref_1D_walk(S,r);
	    }
     	    else if(RANDOM_S)
 	    {
                s = gsl_rng_uniform(r);
            }
     	    else if(SCALE_FREE_S)
 	    {
		if(i < n)
		    s = (double)s_min;
		else
         	    s = ran_scale_free(r,s_min,alpha,N);
            }
	    else if(TOP_HAT_S)
	    {
		
                s = gsl_rng_uniform(r)/(beta/2);
	    }
            G->nodes[i].s = s;
        }
        if(NUM_STUBBORN_NODES > 0)
        {
	    int nodes [NUM_STUBBORN_NODES];
	
	    // We have equal shares of stubborns
	    // ONLY if Q = 2
	    for(i = 0; i < NUM_STUBBORN_NODES; ++i)
	    {
                do{
                    i_ran = gsl_rng_uniform_int(r,N);
                }while(G->nodes[i_ran].infection[0] == strain
		       || node_already_picked(i_ran,nodes,i));
                G->nodes[i_ran].s = S + (i-1) * DELTA_S;
                strain = G->nodes[i_ran].infection[0];

	        if(PREFER_IC)
	        {
		    G->nodes[i_ran].initial_preference = strain;
	        }

	        nodes[i] = i_ran;
	    }
        }
    }
}

int node_already_picked (
			int 			i_ran,
			int * __restrict__	nodes,
			int 			i) 
{
    int j;
    for(j = 0; j<i; ++j)
    {
	if(nodes[j] == i_ran)
	    return 1;
    }
    return 0;
}

void make_flip ( 	node * __restrict__ 	n,
			int * __restrict__	NUMBER_INFECTED,
			int 			Q,
			double 			t,
			gsl_rng * 		r) 
{
    int new_spin, current_spin = n->infection[0];
    if(Q > 2)
    {
    	do{
            //pick random spin
            new_spin = gsl_rng_uniform_int(r,Q) + 1;
    	}while(new_spin == current_spin);
    }
    else
	new_spin = (current_spin % 2) + 1;
    //updating number infected
    NUMBER_INFECTED[current_spin - 1] --;
    NUMBER_INFECTED[new_spin - 1] ++;
    infect(n,new_spin,t);
}


void swap( 	int * __restrict__ 	list,
		int 			length,
		int 			pos,
		int 			new_val) 
{
    //record original value
    int old_val = list[pos];
    //find original position of value to swap
    int previous_pos = find(list,length,new_val);
  
    //finally swap positions
    list[previous_pos] = old_val;
    list[pos] = new_val;
}

int find( 	int * __restrict__ 	list,
		int 			length,
		int 			value) 
{
    /*Returns the index of when a value is first found in a given integer array*/

    int i;
    for(i=0; i<length; ++i)
    {
	if(list[i] == value)
	    return i;
    }
    printf("\n\nERROR: No value found: find() \n\n");
    return -1;
}

void create_random_strain_list( int * __restrict__ 	random_strain_list,
				int 			InfNum,
				double 			delta_inf,
				double 			binom,
				int 			Q,
				gsl_rng * 		r) 
{

    int i,di,num_guilty;
    if (binom > eps)
    {
	num_guilty = gsl_ran_binomial(r,binom,InfNum);
	if(VERBOSE)
  	    printf("Number Guilty: %i\n",num_guilty); 
    	for(i=0; i < num_guilty; ++i)
    	{
            random_strain_list[i]=2;
    	}
	
    	for(i=num_guilty; i < InfNum; ++i)
    	{
            random_strain_list[i]=1;
    	}
    }
    else if (REAL_DATA_INIT_CONDIT)
    {

	
	// make hist of vote data
	/*
        // OR6
	int hist [MIN_NUMBER_OF_NODES+1] = {49,12,0,0,0,49,94};
	// OR12
        int hist [MIN_NUMBER_OF_NODES+1] = {137,66,99,52,0,0,0,0,0,92,160,121,206};
        // CA6
        int hist [MIN_NUMBER_OF_NODES+1] = {28,1,1,0,0,2,19};
	*/

        // CA8_1
        int hist [MIN_NUMBER_OF_NODES+1] = {13,16,29,0,3,1,49,44,16};

        /*
	// CA8_2
        int hist [MIN_NUMBER_OF_NODES+1] = {6,12,15,1,1,2,35,32,17};
        // CA8_3   
        int hist [MIN_NUMBER_OF_NODES+1] = {2,3,7,0,0,0,9,8,3};
        // CA12_1
        int hist [MIN_NUMBER_OF_NODES+1] = {5,56,59,71,2,3,4,3,4,118,95,62,20};
        // CA12_2
        int hist [MIN_NUMBER_OF_NODES+1] = {13,56,77,97,7,8,12,9,4,132,122,93,26};
        // CA12_3
        int hist [MIN_NUMBER_OF_NODES+1] = {5,35,46,56,7,2,3,6,4,96,72,60,10};
        // CA12_4
        int hist [MIN_NUMBER_OF_NODES+1] = {0,10,11,15,3,1,2,1,0,31,20,16,1};
        // CA12_5
        int hist [MIN_NUMBER_OF_NODES+1] = {0,6,4,8,0,0,0,0,0,10,7,3};
        */
        //TEST
        //int hist[MIN_NUMBER_OF_NODES+1] = {2,0,0,0,0,0,0,0,0,0,0,0,100};
        double cdf[MIN_NUMBER_OF_NODES+1] = {};
	hist2cdf(hist,MIN_NUMBER_OF_NODES+1,cdf);
	int j;

        num_guilty = emp_rnd(cdf,InfNum,r);//CHECK
        if(VERBOSE)
            printf("Number Guilty: %i\n",num_guilty);
    	for(i=0; i < num_guilty; ++i)
    	{
            random_strain_list[i]=2;
    	}
	
    	for(i=num_guilty; i < InfNum; ++i)
    	{
            random_strain_list[i]=1;
    	}

    }
    else
    {
    	for(i=0; i < InfNum; ++i)
    	{
            random_strain_list[i]=i % Q + 1;
    	}
    }

    //we create a perturbation betweet two infections
    //order is {1,2,1,2,1,2,1,...}
    di = InfNum * delta_inf/2;
    for(i=0; i<di; ++i)
    {
	random_strain_list[2*i + 1] = 1;
    }
    
    // randomly shuffling infection strains
    gsl_ran_shuffle(r,random_strain_list,InfNum,sizeof(int));    
}

void hist2cdf ( 	int * __restrict__	hist,
		  	int 			length,
		  	double * __restrict__	cdf)
{
    int i,sum = 0;
    int * cumsum=malloc(length*sizeof(int));
    for (i = 0; i < length; ++i)
    {
	sum += hist[i]; // total number in histogram
	cumsum[i] = sum; // cumulative sum up until now
    }
    for (i = length - 1; i >= 0; --i)
    {
	cdf[i] = 1 - cumsum[i]/(double)sum; 
    }
    free(cumsum);
}

int emp_rnd (	    	double * __restrict__	cdf,
			double 			InfNum,
			gsl_rng * 		r)
{
    double rand = gsl_rng_uniform(r);
    int i,cdf_pos = -1;
    for(i = 0; i < InfNum; ++i){
	if(cdf[i] >= rand){
	     cdf_pos = i; 
	     break;
	}
    }
    if(cdf_pos == -1) cdf_pos = InfNum;
    return cdf_pos; 
} 

void frac_of_total ( double * __restrict__ 	FRACTION_INFECTED,
			int * __restrict__ 	NUMBER_INFECTED,
			int 			N,
			int 			Q) 
{

    int i;
    for(i=0; i < Q; ++i)
    {
	if(SAMPLE)
            FRACTION_INFECTED[i] = NUMBER_INFECTED[i]/(double)(N*SAMPLE_FRACTION);
	
	else
            FRACTION_INFECTED[i] = NUMBER_INFECTED[i]/(double)(N);
    }
}

void infect_and_recover (
		graph * __restrict__ 	G,
		int * __restrict__	NUMBER_INFECTED,
             	double 			beta,
	     	double 			gamma,
	     	double 			temp,
             	double 			mvm_p,
             	double 			tau,
             	double 			hung,
             	double 			stop_rate,
	     	double 			mu,
		double 			average_wait_tau,
		double 			wait_tau_alpha,
		double			hung_p_innocent,
		double			hung_p_guilty,
		double			S,
             	gsl_rng * 		r) 
{
    
    if (INFECT_NEUTRAL)
    {
	infect_and_recover_neutral(G,NUMBER_INFECTED,beta,gamma,temp,tau,r);
    }
    else
    {
	infect_and_recover_inout(G,NUMBER_INFECTED,beta,gamma,temp,mvm_p,tau,hung,stop_rate,mu,average_wait_tau,wait_tau_alpha,hung_p_innocent,hung_p_guilty,S,r);
    }
    //REMOVE OLD POINTER
    if(INFECT_NEUTRAL)
    {
    	int node;
    	for (node = 0; node < G->NumNodes; ++node)
            update_infection_states_recover(G,NUMBER_INFECTED,node,tau,RUN_TIME);
    }
    free(G->InfectiousNodes);

}

void infect_and_recover_inout (
		graph * __restrict__ 	G,
		int * __restrict__	NUMBER_INFECTED,
             	double 			beta,
	     	double 			gamma,
	     	double 			temp,
             	double 			mvm_p,
             	double 			tau,
             	double 			hung,
             	double 			stop_rate,
	     	double 			mu,
		double 			average_wait_tau,
		double 			wait_tau_alpha,
		double			hung_p_innocent,
		double			hung_p_guilty,
		double			S,
             	gsl_rng * 		r) 
{
    int Ni = G->Ni,N=G->NumNodes,Q = G->Q,surviving_infection,not_hung,
    //if(RECORD_EQUILIBRIUM_TIME && tau > eps)
    i,youngest_i=0,total_vote,vote;

    double rand,t0 = 0.0, delta_t,p,persuation_rate,rec_gamma = gamma,t_gamma = gamma,
    //if(RECORD_EQUILIBRIUM_TIME && tau > eps)
    node_t=0.0,youngest_t=tau,wait_t = 0,stop,tau_effective;
    if (JURY_TRIAL_TIME)
    {
	persuation_rate = beta/sqrt(S);
    }else{
	persuation_rate = beta;
    }
    double FRACTION_INFECTED [Q];
    long double t = 0.0, dt = 0.0;

    if (STUBBORN_GAMMA)
    {
     	t_gamma = 1.0;
        rec_gamma = 0.0;
    }
    if (RANDOM_WALK_ALPHA)
        stop_rate = 1/ran_ref_1D_walk(1/(stop_rate*stop_rate)*1.27324,r);//mean of stop is EXACTLY "stop"

    while(t < RUN_TIME)
    {
     	//define number of infected nodes
        Ni = G->Ni;
        //define time-step (either infecting in- or out-ward)
        if(Ni > 0)
        {
	
	    total_vote = 0.0;
	    stop = stop_rate * Q_0;//start out with a fixed fraction of	the total stopping rate     	    
	    for (i = 0; i < N; ++i)
            {
            	vote = G->nodes[i].infection[0] - 1;// vote = "0" or "1"
            	total_vote += vote;
	    }
            if(!REMOVE_Q_TIME_DEPENDENCE)
	    {
	    	// dynamics when not hung (these are always the dynamics if we remove hung conditions)	
	    	if ((total_vote/(double)N <= 1-VOTE_THRESHOLD || total_vote/(double)N >= VOTE_THRESHOLD) || (REMOVE_HUNG_CONDITIONS||REMOVE_STOP_HUNG_CONDITION))
		    stop += stop_rate * fabs(0.5 - total_vote/(double)(N));
	    }
	    //printf("Initial Stop: %f\n",stop);

            if(INFECT_INWARD)
                dt = 1.0/((t_gamma + temp + stop) * Ni + N);
            else if (INFECT_OUTWARD)
                dt = 1.0/((1 + t_gamma + temp) * Ni);

            if(t > IGNORE_TIME && !RECORD_LONG_CONSENSUS_TIME)
            {

             	if(t0 < eps)
                    t0 = t - dt;
		
                frac_of_total(FRACTION_INFECTED, NUMBER_INFECTED,N,Q);
                G->m += fabs(dt * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]));
                G->m2 += dt * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]);
                G->m4 += dt * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]);
                calc_mean_fraction_infected(G,FRACTION_INFECTED,dt);
            }

            rand = gsl_rng_uniform(r);
            //recovery
            if(rand < dt * rec_gamma * Ni)
            {
             	recover(G,NUMBER_INFECTED,t,tau,r);
            }
            else if(INFECT_INWARD && rand < dt * (stop + rec_gamma) * Ni)
            {
             	//stop(G,NUMBER_INFECTED,t,tau,r);
 	        t += dt;

	        if(VERBOSE)
		{
	            printf("Consensus at t = %lf\ngamma  = %lf\tbeta = %lf\ttau= %lf\n",(double)t,gamma,beta,tau);
		    printf("# Guilty = %f\n\n",NUMBER_INFECTED[0]/(double)(N));
		}
        	G->ConsensusTime = t;
		if (RECORD_DIST)
         	    G->NumGuilty = total_vote;//NUMBER_INFECTED[0];
		break;
            }
            //flip spin
            else if(rand < dt * (stop + rec_gamma + temp) * Ni)
            {
             	flip_spin(G,NUMBER_INFECTED,tau,t,r);
            }

            else //infect
            {
                if(INFECT_OUTWARD)
                {
                    infect_outward(G,NUMBER_INFECTED,gamma,beta,tau,mu,t,r);
                }
                else if(INFECT_INWARD)
                {
                    infect_inward(G,NUMBER_INFECTED,beta,mvm_p,tau,hung,stop_rate,mu,t,total_vote,r);
                }
	    }
	    if(RECORD_EQUILIBRIUM_TIME && tau > eps)
	    {
	    	//if all nodes stable
	    	if(t > wait_t)
	    	{
		    //NOTE: unexpected behavior if everyone isn't infected at the outset
		    node_t = G->nodes[youngest_i].time_infected;
		    if(t - node_t >= tau)
		    {
	            	for(i = 0; i < N; ++i)
	            	{
		    	    node_t = G->nodes[i].time_infected;
		    	    if(node_t < youngest_t)
		    	    {
		            	youngest_t = node_t;
		            	youngest_i = i;
		    	    }
	    	        }
	    
	    	    	if(t - youngest_t >= tau)
	    	    	{
			    if(VERBOSE)
			        printf("Equilibrium at t = %lf\n",(double)t);
	            	    G->ConsensusTime = t;
		    	    break;
	    	    	}
	    	    	else
	    	    	{
		            wait_t = youngest_t + tau;
	    	    	}
		    }
		    else
		    {
		    	youngest_t = node_t;
		    	wait_t = youngest_t + tau;
		    }
	    	}
	    }
	    else if(RECORD_LONG_CONSENSUS_TIME && !(MVM_P_MIN > eps))
            {
		if( one_infection_survives(NUMBER_INFECTED,Q,N,&surviving_infection)
		    ||(NUMBER_INFECTED[0]/(double)(N) >= VOTE_THRESHOLD && hung_p_guilty > 0.0) 
		    ||(NUMBER_INFECTED[0]/(double)(N) <= 1-VOTE_THRESHOLD && hung_p_innocent > 0.0))
		{
		    p = gsl_rng_uniform(r);
		    not_hung = 0;

		    if (NUMBER_INFECTED[0]/(double)(N) <=1-VOTE_THRESHOLD && p < persuation_rate * hung_p_innocent) not_hung = 1;
		    if (NUMBER_INFECTED[0]/(double)(N) >=VOTE_THRESHOLD && p < persuation_rate * hung_p_guilty)   not_hung = 1;
		    if (not_hung || one_infection_survives(NUMBER_INFECTED,Q,N,&surviving_infection))
		    {
	            	if(VERBOSE)
			{
	            	    printf("Consensus at t = %lf\ngamma  = %lf\tbeta = %lf\ttau= %lf\n",(double)t,gamma,beta,tau);
			    printf("# Guilty = %f\n\n",NUMBER_INFECTED[0]/(double)(N));
			}
        	        G->ConsensusTime = t;
			if (RECORD_DIST)
         	            G->NumGuilty = NUMBER_INFECTED[0];

		        break;
		    }
		}
            }
	}
	else if (!(MVM_P_MIN > eps))
	{
	    if(one_infection_survives(NUMBER_INFECTED,Q,N,&surviving_infection) && RECORD_LONG_CONSENSUS_TIME)
            {
                G->ConsensusTime = t;
            }

	    if(!RECORD_LONG_CONSENSUS_TIME)
	    {
	        t = RUN_TIME;

	        t0 = IGNORE_TIME;
                //frac_of_total(FRACTION_INFECTED, NUMBER_INFECTED,N,Q);
	        G->m  = 0.0;
 	        G->m2 = 0.0;
	        G->m4 = 0.0;
		int i;
		for (i = 0; i < Q; ++i) G->MeanFractionInfected[i] = 0.0;
	    }

	    break;
	}
        t += dt;
	tau_effective = tau;
        if (total_vote/(double)N > 1-VOTE_THRESHOLD && total_vote/(double)N < VOTE_THRESHOLD && (!REMOVE_HUNG_CONDITIONS && !REMOVE_STUBBORNNESS_HUNG_CONDITION ))
        {
            tau_effective /= hung; //hung is the fraction of tau
        }

	for (i = 0; i < N; ++i)
	{

            delta_t = t - G->nodes[i].time_last_checked;
            if (!REMOVE_STUBBORN_TIME_DEPENDENCE && tau_effective > 0.0 && (G->nodes[i].current_stubbornness_p) > 0.0)
            {
            	G->nodes[i].current_stubbornness_p += - delta_t/tau_effective;
            }
       	    else // this avoids NaN
           	G->nodes[i].current_stubbornness_p = 0.0;

       	    // update time last checked, as we update stubbornness
            G->nodes[i].time_last_checked = t;
        }	    		
	if(RECORD_VOTE_TIME)
	{
	    double t0 = 2.5,epsilon = 1/(2*((double)N));
	    for (t0=1.0;t0<6.0;t0++)
	    {
	    	if (fabs(t-t0)<epsilon)
	    	{
			G->v0[(int)t0-1] = total_vote/(double)N;
	    		//printf("%.10f\t%.10f\n",total_vote/(double)N,t);
		}
	    }
	}
    }

    //values are 0 by default
    if(t > IGNORE_TIME && !RECORD_LONG_CONSENSUS_TIME)
    {
        //normalizing m values
	delta_t = t - t0;
	//printf("delta_t = (%lf - %lf) = %lf\n",(double)t,(double)t0,delta_t);
        G->m  = (G->m )/delta_t;
        G->m2 = (G->m2)/delta_t;
        G->m4 = (G->m4)/delta_t;
        normalize_mean_fraction_infected(G,delta_t);
	if(VERBOSE)
	{
	    printf("rho1 = %lf\trho2 = %lf\n",FRACTION_INFECTED[0],FRACTION_INFECTED[1]);
	    printf("m = %lf\n",(double)G->m);
	}
    } 
}
void recover ( 		graph * __restrict__ 	G,
			int * __restrict__	NUMBER_INFECTED,
			double 			tau,
			double 			t,
			gsl_rng * 		r) 
{
    int node = -1, Q = G->Q,pos,strain;

    //random node ends up in the 0 state
    pick_infectious_node(G,&node,r);
    if(node >= 0)//infectious node exists
    {
	if(INFECT_NEUTRAL)
	    update_infection_states_recover(G,NUMBER_INFECTED,node,tau,t);
	else
	{
	    remove_element_int(&(G->InfectiousNodes),&(G->Ni),node,&pos);
	    strain = G->nodes[node].infection[0];
	    NUMBER_INFECTED[strain - 1]--;
	    make_susceptible(&(G->nodes[node]),Q);
	}
    }
}

void pick_infectious_node (
			graph * __restrict__	G,
			int * __restrict__	node,
			gsl_rng * 		r) 
{
    //pick a node that is infectious
    int pos;
    //printf("Ni = %i\n",G->Ni);
    //if all nodes recover, don't pick any node
    if(G->Ni == 0)
	(*node) = -1;
    else
    {
	pos = gsl_rng_uniform_int(r,G->Ni);
	(*node) = G->InfectiousNodes[pos];
    }
}
void update_infection_states_recover (
			graph * __restrict__	G,
			int * __restrict__	NUMBER_INFECTED,
			int 			node,
			double 			tau,
			double 			t) 
{
    /*
	Update states given change in the status of "node":
		- G->Ni
		- G->InfectiousNodes
		- G->nodes[i].num_active_links
		- G->nodes[i].active_neighbors
		- G->nodes[i].num_contrary_links
		- G->nodes[i].contrary_neighbors
		- G->nodes[i].Bj
    */
    int pos, Q = G->Q, N = G->NumNodes, strain = G->nodes[node].infection[0];

    /////////////////////// MAKE THE NODE SUSCEPTIBLE /////////////////////////
    make_susceptible(&(G->nodes[node]),Q);

    /////////////////////// REMOVE NODE FROM INFECTIOUS LIST /////////////////////////
    //update infectious nodes
    remove_element_int(&(G->InfectiousNodes),&(G->Ni),node,&pos);
    NUMBER_INFECTED[strain - 1] --;

    /////////////////////// 	    UPDATE EDGE STATES	      /////////////////////////
    /*
	Update:
		- node -> active_edges //remove active links to susceptible neighbors
		- active_edges <- neighbors //add active links from infected neighbors
		- node <--contrary edges-->  neighbors //remove!
    */
    //remove active links:
    remove_list_int(&(G->nodes[node].active_neighbors),&(G->nodes[node].num_active_edges));

    //remove contrary links
    remove_contrary_links(G,node);

    //make links from infected neighbor to current node  active
    add_active_links_to_neighbors(G,node);

    ///////////////////////     UPDATE CUMULATIVE SUMS     /////////////////////////
	
    //calculate cumulative sum for SumBi, and Active Edges
    CalcCumSumBi(G,tau,t);
    CalcCumSumActiveEdges(G);

    ///////////// 	UPDATE NUMBER OF ACTIVE EDGES, WEIGHTED SUM OF CUMULATIVE EDGES	 ////////////
    //calculate number of active edges
    G->NumActiveEdges = G->CumSumActiveEdges[N - 1];
    //calculate SumBij
    G->SumBij = G->CumSumBi[N - 1];
}
void update_infection_states_infect (
			graph * __restrict__	G,
			int * __restrict__	NUMBER_INFECTED,
			int 			node,
			int	 		strain,
			double 			tau,
			double 			t) 
{
    /*
	Update states given change in the status of "node":
		- G->Ni
		- G->InfectiousNodes
		- G->nodes[i].num_active_links
		- G->nodes[i].active_neighbors
		- G->nodes[i].num_contrary_links
		- G->nodes[i].contrary_neighbors
		- G->nodes[i].Bj
    */
    int N = G->NumNodes, neighbor;

    /////////////////////// MAKE THE NODE INFECTED /////////////////////////
    infect(&(G->nodes[node]),strain,t);

    /////////////////////// ADD NODE TO INFECTIOUS LIST /////////////////////////
    //update infectious nodes
    add_element_int(&(G->InfectiousNodes),&(G->Ni),node);    
    NUMBER_INFECTED[strain - 1] ++;

    if (INFECT_NEUTRAL)
    {
    /////////////////////// 	    UPDATE EDGE STATES	      /////////////////////////
    /*
	Update:
		- node -> active_edges
		- active_edges <- neighbors
		- node <--contrary edges-->  neighbors
    */
    // remove active edges from neighbors
    remove_active_links(G,node);
    //find contrary/active edges
    add_links(G,node,neighbor,strain);

    /////////////////////// 	UPDATE CUMULATIVE SUMS	 /////////////////////////
    //calculate cumulative sum for SumBi, and Active Edges
    CalcCumSumBi(G,tau,t);
    CalcCumSumActiveEdges(G);

    ///////////// 	UPDATE NUMBER OF ACTIVE EDGES, WEIGHTED SUM OF CUMULATIVE EDGES	 ////////////
    //calculate number of active edges
    G->NumActiveEdges = G->CumSumActiveEdges[N - 1];
    //calculate SumBik
    G->SumBij = G->CumSumBi[N - 1];
    }
}

void update_infection_states_contrary (
                        graph * __restrict__ 	G,
			int * __restrict__	NUMBER_INFECTED,
                        int 			node,
                        int 			strain,
                        double 			tau,
                        double 			t) 
{
    /*
       	Update states given change in the status of "node":
                - G->nodes[i].num_contrary_links
                - G->nodes[i].contrary_neighbors
                - G->nodes[i].Bj
    */
    int i, neighbor, N = G->NumNodes, old_strain = G->nodes[node].infection[0];


    //////////////	GIVE NODE CONTRARY OPINION	//////////////
    infect(&(G->nodes[node]),strain,t);

    //////////////	UPDATE NUMBER OF INFECTED NODES     //////////////
    NUMBER_INFECTED[old_strain- 1]--;
    NUMBER_INFECTED[strain - 1]++;


    /////////////////////// UPDATE EDGE STATES /////////////////////////
    /*
      	Update:
               	- node -> active_edges
                - active_edges <- neighbors
                - node <--contrary edges--> neighbors
    */
    //remove old contrary links
    remove_contrary_links(G,node);

    //add contrary links

    for (i = 0; i < G->nodes[node].degree; i++)
    {
        neighbor = G->nodes[node].neighbors[i];
	add_contrary_links(G,node,neighbor,strain);
    }
    //active links STAYS THE SAME

    /////////////////////// UPDATE CUMULATIVE SUMS /////////////////////////
    //calculate cumulative sum for SumBi, and Active Edges
    if(G->Ni < N)
    {
    	CalcCumSumActiveEdges(G);
    }
    if(tau > eps || G->Ni < N)
	CalcCumSumBi(G,tau,t);

    ///////////// UPDATE NUMBER OF ACTIVE EDGES, WEIGHTED SUM OF CUMULATIVE EDGES ////////////
    //calculate number of active edges
    G->NumActiveEdges = G->CumSumActiveEdges[N - 1];
    //calculate SumBij
    G->SumBij = G->CumSumBi[N - 1];

}


void add_links (	graph * __restrict__	G,
			int 			node,
			int 			neighbor,
			int 			strain) 
{
    int i;
    //add contrary links and active links
    int degree = G->nodes[node].degree;
    for (i = 0; i < degree; i++)
    {
	neighbor = G->nodes[node].neighbors[i];
	add_contrary_links(G,node,neighbor,strain);
        add_active_links(G,node,neighbor);
    }
}

void add_active_links (	graph * __restrict__	G,
			int 			node,
			int 			neighbor) 
{
    
    if(node >= 0 && neighbor >= 0)
    {
	//assume neighbor is susceptible
	if(is_susceptible(&(G->nodes[neighbor])))
	{
    	    add_element_int(&(G->nodes[node].active_neighbors),&(G->nodes[node].num_active_edges),neighbor);
	    //printf("num active neighbors: %i\n",G->nodes[node].num_active_edges);
	}
    }

}
void add_active_links_to_neighbors (
			graph * __restrict__	G,
			int			node)
{
    /*
	we turn ALL infected links into active links.
	We assume node has just recovered.
    */
    int i,neighbor;
    for(i = 0; i < G->nodes[node].degree; ++i)
    {
	neighbor = G->nodes[node].neighbors[i];
	if(neighbor >= 0)
	{
	    if(!is_susceptible(&(G->nodes[neighbor])))
	    {
	    	add_active_links(G,neighbor,node);
	    }
	}
    }
}



void add_contrary_links (
			graph * __restrict__	G,
			int 			node,
			int 			neighbor,
			int 			strain) 
{
    int neighbor_strain,list_size;
    double place_holder_val = 1.0;
    if(node >=0 && neighbor >= 0)
    {
	neighbor_strain = G->nodes[neighbor].infection[0];
        //contrary neighbor
        if(neighbor_strain != strain && !is_susceptible(&(G->nodes[neighbor])))
        {
	    //update node's contrary neighbors, Bj
            list_size = G->nodes[node].num_contrary_edges;
		
            add_element_int(&(G->nodes[node].contrary_neighbors),&(G->nodes[node].num_contrary_edges),neighbor);
      	    add_element_double(&(G->nodes[node].Bj),&list_size,place_holder_val);
		
	    //update node's contrary neighbors, Bj
            list_size = G->nodes[neighbor].num_contrary_edges;
		
            add_element_int(&(G->nodes[neighbor].contrary_neighbors),&(G->nodes[neighbor].num_contrary_edges),node);
            add_element_double(&(G->nodes[neighbor].Bj),&list_size,place_holder_val);
        }
    }
}

void remove_contrary_links (
			graph * __restrict__	G,
			int 			node) 
{
    int list_size, pos = 0, i, contrary_neighbor;

    //remove contrary links for neighbors
    for (i = 0; i < G->nodes[node].num_contrary_edges; i++)
    {
        contrary_neighbor = G->nodes[node].contrary_neighbors[i];

        //remove active links from neighbors
        list_size = G->nodes[contrary_neighbor].num_contrary_edges;

        remove_element_int( &(G->nodes[contrary_neighbor].contrary_neighbors),
                            &(G->nodes[contrary_neighbor].num_contrary_edges),
                            node,
                            &pos);

      	//remove Bij from neighbors
        remove_element_double(&(G->nodes[contrary_neighbor].Bj),
                              &list_size,
                              pos);
    }
    //remove Bj, and contrary edges for "node"
    list_size = G->nodes[node].num_contrary_edges;
    remove_list_double(&(G->nodes[node].Bj),&(G->nodes[node].num_contrary_edges));
    remove_list_int(&(G->nodes[node].contrary_neighbors),&list_size);
}


void remove_active_links (
			graph * __restrict__	G,
			int 			node) 
{
    int pos, i, neighbor;

    //remove active links for neighbors
    for (i = 0; i < G->nodes[node].degree; i++)
    {
        neighbor = G->nodes[node].neighbors[i];
	if(neighbor >=0)
	{
            remove_element_int( &(G->nodes[neighbor].active_neighbors),
                            	&(G->nodes[neighbor].num_active_edges),
                            	node,
                            	&pos);
	}
    }
}

void remove_list_int ( 	int * __restrict__ *list,
			int * __restrict__ list_size) 
{
    (*list_size) = 0;
    free(*list);
    (*list) = malloc(0);
}

void remove_list_double (
			double * __restrict__ 	*list,
			int * __restrict__ 	list_size) 
{
    
    (*list_size) = 0;
    free(*list);
    (*list) = malloc(0);
}

void add_element_int (	int * __restrict__	*list,
			int * __restrict__	list_size,
			int 			val) 
{
    //printf("THEN: %i\tval = %i\n",(*list_size),(*list));
    (*list_size) ++;
    (*list) = realloc(*list,sizeof(int) * (*list_size));
    //printf("NOW: %i\tval = %i\n",(*list_size),(*list));

    (*list)[(*list_size) - 1] = val;
}

void add_element_double (
			double * __restrict__	*list,
                        int * __restrict__ 	list_size,
                        double 			val) 
{
    (*list_size) ++;
    (*list) = realloc(*list,sizeof(double) * (*list_size));
    (*list)[(*list_size) - 1] = val;
}

void remove_element_int (
			int * __restrict__ 	*list,
			int * __restrict__ 	list_size,
			int 			value,
			int * __restrict__	pos) 
{
    int i;
    if ((*list_size) > 0)
    {
        (*pos) = find_pos_int(*list,(*list_size),value);
	if((*pos) >= 0)
	{
	    if((*pos) < (*list_size) - 1)
	    {
	        if((*pos) >= (*list_size) - 1)
    	            printf("%i\t%i\n",(*pos),(*list_size));

                for(i = (*pos); i < (*list_size) - 1; i++)
                {
                    (*list)[i] = (*list)[i+1];
		}
	    }

            (*list_size)--;
            (*list) = realloc(*list,sizeof(int) * (*list_size));
	}
    }
}

void remove_element_double (
			double * __restrict__ 	*list,
			int * __restrict__ 	list_size,
			int 			pos) 
{
    int i;
    //int pos = find_pos_double(*list,(*list_size),value);

    for(i = pos; i < (*list_size) - 1; i++)
    {
        (*list)[i] = (*list)[i+1];
    }
    (*list_size)--;
    (*list) = realloc(*list,sizeof(double) * (*list_size));
}

int find_pos_int (
		int * __restrict__	list,
		int 			list_size,
		int 			value) 
{
    int i;
    for (i = 0; i < list_size; i++)
    {
	if(list[i] == value)
	    return i;
    }
    //no value found
    /*if(VERBOSE)
        printf("\nWarning: No value found for %i: find_pos()",value);*/
    return -1;
}

int find_pos_double (
		double * __restrict__	list,
		int 			list_size,
		double 			value) 
{
    int i;
    for (i = 0; i < list_size; i++)
    {
	if(list[i] == value)
	    return i;
    }
    //no value found
    printf("\nWarning: No value found for %f: find_pos()",value);
    return -1.0;
}

void CalcCumSumBi ( 	graph * __restrict__ 	G,
			double 			tau,
			double 			t) 
{
    int i,j,neighbor,N = G->NumNodes;
    double Bi,Bij,delta_t;
    for (i = 0; i < N; ++i)
    {
	Bi = 0;
	for (j = 0; j < G->nodes[i].num_contrary_edges; ++j)
	{
	    neighbor = G->nodes[i].contrary_neighbors[j];
	    delta_t = t - G->nodes[neighbor].time_infected;
	    if(tau > 0)
	        Bij = 1.0 - delta_t / tau;
	    else //diffusive
		Bij = 1.0;
	    if(Bij <= 0) Bij = 0;
	    G->nodes[i].Bj[j] = Bij;
	    Bi += Bij;
	}

	if(i == 0)
	{
	    G->CumSumBi[0] = Bi;
	}
	else
	{
	    G->CumSumBi[i] = (G->CumSumBi[i-1]) + Bi;
	}
    }
}


void CalcCumSumActiveEdges(	graph * __restrict__ 	G) 
{
    int i, num_active_edges, N = G->NumNodes;
    for (i = 0; i < N; ++i)
    {
	num_active_edges = G->nodes[i].num_active_edges;
	if(i == 0)
	{
	    G->CumSumActiveEdges[0] = num_active_edges;
	}
	else
	{
	    G->CumSumActiveEdges[i] = (G->CumSumActiveEdges[i-1]) + num_active_edges;
	}
    }
}

void flip_spin ( 	graph * __restrict__ 	G,
			int * __restrict__	NUMBER_INFECTED,
			double 			tau,
			double 			t,
			gsl_rng * 		r) 
{
    //random node has it's spin flipped

    int node = -1, Q = G->Q, current_strain, new_strain;
    if(Q > 1)
    {
        pick_infectious_node(G,&node,r);
        if(node >= 0)//found infectious node
        {
	    if(INFECT_NEUTRAL)
	    {
	    	current_strain = G->nodes[node].infection[0];
	    	//find new strain NOT equal to the old one
    	    	if(Q > 2)
       	    	{
    		    do{
            	        //pick random spin
            	        new_strain = gsl_rng_uniform_int(r,Q) + 1;
    		    }while(new_strain == current_strain);
    	    	}
	    	else //Q==2, only 1 other choice
		    new_strain = (current_strain % 2) + 1;

            	update_infection_states_contrary(G,NUMBER_INFECTED,node,new_strain,tau,t);
	    }
	    else
            	make_flip(&(G->nodes[node]),NUMBER_INFECTED,Q,t,r);
        }
    }
}

void infect_and_recover_neutral (
			graph * __restrict__ 	G,
			int * __restrict__	NUMBER_INFECTED,
			double 			beta,
			double 			gamma,
			double 			temp,
			double 			tau,
			gsl_rng * 		r) 
{
    int Ni=0, Ea = 0, N = G->NumNodes, Q = G->Q, surviving_infection;
 
    double rand, t0 = 0.0, SumBij,delta_t;
    long double t = 0.0, dt = 0.0;
    double FRACTION_INFECTED [Q];

    /*
	Based on Baguna, Castellano, & Pastor-Satorras (2013):
  	        
		- dt = 1/(gamma * Ni + beta * (Ea + SumBij))
		    - Ni = Number of infected nodes
		    - Ea = Number of active edges (infected - susceptible link)
		    - SumBij = Weighted sum of contrary edges (i_1 - i_2 link)
		- Let rand be in {0 - 1}
 
		    - Case 1: rand < dt * (gamma * Ni):
			 -pick random node (see function: pick_infectious_node)
			 -recover
			 -update number of infected nodes, active links, etc
		    - Case 1: rand < dt * ((gamma + temp) * Ni):
			 -pick random node (see function: pick_infectious_node)
			 -flip
			 -update number of infected nodes, active links, etc
    		    - Case 2: rand < dt * ((gamma + temp) * Ni + beta * Ea)
			 -pick random active link (see function: pick_active_edge)
			 -infect susceptible node
			 -update number of infected nodes, active links, etc

		    - Case 3 (default):
			 -pick contrary link (see function: pick_contrary_edge)
			 -change infection type
			 -update contrary link info (Bi, SumBij, contrary edges,etc) ONLY
     */
    /*printf("N1(dt) = %i\n",NUMBER_INFECTED[0]);
    printf("N2(dt) = %i\n",NUMBER_INFECTED[1]);
    printf("Ni(dt) = %i\n",G->Ni);
    printf("NumContraryEdges[0](dt) = %i\n",G->nodes[0].num_contrary_edges);
    printf("NumActiveEdges[0](dt) = %i\n",G->nodes[0].num_active_edges);
    printf("NumActiveEdges(dt) = %i\n",G->NumActiveEdges);*/
    
    while (t < RUN_TIME)
    {
	
	Ni = G->Ni;
	Ea = G->NumActiveEdges;
	SumBij = G->SumBij;
	if((gamma + temp) * Ni + beta * Ea + SumBij > eps)
	{
	    dt = 1/((gamma + temp) * Ni + beta * (Ea + SumBij));

	    t += dt;
	
	    if(t > IGNORE_TIME && !RECORD_LONG_CONSENSUS_TIME)
	    {
	        if(t0 < eps)
	            t0 = t - dt;

                frac_of_total(FRACTION_INFECTED, NUMBER_INFECTED,N,Q);

	        G->m += fabs(dt * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]));
 	        G->m2 += dt * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]);
	        G->m4 += dt * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]);
	        calc_mean_fraction_infected(G,FRACTION_INFECTED,dt);
	    }
	
	    //find fraction of nodes infected

	    rand = gsl_rng_uniform(r);
	
	    if(rand < dt * gamma * Ni)
	    {
	        recover(G,NUMBER_INFECTED,t,tau,r);
	    }
	    else if(rand < dt * (temp + gamma) * Ni)
	    {
	        flip_spin(G,NUMBER_INFECTED,tau,t,r);
	    }
	    else if(rand < dt * ((gamma + temp) * Ni + beta * Ea))
	    {
	        infect_active_edge(G,NUMBER_INFECTED,tau,t,r);
	    }
	    else if(SumBij > eps && Q > 1)//infect contrary neighbors
	    {
	        infect_contrary_edge(G,NUMBER_INFECTED,tau,t,r);
	    }

	    if(Q > 1 && one_infection_survives(NUMBER_INFECTED,Q,N,&surviving_infection) && RECORD_LONG_CONSENSUS_TIME)
            {
                G->ConsensusTime = t;
	        break;
            }
	}
	else
	{
	    if(one_infection_survives(NUMBER_INFECTED,Q,N,&surviving_infection) && RECORD_LONG_CONSENSUS_TIME)
            {
                G->ConsensusTime = t;
            }

	    if(!RECORD_LONG_CONSENSUS_TIME)
	    {
	        t = RUN_TIME;

	        t0 = IGNORE_TIME;
		delta_t = t - t0;
                frac_of_total(FRACTION_INFECTED, NUMBER_INFECTED,N,Q);
	        G->m = fabs(delta_t * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]));
 	        G->m2 = delta_t * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]);
	        G->m4 = delta_t * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]) * (FRACTION_INFECTED[1] - FRACTION_INFECTED[0]);
		int i;
		for (i = 0; i < Q; ++i) G->MeanFractionInfected[i] = delta_t * FRACTION_INFECTED[i];
		
	    }

	    break;
	}
    }
    printf("|%lf - %lf| = %lf\n",FRACTION_INFECTED[0],FRACTION_INFECTED[1],fabs(FRACTION_INFECTED[0]-FRACTION_INFECTED[1]));
    printf("|%i - %i| = %i\n",NUMBER_INFECTED[0],NUMBER_INFECTED[1],(NUMBER_INFECTED[0]-NUMBER_INFECTED[1]));

    if(t > IGNORE_TIME && !RECORD_LONG_CONSENSUS_TIME)
    {
    	//normalizing m values
    	delta_t = t - t0;
    	G->m = (G->m )/delta_t;
    	G->m2 = (G->m2)/delta_t;
    	G->m4 = (G->m4)/delta_t;
    	normalize_mean_fraction_infected(G,delta_t);
	int i,count = 0;
	for(i = 0; i < N; ++i)
	{
	    if(G->nodes[i].infection[0] == 1)
		count++;
	}
    }
}

void infect_active_edge (
			graph * __restrict__	G,
			int * __restrict__	NUMBER_INFECTED,
			double 			tau,
			double 			t,
			gsl_rng * 		r) 
{
    /*
	Find active edges from list:
	   - pick random value from 0 -> Ea - 1, x_i
	   - find sup {node, n, with ea > 0, and sum_active_edges[n] < x_i}
	   - neighbor = G->nodes[sup{}].neighbors[x_i - sup{}]
	   - infect, update states
    */

    int infectious_node,susceptible_node,strain;
    pick_active_edge(G,&infectious_node,&susceptible_node,r);
    strain = G->nodes[infectious_node].infection[0];
    update_infection_states_infect(G,NUMBER_INFECTED,susceptible_node,strain,tau,t);
}

void infect_contrary_edge (
	                graph * __restrict__ 	G,
			int * __restrict__	NUMBER_INFECTED,
			double 			tau,
			double 			t,
                        gsl_rng * 		r) 
{
    /*
      	Find active edges from list:
           - pick random value from 0 -> Ea - 1, x_i
           - find sup {node, n, with ea > 0, and sum_active_edges[n] < x_i}
           - neighbor = G->nodes[sup{}].neighbors[x_i - sup{}]
           - infect, update states
    */
    int contrary_node1,contrary_node2,strain;
    pick_contrary_edge(G,&contrary_node1,&contrary_node2,r);
    strain = G->nodes[contrary_node1].infection[0];
    update_infection_states_contrary(G,NUMBER_INFECTED,contrary_node2,strain,tau,t);
}


void pick_active_edge (
			graph * __restrict__	G,
			int * __restrict__	infectious_node,
			int * __restrict__	susceptible_node,
			gsl_rng * 		r) 
{

    int x_i = gsl_rng_uniform_int (r,G->NumActiveEdges) + 1;
    int n1 = 0, n2 = 1,Enp,En2 = G->CumSumActiveEdges[n2], np, pos, N = G->NumNodes;
    //printf("Ea = %i\tx_i = %i\n",G->NumActiveEdges,x_i);
    /*
	bit search for node n, s.t. CumSumActiveEdges[n] > x_i
    */

    //double up
    while (En2 < x_i)
    {
	n1 = n2;
        n2 = n1 * 2;
	if(n2 >= N)
	{
	    n2 = N-1;
	    break;
	}
	En2 = G->CumSumActiveEdges[n2];
    }

    //binary search: find nodes n1 != n2, such that:
    //	- CumSum[n2] >= x_i
    //	- CumSum[n1] < x_i
    // - n2 - n1 = 1 (i.e. n2 = n1 + 1)
    while (n2-n1 > 1)
    {
	//printf("Ea HALVE, n1 = %i, n2 = %i\n",n1,n2);

        np = n1 + ((n2-n1)/2);
 	Enp = G->CumSumActiveEdges[np];
        if (Enp >= x_i) n2=np;
        else n1 = np;
    }

    if(n1 < 0){ printf("ERROR:\npick_active_edge(): n1 < 0\n"); exit(1);}
    if(n1==n2 && n1 > 0){ printf("ERROR:\npick_active_edge(): n1 == n2\n"); exit(1);}

    if(x_i <= G->CumSumActiveEdges[n1])
    {
        (*infectious_node) = n1;
	pos = x_i - 1;
        (*susceptible_node) = G->nodes[n1].active_neighbors[pos];
    }
    else
    {
        (*infectious_node) = n2;
        pos = x_i - G->CumSumActiveEdges[n1] - 1;
	
        if(pos < 0 || pos > G->nodes[n2].num_active_edges - 1)
        {
            printf("pos = %i\tlist_size = %i\nn1 = %i\tn2 = %i\tx_i = %i\n",pos,G->nodes[n2].num_active_edges,n1,n2,x_i);
            printf("CumSum[n1] = %i\nCumSum[n2] = %i\n",G->CumSumActiveEdges[n1],G->CumSumActiveEdges[n2]);
	}
        (*susceptible_node) = G->nodes[n2].active_neighbors[pos];
    }

}

void pick_contrary_edge (
	                graph * __restrict__ 	G,
                        int * __restrict__	contrary_node1,
                        int * __restrict__	contrary_node2,
                        gsl_rng * 		r) 
{
    double x_i = gsl_rng_uniform (r) * (G->SumBij),value,/*value1,value2,*/sumBj = 0.0;
    int n1 = 0,n2 = 1,np,pos=0,j,N = G->NumNodes;
    //double up
    if(x_i < eps)
    {
	printf("ERROR: x_i == 0\nWe are picking a contrary edge even though dynamics are frozen");
	exit(1);
    }
    while (G->CumSumBi[n2] < x_i)
    {
	//printf("Bj DOUBLE, n1 = %i, n2 = %i",n1,n2);

        n1 = n2;
        n2 = n1 * 2;
	if(n2 >= N)
	{
	    n2 = N-1;
	    break;
	}
    }
    //binary search
    while (n2-n1>1)
    {
	//printf("Bj HALVE, n1 = %i, n2 = %i",n1,n2);
        np = n1 + (n2-n1)/2;
        if (G->CumSumBi[np] >= x_i) n2=np;
        else n1 = np;
    }
    if(n1 < 0){ printf("ERROR:\npick_contrary_edge(): n1 < 0\n"); exit(1);}
    if(n1==n2&&n1 > 0){ printf("ERROR:\npick_contrary_edge(): n1 == n2\n"); exit(1);}

    if(x_i <= G->CumSumBi[0])
    {
        (*contrary_node1) = 0;
        value = x_i;
	//printf("SumBij: %f\n",sumBj);
        for (j = 0; j < G->nodes[0].num_contrary_edges; ++j)
        {
	    sumBj += G->nodes[0].Bj[j];
	    if(sumBj >= value)
	    {
	        pos = j;
	        break;
	    }
        }
        (*contrary_node2) = G->nodes[0].contrary_neighbors[pos];
    }
    else
    {
        (*contrary_node1) = n2;
	
        value = x_i - G->CumSumBi[n1];
        for (j = 0; j < G->nodes[n2].num_contrary_edges; ++j)
        {
	    sumBj += G->nodes[n2].Bj[j];
	    if(sumBj >= value)
	    {
	        pos = j;
	        break;
	    }
        }
        (*contrary_node2) = G->nodes[n2].contrary_neighbors[pos];
    }
}

void calc_mean_fraction_infected (
			graph * __restrict__	G,
			double * __restrict__	FRACTION_INFECTED,
			double 			dt) 
{
    int i, Q = G-> Q;
    for(i = 0; i < Q; ++ i)
    {
	G->MeanFractionInfected[i] += dt * FRACTION_INFECTED[i];
    }
}

void normalize_mean_fraction_infected (
			graph * __restrict__	G,
			double 			delta_t) 
{
    int i, Q = G-> Q;
    for(i = 0; i < Q; ++ i)
    {
	G->MeanFractionInfected[i] = G->MeanFractionInfected[i]/delta_t;
    }
}


int one_infection_survives (
			int * __restrict__	NUMBER_INFECTED,
			int 			Q,
			int 			N,
			int * __restrict__	surviving_infection) 
{
    int i, num_surviving_infections = Q;
    (*surviving_infection) = -1;

    for(i = 0; i < Q; ++i)
    {
	//if infection survived
	if(NUMBER_INFECTED[i] == 0)
	{
	    num_surviving_infections --;
	}   
	else
	    (*surviving_infection) = i;
    }
    //printf("%i\n",num_surviving_infections);
    //if there is only one surviving infection
    if(num_surviving_infections == 1)
    {
	return 1;
    }
    //multiple infections survive
    else
	return 0;
}

void infect_outward (   graph * __restrict__    G,
                        int * __restrict__	NUMBER_INFECTED,
                        double                  gamma,
                        double                  beta,
                        double                  tau,
                        double                  mu,
                        double                  t,
                        gsl_rng *               r)
{

    int node,degree,j,j_ran,neighbor,wait_tau=0,wait_time=0,Q = G->Q, N = G->NumNodes;
    int * __restrict__ random_order_j;
    double beta_i;

    pick_infectious_node(G,&node,r);
    double rec_prob = 0.0;
    if (STUBBORN_GAMMA)
    {
     	double delta_t = t - G->nodes[node].time_infected;
        if (delta_t < tau)
            rec_prob = gamma * (1 - delta_t/tau);
    }

    if(node >= 0)//infectious node exists
    {
	// Choose to infect or recover with equal probability
	
     	double p = gsl_rng_uniform(r);
	if (p < 0.5 && STUBBORN_GAMMA)
	{
    	    p = gsl_rng_uniform(r);

            int pos,strain;
            if (p < rec_prob)
            {
                remove_element_int(&(G->InfectiousNodes),&(G->Ni),node,&pos);
                strain = G->nodes[node].infection[0];
                NUMBER_INFECTED[strain - 1]--;
                make_susceptible(&(G->nodes[node]),Q);
            }
	}
        else
	{
            degree = G->nodes[node].degree;
            if (degree > 0)
            {

	        if(!VOTER_MODEL)//beta > eps)
        	{

	            mutate(&(G->nodes[node]),mu,Q,t,r);

	            if(degree != N)
        	    {
            		random_order_j = malloc(sizeof(int)*degree);

	            	randomize(random_order_j,degree,r);
		    }
           	    //change B_2 -> B_2 + dB_2

  	            beta_i = beta + (G->nodes[node].infection[0]-1) * DELTA_BETA;

	            for(j=0; j<degree;++j)
        	    { 
               		if(degree != N)
                    	    j_ran = random_order_j[j];
		    	    //look at its neighbors in a random order

 	               if(K_REGULAR_ANNEALED || degree == N)
			    neighbor = gsl_rng_uniform_int(r,N);
			else
            	    	    neighbor=G->nodes[node].neighbors[j_ran];

		    	if(AVERAGE_WAIT_TAU_MIN > eps)
		    	{
                	    wait_tau=G->nodes[node].wait_tau[j_ran];
                    	    wait_time=G->nodes[node].wait_time[j_ran];
		    	}
		
            		//make sure neighbor is not a self-link or multiple edge
	            	if(neighbor >= 0 && (wait_time==wait_tau||wait_tau < eps))
        	    	{
                	    //if neighbor is susceptible
                    	    if(is_susceptible(&(G->nodes[neighbor])))
                    	    {
	                    	infect_susceptible(G,NUMBER_INFECTED,node,neighbor,beta_i,t,r);
        	            }

	                    //else neighbor is infected with different infection
        	            else
	                    {
        	                become_turncoat(G,NUMBER_INFECTED,node,neighbor,beta_i,tau,t,r);
                	    }
  	    		}
            	    }
	    	    if(degree !=N)
                	free(random_order_j);
	    	}
    		// Invasion Process dynamics
	    	else
    		{
		    //pick random neighbor
		    j_ran = gsl_rng_uniform_int(r,degree);

		    if(degree != N)
        	    	neighbor=G->nodes[node].neighbors[j_ran];
		    else neighbor = j_ran;

	            if(AVERAGE_WAIT_TAU_MIN > eps)
        	    {
	            	wait_tau=G->nodes[node].wait_tau[j_ran];
        	    	wait_time=G->nodes[node].wait_time[j_ran];
            	    }

	            //make sure neighbor is not a self-link or multiple edge
        	    if(neighbor >= 0 && (wait_time==wait_tau||wait_tau < eps))
	            {
		    	//infect with probability beta
	    		become_turncoat(G,NUMBER_INFECTED,node,neighbor,beta,tau,t,r);
		    }
	        }
            }
        }
    }
}

void infect_inward ( 	graph * 		G,
			int * __restrict__	NUMBER_INFECTED,
                        double 			beta,
                        double 			mvm_p,
                        double 			tau,
                        double 			hung,
                        double 			stop_rate,
			double 			mu,
			double 			t,
			int                     total_vote,
			gsl_rng * 		r) 
{
    int j,j_ran,node,degree,neighbor,wait_tau=0,wait_time=0,N=G->NumNodes,Q = G->Q, innocent = 0, guilty = 1;
    double beta_j,t0=0.0;
    
    node = gsl_rng_uniform_int(r,N);
    degree=G->nodes[node].degree;
    
    mutate(&(G->nodes[node]),mu,Q,t,r);
    // record if votes are early (before 10th timestep = 60 minutes)
    /*if(RECORD_VOTE_TIME && (t <10.0))
    {
	for(t0=1.0; t0<5.0;t0+=1.0)
	{
	    // record votes at t0 = 1.0,2.0,3.0,...
	    if (fabs(t-t0)<epsilon)
        	printf("%.10f\t%.10f\n",total_vote/(double)N,t);
	}
    }*/

    // MVM dynamics
    if (mvm_p > eps)
    {
	//WARNING: WE ONLY ASSUME COMPLETE GRAPH. CHECK "NEIGHBORS" AND CHANGE APPROPRIATELY

	// find vote of neighbors, assuming none are recovered
	// default: votes are "1" or "2",
	// we convert votes to "0" or "1" and add numbers up
	int vote, total_neighbor_vote;//, innocent_neighbor,guilty_neighbor;
	double tau_effective,dynamics_p,rand;



	dynamics_p = G->nodes[node].current_stubbornness_p;
	// if no dynamics, then stop
	if (dynamics_p <= 0.0) return;

	rand = gsl_rng_uniform(r);
	if (rand < dynamics_p)
	{ 
	    vote = G->nodes[node].infection[0] - 1;
	    total_neighbor_vote = total_vote - vote;
	    rand = gsl_rng_uniform(r);

	    if (total_neighbor_vote < 0.5*degree)
	    {
		//vote innocent with probability p (majority opinion)
		if (rand < mvm_p)
                    become_turncoat(G,NUMBER_INFECTED,-3+innocent,node,beta,0.0,t,r);
		else // else vote guilty (minority opinion)
                    become_turncoat(G,NUMBER_INFECTED,-3+guilty,node,beta,0.0,t,r);
		//return;
	    }
	    else if (total_neighbor_vote > 0.5*degree){ 
		//vote guilty with probability p (majority opinion)
		if (rand < mvm_p)
                    become_turncoat(G,NUMBER_INFECTED,-3+guilty,node,beta,0.0,t,r);
		else // else vote innocent (minority opinion)
                    become_turncoat(G,NUMBER_INFECTED,-3+innocent,node,beta,0.0,t,r);
		//return;
	    }
	    else{// choose either with equal probability if vote is split
		if (rand < 0.5)
                    become_turncoat(G,NUMBER_INFECTED,-3 + innocent,node,beta,0.0,t,r);
		else
                    become_turncoat(G,NUMBER_INFECTED,-3 + guilty,node,beta,0.0,t,r);
		//return;
	    }

        }

    }

    else if(!VOTER_MODEL)//beta > eps)
    {
        int * __restrict__ random_order_j = malloc(sizeof(int)*degree);
        randomize(random_order_j,degree,r);

	printf("WRONG \n");

        for(j=0; j<degree;++j)
        { //look at its neighbors in a random order
            j_ran = random_order_j[j];

	    if(K_REGULAR_ANNEALED)
	        neighbor = gsl_rng_uniform_int(r,N);

	    else
            	neighbor=G->nodes[node].neighbors[j_ran];
	    
	    beta_j = beta + (G->nodes[neighbor].infection[0]-1) * DELTA_BETA;

	    if(AVERAGE_WAIT_TAU_MIN > eps)
	    {
                wait_tau=G->nodes[node].wait_tau[j_ran];
                wait_time=G->nodes[node].wait_time[j_ran];
	    }

            //make sure neighbor is not a self-link or multiple edge
            if(neighbor >= 0 && (wait_time==wait_tau||wait_tau < eps))
            {
                //if neighbor is susceptible
                if(is_susceptible(&(G->nodes[node])) && !is_susceptible(&(G->nodes[neighbor])))
                {
             	    infect_susceptible(G,NUMBER_INFECTED,neighbor,node,beta_j,t,r);
                }
 	        //else neighbor is infected
                else if (!is_susceptible(&(G->nodes[neighbor])))
                {
                    //if neighbor is infected with a distinct infection
                   become_turncoat(G,NUMBER_INFECTED,neighbor,node,beta_j,tau,t,r);
	        }
	    }
        }
        free(random_order_j);
    }
    // Voter Model dynamics
    else
    {
     //   int * __restrict__ random_order_j = malloc(sizeof(int)*degree);

	printf("WRONG \n");
	//pick a random neighbor
        j_ran = gsl_rng_uniform_int(r,degree);
        neighbor=G->nodes[node].neighbors[j_ran];

        if(AVERAGE_WAIT_TAU_MIN > eps)
        {
            wait_tau=G->nodes[node].wait_tau[j_ran];
            wait_time=G->nodes[node].wait_time[j_ran];
        }

	//make sure neighbor is not a self-link or multiple edge
        if(neighbor >= 0 && (wait_time==wait_tau||wait_tau < eps))
        {
            become_turncoat(G,NUMBER_INFECTED,neighbor,node,beta,tau,t,r);
        }
    }
}

void infect_susceptible (
			graph * __restrict__	G,
			int * __restrict__	NUMBER_INFECTED,
			int 			n,
                        int 			neighbor,
                        double 			beta,
                        double 			t,
                        gsl_rng * 		r) 
{
    int new_strain;
    //infect with probability I
    double rand = gsl_rng_uniform (r);
    if(rand<beta)
    {
	new_strain = G->nodes[n].infection[0];

        //copy_infection(G,neighbor,n,SUSCEPTIBLE,Q,t);
	infect(&(G->nodes[neighbor]),new_strain,t);
	add_element_int(&(G->InfectiousNodes),&(G->Ni),neighbor);
	NUMBER_INFECTED[new_strain - 1]++;
    }
    //record the number of exposures for unopinionated individuals
    else
        G->nodes[neighbor].frequency_exposed ++;
    
}

void become_turncoat (	graph * __restrict__	G,
			int * __restrict__	NUMBER_INFECTED,
			int 			n,
			int 			neighbor,
			double 			beta,
			double 			tau,
			double	 		t,
			gsl_rng * 		r) 
{
    int old_strain,new_strain;
    double rand;
    if (n < -1)
    {
	old_strain = G->nodes[neighbor].infection[0];
	new_strain = n+4;// -3 ->+1 (innocent), -2 -> +2 (guilty)
	NUMBER_INFECTED[old_strain - 1]--;
	NUMBER_INFECTED[new_strain - 1]++;
	infect(&(G->nodes[neighbor]),new_strain,t);

    }
    //check if both nodes have an opinion, and those opinions are distinct
    else if(distinct_infection(&(G->nodes[n]),&(G->nodes[neighbor])))
    {
        double 	turncoat_p = 0.0, 
		delta_t = t - G->nodes[neighbor].time_infected;

	find_turncoat_p(G,&turncoat_p,n,neighbor,beta,tau,delta_t);

	if(turncoat_p > eps)
	{
	    //printf("gsl_rng_uniform(): become_turncoat()\n");

            if (turncoat_p < 1.0) rand = gsl_rng_uniform (r);
            else rand = 0.0;//avoid using rng if dynamics deterministic

            if(rand<turncoat_p)
            {
		old_strain = G->nodes[neighbor].infection[0];
		new_strain = G->nodes[n].infection[0];
		NUMBER_INFECTED[old_strain - 1]--;
		NUMBER_INFECTED[new_strain - 1]++;
		infect(&(G->nodes[neighbor]),new_strain,t);
            }
	    //Record the number of exposures for infected individuals
	    else
	        G->nodes[neighbor].frequency_exposed ++;
	}
    }
}

void find_turncoat_p (	graph * __restrict__	G,
			double * __restrict__	turncoat_p,
			int 			n,
			int 			neighbor,
			double 			beta,
			double 			tau,
			double 			delta_t) 
{
    if (delta_t < tau && tau > 0.0)
    {
        (*turncoat_p) = beta * (1 - delta_t / tau);
	//with this function we freeze the dynamics if theta > tau,
	// else dynamics are VM-like
	if(HEAVISIDE_FREEZE)
	    (*turncoat_p) = beta;
    }	
    else if(tau == 0)
    {
        (*turncoat_p) = beta;
    }
	
    if(RANDOM_S || NUM_STUBBORN_NODES > 0 || JURY_TRIAL_TIME)
    {
	if(RANDOM_S || JURY_TRIAL_TIME)
	{
	    (*turncoat_p) = beta / G->nodes[neighbor].s;
	}
	else
	{
	    // stubborn nodes are more difficult to infect
	    // by a factor (1 - s)
	    (*turncoat_p) = beta * (1 - G->nodes[neighbor].s);

	    // If node prefers initial condition (is biased)
	    if(PREFER_IC
		    && G->nodes[n].infection[0]==G->nodes[neighbor].initial_preference)
	    {
	        (*turncoat_p) = beta * (1 + G->nodes[neighbor].s);
	    }
	}
    }
}


void mutate ( 	node * __restrict__	n,
		double 			mu,
		int 			Q,
		double 			t,
		gsl_rng * 		r) 
{
    if(mu > eps)
    {
        int strain,mutation;
        double rand;

        rand = gsl_rng_uniform(r);

        if (rand < mu && !is_susceptible(n))
        {
	    strain = n->infection[0];//find_strain(n,Q);

	    //clear previous infection
	    make_susceptible(n,Q);

	    if(strain==1)
	    {
	        infect(n,strain + 1,t);
	    }
	    else if(strain==Q)
	    {
	        infect(n,strain - 1,t);
	    }
	    else
	    {
	        // mutation = +/- 1
	        mutation = gsl_rng_uniform_int(r,2)*2 - 1;
	        infect(n,strain + mutation,t);
	    }
	}
    }
}

int find_strain (	node * __restrict__ 	n,
			int 			Q) 
{
    //only applicable if just 1 strain per node
    int i,j;
    for (i=0; i<= Q / STRAINS_PER_SITE; ++i)
    {
	if(n->infection[i] > SUSCEPTIBLE)
	{
	    for (j=0; j<STRAINS_PER_SITE; ++j)
	    {
		if(BIT_POS(n->infection[i],j) > SUSCEPTIBLE)
		{
		    return (i * STRAINS_PER_SITE + j);
		}
	    }
	}
    }
    if(!RECORD_INFECTIONS_ON_GRAPH)
    {
        printf("No strain found\n");
        for (i=0; i<= Q / STRAINS_PER_SITE; ++i)
        {
	    printf("%i\t",n->infection[i]);
        }
	
        printf("\n\n");
    }
    return(-1);
}

unsigned int BIT_POS (	int 	var, 
			int 	pos) 
{
/* position of strain in infection [i] */

    return ((var)&(1<<(pos)));
}

void randomize ( 	int * __restrict__ 	list,
			int 			N,
			const gsl_rng * 	r) 
{
    int i;
    for(i=0; i<N; ++i) list[i] = i;
    gsl_ran_shuffle(r,list,N,sizeof(int));
}

int distinct_infection (node * __restrict__ n,
			node * __restrict__ neighbor) 
{
    /*
	Determine whether two nodes have distinct infections
    */

    int distinct=0;
    //ignore if node is susceptible, neighbor can be susceptible
    if(!is_susceptible(n))
    {
	distinct = !(n->infection[0] == neighbor->infection[0]);
    }
    return distinct;
}

int is_susceptible ( node * __restrict__ n) 
{
    /*
	Determine whether a node is susceptible
    */

    int susceptible = (n->infection[0] == 0);
    return susceptible;
}

void copy_infection ( 	graph * __restrict__ 	G,
			int 			neighbor,
			int 			n,
			int 			is_infected,
			int 			Q,
			double 			t) 
{
    int i;
    for(i=0; i<=Q/STRAINS_PER_SITE; ++i){
        G->nodes[neighbor].infection[i]=G->nodes[n].infection[i];
    }

    if(is_infected && RECORD_EXPOSURE)
    {
        int K = G->nodes[neighbor].frequency_exposed;
        int K_max = 10 * TAU_MAX * AVERAGE_DEGREE_MAX;
 	if(K < K_max)
            G->FreqExposureUntilInfectedHist[K]++;
    }

    G->nodes[neighbor].time_infected = t;
    G->nodes[neighbor].frequency_exposed = 0;
}

void infect ( 	node * __restrict__ 	n,
		int 			strain,
		double 			t) 
{

    //STRAINS_PER_SITE = 2^32 - 1
    n->infection[0] = strain;
    n->time_infected = t;
    n->time_last_checked = t;
    n->frequency_exposed = 0;
}

void count_infections ( graph * __restrict__ 	G,
              		int * __restrict__	NUMBER_INFECTED,
			double 			average_wait_tau,
			double 			wait_tau_alpha,
			gsl_rng * 		r)

{
    //initialization
    int n,strain,N=G->NumNodes,Q = G->Q;

    // Set NUMBER_INFECTED to all 0's
    memset(NUMBER_INFECTED, 0, Q * sizeof(NUMBER_INFECTED[0]));
    
    for(n=0; n<N; ++n)
    {
	if(!is_susceptible(&(G->nodes[n])))
	{
	    strain = G->nodes[n].infection[0];//find_strain(&(G->nodes[n]),Q);

	    if(SAMPLE)
	    {
		double rand = gsl_rng_uniform(r);
		if(rand < SAMPLE_FRACTION)
		{
		    NUMBER_INFECTED[strain - 1]++;//else this is an infected node, so record
		}
	    }
	    else
		NUMBER_INFECTED[strain - 1]+=1;//else this is an infected node, so record

	    //record time so far infected

	    //t0 = G->nodes[n].time_infected;
	    //G->nodes[n].time_infected ++;
		
	}
	// for all nodes, we determine whether links remain on or not
	if(AVERAGE_WAIT_TAU_MAX > eps)
	    link_dynamics(G,n,average_wait_tau,wait_tau_alpha,r);
    }
}


void link_dynamics( 	graph * __restrict__ 	G,
			int 			node,
			double 			average_wait_tau,
			double 			wait_tau_alpha,
			gsl_rng * 		r) 
{
    if(average_wait_tau > eps)
    {
        int j,neighbor;
        for(j=0; j<G->nodes[node].degree; ++j)
        {
	    neighbor = G->nodes[node].neighbors[j];
	    if(neighbor >=0)
	    {
	        if(G->nodes[node].wait_time[j]==G->nodes[node].wait_tau[j])
	        {
	            resample_wait_tau(G,node,j,average_wait_tau,wait_tau_alpha,r);
	        }
	        //time += 1, and wait_time starts at 1,
	        //hence when wait_time == wait_tau -> link active
	        G->nodes[node].wait_time[j]++;
	    }
	}
    }
}

void resample_wait_tau (graph * __restrict__ 	G,
			int 			node,
			int 			neighbor_pos,
			double 			tau_av,
			double 			alpha,
			gsl_rng * 		r) 
{
    /*
	We sample from:
		-Geometric Distribution = p (1-p)^k
		     	where p = 1/tau_av

		-Power Law Distribution = C k^-alpha,
			where we adjust the k=1,2,3... distribution to conform to the average
		     
    */
    int neighbor = G->nodes[node].neighbors[neighbor_pos],N=G->NumNodes;

    if(neighbor >= 0)
    {
	//printf("node: %i, neighbor: %i\n",node,neighbor);
        if(tau_av > eps)
        {
            if(alpha > eps)
            {
    	        G->nodes[node].wait_tau[neighbor_pos] = rand_SFtau(tau_av,alpha,N,r);
	    }
        
            else
 	    {
	        int p = 1/tau_av;
    	        G->nodes[node].wait_tau[neighbor_pos] = gsl_ran_geometric(r,p);
	    }
        }
        else
        {
            G->nodes[node].wait_tau[neighbor_pos] = 0;
        }

        G->nodes[node].wait_time[neighbor_pos] = 0;
    
        synchronize_link(G,node,neighbor_pos);
    }
}

int rand_SFtau (double 		tau_av, 
		double 		alpha, 
		int 		N,
		gsl_rng * 	r) 
{
    /*
	We create a power-law tail distribution,
	with a transient to create a consistant <tau> = tau_av
    */

    int tau_min,tau;
    double rand,P_tau_min_minus_1;
    find_distribution_fast(alpha,tau_av,&tau_min,&P_tau_min_minus_1);

    rand = gsl_rng_uniform(r);
    //this is to ensure that <k> remains constant
    if (rand < P_tau_min_minus_1) tau = tau_min - 1;
    else
    {
	 //continuous approximation. not accurate, but still has power law tail
	 //rand = gsl_rng_uniform(r);
 	 //tau = (int)((tau_min-0.5) * pow((1-rand),-1/(alpha-1)) + 0.5);
         tau = ran_scale_free(r,tau_min,alpha,N);
    }

    return tau;
}

void synchronize_link ( graph * __restrict__ 	G,
			int 			node,
			int 			neighbor_pos) 
{
    /*
	Make sure "wait_time", "wait_tau"
	are the same on both ends of a link
    */
    int neighbor = G->nodes[node].neighbors[neighbor_pos];
    if(neighbor >= 0)
    {
    	int node_pos = find_neighbor_pos(node,&(G->nodes[neighbor]));

    	int wait_time = G->nodes[node].wait_time[neighbor_pos];
    	int wait_tau = G->nodes[node].wait_tau[neighbor_pos];
    	if(node_pos >= 0)
    	{
            G->nodes[neighbor].wait_time[node_pos] = wait_time;
            G->nodes[neighbor].wait_tau[node_pos] = wait_tau;
	}
    }
}

void make_thread_args ( thread_args * 	args,
                        int 		thread,
                        int 		N,
			FILE **		s_deg,
			FILE **	    	s_graph,
                        double 		alpha,
                        double 		k_av,
                        double 		fract_tri_vertices,
                        double 		assort_p,
			int 		num_candidates) 
{
    if(av_deg <= 0 || num_candidates <= 0 || N <= 0)
    {
	printf("\nERROR: Average degree, Q, or N <= 0\n");
	exit(1);
    }
    if(alpha < 0 || fract_tri_vertices < 0)
    {
	printf("\nERROR: Alpha or Clustering coefficient < 0\n");
	exit(1);
    }

    args->thread = thread;
    //buffer for mercenne twister
    args->r = gsl_rng_alloc (gsl_rng_mt19937);
    seed_rand(thread,args->r);
    args->Q = num_candidates;
    make_graph(&(args->G),N,alpha,k_av,assort_p,fract_tri_vertices,args->Q,s_deg,s_graph,args->r);

    int beta_steps = ((int)((BETA_MAX-BETA_MIN)/BETA_STEP_SIZE+1+eps));
    int s_steps = ((int)((S_MAX - S_MIN)/S_STEP_SIZE+1+eps));
    int delta_infected_steps = ((int)((DELTA_INFECTED_MAX - DELTA_INFECTED_MIN)/DELTA_INFECTED_STEP_SIZE+1+eps));
    int hung_p_guilty_steps = ((int)((HUNG_P_GUILTY_MAX - HUNG_P_GUILTY_MIN)/HUNG_P_GUILTY_STEP_SIZE+1+eps));
    int hung_p_innocent_steps = ((int)((HUNG_P_INNOCENT_MAX - HUNG_P_INNOCENT_MIN)/HUNG_P_INNOCENT_STEP_SIZE+1+eps));
    int gamma_steps = ((int)((GAMMA_MAX-GAMMA_MIN)/GAMMA_STEP_SIZE+1+eps));
    int stop_steps = (STOP_MAX - STOP_MIN)/STOP_STEP_SIZE + 1 + eps;
    int binom_steps = (BINOM_INIT_CONDIT_MAX - BINOM_INIT_CONDIT_MIN)/BINOM_INIT_CONDIT_STEP_SIZE + 1 + eps;
    int mvm_p_steps = (MVM_P_MAX - MVM_P_MIN)/MVM_P_STEP_SIZE + 1 + eps;
    int tau_steps = (TAU_MAX - TAU_MIN)/TAU_STEP_SIZE + 1 + eps;
    int hung_steps = (HUNG_RATIO_MAX - HUNG_RATIO_MIN)/HUNG_RATIO_STEP_SIZE + 1 + eps;
    int wait_tau_alpha_steps = ((int)((WAIT_TAU_ALPHA_MAX-WAIT_TAU_ALPHA_MIN)/WAIT_TAU_ALPHA_STEP_SIZE+1+eps));
    int average_wait_tau_steps = ((int)((AVERAGE_WAIT_TAU_MAX-AVERAGE_WAIT_TAU_MIN)/AVERAGE_WAIT_TAU_STEP_SIZE+1+eps));
    int mu_steps = ((int)((MU_MAX-MU_MIN)/MU_STEP_SIZE+1+eps));
    int temp_steps = fabs(T_MAX-T_MIN)/T_STEP_SIZE+1+eps;

    int behavior_space = NUM_RUNS * num_candidates * delta_infected_steps * hung_p_innocent_steps * hung_p_guilty_steps * beta_steps * s_steps * gamma_steps * temp_steps * mvm_p_steps * tau_steps * hung_steps * stop_steps * binom_steps * wait_tau_alpha_steps * average_wait_tau_steps * mu_steps;
    if (VERBOSE)
        printf("BEHAVIOR SPACE: %i\n",behavior_space);

    //we allocate memory for the entire behavior space
    
    /*Allocate fraction of nodes infected with each strain
      for all time, for the entire behavior space = RUN_TIME * BEHAVIOR_SPACE*/

    if(!RECORD_LONG_CONSENSUS_TIME)
    {
        args->FractionInfectedVersusTime	= malloc(sizeof(double)*(RUN_TIME - IGNORE_TIME)*behavior_space);
        args->UVersusTime 			= malloc(sizeof(double)*(RUN_TIME - IGNORE_TIME));

        //mean, variance of fraction infected versus time for each strain in equilibrium
        //args->Variance = malloc(sizeof(double)*behavior_space);
        args->Mean = malloc(sizeof(double)*behavior_space);
        args->m = malloc(sizeof(double)*behavior_space / num_candidates);
        args->m2 = malloc(sizeof(double)*behavior_space / num_candidates);
        args->m4 = malloc(sizeof(double)*behavior_space / num_candidates);

        /*Allocate time_infected distribution for initial 
	transient time, from t=0:TAU_MAX,
           for the entire behavior space*/

	if(RECORD_TAU)
            args->T_distribution = malloc(sizeof(int)*behavior_space*(TAU_MAX+1)*(RUN_TIME - IGNORE_TIME));
    }

    else
    {
	args->ConsensusTimes = malloc(sizeof(double) * behavior_space/num_candidates);
        if(RECORD_DIST)
 	    args->NumsGuilty = malloc(sizeof(int) * behavior_space/num_candidates);
	if(RECORD_FREEZE_P)
        {
            args->FreezeTime = malloc(sizeof(double) * behavior_space/num_candidates);
            args->DynamicsP = malloc(sizeof(double) * behavior_space/num_candidates);
        }
    }

    if(RECORD_INFECTIONS_ON_GRAPH)
    {
        args->whos_infected = malloc(sizeof(double)*behavior_space/num_candidates * N);
    }
    if (RECORD_VOTE_TIME)
    {
	int num_times = 5;//timesteps: 1.0,2.0,3.0,4.0,5.0
	args->v0 = malloc(sizeof(double) * num_times * behavior_space/num_candidates);
    }


    if(RECORD_EXPOSURE)
    {
        //record
        args->TOTAL_FreqExposureUntilInfectedHist = malloc(sizeof(int)*behavior_space*(AVERAGE_DEGREE_MAX * TAU_MAX)*10);
        args->TOTAL_TimeUntilInfectedHist = malloc(sizeof(int)*behavior_space*(TAU_MAX*10+1));//record in steps of 0.1
    }


    if(CORRELATION)
    {
	double N = MIN_NUMBER_OF_NODES;
	int L = sqrt(N) + eps;
	int l = REGION_DIMENSION;
	int r_max = (L/l)/sqrt(2) + 1;
	int num_regions = (L/l) * (L/l);

	args -> Correlation = malloc(sizeof(double) * r_max * behavior_space);
	args -> Rho_i = malloc(sizeof(double) * num_regions * behavior_space);
    }
    if (VERBOSE)
        printf("Thread args made.\n");
}

void make_node (node * __restrict__ 	n,
		int 			ID,
		int 			deg,
		double 			fract_tri_vertices,
		int 			num_candidates) 
{
    //ID is assumed from array position in our code
    //but, in case this is important in the future,
    //we leave this here.
    n->ID = ID;
    n->degree=deg;
    n->tri_vertex_edges = 2 * (int)(deg*fract_tri_vertices/2);//take floor of value
    n->current_number_of_neighbors = 0;
    if(!COMPLETE) n->neighbors = malloc(sizeof(int)*deg);
    if(AVERAGE_WAIT_TAU_MIN > eps)
    {
        n->wait_tau = malloc(sizeof(int)*deg);
        n->wait_time = malloc(sizeof(int)*deg);
    }

    n->infection = malloc(sizeof(unsigned int)*(num_candidates / STRAINS_PER_SITE + 1));
    n->active_neighbors = malloc(0);//NULL;
    n->num_active_edges = 0;
    n->contrary_neighbors = malloc(0);//NULL;
    n->Bj = malloc(0);//NULL;
    n->num_contrary_edges = 0;
}

void make_graph (graph * __restrict__ 	G,
                 int 			N,
                 double 		alpha,
                 double 		k_av,
                 double 		assort_p,
                 double 		fract_tri_vertices,
		 int 			Q,
		 FILE ** 		s_deg,
		 FILE ** 		s_graph,
                 gsl_rng * 		r) 
{
    // graph properties
    G->NumNodes = N;
    G->NumEdges = 0;
    G->alpha = 0.0;
    G->average_degree = k_av;//we update this for complete/grid graph
    G->assort_p = assort_p;
    G->nodes = malloc(sizeof(node)*N);
    if(fabs(assort_p) > eps)
        G->degree_sum = malloc(sizeof(int)*N);

    // dynamics
    G->Q = Q;
    G->ConsensusTime = 0.0;
    G->NumGuilty = 0;
    G->MeanFractionInfected = malloc(sizeof(long double) * Q);
    if(INFECT_NEUTRAL)
    {
        G->CumSumBi = malloc(sizeof(double)*N);
        G->CumSumActiveEdges = malloc(sizeof(double)*N);
    }
    if(RECORD_EXPOSURE)
    {
        G->FreqExposureUntilInfectedHist = malloc(sizeof(int) * AVERAGE_DEGREE_MAX * TAU_MAX * 10);
        G->TimeUntilInfectedHist = malloc(sizeof(int) * (TAU_MAX * 10 + 1));//record in steps of 0.1
    }

    if(VERBOSE)
	printf("Making nodes...\n");
    if(REAL_DATA_DISTRIBUTION && alpha < eps && K_REGULAR == 0)
    {
	make_file_nodes(G,s_deg,s_graph);
    }
    else
    {
    if(GRID_GRAPH)
    {
	make_grid(G,N,assort_p);
    }
    else
    {
        //make SF degree network IF alpha > 0
        if(alpha > eps)
        {
            int k_max = N*K_CUT_OFF;
            make_SFdegree_nodes (G,N,alpha,k_max,k_av,assort_p,fract_tri_vertices,r);
        }
        else
        {
	    if(K_REGULAR || K_REGULAR_ANNEALED)
	    {
		make_Kregular_nodes(G,N,k_av,fract_tri_vertices);
	    }

	    else if(COMPLETE)
	    {
		make_complete_nodes(G,N);
	    }
	    else
	    {
                make_poisson_degree_nodes(G,N,k_av,fract_tri_vertices,assort_p,r);
	    }
        }
    }
    //whether a grid or random graph, we connect nodes with a single function
    if(R_C > eps)
	spatial_connect_nodes(G,r);
    else if (!COMPLETE)
	connect_nodes(G,r);
    }
}

void make_file_nodes (	graph * __restrict__ 	G,
			FILE ** 		s_deg,
			FILE ** 		s_graph)
{
/*
    Read in a graph file, and return the completely connected graph
	s_deg: stream of degree
	s_graph: stream of which nodes connect where
*/
    int i, j, degree, neighbor, N = G->NumNodes,Ntest = 0,Q = 0;
    double skip, fract_tri_vertices = 0.0;
    char ignore1[6], ignore2[2], ignore3[6], ignore4[2];

   ////////////////////////////  NUMBER OF EDGES (ASSUME <K> = 10) ////////////////////////
    G->NumEdges = N * 10/2.0;

   ////////////////////////////  	SKIP SIMULATION DATA	   ////////////////////////

    fscanf(*s_graph," %s %s %i %s %s %i\n",ignore1,ignore2,&Q,ignore3,ignore4,&Ntest);
    if(Q == 0)
    {
	printf("ERROR: Q = 0\n");
        printf("N: %i\tNtest: %i\nQ: %i\n",N,Ntest,Q);

	exit(1);
    }
    if(Ntest != N)
    {
	printf("%s %s %s %s\n",ignore1,ignore2,ignore3,ignore4);
	printf("ERROR: Files read incompatible N: %i != %i\n",Ntest,N);
	exit(1);
    }

    //ignore simulation data
    for(i = 0; i < Q; i++)
        fscanf(*s_graph," %lf",&skip);

   ////////////////////////////  	MAKE/CONNECT NODES	   ////////////////////////
    for(i = 0; i < N; ++i)
    {
	//Find degree
	fscanf(*s_deg," %i",&degree);
	G->nodes[i].degree = degree;
        //printf("k: %i\n",degree);
	
	//make nodes with known degree, clustering coefficient (ignore), and Q
        make_node(&(G->nodes[i]),i,degree,fract_tri_vertices,Q);

	//read connections from file
	for(j = 0; j<degree; ++j)
	{
	    fscanf(*s_graph," %i",&neighbor);
	    G->nodes[i].neighbors[j] = neighbor;
	}
    }
}



void make_grid (	graph * __restrict__	G,
			int 			N,
			double 			assort_p) 
{
   /*
	We are creating a grid with
		- K = 4
		- Periodic BC
        Periodic BC:
		-{i-1,i+L,i-L,i+1} (mod N = L^2)
   */

    int L = sqrt((double)(N))+eps;
    int l = REGION_DIMENSION;
    int num_regions = (L/l) * (L/l);
    int region = 0;
    int K = 4;
    int i, j, pos,Q = G->Q;

    //Updating degree and number of edges
    G->NumEdges = (K * N)/2; // property of grid graphs
    G->average_degree = K;

    //Allocating in what region each node is
    G->Regions = malloc(sizeof(int)*N);
    //allocating average infection (rho_i) for each region
    G->Rho_i = malloc(sizeof(double)* num_regions);

    //////////////////////////////////////////
    // make nodes with degree 4 //
    //////////////////////////////////////////
    for(i=0; i<L; i++)
    {
        for(j=0; j<L; j++)
        {
	    make_node(&(G->nodes[i* L + j]), i * L + j, K, 0.0,Q);
            if(fabs(assort_p) > eps)
            {
                if(i==0 && j == 0)
                {
                    G->degree_sum[0]=K;
                }
                else
                {
             	    G->degree_sum[i * L + j] = G->degree_sum[i * L + j - 1] + K;
                }
	    }
	}
    }
    if(VERBOSE)
	printf("Nodes mode.\n");
    if(VERBOSE)
	printf("Making regions...\n");


    ////////////////////////////
    //     define regions     //
    ////////////////////////////
    /*
    Algorithm:
	-l = region dimension
	-if i % l == 0, region = j -> j+1, for j = 1 ... L % l (assume its divisable)
    */
    if(l*l > N) 
    {
	printf("ERROR: region is larger than grid: make_grid()\n"); exit(1);
    }
    if(L % l != 0) 
    {
	printf("ERROR: grid dimenion is not divisable by region dimension: make_grid()\n"); 
	exit(1);
    }

    for(i=0; i<N; ++i)
    {
	pos = i % l + l * ((i/L) % l);
 
	
 	if(i > 0)
	{
  	    if((i % l) == 0) region++;
  	    if((region % (L/l)) == 0) region = L/l * ((i + 1)/(L* l));
  	}
 
	G->Regions[l * l * region + pos] = i;
    }
    if(VERBOSE)
	printf("Regions made.\n");

}

void make_SFdegree_nodes ( graph * __restrict__ G,
                           int 			N,
                           double 		alpha,
                           int 			k_max,
                           double 		k_av,
                           double 		assort_p,
                           double 		fract_tri_vertices,
                           const gsl_rng *	r) 
{
    int i,n,k_min,node_degree,Q = G->Q;
    double p;
    G->alpha = alpha;
    find_distribution(alpha,N,k_av,&k_min,k_max,&n);

    for(i=0; i<N; i++)
    {
        //node_degree = degrees[i];
	p = gsl_rng_uniform(r);
	if(p < (double)n/N)
	{
	    node_degree = k_min - 1;
	}
	else
	{
     	    do{
            	node_degree = ran_scale_free(r,k_min,alpha,N);
            }while(node_degree > k_max);
	}

        if(fabs(assort_p) > eps)
        {
            if(i==0)
            {
                G->degree_sum[0]=node_degree;
            }
            else
            {
                G->degree_sum[i] = G->degree_sum[i-1] + node_degree;
            }
        }
        G->NumEdges += node_degree;

        make_node(&(G->nodes[i]),i,node_degree,fract_tri_vertices,Q);
    }
 
    G->NumEdges = (G->NumEdges)/2;
}

void find_distribution_fast (	double 			alpha,
                                double 			k_av,
                                int * __restrict__ 	k_min,
                                double * __restrict__ 	P_k_min_minus_1) 
{
    double k;
    double kprime_av;
    if(k_av <= 0.0)
    {
	printf("\nERROR: Average degree <= 0\n");
	exit(1);
    }
    //assume k_av > 0
    if(k_av < 2.0)
    {
     	(*k_min) = 2;
        kprime_av = av_deg_fast(alpha,(double)(*k_min));
    }
    else{
        for(k=2; k<=k_av+eps; ++k)
        {
     	    kprime_av = av_deg_fast(alpha,k);

            if(kprime_av >= k_av)
            {
                (*k_min) = k;
                break;
	    }
        }
    }

    //P_k_min_minus_1 is the fraction of degrees (wait_tau's, etc) with k_min-1

    (*P_k_min_minus_1) = (kprime_av - k_av)/(kprime_av - ((*k_min) - 1));
}

void find_distribution (	double 			alpha,
                                int 			N,
                                double 			k_av,
                                int * __restrict__ 	k_min,
                                int 			k_max,
                                int * __restrict__ 	n) 
{

    double kprime_av;
    //find k_min s.t. avg(k') > k_av
    int k;
    for(k=2; k<=k_av; ++k)
    {
        kprime_av = av_deg(alpha,k,k_max);
        if(kprime_av > k_av)
        {
            (*k_min) = k;
            break;
        }
    }

    if(k_av < 2)
    {
     	(*k_min) = 2;
        kprime_av = av_deg(alpha,*k_min,k_max);
    }


    (*n) = (int)(N*(kprime_av - k_av)/(kprime_av - (*k_min - 1)));
}

double av_deg_fast ( double 	alpha,
                     double 	k_min) 
{
    double c,deg_av;
    c = 1/gsl_sf_hzeta(alpha,k_min);
    deg_av = c * gsl_sf_hzeta(alpha-1,k_min);
    return deg_av;
}

double av_deg ( double 		alpha,
                int 		k_min,
                int 		k_max) 
{
    double c,deg_av;

    c = 1/finite_size_poly_sum(alpha,k_min,k_max);

    deg_av = c * finite_size_poly_sum(alpha-1,k_min,k_max);

    return deg_av;

}

double finite_size_poly_sum( 	double 	alpha,
                                int 	k_min,
                                int 	k_max) 
{
    //clarified code such that we sum from k_min -> k_max, not N/2

    double s = gsl_sf_hzeta(alpha,k_min) - gsl_sf_hzeta(alpha,k_max+1);

    return s;
}

int ran_scale_free ( 	const gsl_rng * r,
                        int 		k_min,
                        double 		alpha,
			int 		N) 
{
    //Below method based on:
    // arXiv: 0706.1062 (By: Clauset, Shalizi, & Newman)


    double rand = gsl_rng_uniform(r);
    if (k_min>5)
        return (int)((k_min-0.5) * pow((1-rand),-1/(alpha-1)) + 0.5);
    else
    {
     	double x1 = 0.0,x2 = k_min,xp;
        double Px=2;
        //double up
        while (Px>1-rand)
        {
            x1=x2;
            x2=x1 * 2;

            Px = gsl_sf_hzeta(alpha,x2)/gsl_sf_hzeta(alpha,k_min);
        }
        //binary search
        while (x2-x1>1)
        {
            xp = x1 + ((x2-x1)/2);
            Px = gsl_sf_hzeta(alpha,xp)/gsl_sf_hzeta(alpha,k_min);
            if (Px<1-rand) x2=xp;

            else x1=xp;
        }
    
        if ((int) x1 <= N) return (int) x1;
        else return N;
    }

}

double ran_ref_1D_walk (
                        double 		T,
			const gsl_rng * r) 
{
     //T = T * 3600/3;//convert hours to timesteps (3 seconds)

     double delta = 0.0001;//error threshold
     double rand = gsl_rng_uniform(r);
     double x1=0,x2=sqrt(T),xp;
     double Px = 1 - gsl_sf_exp(-x2 * x2 / T);

     //double up
     while(Px < rand)
     {
         x1=x2;
         x2=x1 * 2;

         Px = 1 - gsl_sf_exp(-x2 * x2 / T);
     }

     //binary search
     while (x2-x1>delta)
     {
         xp = x1 + ((x2-x1)/2);
         Px = 1 - gsl_sf_exp(-xp * xp / T);

         if (Px<rand) x1 = xp;

         else x2 = xp;
     }

     //printf("rand: %f\nPx: %f\n s: %f\n",rand,Px,x1);
     if (x2 <= T) return (x1 + x2)/2;
     else return T;    
}


void make_Kregular_nodes (
			graph * __restrict__ 	G,
                        int 			N,
                        double 			av_deg,
			double 			fract_tri_vertices) 
{
    int Q = G->Q, node_degree = av_deg, i;

    G->average_degree = node_degree;
    G->NumEdges = N * node_degree/2.0;

    for (i=0; i<N; ++i)
    {
        make_node(&(G->nodes[i]),i,node_degree,fract_tri_vertices,Q);
    }
}

void make_complete_nodes (
                        graph * __restrict__ 	G,
                        int 			N) 
{
    int Q = G->Q;
    // links to all nodes but self
    // degree is therefore N-1
    G->average_degree = N-1;
    G->NumEdges = N * (N-1)/2;

    int i, node_degree = N-1;
    for (i=0; i<N; ++i)
    {
     	make_node(&(G->nodes[i]),i,node_degree,0.0,Q);
    }
}

void make_poisson_degree_nodes (graph * __restrict__ 	G,
				int 			N,
				double 			av_deg,
				double 			fract_tri_vertices,
				double 			assort_p,
				gsl_rng * 		r) 
{
    int Q = G->Q;

    int i, node_degree;
    for (i=0; i<N; ++i)
    {
	
        node_degree = gsl_ran_poisson(r,av_deg);
	if(fabs(assort_p) > eps)
	{
	    if(i==0)
	    {
	        G->degree_sum[0] = node_degree;
	    }
	    else
	    {
	        G->degree_sum[i] = G->degree_sum[i-1] + node_degree;
	    }
	}
	G->NumEdges += node_degree;

        make_node(&(G->nodes[i]),i,node_degree,fract_tri_vertices,Q);
    }
    G->NumEdges = (G->NumEdges)/2;
}

void connect_nodes ( 	graph * __restrict__ 	G,
                        const gsl_rng *		r) 
{
    if(VERBOSE)
        printf("\n Connecting nodes...\n");

    int i,
	i_ran,
	j,
        N = G->NumNodes,
	triangle_edges,
	current_num_neighbors,
	degree,
	neighbors_number_of_neighbors,
	neighbor,
	neighbor1,
	neighbor2,
	multi_edge,
	multi_edge1,
	multi_edge2,
	node_number_of_neighbors,
	neighbor1_number_of_neighbors,
	neighbor2_number_of_neighbors,
	x_pos,y_pos;

    if(GRID_GRAPH)
    {
        int L = sqrt((double)N) + eps;
        for(i = 0; i< L; ++i)
        {
	    for(j = 0; j < L; ++j)
	    {
		x_pos = (j - 1 + L) % L;
		if(x_pos < 0) x_pos += L;
		y_pos = i * L;
		neighbor = (x_pos + y_pos) % N;
		if(neighbor < 0) neighbor += N;
     	        G->nodes[i * L + j].neighbors[0] = neighbor;
		//printf("%i \n",neighbor);

		x_pos = (j + 1 + L) % L;
		if(x_pos < 0) x_pos += L;
		y_pos = i * L;
		neighbor = (x_pos + y_pos) % N;
		if(neighbor < 0) neighbor += N;
                G->nodes[i * L + j].neighbors[1] = neighbor;
		//printf("%i \n",neighbor);

		x_pos = j;
		y_pos = (i - 1 + L) % L;
		if(y_pos < 0) y_pos += L;
		y_pos = y_pos * L;
		neighbor = (x_pos + y_pos) % N;
		if(neighbor < 0) neighbor += N;
                G->nodes[i * L + j].neighbors[2] = neighbor;
		//printf("%i \n",neighbor);

		x_pos = j;
		y_pos = (i + 1 + L) % L;
		if(y_pos < 0) y_pos += L;
		y_pos = y_pos * L;
		neighbor = (x_pos + y_pos) % N;
		if(neighbor < 0) neighbor += N;
                G->nodes[i * L + j].neighbors[3] = neighbor;
		//printf("%i \n",neighbor);
	    }
	}
	if(VERBOSE)
            printf("GRID CONNECTED \n");
    }
    else if(SPATIALLY_CLOSE)
    {
	/*
	    Organize nodes on a square lattice:
		.   .   .   .   .   .
		.   .   .   .   .   .
		.   .   .   .   .   . 
		.   .   .   .   .   .
		.   .   .   .   .   .
	     Note: nodes have degree heterogeneity (e.g., scale-free)
	     Connect from closest outward:
		e.g., degree k = 8
		.   .   .
		 \  |  /
		.__\./__.
		   /|\
		./  .  \.
	     Edges are *directed* to avoid double edges
	     
	     Implication: neutral infections cannot work! Outward is OK!
	*/
	/*
	     Algorithm:
		dist = 0;
		node_k = ...;
		current_k = 0;
		While(node isn't fully connected)
		    1) dist+= D
		    2) connect to all nodes with d(i,j) < dist:
			- (x_i,y_i) is current coordinate
			- dx, dy is delta x, y pos
			- y_min, y_max are respective min, max values for y to reduce search space
			- x_min, x_max are respective min, max values for x
			- move outward by shell:
			    -- dx = 0, dy = -floor(dist), -floor(dist) + 1, ... , floor(dist) - 1, floor(dist)
			    -- dx = 1, dy = -floor(dist), -floor(dist) + 1, ... , floor(dist) - 1, floor(dist)
				...
			    -- if d = sqrt(dx^2 + dy^2) > dist, y < 0 dy++, 
			    -- if d = sqrt(dx^2 + dy^2) < dist, >dist - D, y_max = dy 
				    -current_k ++;	
			    -- if d = sqrt(dx^2 + dy^2) < dist, <dist - D, y_min = dy, skip to y = -y
			- STOP when y > 0 reaches y_max (d > dist) OR current_k == node_k
	*/
	double dist, shell_width=0.5,dx_j,dy_j,d;
	int y_max, current_k, node_k,x_i,y_i,outer_shell = 0;
	int L = sqrt(N) + eps, node1, node2, old_neighbor_pos1, old_neighbor_pos2, num_times;
	int l = REGION_DIMENSION;
        int num_regions = (L/l) * (L/l);
        int region = 0, pos;

        int edges [4], reordered_edges [2];

	for(i = 0; i < N; ++i)
	{
	    //(x,y) position
	    x_i = i % L;
	    y_i = i/L;
	    //node degree
	    node_k = G->nodes[i].degree;
	    //current degree (i.e., node not connected)
	    current_k = 0;
	    //move outward from a shells of radius "dist"
	    //and width D
	    dist = 0;
	    //while node isn't connected
	    while(current_k < node_k)
	    {
		//look at a larger shell
		dist+= shell_width;
		
		//stop if you risk double edges
		if(dist > sqrt(N)) break;
		
		//look no further than the edge of the shell
		y_max = dist;
		//for x values around x_i
		for(dx_j = -dist; dx_j <= dist; dx_j++)
		{
		    //assume we haven't reached the outer shell: dist = sqrt(dx^2+dy^2)
  		    outer_shell = 0;
		    //for y values around y_i
		    for(dy_j = -dist; dy_j <= y_max; dy_j++)
		    {
			
			//current distance from x_i, y_i
			d = sqrt(dx_j*dx_j + dy_j*dy_j);
			//if d is within the shell... (if outside the shell, we ignore)
			if(d < dist && d > dist - shell_width)
			{
			    
			    //if we hadn't reached the outer shell yet
			    if(outer_shell == 0)
			    {
				//y_max is the outer shell y value
				y_max = dy_j;
				outer_shell = 1;
			    }
			    
			    //associated neighbor we will connect with
			    //if we don't reach edge
			    y_pos = (y_i + (int)dy_j) * L;
			    x_pos = (x_i + (int)dx_j) % L;
			    //make operation modulo
			    if(x_pos < 0) x_pos += L;

			    neighbor = (x_pos + y_pos) % N;
			    //make operation modulo
			    if(neighbor < 0) neighbor += N;

			    if(current_k < node_k)
			    {
				//if we haven't completely connected the node, attach to neighbors
			    	G->nodes[i].neighbors[current_k]=neighbor;
			    	current_k++;
			    }
			    else//else node is filled, so we break
				break;
			}
			//if we are inside the shell, we skip to the other side
			if(d < dist - shell_width)
			{
			    dy_j = -dy_j;
			}
		    }
		}
	    }
	}

	//once our spatial scalefree network is created, we randomize links...
        if(fabs(G->assort_p) > eps)
        {
	    if(VERBOSE)
	        printf("BEGINNING RE-WIRING\n");
	    //list of two random edges, and outer nodes of re-ordered edges
	    //rewire this number of times (fraction of "assort_p")
	    num_times = (G->NumEdges) * G->assort_p;

            for (i=0; i<num_times; ++i)
            {
		//find random (directed) edges
	        find_random_edges(G,edges,r);
		
		// initialize values of "reordered_edges"
   	        for(j=1; j<4; j+= 2) reordered_edges[(j-1)/2] = edges[j];
	 	// reorder edges
	        gsl_ran_shuffle(r,reordered_edges,2,sizeof(int));
		//printf("edges:\n %i\t%i\n%i\t%i\n",edges[0],edges[1],edges[2],edges[3]);
		//printf("reordered edges:\n %i\t%i\n%i\t%i\n",edges[0],reordered_edges[0],edges[2],reordered_edges[1]);
		
		// root node
		node1 = edges[0];
		// ordered link root node connects to
		neighbor1 = reordered_edges[0];
		node2 = edges[2];
		neighbor2 = reordered_edges[1];

		//position of old neighbor
		old_neighbor_pos1 = find_neighbor_pos(edges[1],&(G->nodes[node1]));
		//new neighbor value...
		G->nodes[node1].neighbors[old_neighbor_pos1] = neighbor1;
		old_neighbor_pos2 = find_neighbor_pos(edges[3],&(G->nodes[node2]));
		G->nodes[node2].neighbors[old_neighbor_pos2] = neighbor2;
            }
        }

	if(CORRELATION)
	{
            //Allocating in what region each node is
            G->Regions = malloc(sizeof(int)*N);
            //allocating average infection (rho_i) for each region
            G->Rho_i = malloc(sizeof(double)* num_regions);

            if(VERBOSE)
		printf("Making regions...\n");

            ////////////////////////////
            //     define regions     //
            ////////////////////////////
            /*
            Algorithm:
                -l = region dimension
                -if i % l == 0, region = j -> j+1, for j = 1 ... L % l (assume its divisable)
    	    */
    	    if(l*l > N)
    	    {
                printf("ERROR: region is larger than grid: spatial_connect_nodes()\nl^2 = %i\tN = %i",l*l,N); exit(1);
            }
            if(L % l != 0)
            {
     	        printf("ERROR: grid dimenion is not divisable by region dimension: spatial_connect_nodes()\n");
                exit(1);
    	    }
    	    for(i=0; i<N; ++i)
    	    {
     	        pos = i % l + l * ((i/L) % l);

                if(i > 0)
                {
            	    if((i % l) == 0) region++;
            	    if((region % (L/l)) == 0) region = L/l * ((i + 1)/(L* l));
                }
                G->Regions[l * l * region + pos] = i;
    	    }
    	    if(VERBOSE)
                printf("Regions made.\n");
	}
    }

    else if(COMPLETE)
    {
	/*
        for(i = 0; i< N; ++i)
        {
	    //if(VERBOSE)
	    //    printf("\n");

	    for(j = 0; j < N; ++j)
	    {
     	        G->nodes[i].neighbors[j] = j;
		//if(VERBOSE)
		//    printf("%i ",j);
	    }
	}
	*/
	if(VERBOSE)
            printf("COMPLETE GRAPH CONNECTED \n");

    }

    //RANDOM
    else
    {

        /*********************************************************************************/
        ///// Randomly Connect Nodes (no clustering or assortativity) /////
        /*********************************************************************************/

        //randomize(random_order_i,N,r);

       for (i=0; i<N; ++i)
        {
            i_ran = i;//gsl_rng_uniform_int(r,N);//random_order_i[i];
            degree=G->nodes[i_ran].degree;
            current_num_neighbors = G->nodes[i_ran].current_number_of_neighbors;
            triangle_edges = G->nodes[i_ran].tri_vertex_edges;
            //create random neighbors
            for (j = current_num_neighbors;
                  j < (degree - triangle_edges);
                  ++j)
            {

             	neighbor = find_neighbor(G,i_ran,j,&multi_edge,r);

                G->nodes[i_ran].current_number_of_neighbors+=1;

                if (!multi_edge && i_ran!=neighbor)
                {
                    G->nodes[i_ran].neighbors[j]=neighbor;

                    neighbors_number_of_neighbors = G->nodes[neighbor].current_number_of_neighbors;
                    G->nodes[neighbor].neighbors[neighbors_number_of_neighbors] = i_ran;
                    G->nodes[neighbor].current_number_of_neighbors+=1;
                }
                else
                //link is a multiple of some other edge or is a self-loop
                {
                    //set to a negative number (meaning we ignore this link)
                    G->nodes[i_ran].neighbors[j]=-1;
                }
            }
        }
    }

    /*********************************************************************************/
    ///////////////// Assortative Degree Mixing /////////////////
    /*********************************************************************************/
    /*
	Method here is based on: R. Xulvi-Brunet, I. Sokolov,
	"Reshuffling scale-free networks: From random to assortative", Phys Rev. E (2004)

	Link: http://journals.aps.org/pre/pdf/10.1103/PhysRevE.70.066102

	Algorithm:
		Pick two edges at random
		Sort nodes with respect to degrees
		With probability p,
			Assortative: rewire larger (smaller) degrees together
			Disassortative: rewire largest and smallest degrees together
		Otherwise, rewire randomly
		If edges already exist between nodes, discard choice, find new pair
    */

    int num_times;
    
    if(fabs(G->assort_p) > eps && !SPATIALLY_CLOSE)
    {
	if(VERBOSE)
	    printf("BEGINNING RE-WIRING\n");

        int random_nodes [4], ordered_nodes [4];

	num_times = (G->NumEdges)/2;

	////////////////////////////////////////////////////////
	// If a grid, we are using this algorithm to create a //
	// small-world property instead of assortative mixing //
	////////////////////////////////////////////////////////
	if(GRID_GRAPH)
	    num_times = num_times * G->assort_p;
	////////////////////////////////////////////////////////
    
        for (i=0; i<num_times; ++i)
        {
	
	    find_random_edges(G,random_nodes,r);
	    
   	    for(j=0; j<4; ++j) ordered_nodes[j] = random_nodes[j];
	
	    if(GRID_GRAPH)
	        gsl_ran_shuffle(r,ordered_nodes,4,sizeof(int));
	    else
	        //sort random nodes from smallest to largest
	        sort(G,ordered_nodes);

	    rewire(G,random_nodes,ordered_nodes,G->assort_p,r);
	
        }
    }
    
    /*********************************************************************************/
    /////////////////// Randomly Cluster Nodes /////////////////
    /*********************************************************************************/

    /*
	Method here is based on: M. Newman, "Random graphs with clustering", PRL (2009)

	Link: http://arxiv.org/pdf/0903.4009v1.pdf

    */
    if(FRACT_TRI_VERTEX_MAX > eps)
    {
	int random_order_i[N];
    	randomize(random_order_i,N,r);

    	for (i=0; i<N; ++i)
    	{
	    i_ran = random_order_i[i];
     	    degree=G->nodes[i_ran].degree;
            current_num_neighbors = G->nodes[i_ran].current_number_of_neighbors;
	
            //for every pair of random neighbors, connect to each other
            for (j = current_num_neighbors;
            	j < degree;
             	j+=2)
            {
            	find_tri_neighbors(G,i_ran,j,&multi_edge1,&multi_edge2,&neighbor1,&neighbor2,r);

	    	node_number_of_neighbors = G->nodes[i_ran].current_number_of_neighbors;

            	if (!multi_edge1 && !multi_edge2 && i_ran!=neighbor1 && i_ran!=neighbor2 && neighbor1 != neighbor2)
            	{
             	    neighbor1_number_of_neighbors = G->nodes[neighbor1].current_number_of_neighbors;
                    neighbor2_number_of_neighbors = G->nodes[neighbor2].current_number_of_neighbors;

                    G->nodes[i_ran].neighbors[node_number_of_neighbors ] = neighbor1;
                    G->nodes[i_ran].neighbors[node_number_of_neighbors+1] = neighbor2;

	 	    G->nodes[neighbor1].neighbors[neighbor1_number_of_neighbors ] = i_ran;
                    G->nodes[neighbor1].neighbors[neighbor1_number_of_neighbors+1] = neighbor2;

                    G->nodes[neighbor2].neighbors[neighbor2_number_of_neighbors ] = i_ran;
                    G->nodes[neighbor2].neighbors[neighbor2_number_of_neighbors+1] = neighbor1;

                    G->nodes[neighbor1].current_number_of_neighbors += 2;
                    G->nodes[neighbor2].current_number_of_neighbors += 2;
                    G->nodes[i_ran].current_number_of_neighbors += 2;
            	}
            else
     	    	{
	            G->nodes[i_ran].neighbors[node_number_of_neighbors ]=-1;
                    G->nodes[i_ran].neighbors[node_number_of_neighbors+1]=-1;
                    G->nodes[i_ran].current_number_of_neighbors+=2;
            	}  
	    }
    	}
    }

    
    if(VERBOSE)
        printf("Graph Created\n");
}


void spatial_connect_nodes ( 	
			graph * __restrict__ 	G,
                        const gsl_rng *		r) 
{
    int i,node,N = G->NumNodes,num_neighbors_i=0,k_i;
    int L = sqrt((double)(N))+eps;
    int l = REGION_DIMENSION;
    int num_regions = (L/l) * (L/l);
    int region = 0, pos;

    //Allocating in what region each node is
    G->Regions = malloc(sizeof(int)*N);
    //allocating average infection (rho_i) for each region
    G->Rho_i = malloc(sizeof(double)* num_regions);

    for (i=0; i<N; ++i)
    {
	k_i = G->nodes[i].degree;
	while(num_neighbors_i < k_i)
	{
	    node = gsl_rng_uniform_int(r,N);
	    connect_node(G,i,node,r);
	    num_neighbors_i = G->nodes[i].current_number_of_neighbors;
	}
    }

   if(CORRELATION)
   {
       if(VERBOSE)
            printf("Making regions...\n");

        ////////////////////////////
        //     define regions     //
        ////////////////////////////
        /*
    	Algorithm:
            -l = region dimension
            -if i % l == 0, region = j -> j+1, for j = 1 ... L % l (assume its divisable)
    	*/
    	if(l*l > N)
    	{
            printf("ERROR: region is larger than grid: spatial_connect_nodes()\nl^2 = %i\tN = %i",l*l,N); exit(1);
    	}
    	if(L % l != 0)
    	{
     	    printf("ERROR: grid dimenion is not divisable by region dimension: spatial_connect_nodes()\n");
            exit(1);
    	}
    	for(i=0; i<N; ++i)
    	{
     	    pos = i % l + l * ((i/L) % l);


            if(i > 0)
            {
            	if((i % l) == 0) region++;
            	if((region % (L/l)) == 0) region = L/l * ((i + 1)/(L* l));
            }

            G->Regions[l * l * region + pos] = i;
    	}
    	if(VERBOSE)
            printf("Regions made.\n");
    }
}

double dist (int i, int j, int L)
{
    double x_i,y_i,x_j,y_j,dx,dy,d;
    x_i = i % L;
    y_i = i/L;
    x_j = j % L;
    y_j = j/L;

    dx = fabs(x_i - x_j);
    if(L - dx < dx)
	dx = L-dx;
    dy = fabs(y_i - y_j);
    if(L - dy < dy)
	dy = L-dy;
    d = sqrt(dx*dx + dy*dy);
    return d;
}

void connect_node (	graph * __restrict__ 	G,
			int 			x,
			int			y,
			const gsl_rng *		r)
{
    int num_neighbors_x, k_x, num_neighbors_y, k_y,N,L;
    N = G->NumNodes;
    L = sqrt(N) + eps;
    double p,d, spatial_p;

    num_neighbors_x = G->nodes[x].current_number_of_neighbors;
    k_x 	    = G->nodes[x].degree;
    num_neighbors_y = G->nodes[y].current_number_of_neighbors;
    k_y 	    = G->nodes[y].degree;

    if(num_neighbors_x < k_x && num_neighbors_y < k_y)
    {
        d = dist(x,y,L);	
	p = gsl_rng_uniform(r);
	spatial_p = 1/R_C * gsl_sf_exp(-d/R_C);
	if(p < spatial_p)
	{
	    //connect nodes
	    G->nodes[x].neighbors[num_neighbors_x] = y;
	    G->nodes[x].current_number_of_neighbors ++;

	    G->nodes[y].neighbors[num_neighbors_y] = x;
	    G->nodes[y].current_number_of_neighbors ++;
	}
    }
}

void find_random_edges (graph * __restrict__ 	G,
			int 			random_nodes [],
			const gsl_rng * 	r) 
{
    int first_edge_pos = 0;
    pick_edge(G,random_nodes,first_edge_pos,r);

    int second_edge_pos = 2;
    pick_edge(G,random_nodes,second_edge_pos,r);

}


void pick_edge (graph * __restrict__ 	G,
		int 			edge_array [],
		int 			index,
		const gsl_rng * 	r) 
{
    int random_edge,
	index1,
	index2,
	N=G->NumNodes,
	E=G->NumEdges;
    do{
        do{
            //pick a random edge
            random_edge = gsl_rng_uniform_int(r,2*E);
	   	
            //find a node attached to that random edge
            index1 = node_search_int(G->degree_sum,random_edge,N);
	
        }while(G->nodes[index1].degree == 0);
	    

        if (index1 > 0)
        {
            index2 = G->nodes[index1].neighbors[random_edge % G->degree_sum[index1-1]];
        }
        else
        {
     	    index2 = G->nodes[index1].neighbors[random_edge];
        }
    }while(index2 < 0);
    //printf("%i - %i\n",index1,index2);
	
    edge_array[index] = index1;
    edge_array[index + 1] = index2;
}

int node_search_int(	int * __restrict__	degree_sum,
			int 			random_edge, 
			int 			N) 
{
    int x1=0,x2=0,xp;

    //double up
    while (degree_sum[x2] <= random_edge)
    {
        x1=x2;
	if(x2==0) x2++;
	else x2=x1*2;

	if(x2 >= N)
	{
	    x2 = N - 1;
	    break;
	}
    }

    //binary search
    while (x2-x1>1)
    {
        xp = x1 + ((x2-x1)/2);
        if (degree_sum[xp]>random_edge) x2=xp;

        else x1=xp;
    }

    if (x2 < N) return x2;
    else return N - 1;
}


void sort ( 	graph * __restrict__ 	G,
		int * 			node_array) 
{
    //sort DEGREE
    int i,degree_array [4];

    const size_t STRIDE = 1;
    size_t n = 4;
    for (i=0; i<4; ++i)
    {
	degree_array[i]=G->nodes[node_array[i]].degree;

    }

    gsl_sort2_int(degree_array,STRIDE,node_array,STRIDE,n);

}

void rewire ( 	graph * __restrict__ 	G,
		int * 			random_nodes,
		int * 			ordered_nodes,
		double 			p,
		const gsl_rng * 	r) 
{
    // if rand < |p|
    // rewire assortatively if p > 0, else dis-assortatively
    // if links already exist, discard
    // else remove old links and replace with new one

    //else if rand > |p|
    // randomly re-wire
    // remove old links, replace with new ones
    
    double rand = gsl_rng_uniform(r);
    int new_neighbor_0;
    int new_neighbors [4];

    if(rand < fabs(p))
    {
	
	if(p<0)
	{
	
	    //dis-assortative
	    //new links: ordered_nodes[0] - ordered_nodes[3], ordered_nodes[1] - ordered_nodes[2]
	    new_neighbors [0] = ordered_nodes[0];
	    new_neighbors [1] = ordered_nodes[3];
	    new_neighbors [2] = ordered_nodes[1];
	    new_neighbors [3] = ordered_nodes[2];
	    
 	    rewire_links (G,new_neighbors,random_nodes);
	
	}
	
	else
	{
	    //assortative
	    //new links: ordered_nodes[0] - ordered_nodes[1], ordered_nodes[2] - ordered_nodes[3]
            new_neighbors [0] =	ordered_nodes[0];
            new_neighbors [1] =	ordered_nodes[1];
            new_neighbors [2] =	ordered_nodes[2];
            new_neighbors [3] =	ordered_nodes[3];

	    rewire_links (G,new_neighbors,random_nodes);
	}
    }

    else
    {
	
	//rewire randomly

	//random_node[0] neighbor = random_node[neighbor_1]
	new_neighbor_0 = gsl_rng_uniform_int(r,2) + 1;
	
	if(new_neighbor_0 == 1)
	{
	    // links: 0 - 1, 2 - 3

            new_neighbors [0] =	random_nodes[0];
            new_neighbors [1] = random_nodes[1];
            new_neighbors [2] =	random_nodes[2];
            new_neighbors [3] =	random_nodes[3];

	}
	
	else if(new_neighbor_0 == 2)
        {
	    // links: 0 - 2, 1 - 3

            new_neighbors [0] = random_nodes[0];
            new_neighbors [1] = random_nodes[2];
            new_neighbors [2] = random_nodes[1];
            new_neighbors [3] = random_nodes[3];

        }
	else if(new_neighbor_0 == 3)
        {
	    // links: 0 - 3, 1 - 2

            new_neighbors [0] = random_nodes[0];
            new_neighbors [1] = random_nodes[3];
            new_neighbors [2] = random_nodes[1];
            new_neighbors [3] = random_nodes[2];
	}
	
	rewire_links (G,new_neighbors,random_nodes);
	
    }

}

void rewire_links ( 	graph * __restrict__ 	G,
			int * __restrict__ 	new_neighbors,
			int * __restrict__	old_neighbors) 
{

    //if new_neighbor links don't already exist
    if(!(already_neighbors(G,new_neighbors)))
    {
	replace_links(G,old_neighbors,new_neighbors);
    }
}

int already_neighbors ( graph * __restrict__ 	G,
			int * __restrict__	new_neighbors) 
{
    int i,j,node=-1,neighbor=-1;
    for(i=0; i<4; i+=2)
    {
	node = new_neighbors[i];
	neighbor = new_neighbors[i+1];

	//if node has degree 0
	if(G->nodes[node].degree == 0) return 1;
	
	for(j=0; j<G->nodes[node].degree; ++j)
	{
	    if(neighbor == G->nodes[node].neighbors[j] || neighbor == node)
	    {
		return 1;
	    }
	    
	}
    }
    //else not a neighbor
    return 0;
}

void replace_links ( 	graph * __restrict__ 	G,
			int * __restrict__ 	old_neighbors,
			int * __restrict__	new_neighbors) 
{
    /*
	replace old links with new links:

    		-find position of old neighbor for each node
    		-replace old neighbors with new neighbors
    */
    int i,n1,old_neighbor,old_neighbor_pos,new_neighbor;

    for (i=0; i<4; ++i)
    {
	n1 = new_neighbors[i];
	old_neighbor = find_neighbor_from_list(n1,old_neighbors);

	new_neighbor = find_neighbor_from_list(n1,new_neighbors);

	old_neighbor_pos = find_neighbor_pos (old_neighbor,&(G->nodes[n1]));

	if (old_neighbor_pos >=0)
	{
	    G->nodes[n1].neighbors[old_neighbor_pos] = new_neighbor;
	}
	else
	{
	    printf("\nERROR: old_neighbor position not found: replace_links()\n");
	}
    }
}

int find_neighbor_from_list (	int n1, 
				int neighbor_list []) 
{
    int i;
    for (i=0; i<4; ++i)
    {
	if(n1 == neighbor_list[i])
	{
	    if(i%2 ==0)
	    {
		return neighbor_list[i+1];
	    }
	    else
	    {
		return neighbor_list[i-1];
	    }
	}
    }
    //if n1 is not in neighbor_list
    printf("ERROR: n1 is not in neighbor list: find_neighbor_from_list\n");
    return -1;
}

int find_neighbor_pos (	int 			neighbor, 
			node * __restrict__ 	n) 
{

    int i;
    for (i = 0; i<n->degree; ++i)
    {
	if (n->neighbors[i] == neighbor)
	{
	    //return position of neighbor
	    return i;
	}
    }
    //if neighbor not found
    printf("\n %i ERROR: neighbor not found: find_neighbor_pos()\n",neighbor);
    return -1;
}

int find_neighbor( 	graph * __restrict__ 	G,
			int 			node,
                        int 			j,
                        int * __restrict__ 	multi_edge,
                        const gsl_rng * 	r) 
{
    double rand = 0;
    int neighbor,k,N;

    N = G->NumNodes;

    do{

          neighbor = gsl_rng_uniform_int (r, N);
          (*multi_edge) = 0;
          for(k=0; k<j;++k)
          {
              if(G->nodes[node].neighbors[k]==neighbor)
              {
               	  (*multi_edge) = 1;
                  break;
              }
          }
     	  rand=0;
	  if(*multi_edge || node==neighbor)
          {
              rand = gsl_rng_uniform(r);
	  }

    }while (G->nodes[neighbor].current_number_of_neighbors
            ==G->nodes[neighbor].degree-G->nodes[neighbor].tri_vertex_edges
            || rand > RETRY_P);
    return neighbor;
}

void find_tri_neighbors(	graph * __restrict__ 	G,
                                int 			node,
                                int 			current_number_of_neighbors,
                                int * __restrict__ 	multi_edge1,
                                int * __restrict__ 	multi_edge2,
                                int * __restrict__ 	neighbor1,
                                int * __restrict__ 	neighbor2,
                                const gsl_rng * 	r) 
{
    //returns 2 neighbors and information on whether either is a multi-edge

    double rand = 0;
    int k,N;

    N = G->NumNodes;

    do{
          (*neighbor1) = gsl_rng_uniform_int (r, N);

          (*multi_edge1) = 0;
          for(k=0; k<current_number_of_neighbors;++k)
          {
              if(G->nodes[node].neighbors[k]==(*neighbor1))
              {
               	  (*multi_edge1) = 1;
                  break;
              }
          }
          rand=0;
          if((*multi_edge1) || node==(*neighbor1))
          {
              rand = gsl_rng_uniform(r);
          }

    }while (G->nodes[*neighbor1].current_number_of_neighbors
            ==G->nodes[*neighbor1].degree
            || rand > RETRY_P);
    do{
          (*neighbor2) = gsl_rng_uniform_int (r, N);

          (*multi_edge2) = 0;
          for(k=0; k<current_number_of_neighbors;++k)
          {
              if(G->nodes[node].neighbors[k]==*neighbor2)
              {
               	  (*multi_edge2) = 1;
                  break;
              }
          }

          rand=0;
          if((*multi_edge2) || node==(*neighbor2) || (*neighbor1)==(*neighbor2))
          {
              rand = gsl_rng_uniform(r);
          }

    }while (G->nodes[*neighbor2].current_number_of_neighbors
            ==G->nodes[*neighbor2].degree
            || rand > RETRY_P);
}

void free_thread_args (thread_args * args) 
{
    if(!RECORD_LONG_CONSENSUS_TIME)
    {
        free(args->FractionInfectedVersusTime);
        free(args->UVersusTime);
        free(args->Mean);
        free(args->m);
        free(args->m2);
        free(args->m4);
	if(RECORD_TAU)
            free(args->T_distribution);
    }
    
    else
    {
	free(args->ConsensusTimes);
	if(RECORD_DIST)
  	    free(args->NumsGuilty);
	if(RECORD_FREEZE_P)
        {
            free(args->FreezeTime);
            free(args->DynamicsP);
        }
    }
    if(RECORD_VOTE_TIME)
    {
        free(args -> v0);
    }

    if(RECORD_EXPOSURE)
    {
        free(args->TOTAL_FreqExposureUntilInfectedHist);
        free(args->TOTAL_TimeUntilInfectedHist);
    }

    if(CORRELATION)
    {
	free(args->Correlation);
	free(args->Rho_i);
    }

    gsl_rng_free (args->r);
    free_graph(&(args->G));
}

void free_node (node * __restrict__ n) 
{
    if (!COMPLETE)
        free(n->neighbors);
    free(n->infection);
    free(n->active_neighbors);
    free(n->contrary_neighbors);
    free(n->Bj);
    if(AVERAGE_WAIT_TAU_MIN > eps)
    {
        free(n->wait_tau);
        free(n->wait_time);
    }
}


void free_graph (graph * __restrict__ G) 
{
    //we free this just after each simulation
    //free(G->InfectiousNodes);
    free(G->MeanFractionInfected);

    if(fabs(G->assort_p) > eps)
        free(G->degree_sum);
    if(INFECT_NEUTRAL)
    {
        free(G->CumSumBi);
        free(G->CumSumActiveEdges);
    }
    if(RECORD_EXPOSURE)
    {
        free(G->FreqExposureUntilInfectedHist);
        free(G->TimeUntilInfectedHist);
    }
    if(CORRELATION)
    {
        free(G->Regions);
        free(G->Rho_i);
    }

    if (RECORD_VOTE_TIME)
    {
        free(G->v0);
    }

    //before we free the "nodes" pointer
    //we must first free each node's associated pointer
    int i,N=G->NumNodes;
    for (i=0; i<N; ++i)
        free_node(&(G->nodes[i]));
    //now each node's pointers are free, we are done
    free(G->nodes);
}
