//GA parameters
#define PI 3.1415926
#define POPSIZE 25 //control population size, excluding master
#define MAXGEN 200 //maximum generations for control GA
#define NCTRL 12 //total no. of controls (10 PDF, 1 width, 1 amplitude)
#define NSAVED 3 //best gene saved in each generation (elitism)
#define NUSED 20 //no. of gene templates used for crossover and mutation
#define NMUT 10 //number of mutations
#define PMUT 2 //no. of mutation points

//current input boundaries
#define CTRL_L 0 //lower control limit
#define CTRL_H 50 //high control limit
#define CTRL_S 1 //control step size
#define DUR_L 1 //duration low
#define DUR_H 100 //duration high
#define DUR_S 1 //duration step
#define AMP_L 10 //amplitude low
#define AMP_H 100 //amplitude high
#define AMP_S 10 //amplitude step

void popInitializer( int [][NCTRL] ); //random current generator
int randGen( int, int, int ); //random number generator
void piksr2( int, double [], int [][NCTRL] ); //sorting
void crossOver( int [][NCTRL] ); //crossover
void muTation( int [][NCTRL] ); //mutation
