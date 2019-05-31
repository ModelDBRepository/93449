//Note: 131 Ism spikes,
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>
#include "control.h"
#include "integrator.h"

using namespace std;

int main( int argc, char** argv )
{

int myRank, poolSize;
MPI_Init( &argc, &argv );
MPI_Comm_size( MPI_COMM_WORLD, &poolSize );
MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

srand( time( NULL ) );

if( myRank == 0 ) { //if I am the master
   MPI_Status statusMaster;
   int i, j, destination, tagSent, tagReceived, counter;
   int ctrlGenome[POPSIZE][NCTRL]={0};
   int geneSent[NCTRL]={0};
   double resultReceived[1]={0};
   double score[POPSIZE]={0};
   double temp=0;

   //input and output files
   ofstream ctrlFile( "process.dat", ios::out );
   ofstream scoreFile( "score.dat", ios::out );
   if( !scoreFile || !ctrlFile ) {
      cerr << "Output files can't be created"; exit(1);
      }

   //genome evolution 
   popInitializer( ctrlGenome ); //genome initialization
   scoreFile << "#Generation & minimum score & mean score" << endl;

   counter=1;
   while( counter<=MAXGEN ) {
      //send out one genome to each slave
      for( destination=1; destination<=POPSIZE; destination++ ) {
         tagSent=destination-1;
         for( i=0; i<NCTRL; i++ ) geneSent[i] = ctrlGenome[tagSent][i];
         MPI_Send( geneSent, NCTRL, MPI_INT, destination, tagSent,
                   MPI_COMM_WORLD );
         }

      //get the result from each slave 
      for( i=1; i<=POPSIZE; i++ ) {
         MPI_Recv( resultReceived, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 
		   MPI_ANY_TAG, MPI_COMM_WORLD, &statusMaster );
         tagReceived = statusMaster.MPI_TAG;
         score[tagReceived]=resultReceived[0];//NOTE: score minimization
         }

      //sort the scores in ascending order
      piksr2( POPSIZE, score, ctrlGenome ); 

      //output the scores
      ctrlFile << "Generation no. " << counter << endl;
      for( i=0; i<POPSIZE; i++ ) {
         ctrlFile << "Control: ";
         for( j=0; j<NCTRL; j++ ) ctrlFile << ctrlGenome[i][j] << " ";
         ctrlFile << endl << "Score = " << score[i] << endl;
         }
      scoreFile << counter << " " << score[0] << " ";//output the best score
      for( i=0; i<POPSIZE; i++ ) temp+=score[i]/POPSIZE;
      scoreFile << temp << endl; //output the mean score
      temp = 0;

      //crossover and mutation
      crossOver( ctrlGenome );
      muTation( ctrlGenome );

      counter++; //proceed to the next generation
      }
   ctrlFile.close();
   scoreFile.close();
   }

else { //if I am not the master
   MPI_Status statusSlave;
   int i, j, tagReceived;
   int geneReceived[NCTRL]={0}; 
   double resultSent[1]={0};

   for( i=1; i<=MAXGEN; i++ ) {
      MPI_Recv( geneReceived, NCTRL, MPI_INT, 0, MPI_ANY_TAG,
                MPI_COMM_WORLD, &statusSlave ); //receive a genome
      tagReceived = statusSlave.MPI_TAG;

      //genome evaluation
      if( integrator( geneReceived, resultSent ) == 0 ) 
         MPI_Send( resultSent, 1, MPI_DOUBLE, 0, tagReceived, MPI_COMM_WORLD ); 
      else {
         cerr << "Integration failed with control" << endl;
	 for( j=0; j<NCTRL; j++ ) cerr << geneReceived[j] << " ";
	 cerr << endl << endl;
         resultSent[0]=0;
         MPI_Send( resultSent, 1, MPI_DOUBLE, 0, tagReceived, MPI_COMM_WORLD ); 
	 }
      }
   }

MPI_Finalize( );

return 0;
}


//random control current generator
void popInitializer( int ctrl[][NCTRL] )
{
int i, j;
                                                                                
for( i=0; i<POPSIZE; i++ ) {
   for( j=0; j<NCTRL-2; j++ ) //set distribution function
      ctrl[i][j] = randGen( CTRL_L, CTRL_H, CTRL_S );
   ctrl[i][NCTRL-2]=randGen( DUR_L, DUR_H, DUR_S ); //set pulse duration
   ctrl[i][NCTRL-1]=randGen( AMP_L, AMP_H, AMP_S ); //set pulse amplitude
   }
}
                                                                                
                                                                                
//random number generator
int randGen( int low, int high, int step )
{
int nstep, rstep;
                                                                                
nstep = int( ( high - low ) / step ) + 1;
                                                                                
rstep = int( double( nstep ) * rand( ) / ( RAND_MAX + 1.0 ) );
                                                                                
return ( low + rstep * step );
}

//sorting
void piksr2( int n, double arr[ ], int brr[ ][ NCTRL ] )
{
int i, j, k, b[ NCTRL ] = { 0 };
double a;
                                                                                
for( j = 1; j < n; j++ ) {
   a = arr[ j ];
   for( k = 0; k < NCTRL; k++ )
      b[ k ] = brr[ j ][ k ];
   i = j - 1;
   while( i >= 0 && arr[ i ] > a ) {
      arr[ i + 1 ] = arr[ i ];
      for( k = 0; k < NCTRL; k++ )
         brr[ i + 1 ][ k ] = brr[ i ][ k ];
      i--;
      }
   arr[ i + 1 ] = a;
   for( k = 0; k < NCTRL; k++ )
      brr[ i + 1 ][ k ] = b[ k ];
   }
}


//crossover
void crossOver( int arr[ ][ NCTRL ] )
{
int a, b, c, d, i, j, temp;
int brr[POPSIZE][NCTRL];
                                                                                
for( i=0; i<POPSIZE; i++ )
   for( j=0; j<NCTRL; j++ )
      brr[i][j]=0;
                                                                                
//crossover
for( i=NSAVED; i<POPSIZE-NMUT; i++ ) {
   a = int( double( NCTRL ) * rand( ) / ( RAND_MAX + 1.0 ) );
   b = int( double( NCTRL ) * rand( ) / ( RAND_MAX + 1.0 ) );
   c = int( double( NUSED ) * rand( ) / ( RAND_MAX + 1.0 ) );
   d = int( double( NUSED ) * rand( ) / ( RAND_MAX + 1.0 ) );
                                                                                
   while( a==b ) a = int( double( NCTRL ) * rand( ) / ( RAND_MAX + 1.0 ) );
   if( a > b ) {
      temp = a;
      a = b;
      b = temp;
      }
   while( c==d ) c = int( double( NUSED ) * rand( ) / (RAND_MAX+1.0) );
                                                                                
   for( j=0; j<NCTRL; j++ ) brr[i][j] = arr[c][j];
   for( j=a; j<b; j++ ) brr[i][j] = arr[d][j];
   }
                                                                                
//put crossover result back to arr[][]
for( i=NSAVED; i<POPSIZE-NMUT; i++ )
   for( j=0; j<NCTRL; j++ )
      arr[i][j]=brr[i][j];
}


//mutation
void muTation( int arr[ ][ NCTRL ] )
{
int a=0, b=0, i, j;
                                                                                
//two point mutation
for( i=POPSIZE-NMUT; i<POPSIZE; i++ ) {
   a = int( double( NUSED ) * rand( ) / ( RAND_MAX + 1.0 ) );
   for( j=0; j<NCTRL; j++ ) arr[i][j]=arr[a][j];
   for( j=0; j<PMUT; j++ ) {
      b = int( double( NCTRL ) * rand( ) / ( RAND_MAX + 1.0 ) );
      if( b==(NCTRL-2) ) arr[i][b]=randGen( DUR_L, DUR_H, DUR_S ); //width
      else if( b==(NCTRL-1) ) arr[i][b]=randGen( AMP_L, AMP_H, AMP_S ); //amp
      else arr[i][b] = randGen( CTRL_L, CTRL_H, CTRL_S ); //PDF
      }
   }
}
