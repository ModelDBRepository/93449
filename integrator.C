#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "control.h"
#include "integrator.h"

using namespace std;

static int ii=0;
static double dbsInput[NCELL][NSTEP]={0};
static double stngpe[STNGPE]={0};

int integrator( int ctrl[], double cost[] )
{
int i, j, iopt, istate, itask, itol, iwork[ LIW ], jt, liw, lrw, neq, stopTag=0;
double atol[ NEQ ], rtol, rwork[ LRW ], t, tout, y[ NEQ ];
double volt[NCELL][NSTEP]={0}, ctrlInput[NSTEP]={0}, control[NCTRL]={0};

//integrator parameter setting
neq = NEQ;
t = 0.0;
tout = TOUT;
itol = ITOL;
rtol = RTOL;
itask = ITASK;
istate = ISTATE;
iopt = IOPT;
lrw = LRW;
liw = LIW;
jt = JT;
for( i = 0; i < NEQ; i++ ) {
   y[i]=0;
   atol[i]=ATOL;
   if( i%NEQN==0 ) atol[i]=ATOLV;
   }

//STN initial concentrations for v,n,h,r,s,ca
y[0]=-61.37;y[1]=0.42;y[2]=0.098;y[3]=0.086;y[4]=0.6;y[5]=0.16;
y[6]=-56.61;y[7]=0.36;y[8]=0.11;y[9]=0.09;y[10]=0.42;y[11]=0.15;
y[12]=-51.99;y[13]=0.3;y[14]=0.15;y[15]=0.09;y[16]=0.3;y[17]=0.15;
y[18]=-60.72;y[19]=0.46;y[20]=0.09;y[21]=0.099;y[22]=0.75;y[23]=0.15;
y[24]=-59.69;y[25]=0.38;y[26]=0.099;y[27]=0.11;y[28]=0.57;y[29]=0.15;
y[30]=-60.86;y[31]=0.066;y[32]=0.33;y[33]=0.11;y[34]=0.036;y[35]=0.15;
y[36]=-63.57;y[37]=0.084;y[38]=0.29;y[39]=0.11;y[40]=0.036;y[41]=0.15;
y[42]=-55.28;y[43]=0.11;y[44]=0.25;y[45]=0.1;y[46]=0.049;y[47]=0.15;

//GPe initial conditions for v,n,h,r,s,ca
y[48]=-73.58;y[49]=0.19;y[50]=0.69;y[51]=0.58;y[52]=0.14;y[53]=0.16;
y[54]=-76.1;y[55]=0.077;y[56]=0.88;y[57]=0.75;y[58]=0.028;y[59]=0.15;
y[60]=-69.24;y[61]=0.23;y[62]=0.63;y[63]=0.57;y[64]=0.2;y[65]=0.16;
y[66]=-86.2;y[67]=0.088;y[68]=0.85;y[69]=0.71;y[70]=0.04;y[71]=0.15;
y[72]=-65.58;y[73]=0.27;y[74]=0.58;y[75]=0.62;y[76]=0.31;y[77]=0.16;
y[78]=-82.18;y[79]=0.12;y[80]=0.80;y[81]=0.67;y[82]=0.063;y[83]=0.16;
y[84]=-66.91;y[85]=0.33;y[86]=0.51;y[87]=0.69;y[88]=0.46;y[89]=0.16;
y[90]=-77.61;y[91]=0.15;y[92]=0.75;y[93]=0.62;y[94]=0.095;y[95]=0.16;

//GPi initial conditions for v,n,h,r,s,ca
y[96]=-70.96;y[97]=0.16;y[98]=0.79;y[99]=0.72;y[100]=0.13;y[101]=0.16;
y[102]=-68.51;y[103]=0.31;y[104]=0.52;y[105]=0.59;y[106]=0.63;y[107]=0.17;
y[108]=-74.81;y[109]=0.13;y[110]=0.81;y[111]=0.71;y[112]=0.18;y[113]=0.17;
y[114]=-79.29;y[115]=0.85;y[116]=0.22;y[117]=0.53;y[118]=0.89;y[119]=0.17;
y[120]=-78.54;y[121]=0.15;y[122]=0.7;y[123]=0.65;y[124]=0.46;y[125]=0.17;
y[126]=-62.81;y[127]=0.24;y[128]=0.7;y[129]=0.61;y[130]=0.07;y[131]=0.16;
y[132]=-82.04;y[133]=0.12;y[134]=0.75;y[135]=0.63;y[136]=0.34;y[137]=0.17;
y[138]=-67.61;y[139]=0.19;y[140]=0.74;y[141]=0.68;y[142]=0.18;y[143]=0.16;

//TC initial conditions for v, h, r
y[144] = -69.68; y[145] = 1.0; y[146] = 0.00149; 

//current input
for( i=0; i<NCTRL; i++ ) control[i]=ctrl[i];
control[NCTRL-2]=ctrl[NCTRL-2]*TSTEP; //pulse duration
randInput( control, ctrlInput );
                                                                                
//heterogeneity in DBS input
for( i=0; i<NCELL; i++ )
   for( j=0; j<NSTEP; j++ ) dbsInput[i][j]=ctrlInput[j]*hetero( i, DEV_DBS );
//heterogeneity in the network itself
for( i=0; i<NCELL*2; i++ ) stngpe[i]=hetero( i+1, DEV_NET );
srand( time( NULL ) );

//integration
for( i=0; i<NSTEP+NDELAY && stopTag==0; i++ ) {
   ii=i;//parameter transfer
   lsoda_( fex, &neq, y, &t, &tout, &itol, &rtol, atol, &itask, &istate,
         &iopt, rwork, &lrw, iwork, &liw, jdum, &jt ); //integrate
   if( istate==-1 ) istate = 2; //reset istate
   else if( istate!=2 ) {
      stopTag=-1; cerr<<"istate = "<<istate<<", t = "<<t<<endl<<endl;
      }
   else { 
      if( i>=NDELAY ) {
         for( j=0; j<NCELL; j++ )
            volt[j][i-NDELAY]=y[NCELL*NEQN*2+4+NEQN*j];
         }
      tout+=TSTEP;
      }
   }

//cost function calculation
if( stopTag==0 ) cost[0]=costCal( volt );
else cost[0]=1.0e8;

return stopTag;
}


//ODEs
void fex( int *neq, double *t, double *y, double *ydot ) //ODEs
{
int i;

//STN equations
for( i = 0; i < NCELL*NEQN; i += NEQN )
   {
   ydot[i] = (-Il(y[i]) - Ik(y[i], y[i+1]) - Ina(y[i], y[i+2])
             - It(y[i], y[i+3]) - Ica(y[i]) - Iahp(y[i], y[i+5])
             - ( IsynGS(y[i], i, y) - ISTN ) * stngpe[i/NEQN]
             + Idbs(ii,i))/P_cm;//v
   ydot[i+1] = P_phin * (ninf(y[i]) - y[i+1]) / taun(y[i]);//n
   ydot[i+2] = P_phih * (hinf(y[i]) - y[i+2]) / tauh(y[i]);//h
   ydot[i+3] = P_phir * (rinf(y[i]) - y[i+3]) / taur(y[i]);//r
   ydot[i+4] = P_alpha * (1.0-y[i+4]) * Hinf(y[i]-P_thetag) - P_beta*y[i+4];//s
   ydot[i+5] = P_phica*P_eps*(-Ica(y[i])-It(y[i],y[i+3])-P_kca*y[i+5]);//ca
   }

//GPe equations
for( i = NCELL*NEQN; i < NCELL*NEQN*2; i += NEQN )
   {
   ydot[i] = (-Ilg(y[i]) - Ikg(y[i], y[i+1]) - Inag(y[i], y[i+2])
             - Itg(y[i], y[i+3]) - Icag(y[i])
             + ( IGPE - IsynGG(y[i], i, y) ) * stngpe[i/NEQN]
             - Iahpg(y[i],y[i+5])
             -P_gSG*(y[i]-P_vSG)*y[i-NCELL*NEQN+4])/P_cm*stngpe[i/NEQN];//v
   ydot[i+1] = P_phing * (ninfg(y[i]) - y[i+1]) / taung(y[i]);//n
   ydot[i+2] = P_phihg * (hinfg(y[i]) - y[i+2]) / tauhg(y[i]);//h
   ydot[i+3] = P_phirg * (rinfg(y[i]) - y[i+3]) / P_taurg;//r
   ydot[i+4] = P_alphag * (1-y[i+4])*Hinfg(y[i]-P_thetagg) - P_betag*y[i+4];//s
   ydot[i+5] = P_epsg*(-Icag(y[i])-Itg(y[i],y[i+3])-P_kcag*y[i+5]);//ca
   }

//GPi equations
for( i=NCELL*NEQN*2; i<NCELL*NEQN*3; i+=NEQN )
   {
   ydot[i] = ( -Ili(y[i]) - Iki(y[i],y[i+1]) - Inai(y[i],y[i+2])
             - Iti(y[i],y[i+3]) - Icai(y[i]) + IGPI - Iahpi(y[i],y[i+5])
	     - P_gGI*(y[i]-P_vGI)*y[i-NCELL*NEQN+4]
	     - P_gSI*(y[i]-P_vSI)*y[i-NCELL*NEQN*2+4] )/P_cm;//v
   ydot[i+1] = P_phini * (ninfi(y[i]) - y[i+1]) / tauni(y[i]);//n
   ydot[i+2] = P_phihi * (hinfi(y[i]) - y[i+2]) / tauhi(y[i]);//h
   ydot[i+3] = P_phiri * (rinfi(y[i]) - y[i+3]) / P_tauri;//r
   ydot[i+4] = P_alphai * (1-y[i+4])*Hinfi(y[i]-P_thetagi) - P_betai*y[i+4];//s
   ydot[i+5] = P_epsi*(-Icai(y[i])-Iti(y[i],y[i+3])-P_kcai*y[i+5]);//ca
   }

//TC1 equations
for( i=NCELL*NEQN*3; i<NCELL*NEQN*3+NTCEQN; i+=NTCEQN )
   {
   ydot[i] = ( -Ilt(y[i]) - Ikt(y[i],y[i+1]) - Inat(y[i],y[i+1])
               - Itt(y[i],y[i+2]) - IsynIT1(y,y[i]) + Ism(*t) )/P_cm; //v
   ydot[i+1] = (hinft(y[i]) - y[i+1]) / tauht(y[i]); //h
   ydot[i+2] = (rinft(y[i]) - y[i+2]) / taurt(y[i]); //r
   }
                                                                                
//TC2 equations
for( i=NCELL*NEQN*3+NTCEQN; i<NEQ; i+=NTCEQN )
   {
   ydot[i] = ( -Ilt(y[i]) - Ikt(y[i],y[i+1]) - Inat(y[i],y[i+1])
               - Itt(y[i],y[i+2]) - IsynIT2(y,y[i]) + Ism(*t) )/P_cm; //v
   ydot[i+1] = (hinft(y[i]) - y[i+1]) / tauht(y[i]); //h
   ydot[i+2] = (rinft(y[i]) - y[i+2]) / taurt(y[i]); //r
   }
}


//random inhibition of GPe to STN
double IsynGS(const double volt, int i, const double yval[])
{
int a = 0, b = 0;

if( i == 0 ) {a = 3; b = 7;}
else if( i == NEQN ) {a = 4; b = 8;}
else if( i == NEQN*2 ) {a = 5; b = 1;}
else if( i == NEQN*3 ) {a = 6; b = 2;}
else if( i == NEQN*4 ) {a = 7; b = 3;}
else if( i == NEQN*5 ) {a = 8; b = 4;}
else if( i == NEQN*6 ) {a = 1; b = 5;}
else if( i == NEQN*7 ) {a = 2; b = 6;}
else cerr << "i = " << i << endl;
a -= 1; b -= 1;

return P_gGS*(volt-P_vGS)*(yval[(NCELL+a)*NEQN+4] + yval[(NCELL+b)*NEQN+4]);
}


//neighboring inhibition of GPe to GPe
double IsynGG( const double volt, int index, const double yval[ ] )
{
int a = 0, b = 0;
//cerr << "i = " << index << endl;
if( index == NCELL*NEQN ) {a = 2; b = 8;}
else if( index == (NCELL+1)*NEQN ) {a = 1; b = 3;}
else if( index == (NCELL+2)*NEQN ) {a = 2; b = 4;}
else if( index == (NCELL+3)*NEQN ) {a = 3; b = 5;}
else if( index == (NCELL+4)*NEQN ) {a = 4; b = 6;}
else if( index == (NCELL+5)*NEQN ) {a = 5; b = 7;}
else if( index == (NCELL+6)*NEQN ) {a = 6; b = 8;}
else if( index == (NCELL+7)*NEQN ) {a = 7; b = 1;}
else cerr << "index = " << index << endl;
a -= 1; b -= 1;
return P_gGG*(volt-P_vGG)*(yval[(NCELL+a)*NEQN+4]+yval[(NCELL+b)*NEQN+4]);
}

//synaptic current from GPi to TC1
                                                                                
double IsynIT1( const double yval[], const double v )
{
return P_gIT * (v-P_vIT) * (yval[NCELL*NEQN*2+4]+yval[NCELL*NEQN*2+4+NEQN]
                +yval[NCELL*NEQN*2+4+NEQN*4]+yval[NCELL*NEQN*2+4+NEQN*5]);
}
                                                                                
//synaptic current from GPi to TC2
                                                                                
double IsynIT2( const double yval[], const double v )
{
return P_gIT * (v-P_vIT) * (yval[NCELL*NEQN*2+4+NEQN*2]
       +yval[NCELL*NEQN*2+4+NEQN*3]+yval[NCELL*NEQN*2+4+NEQN*6]
       +yval[NCELL*NEQN*2+4+NEQN*7]);
}

//dummy function

void jdum( int *neq, double *t, double *y, int *ml, int *mu, double *pd,
           int *nrowpd )
{
}

//STN currents

double Il( const double v )
{return P_gl * (v - P_vl);}

double Ik( const double v, const double n )
{return P_gk * pow(n, 4) * (v - P_vk);}

double Ina( const double v, const double h )
{return P_gna * pow(minf(v), 3) * h * (v - P_vna);}

double It( const double v, const double r )
{return P_gt * pow(ainf(v), 3) * pow(binf(r), 2) * (v - P_vca);}

double Ica( const double v )
{return P_gca * pow(sinf(v), 2) * (v - P_vca);}

double Iahp( const double v, const double ca )
{return P_gahp * (v - P_vk) * ca / ( ca + P_k1);}

//STN functions

double minf( const double v )
{return 1.0/(1.0 + exp(-(v-P_thetam)/P_sigmam));}

double ainf( const double v )
{return 1.0/(1.0+exp(-(v-P_thetaa)/P_sigmaa));}

double sinf( const double v )
{return 1.0/(1.0+exp(-(v-P_thetas)/P_sigmas));}

double ninf( const double v )
{return 1.0/(1.0+exp(-(v-P_thetan)/P_sigman));}

double hinf( const double v )
{return 1.0/(1.0+exp(-(v-P_thetah)/P_sigmah));}

double rinf( const double v )
{return 1.0/(1.0+exp(-(v-P_thetar)/P_sigmar));}

double binf( const double r )
{return 1.0/(1.0+exp((r-P_thetab)/P_sigmab))-1.0/(1.0+exp(-P_thetab/P_sigmab));}

double Hinf( const double v )
{return 1.0/(1.0+exp(-(v-P_thetaH)/P_sigmaH));}

double taun( const double v )
{return P_taun0 + P_taun1/(1.0+exp(-(v-P_thn)/P_sigmant));}

double tauh( const double v )
{return P_tauh0 + P_tauh1/(1.0+exp(-(v-P_thh)/P_sigmaht));}

double taur( const double v )
{return P_taur0 + P_taur1/(1.0+exp(-(v-P_thr)/P_sigmart));}


//GPe currents

double Ilg( const double v )
{return P_glg * (v - P_vlg);}

double Ikg( const double v, const double n )
{return P_gkg * pow(n, 4) * (v - P_vkg);}
                                                                                
double Inag( const double v, const double h )
{return P_gnag * pow(minfg(v), 3) * h * (v - P_vnag);}
                                                                                
double Itg( const double v, const double r )
{return P_gtg * pow(ainfg(v), 3) * r * (v - P_vcag);}
                                                                                
double Icag( const double v )
{return P_gcag * pow(sinfg(v), 2) * (v - P_vcag);}

double Iahpg( const double v, const double ca )
{return P_gahpg * (v - P_vkg) * ca / ( ca + P_k1g);}

//GPe functions

double minfg( const double v )
{return 1.0/(1.0 + exp(-(v-P_thetamg)/P_sigmamg));}
                                                                                
double ainfg( const double v )
{return 1.0/(1.0+exp(-(v-P_thetaag)/P_sigmaag));}
                                                                                
double sinfg( const double v )
{return 1.0/(1.0+exp(-(v-P_thetasg)/P_sigmasg));}
                                                                                
double ninfg( const double v )
{return 1.0/(1.0+exp(-(v-P_thetang)/P_sigmang));}
                                                                                
double hinfg( const double v )
{return 1.0/(1.0+exp(-(v-P_thetahg)/P_sigmahg));}
                                                                                
double rinfg( const double v )
{return 1.0/(1.0+exp(-(v-P_thetarg)/P_sigmarg));}
                                                                                
double Hinfg( const double v )
{return 1.0/(1.0+exp(-(v-P_thetaHg)/P_sigmaHg));}

double taung( const double v )
{return P_taun0g + P_taun1g/(1.0+exp(-(v-P_thng)/P_sigmantg));}
                                                                                
double tauhg( const double v )
{return P_tauh0g + P_tauh1g/(1.0+exp(-(v-P_thhg)/P_sigmahtg));}

//GPi currents

double Ili( const double v )
{return P_gli * (v - P_vli);}

double Iki( const double v, const double n )
{return P_gki * pow(n, 4) * (v - P_vki);}

double Inai( const double v, const double h )
{return P_gnai * pow(minfi(v), 3) * h * (v - P_vnai);}

double Iti( const double v, const double r )
{return P_gti * pow(ainfi(v), 3) * r * (v - P_vcai);}

double Icai( const double v )
{return P_gcai * pow(sinfi(v), 2) * (v - P_vcai);}

double Iahpi( const double v, const double ca )
{return P_gahpi * (v - P_vki) * ca / ( ca + P_k1i);}

//GPi functions

double minfi( const double v )
{return 1.0/(1.0 + exp(-(v-P_thetami)/P_sigmami));}

double ainfi( const double v )
{return 1.0/(1.0+exp(-(v-P_thetaai)/P_sigmaai));}

double sinfi( const double v )
{return 1.0/(1.0+exp(-(v-P_thetasi)/P_sigmasi));}

double ninfi( const double v )
{return 1.0/(1.0+exp(-(v-P_thetani)/P_sigmani));}

double hinfi( const double v )
{return 1.0/(1.0+exp(-(v-P_thetahi)/P_sigmahi));}

double rinfi( const double v )
{return 1.0/(1.0+exp(-(v-P_thetari)/P_sigmari));}

double Hinfi( const double v )
{return 1.0/(1.0+exp(-(v-P_thetaHi)/P_sigmaHi));}

double tauni( const double v )
{return P_taun0i + P_taun1i/(1.0+exp(-(v-P_thni)/P_sigmanti));}

double tauhi( const double v )
{return P_tauh0i + P_tauh1i/(1.0+exp(-(v-P_thhi)/P_sigmahti));}

//TC currents
                                                                                
double Ilt( const double v )
{return P_glt * (v - P_vlt);}
                                                                                
double Ikt( const double v, const double h )
{return P_gkt * pow((0.75*(1.0-h)), 4) * (v - P_vkt);}
                                                                                
double Inat( const double v, const double h )
{return P_gnat * pow(minft(v), 3) * h * (v - P_vnat);}
                                                                                
double Itt( const double v, const double r )
{return P_gtt * pow(ainft(v), 2) * r * (v - P_vcat);}

//TC functions
                                                                                
double minft( const double v )
{return 1.0/(1.0 + exp(-(v-P_thetamt)/P_sigmamt));}
                                                                                
double ainft( const double v )
{return 1.0/(1.0 + exp(-(v-P_thetaat)/P_sigmaat));}
                                                                                
double hinft( const double v )
{return 1.0/(1.0 + exp(-(v-P_thetaht)/P_sigmaht));}
                                                                                
double rinft( const double v )
{return 1.0/(1.0 + exp(-(v-P_thetart)/P_sigmart));}
                                                                                
double tauht( const double v )
{return 1.0/(0.128*exp(-(v+46.0)/18.0)+4.0/(1.0+exp(-(v+23.0)/5.0)));}
                                                                                
double taurt( const double v )
{return 28.0+exp(-(v+25.0)/10.5);}
                                                                                
double Heaviside( const double var )
{
if( var < 0 ) return 0.0;
else return 1.0;
}

//other current functions
double Idbs( const unsigned long int cnt, const int index ) //DBS current to STN
{
if( cnt<NDELAY ) 
   return 0.0;
else
   return dbsInput[index/NEQN][cnt-NDELAY];
}

//generate random DBS current
void randInput( const double pdf[], double current[] )
{
int i, j;
                                                                                
srand( 1 );
for( i=0; i<NSTEP; i++ ) current[i]=0; //re-initialization
i=(int)(randGen2(pdf,NCTRL,ILOW,ISTEP)/TSTEP); //first pulse time
for( ; i<NSTEP; i+=(int)(randGen2(pdf,NCTRL,ILOW,ISTEP)/TSTEP) ) {
   for( j=0; j<(int)(pdf[NCTRL-2]/TSTEP) && i+j<NSTEP; j++ ) //set pulse width
      current[i+j]=pdf[NCTRL-1]; //set pulse amplitude
   }
srand( time( NULL ) );
}

//random number generator for specified probability distribution function
double randGen2( const double dist[], const int size, const double low,
                 const double step )
{
int i;
double total=0, y=0;
                                                                                
for( i=0; i<size-2; i++ ) total+=dist[i]; //NOTE: size-2!
y=total*rand()/(RAND_MAX+1.0);
total=0;
for( i=0; y>total; i++ ) total+=dist[i];
return low+step*i-(total-y)/dist[i-1]*step;
//return low+step*i; //if one PDF parameter only specifies one frequency
}

//cost function value calculation
double costCal( const double arr[][NSTEP] )
{
int i, j, k, m=0, n1, n2, nspk[NCELL]={0};
int autoCnt[NCELL][NBIN]={0}, crossCnt[COMBI][NBIN]={0};
double tspk[NCELL][NSPK]={0}, autoStd=0, crossStd=0, total=0;
float autoConv[NCELL][DURATION2]={0}, crossConv[COMBI][DURATION2]={0};
                                                                                
//count Gpi spike number and time
for( i=0; i<NCELL; i++ ) nspk[i]=spkCount( arr[i], tspk[i] );
                                                                                
//calculate spike time distance autocorrelation and convolution
for( i=0; i<NCELL; i++ ) autoCal( tspk[i], autoCnt[i] ); //autocorrelation
for( i=0; i<NCELL; i++ ) convCal( autoCnt[i], autoConv[i] ); //convolution
for( i=0; i<NCELL; i++ ) autoStd+=stdCal( autoConv[i] ); //standard deviation
                                                                                
//calculate spike time distance crosscorrelation and convolution
for( n1=0; n1<1; n1++ ) //crosscorrelation with Gpi1
   for( n2=1; n2<NCELL; n2++ ) { //for all other Gpi cells
      for( i=0; i<NSPK; i++ ) {
         if( tspk[n1][i]>NSTART*TSTEP && tspk[n1][i]<(NEND)*TSTEP ) {
            for( j=0; j<NSPK; j++ ) {
               if(tspk[n2][j]>0&&fabs(tspk[n1][i]-tspk[n2][j])<DURATION*TSTEP){
                  k=(int)(fabs((tspk[n1][i]-tspk[n2][j])/BIN_RESLN/TSTEP));
                  crossCnt[m][k]++;
                  if(fabs(tspk[n1][i]-tspk[n2][j])<TSTEP) total+=1;
                  }
               }
            }
         }
      m++;
      }
for( i=0; i<COMBI; i++ ) convCal( crossCnt[i], crossConv[i] ); //convolution
for( i=0; i<COMBI; i++ ) crossStd+= stdCal( crossConv[i] ); //standard deviation

//cost function calculation
if( total==0 ) return 1.0e8;
else if( total<100 ) return 1.0e8/total;
else return autoStd*W1+crossStd*W2+total*W3;
}


                                                                         
//count the number of spikes for Gpi synaptic output
int spkCount( const double spk[], double tispk[] )
{
int i, flag=-1, cnt=0;
double temp=spk[0];
                                                                                
for( i=0; i<NSPK; i++ ) tispk[i]=-1.0;
for( i=1; i<NSTEP; i++ ) {
   if( spk[i]<VTHLD ) flag=-1;
   else if( spk[i]>temp ) flag=1;
   else if( flag==1 ) { flag=2; tispk[cnt]=TSTEP*i; cnt++; }
   temp=spk[i];
   }
if( cnt>NSPK ) cerr << "Gpi spike no. > NSPK" << endl;
return cnt;
}


//calculate autocorrelation of Gpi synaptic output spike distance
void autoCal( const double tispk[], int count[] )
{
int i, j, k;
for( i=0; i<NSPK; i++ ) {
   if( tispk[i]>NSTART*TSTEP && tispk[i]<(NEND)*TSTEP ) {
      for( j=0; j<NSPK; j++ ) {
         if( tispk[j]>0 && i!=j && fabs(tispk[j]-tispk[i])<DURATION*TSTEP ) {
            k=(int)(fabs((tispk[j]-tispk[i])/BIN_RESLN/TSTEP));
            count[k]++;
            }
         }
      }
   }
}


//convolution calculation
void convCal( const int before[], float after[] )
{
int i;
float response[DURATION]={0}, signal[DURATION]={0};
                                                                                
//set up the Gaussian response function
for( i=0; i<(MSTEP+1)/2; i++ )
   response[i]=exp(-i*i/2.0/(float)(WIDTH*WIDTH));
for( i=(MSTEP+1)/2; i<MSTEP; i++ ) response[i]=response[MSTEP-i];
                                                                                
//reorganize the autocorrelation data
for( i=0; i<NBIN; i++ ) signal[i]=(float)(before[i]);
//for( i=0; i<DURATION-(MSTEP+1)/2; i++ ) signal[i]=(float)(before[i]);
//for( i=DURATION-(MSTEP+1)/2; i<DURATION; i++ ) signal[i]=0.0;
                                                                                
//convolution and output
convlv( signal-1, DURATION, response-1, MSTEP, 1, after-1 );
for( i=0; i<DURATION; i++ ) //normalization
   after[i]=after[i]/WIDTH/sqrt(2.0*PI);
}
                                                                                
                                                                                
                                                                                
//standard deviation calculation of auto- and cross- correlations
double stdCal( const float data[] )
{
int i;
float total=0, std=0;
for( i=0; i<NBIN; i++ ) total+=data[i];
total/=NBIN;
for( i=0; i<NBIN; i++ )
   std+=(data[i]-total)*(data[i]-total);
return sqrt(std/NBIN)/total;
}


double Ism( const double time ) //Real Ism
{
return ISM * Heaviside( sin(2.0*PI*time/PSM) )
       * ( 1.0 - Heaviside( sin(2.0*PI*(time+DSM)/PSM) ) );
}
                                                                                
double Ism2( const double time ) //Ism with longer duration, used in peak count
{
return ISM * Heaviside( sin(2.0*PI*(time-SHIFT)/PSM) )
       * (1.0 - Heaviside( sin(2.0*PI*(time+DSM)/PSM) ) );
}
                                                                                
//random number generator
double hetero( const int seed, const double dev )
{
srand( seed );
return 1.0-dev+dev*2.0*rand()/(RAND_MAX+1.0);
}
