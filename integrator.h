//general parameters 
#define PI 3.1415926
#define NCELL 8 //number STN, GPe, or GPi cells
#define STNGPE NCELL*2 //sum of STN and GPe cells
#define NEQN 6 //number of ODEs per cell
#define NTCCELL 2 //number of TC cells
#define NTCEQN 3 //number of ODEs per TC cell
#define P_cm 1.0

//STN parameters
#define P_phin 0.75
#define P_phih 0.75
#define P_phica 0.75 //0.75 in Terman's code, not available from his paper.
#define P_eps 3.75e-5 //5.0e-5 in Terman's code
#define P_kca 22.5
#define P_phir 0.2
#define P_alpha 1.0 
#define P_thetag 30.0
#define P_beta 0.05 
#define P_gGS 0.9 
#define P_vGS -100.0 
#define P_gl 2.25
#define P_vl -60.0
#define P_gk 45.0
#define P_vk -80.0
#define P_gna 37.5
#define P_vna 55.0
#define P_gt 0.5
#define P_vca 140.0
#define P_gca 0.5
#define P_gahp 9.0
#define P_k1 15.0
#define P_thetam -30.0
#define P_sigmam 15.0
#define P_thetaa -63.0
#define P_sigmaa 7.8
#define P_thetas -39.0
#define P_sigmas 8.0
#define P_thetan -32.0
#define P_sigman 8.0
#define P_thetah -39.0
#define P_sigmah -3.1
#define P_thetar -67.0
#define P_sigmar -2.0
#define P_thetab 0.4
#define P_sigmab -0.1
#define P_thetaH -39.0
#define P_sigmaH 8.0
#define P_taun0 1.0
#define P_taun1 100.0
#define P_thn -80.0
#define P_sigmant -26.0
#define P_tauh0 1.0
#define P_tauh1 500.0
#define P_thh -57.0
#define P_sigmaht -3.0
#define P_taur0 40.0
#define P_taur1 17.5
#define P_thr 68.0
#define P_sigmart -2.2
#define ISTN 25.0 

//GPe parameters
#define P_phing 0.1
#define P_phihg 0.05
#define P_phirg 1.0
#define P_alphag 1.0 
#define P_thetagg 20.0
#define P_betag 0.1 
#define P_epsg 1.0e-4 
#define P_kcag 15.0 
#define P_gSG 0.6 //0.3 in Terman 2004
#define P_vSG 0.0 
#define P_gGG 0.0005 //normal = 2.0, PD = 0.05 or 0.0005
#define P_vGG -80.0 
#define P_glg 0.1
#define P_vlg -55.0
#define P_gkg 30.0
#define P_vkg -80.0
#define P_gnag 120.0
#define P_vnag 55.0
#define P_gtg 0.5
#define P_vcag 120.0
#define P_gcag 0.15
#define P_gahpg 30.0
#define P_k1g 30.0
#define P_thetamg -37.0
#define P_sigmamg 10.0
#define P_thetaag -57.0
#define P_sigmaag 2.0
#define P_thetasg -35.0
#define P_sigmasg 2.0
#define P_thetang -50.0
#define P_sigmang 14.0
#define P_thetahg -58.0
#define P_sigmahg -12.0
#define P_thetarg -70.0
#define P_sigmarg -2.0
#define P_thetaHg -57.0
#define P_sigmaHg 2.0
#define P_taun0g 0.05
#define P_taun1g 0.27
#define P_thng -40.0
#define P_sigmantg -12.0
#define P_tauh0g 0.05
#define P_tauh1g 0.27
#define P_thhg -40.0
#define P_sigmahtg -12.0
#define P_taurg 30.0
#define IGPE -13.3 //normal = -2.0, PD = -13.0

//GPi parameters
#define P_phini 0.1
#define P_phihi 0.05
#define P_phiri 1.0
#define P_alphai 2.0
#define P_thetagi 20.0
#define P_betai 0.08
#define P_epsi 1.0e-4
#define P_kcai 15.0
#define P_gSI 0.2 //0.3 in Terman 2004
#define P_vSI 0.0 
#define P_gGI 1.7 //1.0 in Terman 2004
#define P_vGI -100.0 
#define P_gli 0.1
#define P_vli -55.0
#define P_gki 30.0
#define P_vki -80.0
#define P_gnai 120.0
#define P_vnai 55.0
#define P_gti 0.5
#define P_vcai 120.0
#define P_gcai 0.15
#define P_gahpi 30.0
#define P_k1i 30.0
#define P_thetami -37.0
#define P_sigmami 10.0
#define P_thetaai -57.0
#define P_sigmaai 2.0
#define P_thetasi -35.0
#define P_sigmasi 2.0
#define P_thetani -50.0
#define P_sigmani 14.0
#define P_thetahi -58.0
#define P_sigmahi -12.0
#define P_thetari -70.0
#define P_sigmari -2.0
#define P_thetaHi -57.0
#define P_sigmaHi 2.0
#define P_taun0i 0.05
#define P_taun1i 0.27
#define P_thni -40.0
#define P_sigmanti -12.0
#define P_tauh0i 0.05
#define P_tauh1i 0.27
#define P_thhi -40.0
#define P_sigmahti -12.0
#define P_tauri 30.0
#define IGPI 3.0 

//TC parameters (obtained from Terman 2004)
#define P_gIT 0.08 //0.06 in Terman 2004
#define P_vIT -85.0
#define P_glt 0.05
#define P_vlt -70.0
#define P_gkt 5.0
#define P_vkt -90.0
#define P_gnat 3.0
#define P_vnat 50.0
#define P_gtt 5.0
#define P_vcat 0.0
#define P_thetamt -37.0
#define P_sigmamt 7.0
#define P_thetaat -60.0
#define P_sigmaat 6.2
#define P_thetaht -41.0
#define P_sigmahtt -4.0
#define P_thetart -84.0
#define P_sigmartt -4.0
#define ISM 10.0 //5.0 in Terman 2004, 7.5-10 may be good
#define PSM 50.0 //25.0 in Terman 2004, 50 in his figures
#define DSM 5.0
#define SHIFT 2.5 //phase shift when counting spike number
#define THLD_TC -30.0 // threshold for TC voltage spikes

//integrator parameters
#define NEQ (NCELL*NEQN*3+NTCCELL*NTCEQN) //number of ODEs in total
#define LRW ( 22 + NEQ * ( NEQ + 9 ) )
#define LIW ( 20 + NEQ )
#define ITOL 2
#define RTOL 1.0e-10 
#define ATOL 1.0e-12 
#define ATOLV 1.0e-10 //absolute tolerance for v
#define ITASK 1
#define ISTATE 1
#define IOPT 0
#define JT 2
#define TOUT 0.1 //first output time
#define TSTEP 0.1 //output step size
#define NSTEP 65536
#define NDELAY 25000 //No. of steps before the control is applied

//current input parameters
#define ILOW 1.0 //lowest current pulse distance
#define ISTEP 5.0 //step size for increasing pulse distance
#define DEV_DBS 0.1 //white noise of I_{DBS}
#define DEV_NET 0.01 //heterogeneity in the network itself (0 to 1)

//parameters for counting spike time
#define NSPK 2000 //maximum no. of Gpi spikes
#define VTHLD 0.5 //Gpi synaptic output threshold for counting spikes
                                                                                
//parameters for calculating auto- and cross-correlations
#define NSTART 10000 //when to start calculating correlations (in steps)
#define NEND 55000 //when to end calculating correlations (in steps)
#define DURATION 8192 //duration for correlation calculation (in steps)
#define BIN_RESLN 100 //binning resolution of Gpi spikes (in steps)
#define NBIN DURATION/BIN_RESLN+1 //maximum no. of Gpi bins
#define COMBI NCELL-1 //no. of comparisons for crosscorrelation
                                                                                
//parameters for convolution
#define DURATION2 DURATION*2
#define MSTEP 21 //width of the response function (in bins)
#define WIDTH 1 //Gaussian width (in bins)

//parameters for cost function calculation
#define W1 200.0 //1.0 //weight for autocorrelation term
#define W2 100.0 //weight for crosscorrelation term
#define W3 0.0002 //weight for Gpi spike number in first bin term
#define W4 0 //0.1 //current cost weight

//function declaration

int integrator( int [], double [] ); //my integrator

extern "C" void lsoda_( void( int *, double *, double *, double * ), int *,
double *, double *, double *, int *, double *, double *, int *, int *, int *,
double *, int *, int *, int *, void( int *, double *, double *, int *, int *,
double *, int * ), int * ); //ODE integrator

void fex( int *, double *, double *, double * ); //ODE container

void jdum( int *, double *, double *, int *, int *, double *, int * ); //dummy

//STN currents
double Il( const double v );
double Ik( const double v, const double n );
double Ina( const double v, const double h );
double It( const double v, const double r );
double Ica( const double v );
double Iahp( const double v, const double ca );

//STN functions
double minf( const double v );
double ainf( const double v );
double sinf( const double v );
double ninf( const double v );
double hinf( const double v );
double rinf( const double v );
double binf( const double r );
double Hinf( const double v );
double taun( const double v );
double tauh( const double v );
double taur( const double v );

//GPe currents
double Ilg( const double v );
double Ikg( const double v, const double n );
double Inag( const double v, const double h );
double Itg( const double v, const double r );
double Icag( const double v );
double Iahpg( const double v, const double ca );

//GPe functions
double minfg( const double v );
double ainfg( const double v );
double sinfg( const double v );
double ninfg( const double v );
double hinfg( const double v );
double rinfg( const double v );
double Hinfg( const double v );
double taung( const double v );
double tauhg( const double v );

//GPi currents
double Ili( const double v );
double Iki( const double v, const double n );
double Inai( const double v, const double h );
double Iti( const double v, const double r );
double Icai( const double v );
double Iahpi( const double v, const double ca );

//GPi functions
double minfi( const double v );
double ainfi( const double v );
double sinfi( const double v );
double ninfi( const double v );
double hinfi( const double v );
double rinfi( const double v );
double Hinfi( const double v );
double tauni( const double v );
double tauhi( const double v );

//TC currents
double Ilt( const double v );
double Ikt( const double v, const double h );
double Inat( const double v, const double h );
double Itt( const double v, const double r );
double taurt( const double v );
double tauht( const double v );

//TC functions
double minft( const double v );
double ainft( const double v );
double hinft( const double v );
double rinft( const double v );
double Heaviside( const double var );

//other functions
//inhibition of GPe to STN
double IsynGS(const double, int, const double [ ]); 
//neighboring inhibition of GPe to GPe
double IsynGG(const double, int, const double [ ] ); 
//synaptic current from GPi to TC
double IsynIT1( const double [], const double );
double IsynIT2( const double [], const double );
//DBS current to STN
double Idbs( const unsigned long int, const int ); 
//cost value calculation
//Ism current to TC
double Ism( const double );
//Ism current to TC with delay, used in peak counting
double Ism2( const double );
//heterogeneity / random number generator
double hetero( const int, const double );
//random square pulse DBS input generated from specified PDF
void randInput( const double [], double [] );
//random number generator for specified probability distribution function
double randGen2( const double [], const int, const double, const double );
//auto- and cross-correlation calculation
double costCal( const double [][NSTEP] );
//count the number of spikes
int spkCount ( const double [], double [] );
//calculate autocorrelation
void autoCal( const double [], int [] );
//convolution calculation
void convCal( const int [], float [] );
//calculate standard deviation
double stdCal ( const float[] );
                                                                                
void four1( float [], unsigned long int, int ); //fft
void realft( float [], unsigned long, int );
void twofft( float [], float [], float [], float [], unsigned long );
void convlv( float [], unsigned long, float [], unsigned long, int, float [] );
