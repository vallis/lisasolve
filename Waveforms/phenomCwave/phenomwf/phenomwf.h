/******************************************************************************/
/* phenomwf.h                                                                 */
/* Emma Robinson 2011                                                         */
/* references are to http://arxiv.org/abs/1005.3306 unless otherwise stated   */
/******************************************************************************/

/******************************************************************************/
/* structure types                                                            */
/******************************************************************************/

/* phenParams. Table II */
typedef struct tagPhenPars{
    double alpha1;
    double alpha2;
    double alpha3;
    double alpha4;
    double alpha5;
    double alpha6;
    double gamma1;
    double delta1;
    double delta2;
} phenPars;

/* amplitude PN params Eq. A5 */
typedef struct tagPNPars{
    double A0;
    double A1;
    double A2;
    double A3;
    double A4;
    double A5re;
    double A5im;
    double A6re;
    double A6im;
} PNPars;

/* taylorParams Eq. A3*/
typedef struct tagTaylorPars{
    double a0;
    double a1;
    double a2;
    double a3;
    double a4;
    double a5;
    double a6;
    double a7;
} taylorPars;

/* phase PN params Eq. A4*/
typedef struct tagphasePNPars{
    double alphaPN0;
    double alphaPN1;
    double alphaPN2;
    double alphaPN3;
    double alphaPN4;
    double alphaPN5;
    double alphaPN6;
    double alphaPN7;
} phasePNPars;

/******************************************************************************/
/* Function definitions                                                       */
/******************************************************************************/

/* waveform function */
void phenomwf( double *wfreal ,
             double *wfimag ,
             double *f ,
             int n ,
             double df ,
             double fMin ,
             double fMax ,
             double eta ,
             double chi ,
             double Mtot ,
             double dL );

/* compute the spin of the final BH Details from 
Astrophys.J.674:L29-L32, arXiv:0710.3345 [gr-qc].*/
double AEIFinalSpin( double eta ,
                     double a );

/* parameter calculation functions */

/* Table II */
void calcPhenParams( phenPars *phenParams ,
                     double eta ,
                     double chi );

/* Eq A5 */
void calcPNParams( PNPars *PNParams ,
                    double eta ,
                    double chi1 ,
                    double chi2 );

/* Eq A3 */
void calcTaylorParams( taylorPars *taylorParams ,
                       double eta,
                       double chi ,
                       double chi1 ,
                       double chi2);

/* Eq A4 */
void calcPhasePNParams( phasePNPars *phasePNParams ,
                         double eta ,
                         double chi ,
                         double chi1 ,
                         double chi2 );

/* amplitude functions */
double AmpPNFn ( double x ,
                  double eta , 
                  PNPars *params ,
                  taylorPars *tparams);

void PNAmplitude22 ( double *Amp ,
                     double x ,
                     PNPars *params );

double XdotT4Fn ( double x ,
                    double eta , 
                    taylorPars *params );

/* phase functions */
double PhasePNFn( double f ,
                    double eta ,
                    phasePNPars *params );

double PMPhaseFn ( double f ,
                    double eta , 
                    phenPars *params );

double RDPhaseFn ( double f ,
                    double eta , 
                    double fTrans ,
                    phenPars *params );

/*** auxilliary functions ***/

/* Lorentzian */
double LorentzianFn ( double freq,
			 double fRing,
			 double sigma);

/* window functions */
double TanhWindowPlus ( double freq,
                           double f0,
			   double sigma);

double TanhWindowMinus ( double freq,
                           double f0,
			   double sigma);

/* cardoso's functions */
double fQNM( double a );

double Qual( double a );

/******************************************************************************/
/* phenomenological parameters                                                */
/* from Table II in http://arxiv.org/abs/1005.3306                            */
/******************************************************************************/

#define BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA01 -2.417e-03
#define BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA02 -1.093e-03
#define BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA11 -1.917e-02
#define BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA10  7.267e-02
#define BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA20 -2.504e-01

#define BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA01  5.962e-01
#define BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA02 -5.600e-02
#define BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA11  1.520e-01
#define BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA10 -2.970e+00
#define BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA20  1.312e+01

#define BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA01 -3.283e+01
#define BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA02  8.859e+00
#define BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA11  2.931e+01
#define BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA10  7.954e+01
#define BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA20 -4.349e+02

#define BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA01  1.619e+02
#define BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA02 -4.702e+01
#define BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA11 -1.751e+02
#define BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA10 -3.225e+02
#define BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA20  1.587e+03

#define BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA01 -6.320e+02
#define BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA02  2.463e+02
#define BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA11  1.048e+03
#define BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA10  3.355e+02
#define BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA20 -5.115e+03

#define BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA01 -4.809e+01
#define BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA02 -3.643e+02
#define BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA11 -5.215e+02
#define BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA10  1.870e+03
#define BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA20  7.354e+02

#define BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA01  4.149e+00
#define BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA02 -4.070e+00
#define BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA11 -8.752e+01
#define BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA10 -4.897e+01
#define BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA20  6.665e+02

#define BBHNEWPHENOMCOEFFSH_DELTA1_ZETA01 -5.472e-02
#define BBHNEWPHENOMCOEFFSH_DELTA1_ZETA02  2.094e-02
#define BBHNEWPHENOMCOEFFSH_DELTA1_ZETA11  3.554e-01
#define BBHNEWPHENOMCOEFFSH_DELTA1_ZETA10  1.151e-01
#define BBHNEWPHENOMCOEFFSH_DELTA1_ZETA20  9.640e-01

#define BBHNEWPHENOMCOEFFSH_DELTA2_ZETA01 -1.235e+00
#define BBHNEWPHENOMCOEFFSH_DELTA2_ZETA02  3.423e-01
#define BBHNEWPHENOMCOEFFSH_DELTA2_ZETA11  6.062e+00
#define BBHNEWPHENOMCOEFFSH_DELTA2_ZETA10  5.949e+00
#define BBHNEWPHENOMCOEFFSH_DELTA2_ZETA20 -1.069e+01

/******************************************************************************/
/* params for the tanh-window functions                                       */
/* Eq. 5.8 in http://arxiv.org/abs/1005.3306                                  */
/******************************************************************************/

#define BBHNEWPHENOMCOEFFSH_D_A 0.015 /* Width of window for amplitude fns    */
#define BBHNEWPHENOMCOEFFSH_D_P 0.005 /* Width of window for phase fns        */
#define BBHNEWPHENOMCOEFFSH_F0_COEFF 0.98 /* transition freq for amplitude fns*/
#define BBHNEWPHENOMCOEFFSH_F1_COEFF 0.1  /* lower transn freq for phase fns  */
#define BBHNEWPHENOMCOEFFSH_F2_COEFF 1.0  /* upper transn freq for phase fns  */


