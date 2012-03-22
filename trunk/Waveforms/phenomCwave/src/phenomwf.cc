/******************************************************************************/
/* phenomwf.c                                                                 */
/* Emma Robinson 2011                                                         */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phenomwf.h"
#include "Constants.hh"
#include <iostream>
#include <fstream>
#include <iomanip>

/******************************************************************************/
/* Constants                                                                  */
/******************************************************************************/

const double eulergamma = 0.5772156649015328606065120900824024L;
const double pi   = 3.1415926535897932384626433832795029L;
const double mtsun = 4.92549095e-6;
const double pct = LISAWP::LISAWP_PC_SI/LISAWP::LISAWP_C_SI; /* parsec in sec */


/******************************************************************************/
/* Function to calculate the phenomenological waveforms from                  */
/* Santamaria et al: http://arxiv.org/abs/1005.3306                           */
/******************************************************************************/

void phenomwf( double *wfreal , double *wfimag , double *f , int n , 
               double fMin , double fMax , double eta , 
               double chi , double Mtot)
{

    double etasq , etacub, m1, m2;
    double pisq = pi*pi;
    double delta, chi1, chi2, finalSpin, fRD, Q, freq, x;
    double AmpPN , AmpPM, AmpRD, AmpPhenom;
    double PsiPN, PsiPM, PsiRD, PhiPhenom;
    double L, wPlus, wMinus, w1Plus, w1Minus, w2Plus, w2Minus;
    phenPars phenParams;
    PNPars PNParams;
    taylorPars taylorParams;
    phasePNPars phasePNParams;
    int i;

    /* dimensionless mass */
    m1=(1-(4.0*eta))/2.0;
    m2=(1+(4.0*eta))/2.0;

    /* asymetry parameter */
    delta = sqrt(1.-4.*eta);

    /* we assume chi1 = chi2 = chi, following the notebook */
    chi1 = chi;
    chi2 = chi;
    
    /* total spin of final black hole */
    finalSpin = AEIFinalSpin(eta,chi);

    /* calculate the phenom parameters (see Table II) */
    calcPhenParams( &phenParams , eta , chi );

    /* fQNM frequency in units of M */
    fRD = fQNM( fabs(finalSpin) );
    /* Q: quality factor of ringdown */
    Q = Qual( fabs(finalSpin) );

    /* params for the PN part of the amplitude Eq A5*/
    calcPNParams( &PNParams , eta , chi1 , chi2 );
    /* Coefficients of the TaylorT4 expansion of dx/dt Eq A3*/
    calcTaylorParams( &taylorParams , eta , chi , chi1 , chi2 );
    /* coefficients of the PN part of the phase Eq A4*/
    calcPhasePNParams( &phasePNParams , eta , chi , chi1 , chi2 );

    for ( i=0; i<n ; i++ )
    {
      freq = f[i] * Mtot * mtsun;

      if ( f[i] < fMin || f[i] > fMax || freq > .15 )
      {
        wfreal[i] = 0;
        wfimag[i] = 0;
      }
      else
      {

        x = pow( freq * pi , 2.0/3.0 );
        /* PN amplitude from Eq 3.10 */
        AmpPN = AmpPNFn ( x , eta , &PNParams , &taylorParams );
        /* pre-merger contribution Eq 5.11 */
        AmpPM = AmpPN + ( phenParams.gamma1 * pow( freq , 5.0/6.0 ) );

        /* ringdown contribution */
        L = LorentzianFn ( freq, fRD, fRD * phenParams.delta2 / Q );
        /* to correct for different defs of Lorentzian */
        L *= 2.0 * pi * fRD * phenParams.delta2 / Q;
  
        /* ringdown amplitude Eq 5.12 */
        AmpRD = phenParams.delta1 * L * pow( freq , -7.0/6.0 );

        /* window functions Eq 5.8*/
        wPlus = TanhWindowPlus( freq , 
                                BBHNEWPHENOMCOEFFSH_F0_COEFF * fRD ,
                                BBHNEWPHENOMCOEFFSH_D_A );
        wMinus = TanhWindowMinus( freq , 
                                BBHNEWPHENOMCOEFFSH_F0_COEFF * fRD ,
                                BBHNEWPHENOMCOEFFSH_D_A );

        /* phenom amplitude Eq 5.13 */
        AmpPhenom = ( (AmpPM*wMinus) + (AmpRD*wPlus) ) ; 
  
        /*** phase ***/

        /* PN contribution Eq 3.13 with t0=phi0=0*/
        PsiPN = PhasePNFn( freq , eta , &phasePNParams );

        /* pre-mreger stage Eq 5.3 */
        PsiPM = PMPhaseFn ( freq , eta , &phenParams );

        /* ringdown stage Eq 5.7 */
        PsiRD = RDPhaseFn ( freq , eta , BBHNEWPHENOMCOEFFSH_F2_COEFF * fRD , 
                            &phenParams );
  
        /* window functions Eq 5.8*/
        w1Plus = TanhWindowPlus(  freq , 
                                      BBHNEWPHENOMCOEFFSH_F1_COEFF * fRD ,
                                      BBHNEWPHENOMCOEFFSH_D_P );
        w1Minus = TanhWindowMinus( freq , 
                                       BBHNEWPHENOMCOEFFSH_F1_COEFF * fRD,
                                       BBHNEWPHENOMCOEFFSH_D_P );
        w2Plus = TanhWindowPlus(  freq , 
                                      BBHNEWPHENOMCOEFFSH_F2_COEFF * fRD ,
                                      BBHNEWPHENOMCOEFFSH_D_P );
        w2Minus = TanhWindowMinus( freq , 
                                       BBHNEWPHENOMCOEFFSH_F2_COEFF * fRD ,
                                       BBHNEWPHENOMCOEFFSH_D_P );

        /* phenom phase Eq 5.9*/
        PhiPhenom = ( PsiPN * w1Minus ) + 
                    ( w1Plus * w2Minus * PsiPM ) + 
                    ( PsiRD * w2Plus );

        /*** generate the waveform ***/
        if (AmpPN == 0.0)
           AmpPhenom = 0.0;

        /* real */
        wfreal[i] = AmpPhenom * cos( PhiPhenom );

        /* imag */
        wfimag[i] = AmpPhenom * sin (PhiPhenom );

      }
    }

    return;
}

/**********************************************************
The same as above but returns amplitude and phase
**********************************************************/

void phenomwf2( double *amp, double *phase, double *f , int n , 
               double fMin , double fMax , double eta , 
               double chi , double Mtot , double dL )
{

    double etasq , etacub, m1, m2;
    double pisq = pi*pi;
    double delta, chi1, chi2, finalSpin, fRD, Q, freq, x;
    double AmpPN , AmpPM, AmpRD, AmpPhenom;
    double PsiPN, PsiPM, PsiRD, PhiPhenom;
    double L, wPlus, wMinus, w1Plus, w1Minus, w2Plus, w2Minus;
    double Msec, DLsec;
    phenPars phenParams;
    PNPars PNParams;
    taylorPars taylorParams;
    phasePNPars phasePNParams;
    int i;

    /* dimensionless mass */
    m1=(1-(4.0*eta))/2.0;
    m2=(1+(4.0*eta))/2.0;
   

    /* asymetry parameter */
    delta = sqrt(1.-4.*eta);

    /* we assume chi1 = chi2 = chi, following the notebook */
    chi1 = chi;
    chi2 = chi;
    
    /* total spin of final black hole */
    finalSpin = AEIFinalSpin(eta,chi);
    if (finalSpin >= 0.99)
       finalSpin = 0.98;

    /* calculate the phenom parameters (see Table II) */
    calcPhenParams( &phenParams , eta , chi );

    /* fQNM frequency in units of M */
    fRD = fQNM( fabs(finalSpin) );
     //std::cout << " frd=    " << fRD  << " fin Spin = " << finalSpin<< std::endl;
    /* Q: quality factor of ringdown */
    Q = Qual( fabs(finalSpin) );

    /* params for the PN part of the amplitude Eq A5*/
    calcPNParams( &PNParams , eta , chi1 , chi2 );
    /* Coefficients of the TaylorT4 expansion of dx/dt Eq A3*/
    calcTaylorParams( &taylorParams , eta , chi , chi1 , chi2 );
    /* coefficients of the PN part of the phase Eq A4*/
    calcPhasePNParams( &phasePNParams , eta , chi , chi1 , chi2 );

    Msec = Mtot*mtsun;
    DLsec = dL*pct;
    //std::cout << "Stas: DL(sec) = "<< DLsec << std::endl;
    
    /*fprintf(stdout, "Ampl factor  =  %e\n",  Msec*Msec/DLsec);*/
    for ( i=0; i<n ; i++ )
    {
      freq = f[i] * Msec;

      if ( f[i] < fMin || f[i] > fMax || freq > .15 )
      {
        amp[i] = 0;
        phase[i] = 0;
      }
      else
      {

        x = pow( freq * pi , 2.0/3.0 );
        /* PN amplitude from Eq 3.10 */
        AmpPN = AmpPNFn ( x , eta , &PNParams , &taylorParams );
        /* pre-merger contribution Eq 5.11 */
        AmpPM = AmpPN + ( phenParams.gamma1 * pow( freq , 5.0/6.0 ) );
 
        //std::cout << f[i] << "   " << AmpPN << "   " << AmpPM  << "   " << x<< std::endl;
        /* ringdown contribution */
        L = LorentzianFn ( freq, fRD, fRD * phenParams.delta2 / Q );
        /* to correct for different defs of Lorentzian */
        L *= 2.0 * pi * fRD * phenParams.delta2 / Q;
  
        /* ringdown amplitude Eq 5.12 */
        AmpRD = phenParams.delta1 * L * pow( freq , -7.0/6.0 );

        /* window functions Eq 5.8*/
        wPlus = TanhWindowPlus( freq , 
                                BBHNEWPHENOMCOEFFSH_F0_COEFF * fRD ,
                                BBHNEWPHENOMCOEFFSH_D_A );
        wMinus = TanhWindowMinus( freq , 
                                BBHNEWPHENOMCOEFFSH_F0_COEFF * fRD ,
                                BBHNEWPHENOMCOEFFSH_D_A );

        /* phenom amplitude Eq 5.13 */
       // AmpPhenom = ( (AmpPM*wMinus) + (AmpRD*wPlus) ) ;
        AmpPhenom = ( (AmpPM*wMinus) + (AmpRD*wPlus) )  * Msec *Msec / (DLsec );
                     // * Msec * Msec  / (DLsec ); 
  
       
        /*** phase ***/

        /* PN contribution Eq 3.13 with t0=phi0=0*/
        PsiPN = PhasePNFn( freq , eta , &phasePNParams );

        /* pre-mreger stage Eq 5.3 */
        PsiPM = PMPhaseFn ( freq , eta , &phenParams );

        /* ringdown stage Eq 5.7 */
        PsiRD = RDPhaseFn ( freq , eta , BBHNEWPHENOMCOEFFSH_F2_COEFF * fRD , 
                            &phenParams );
  
        /* window functions Eq 5.8*/
        w1Plus = TanhWindowPlus(  freq , 
                                      BBHNEWPHENOMCOEFFSH_F1_COEFF * fRD ,
                                      BBHNEWPHENOMCOEFFSH_D_P );
        w1Minus = TanhWindowMinus( freq , 
                                       BBHNEWPHENOMCOEFFSH_F1_COEFF * fRD,
                                       BBHNEWPHENOMCOEFFSH_D_P );
        w2Plus = TanhWindowPlus(  freq , 
                                      BBHNEWPHENOMCOEFFSH_F2_COEFF * fRD ,
                                      BBHNEWPHENOMCOEFFSH_D_P );
        w2Minus = TanhWindowMinus( freq , 
                                       BBHNEWPHENOMCOEFFSH_F2_COEFF * fRD ,
                                       BBHNEWPHENOMCOEFFSH_D_P );

        /* phenom phase Eq 5.9*/
        PhiPhenom = ( PsiPN * w1Minus ) + 
                    ( w1Plus * w2Minus * PsiPM ) + 
                    ( PsiRD * w2Plus );

        /*** generate the waveform ***/

        /* amplitude */
        if (AmpPN == 0.0)
               AmpPhenom = 0.0;
        amp[i] = AmpPhenom;

        /* phase */
        phase[i] = PhiPhenom;

      }
    }

    return;
}



/******************************************************************************/
/* calculate parameters for the various wf functions                          */
/******************************************************************************/

/* the phenomenological parameters Table II*/
void calcPhenParams( phenPars *phenParams , double eta , double chi )
{
    double etasq = eta*eta;
    double etachi = eta*chi;
    double chisq = chi*chi;

    phenParams->alpha1 = BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA01 * chi + 
             BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA02 * chisq +
             BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA11 * etachi +
             BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA10 * eta + 
             BBHNEWPHENOMCOEFFSH_ALPHA1_ZETA20 * etasq;

    phenParams->alpha2 = BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA01 * chi + 
             BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA02 * chisq +
             BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA11 * etachi +
             BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA10 * eta + 
             BBHNEWPHENOMCOEFFSH_ALPHA2_ZETA20 * etasq;

    phenParams->alpha3 = BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA01 * chi + 
             BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA02 * chisq +
             BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA11 * etachi +
             BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA10 * eta + 
             BBHNEWPHENOMCOEFFSH_ALPHA3_ZETA20 * etasq;

    phenParams->alpha4 = BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA01 * chi + 
             BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA02 * chisq +
             BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA11 * etachi +
             BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA10 * eta + 
             BBHNEWPHENOMCOEFFSH_ALPHA4_ZETA20 * etasq;

    phenParams->alpha5 = BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA01 * chi + 
             BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA02 * chisq +
             BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA11 * etachi +
             BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA10 * eta + 
             BBHNEWPHENOMCOEFFSH_ALPHA5_ZETA20 * etasq;

    phenParams->alpha6 = BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA01 * chi + 
             BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA02 * chisq +
             BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA11 * etachi +
             BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA10 * eta + 
             BBHNEWPHENOMCOEFFSH_ALPHA6_ZETA20 * etasq;

    phenParams->gamma1 = BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA01 * chi + 
             BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA02 * chisq +
             BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA11 * etachi +
             BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA10 * eta + 
             BBHNEWPHENOMCOEFFSH_GAMMA1_ZETA20 * etasq;

    phenParams->delta1 = BBHNEWPHENOMCOEFFSH_DELTA1_ZETA01 * chi + 
             BBHNEWPHENOMCOEFFSH_DELTA1_ZETA02 * chisq +
             BBHNEWPHENOMCOEFFSH_DELTA1_ZETA11 * etachi +
             BBHNEWPHENOMCOEFFSH_DELTA1_ZETA10 * eta + 
             BBHNEWPHENOMCOEFFSH_DELTA1_ZETA20 * etasq;

    phenParams->delta2 = BBHNEWPHENOMCOEFFSH_DELTA2_ZETA01 * chi + 
             BBHNEWPHENOMCOEFFSH_DELTA2_ZETA02 * chisq +
             BBHNEWPHENOMCOEFFSH_DELTA2_ZETA11 * etachi +
             BBHNEWPHENOMCOEFFSH_DELTA2_ZETA10 * eta + 
             BBHNEWPHENOMCOEFFSH_DELTA2_ZETA20 * etasq;

    return;
}

/* coefficients for the PN part of the amplitude Eq A5*/
void calcPNParams( PNPars *PNParams , double eta , double chi1 , double chi2 )
{
    double chi12 = chi1*chi2;
    double etasq = eta*eta;
    double etacub = eta*etasq;
    double pisq = pi*pi;

    PNParams->A0 = 1.0;
    PNParams->A1 = 0.0;
    PNParams->A2 = (55.0*eta/42.0) - (107.0/42.0);
    PNParams->A3 = (2.0*pi) + 
         (-2.0*sqrt(1.0-(4.0*eta))*(chi1-chi2)/3.0) - 
         (2.0*(1-eta)*(chi1+chi2)/3.0);
    PNParams->A4 = (-2173.0/1512.0) -
         (eta*((1069.0/216.0)-(2.0*chi12)))+
         (etasq*2047.0/1512.0);
    PNParams->A5re = (-107.0*pi/21.0) +
         (eta*(34.0*pi/21.0));
    PNParams->A5im = -24.0 * eta;
    /* A6 doesn't include the log(16x) part */
    PNParams->A6re = (27027409.0/646800.0) - 
            (856.0*eulergamma/105.0) +
            (2.0*pisq/3.0) + 
            (eta*((41.0*pisq/96.0)-(278185.0/33264.0))) -
            (etasq*20261.0/2772.0) +
            (etacub*114635.0/99792.0);
    PNParams->A6im = (428.0*pi/105.0);
    return;
}

/* coefficients for the TaylorT4 expansion of xdot Eq A3*/
void calcTaylorParams( taylorPars *taylorParams , double eta, 
                       double chi , double chi1 , double chi2)
{
    double etasq = eta*eta;
    double etacub=etasq*eta;
    double chisq = chi*chi;
    double chicub=chi*chisq;
    double chi12 = chi1*chi2;
    double etachi = eta*chi;
    double pisq=pi*pi;

    taylorParams->a0 = 1.0;
    taylorParams->a1 = 0.0;
    taylorParams->a2 = (-743.0/336.0) - (11.0*eta/4.0);
    taylorParams->a3 = (4.0*pi) - 
         (113.0*chi/12.0) + 
         (19.0*eta*(chi1+chi2)/6.0);
    taylorParams->a4 = (34103.0/18144.0) + 
         (5.0*chisq) +
         (eta*((13661.0/2016.0)-(chi12/8.0))) +
         (59.0*etasq/18.0);
    taylorParams->a5 = (-1.0*pi*((4159.0/672.0)+(189.0*eta/8.0))) -
         (chi*((31571.0/1008.0)-(1165.0*eta/24.0))) +
       ((chi1+chi2)*((21863.0*eta/1008.0)-(79.0*etasq/6.0))) -
         (3.0*chicub/4.0) +
         (9.0*etachi*chi12/4.0);
      /* a6 does not include the log(16x) part */
    taylorParams->a6 = (16447322263.0/139708800.0) -
         (1712.0*eulergamma/105.0) +
         (16.0*pisq/3.0) +
         (eta*((451.0*pisq/48.0)
                              -(56198689.0/217728.0))) +
         (541.0*etasq/896.0) -
         (5605.0*etacub/2592.0) -
         (80.0*pi*chi/3.0) +
         (((20.0*pi/3.0)-(1135.0*chi/36.0))*eta*(chi1+chi2)) +
         (((64153.0/1008.0)-(457.0*eta/36.0))*chisq) -
         (((787.0*eta/144.0)-(3037.0*etasq/144.0))*chi12);
    taylorParams->a7 = (-1.0*pi * ( (4415.0/4032.0) - 
                     (358675.0*eta/6048.0) -
                     (91495.0*etasq/1512.0) ) ) -
         (chi*( (2529407.0/27216.0) -
                (845827.0*eta/6048.0) +
                (41551.0*etasq/864.0) )) +
         ((chi1+chi2)*( (1580239.0*eta/54432.0) -
                        (451597.0*etasq/6048.0) +
                        (2045.0*etacub/432.0) +
                        (107.0*etachi*chi/6.0) -
                        (5.0*etasq*chi12/24.0) )) + 
         (12.0*pi*chisq) - 
         (chicub*( (1505.0/24.0) +
                        (eta/8.0) )) +
         (chi*chi12*( (101.0*eta/24.0) +
                          (3.0*etasq/8.0) ));
    return;
}

/* coefficients for the PN part of the phase Eq A4*/
void calcPhasePNParams( phasePNPars *phasePNParams , double eta , 
                        double chi , double chi1 , double chi2 )
{
    double etasq = eta*eta;
    double etacub = eta*etasq;
    double chisq = chi*chi;
    double chicub = chi*chisq;
    double etachi = eta*chi;
    double chi12 = chi1*chi2;
    double pisq = pi*pi;

    phasePNParams->alphaPN0 = 1.0;
    phasePNParams->alphaPN1 = 0.0;
    phasePNParams->alphaPN2 = (3715.0/756.0) + (55.0*eta/9.0);
    phasePNParams->alphaPN3 = (-16.0*pi) + (113.0*chi/3.0) - (38.0*eta*(chi1+chi2)/3.0);
    phasePNParams->alphaPN4 = (15293365.0/508032.0) - 
                (50.0*chisq) + 
                (eta*((27145.0/504.0)+(5.0*chi12/4.0))) +
                (3085.0*etasq/72.0);
      /* alphaPN5 doesn't include the ln(pi*f) part */
    phasePNParams->alphaPN5 = (pi*((38645.0/756.0)-(65.0*eta/9.0))) -
                (chi*((735505.0/2268.0)+(130.0*eta/9.0))) +
                ((chi1+chi2)*((12850.0*eta/81.0)+(170.0*etasq/9.0))) - 
                (10.0*chicub/3.0) +
                (10.0*etachi*chi12);
      /* alphaPN6 doesn't include the ln(64*pi*f) part */
    phasePNParams->alphaPN6 = (11583231236531.0/4694215680.0) - 
                (640.0*pisq/3.0) - 
                (6848.0*eulergamma/21.0) +  
                (eta*((2255.0*pisq/12.0)-
                      (15737765635.0/3048192.0))) +
                (76055.0*etasq/1728.0) - 
                (127825.0*etacub/1296.0) +
                (2920.0*pi*chi/3.0) - 
                ((175.0-(1490.0*eta))*chisq/3.0) -
                ( ( (1120.0*pi/3.0) - 
                    (1085.0*chi/3.0) ) * eta * (chi1+chi2)) +
                ( ( (26945.0*eta/336.0) - 
                    (2365.0*etasq/6.0) )*chi12);
    phasePNParams->alphaPN7 = (pi*((77096675.0/254016.0)+
                        (378515.0*eta/1512.0)-
                        (74045.0*etasq/756.0))) - 
                (chi*( (20373952415.0/3048192.0) +
                       (150935.0*eta/224.0) -
                       (578695.0*etasq/432.0) )) +
                ((chi1+chi2)*( (4862041225.0*eta/1524096.0) +
                               (1189775.0*etasq/1008.0) - 
                               (71705.0*etacub/216.0) -
                               (830.0*eta*chisq/3.0) +
                               (35.0*etasq*chi12/3.0) )) -
                (560.0*pi*chisq) +
                (20.0*pi*eta*chi12) +
                (chicub*( (94555.0/168.0) -
                               (85.0*eta) )) +
                (chi*chi12*( (39665.0*eta/168.0) +
                                 (255.0*etasq) ));
    return;
}

/******************************************************************************/
/* Amplitude functions                                                        */
/******************************************************************************/

/* Calculate the PN amplitude in the freq domain Eq 3.10*/
double AmpPNFn ( double x ,
                  double eta , 
                  PNPars *params ,
                  taylorPars *tparams )
{
  double Amp, Amp0, XdotT4;
  double* Amp22;
  double fractionre; 
  double fractionim;

  //Amp22 = malloc(sizeof(double)*2);
  Amp22 = new double[2];

  Amp0 = sqrt( 2.0 * pi / ( 3.0 * sqrt(x) ) );

  PNAmplitude22 ( Amp22 , x , params );
  /* real part */
  Amp22[0] *= 8.0*eta*x*sqrt(pi/5.0);
  /* imag part */
  Amp22[1] *= 8.0*eta*x*sqrt(pi/5.0);
  XdotT4 = XdotT4Fn ( x , eta , tparams );
 
 // std::cout << "Stas: XdotT4 = " << XdotT4 << std::endl;
 
  if (XdotT4 > 0.0){
     fractionre=Amp22[0]/sqrt(XdotT4);
     fractionim=Amp22[1]/sqrt(XdotT4);
  }else{
     fractionre = 0.0;
     fractionim = 0.0;
  }

  Amp = Amp0 * sqrt((fractionre*fractionre)+(fractionim*fractionim));

  delete [] Amp22;
  return(Amp);
}

/* Function to calculate the 3PN time-domain amplitude for l=m=2, Eq 3.14*/
void PNAmplitude22 ( double *Amp ,
                     double x ,
                     PNPars *params )
{

  double x2=x*x;
  double x3=x2*x;
  double x52=sqrt(x2*x3);
  Amp[0] = params->A0 +
        ( params->A1 * sqrt(x) ) +
        ( params->A2 * x ) +
        ( params->A3 * sqrt(x3) ) +
        ( params->A4 * x2 ) +
        ( params->A5re * x52 ) +
        ( (params->A6re - (428.0*log(16.0*x)/105.0))* x3 );
  Amp[1] = ( params->A5im * x52 ) +
           ( params->A6im* x3 );
  
  return;

}

/* Function to calculate the TaylorT4 expansion of dx/dt Eq 3.15*/
double XdotT4Fn ( double x ,
                    double eta , 
                    taylorPars *params )
{

  double Xdot, Amp0;
  double x2=x*x;
  double x3=x2*x;
  double x52=sqrt(x2*x3);
  double x72=x52*x;
  Amp0 = 64.0 * eta * pow(x, 5.0) / 5.0;

  Xdot = params->a0 +
        ( params->a1 * sqrt(x) ) +
        ( params->a2 * x ) +
        ( params->a3 * sqrt(x3) ) +
        ( params->a4 * x2 ) +
        ( params->a5 * x52 ) +
        ( ( params->a6 - (856.0*log(16.0*x)/105.0) ) * x3 ) +
        ( params->a7 * x72 );

  Xdot *= Amp0;
 
  return(Xdot);

}

/******************************************************************************/
/* Phase functions                                                            */
/******************************************************************************/

/* PN part of the phase Eq 3.13 */
double PhasePNFn( double f ,
                    double eta ,
                    phasePNPars *params )
{

  double PsiPN;
  double pif13 = pow( pi*f , 1.0/3.0 );
  double pif23 = pif13*pif13;
  double pif43 = pif23*pif23;
  double pif53 = pif43*pif13;
  double pif63 = pif43*pif23;
  double pif73 = pif63*pif13;

  PsiPN = params->alphaPN0 +
           params->alphaPN1*pif13 +
           params->alphaPN2*pif23 +
           params->alphaPN3* pi*f  +
           params->alphaPN4*pif43 +
           params->alphaPN5*(1.0+log(pi*f))*pif53 +
           (params->alphaPN6 - 6848.0*log(64.0*pi*f)/63.0)* pif63 +
           params->alphaPN7*pif73;

  PsiPN *= 3.0 / (pif53 * 128.0 * eta );
  PsiPN -= (pi / 4.0);
  
  return PsiPN;
}

/* the pre-merger phase Eq 5.3*/
double PMPhaseFn ( double f ,
                    double eta , 
                    phenPars *params )
{
  double PsiPM;
  double fm13 = pow(f,-1.0/3.0);
  double fm53 = fm13 * fm13 * fm13 * fm13 * fm13;
  double f23 = 1/(fm13 * fm13);
  
  PsiPM = ( (params->alpha1*fm53) + 
            (params->alpha2/f) +
            (params->alpha3*fm13) +
             params->alpha4 +
            (params->alpha5*f23) +
            (params->alpha6*f) ) / eta;

  return(PsiPM);
}

/* the ring-down phase Eq 5.7*/
double RDPhaseFn ( double f ,
                    double eta , 
                    double fTrans ,
                    phenPars *params )
{
  double PsiRD, beta1, beta2;
  double fTransm13 = pow(fTrans,-1.0/3.0);
  double fTransm43 = fTransm13*fTransm13*fTransm13*fTransm13;
  double fTransm63 = fTransm43*fTransm13*fTransm13;
  double fTransm83 = fTransm43*fTransm43;

  beta2 = ( (-5.0*params->alpha1*fTransm83/3.0) -
            (params->alpha2*fTransm63) -
            (params->alpha3*fTransm43/3.0) +
            (2.0*params->alpha5*fTransm13/3.0) +
            params->alpha6 ) / eta;
  beta1 = PMPhaseFn ( fTrans , eta , params ) - (fTrans * beta2);
 

  PsiRD = beta1 + ( beta2 * f );

  return(PsiRD);
}

/******************************************************************************/
/* Auxilliary functions                                                       */
/******************************************************************************/

/******************************************************************************/
/* Final spin of the resulting black hole                                     */
/* compute the spin of the final BH Details from                              */
/*Astrophys.J.674:L29-L32, arXiv:0710.3345 [gr-qc].                           */
/******************************************************************************/
double AEIFinalSpin( double eta ,
                     double a )
{
  double aFin;
  double s4 = -0.129;
  double s5 = -0.384;
  double t0 = -2.686;
  double t2 = -3.454;
  double t3 = 2.353;
  double etasq = eta * eta;
  
  aFin = a + (s4 * a * a * eta) + (s5 * a * etasq) + (t0 * a * eta) + 
         (2.0 * sqrt(3.0) * eta) + (t2 * etasq) + 
         (t3 * etasq * eta);

  return aFin;
}


/******************************************************************************/
/* Lorentzian                                                                 */
/******************************************************************************/

double LorentzianFn ( double freq,
			 double fRing,
			 double sigma)
{
  double out;

  out = sigma / (2 * pi * ((freq - fRing)*(freq - fRing)
			       + sigma*sigma / 4.0));

  return(out);
}

/******************************************************************************/
/* Window functions Eq 5.8                                                    */
/******************************************************************************/

double TanhWindowPlus ( double freq,
                           double f0,
			   double sigma)
{
  double out, fact;
  fact = 4.0 * (freq-f0) / sigma;
  out = (1.0 + tanh(fact)) / 2.0;
  return(out);
}

double TanhWindowMinus ( double freq,
                           double f0,
			   double sigma)
{
  double out, fact;
  fact = 4.0 * (freq-f0) / sigma;
  out = (1 - tanh(fact)) / 2.0;
  return(out);
}

/******************************************************************************/
/* Functions to compute the QNM frequency (fQNM) and the quality factor (Q)   */
/* of the ringdown.                                                           */
/* from Victor Cardoso's webpage                                              */
/* http://gamow.ist.utl.pt/~vitor/ringdown.html for l=2, m=2, n=0 mode        */
/******************************************************************************/

double fQNM( double a )
{

    double out;
    double f1=1.5251;
    double f2=-1.1568;
    double f3=0.1292;

    out = ( f1 + (f2*pow((1.0-a),f3)) );
    out /= (2.0 * pi);

    return(out);

}

double Qual( double a )
{

    double out;
    double q1 = .7;
    double q2 = 1.4187;
    double q3 = -.499;

    out = q1 + ( q2 * pow( (1.0-a) , q3 ) );

    return(out);
}
