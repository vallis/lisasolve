#include <iostream>
#include <fstream>
#include <iomanip>
#include "FreqAK_RA.hh"
#include "Constants.hh"
#include "EMRItemplate.hh"
#include "LISAmotion.hh"
#include "Macros.hh"
#include "FastFT.hh"

using namespace LISAWP;

/**
This test code generates orbital evolution under radiation reaction
and computes AK (analytic kludge) waveform
*/


int main(){
 
  std::cout.precision(10);
  std::string spr = "     ";
  
  double mu = 10.;
  double MBHMass = 5.e6;
  double spin = 0.9;
  double lam = 0.5;
  double thS = 0.5; /* Sky position - colatitude in [0,pi] */
  double phS = 2.; /* Sky position - azimuth in [0,2*pi] */
  double thK = 2.; /* Orientation of spin - colatitude */
  double phK = 3.4;
  double D = 1.e6;
  
  thS = LISAWP_PI/2. + 1.41414031334;
  phS = 4.89367416229;
  thK = 1.70027801062;
  phK = 2.15339027642;
  spin = 0.651010374883;
  mu = 10.1314155287;
  MBHMass = 10397935.923;
  double ph0 = 5.08966915448;
  double nu0 = 0.000167447160358;
  double e0 = 0.252402550324; 
  double gam0 = 2.02625005151;
  double al0 = 4.08747792161;
  lam = 1.20560687892;
  D = 328512268.635;
  double Tpl = 62613029.9325;
  
  // high freq
   thS = LISAWP_PI/2. +  0.118089064616;
   phS = 4.95326609056;
   thK = 1.70965366389;
   phK = 2.19422416585;
   spin = 0.650048619921;
   mu = 9.74784672881;
   MBHMass = 975650.297375;
   ph0 = 3.71212583344;
   nu0 = 0.000999762694898;
   e0 = 0.360970388293; 
   gam0 = 2.11218737343;
   al0 = 0.376758533458;
   lam = 0.510993867146;
   D = 1844211338.18;
   Tpl = 45949211.6256;//45950976
   double al_pl = 4.28437707876;
   double ph_pl = 5.60608267871;
   double gam_pl = 2.16842966408;
   double e_pl = 0.226069543231;
   double nu_pl = pow( (1.0-e_pl*e_pl)/(6.0+2.0*e_pl), 1.5 )/(LISAWP_TWOPI * MBHMass*LISAWP_MTSUN_SI);
   
  
  
  EMRItemplate S(MBHMass, mu, spin, lam);
  
  S.SetPosition(thS, phS, thK, phK, D);
  S.t0 = 0.0;
  //S.e0 = 0.2;
  //S.nu0 = 0.000362112;
  //S.Phi0 = 0.0;
  //S.alpha0 = 0.0;
  //S.gamma0 = 0.0;
  S.e0 = e0;
  S.nu0 = nu0;
  S.Phi0 = ph0;
  S.alpha0 = al0;
  S.gamma0 = gam0;
  
  S.e0 = e_pl;
  S.nu0 = nu_pl;
  S.Phi0 = ph_pl;
  S.alpha0 = al_pl;
  S.gamma0 = gam_pl;
  
  double timestep = 15.0;
  double maxDur =  62914560.;
  double dtPhase = 2048.;
  int Nps = (int)floor(maxDur/timestep);
  S.XAalpha.ExpandEmpty(Nps);
  S.YEbeta.ExpandEmpty(Nps);
  S.ZTgamma.ExpandEmpty(Nps);
  S.XAalpha = 0.0;
  S.YEbeta = 0.0;
  S.ZTgamma = 0.0;
  
  double arm = 16.6782;
  
  
  
  FreqAK_RA frAKRA(timestep, maxDur, dtPhase, arm);
//  for (int i = 0; i<10; i++){
     //frAKRA.PhaseEv_RA(S, 0.0, maxDur, 0.0);
     frAKRA.PhaseEv_RA(S, Tpl, Tpl, 0.0);
     frAKRA.ConstructWaveFreq_RA(S, S.XAalpha, S.YEbeta, S.ZTgamma);
 // }     
  
   Matrix<double> X(Nps);
   Matrix<double> Y(Nps);
   Matrix<double> Z(Nps);
   
   FastFT fft;
   fft.Inverse(S.XAalpha, X);
   fft.Inverse(S.YEbeta, Y);
   fft.Inverse(S.ZTgamma, Z);
   
   std::ofstream fout("Data/FreqAK.dat");
   double T;
   for (int i=1; i<Nps; i++){
       T = (double)i*timestep;
       fout << std::setprecision(10) << T << spr << X(i)/maxDur << spr << Y(i)/maxDur << spr << Z(i)/maxDur << std::endl; 
   }
   fout.close();
     
  
  return(0);
  
}