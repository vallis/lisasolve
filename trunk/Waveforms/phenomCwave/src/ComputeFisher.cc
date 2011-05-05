/*
Copyright (C) 2005  S. Babak

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/      
/******************  CVS info ************************ 
#define CVSTAG "$Name:  $"
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/src/RandomGen.cc,v 1.1 2007/09/18 21:50:43 stba Exp $" 
*/



#include "ComputeFisher.hh"

namespace LISAWP{
   
ComputeFisherC::ComputeFisherC(double arm, double oneyear, std::string config, double Fmax, double delf, double duration){
   
   L = arm;
   year = oneyear;
   conf = config;
   fmax = Fmax;
   df = delf;
   Tobs = duration;
   
}

ComputeFisherC::~ComputeFisherC()
{
   std::cout << "The end \n";   
}

void ComputeFisherC::ComputeRAFisher4links(BBHTemplate S, int n, double* &freq, double* &Sn, Matrix<double> &Fisher){
   
   // Computing derivatives numerically
   std::complex<double>* X;
   std::complex<double>* Y;
   std::complex<double>* Z;
   X = new  std::complex<double>[n];
   Y = new  std::complex<double>[n];
   Z = new  std::complex<double>[n];
   
   std::complex<double>* Xn;
   Xn = new  std::complex<double>[n];
   
   BBHTemplate S2(S);
   
   ComputeWaveXYZ(S, n, freq, X, Y, Z);
   double tc = S.tc;
   double ThetaQ = pow(0.2*S.eta/S.Mt*tc, -0.25);
   double x = 0.25*ThetaQ*( 1. + (743./4032. + 11.*S.eta/48.)*ThetaQ );
   fminS = pow(x, 1.5)/(LISAWP_PI*S.Mt);
   // find max freq of signal
   std::ofstream fout4("Data/CheckW.dat");
   for (int i=0; i<n; i++){
//      fout4 << freq[i] << "    " << abs(X[i]) << std::endl;
      if (freq[i] > fminS ){
         if (abs(X[i]) == 0.0){
            fmaxS = freq[i];
            break;
         }
      }
   }

   std::cout << "Freq. range of the signal: " << fminS << "   " << fmaxS << std::endl; 
   
   
   // compute derivatives w.r.t. total mass
   
   double dM = 1.e-6*S.M; 
   S2.M = S.M + dM;
   S2.m1 = 0.5*S2.M*(1.+ sqrt(1. - 4.*S2.eta));
   S2.m2 = 0.5*S2.M*(1.- sqrt(1. - 4.*S2.eta));
   S2.ComputeAuxParams();
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXm;
   dXm = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXm[i] = (Xn[i] - X[i])/dM;
      fout4 << freq[i] << "     " << dXm[i] << std::endl;
   }
   
   // Derivative w.r.t. eta
   
   S2 = S;
   double deta = 1.e-6;
   if (S2.eta < 0.25 - deta){
      S2.eta = S.eta + deta;
   }else{
      deta *= -1.;
      S2.eta = S.eta + deta;  
   }   
   S2.m1 = 0.5*S2.M*(1.+ sqrt(1. - 4.*S2.eta));
   S2.m2 = 0.5*S2.M*(1.- sqrt(1. - 4.*S2.eta));
   S2.ComputeAuxParams();
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXeta;
   dXeta = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXeta[i] = (Xn[i] - X[i])/deta;
   }
   
   // Derivative w.r.t. chi
   
   S2 = S;
   double dchi = 1.e-5;
   if (S2.chi < 1.-dchi){
      S2.chi = S.chi + dchi;
   }else{
      dchi *= -1.;
      S2.chi = S.chi + dchi;
   }
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXchi;
   dXchi = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXchi[i] = (Xn[i] - X[i])/dchi;
   }
   
   
   // Derivatives w.r.t. sky postition
   
   S2 = S;
   double dthetaS = 1.e-4;
   if (S2.thetaS < LISAWP_PI - dthetaS){
      S2.thetaS = S.thetaS + dthetaS;
   }else{
      dthetaS *= -1.;
      S2.thetaS = S.thetaS + dthetaS;
   }
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXthS;
   dXthS = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXthS[i] = (Xn[i] - X[i])/dthetaS;
   }
   
   S2 = S;
   double dphiS = 1.e-4;
   S2.phiS = S.phiS + dphiS;
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXphS;
   dXphS = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXphS[i] = (Xn[i] - X[i])/dphiS;
   }
   
   // Derivative w.r.t log(Tc)
   
   S2 = S;
   double delTc = 0.1;
   S2.tc = S.tc + delTc;
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXtc;
   dXtc = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXtc[i] = S.tc*(Xn[i] - X[i])/delTc;
   }
   
   // Derivative w.r.t. psi (polariz)
   
   S2 = S;
   double dpsi = 1.e-3;
   S2.psi = S.psi + dpsi;
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXpsi;
   dXpsi = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXpsi[i] = (Xn[i] - X[i])/dpsi;
   }
   
   // Derivative w.r.t. phi0 (initial orbital phase)
   
   S2 = S;
   double dphi0 = 1.e-3;
   S2.phi0 = S.phi0 + dphi0;
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXphi0;
   dXphi0 = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXphi0[i] = (Xn[i] - X[i])/dphi0;
   }
   
   // Derivative w.r.t. inclination
   
   S2 = S;
   double diota = 1.e-3;
   S2.iota = S.iota + diota;
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXiota;
   dXiota = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXiota[i] = (Xn[i] - X[i])/diota;
   }
   
   // Derivative w.r.t log(DL)
   
   S2 = S;
   double delDL = 1.e-6*S.dist;
   S2.dist = S.dist + delDL;
   
   ComputeWaveXYZ(S2, n, freq, Xn, Y, Z);
   
   std::complex<double>* dXdl;
   dXdl = new  std::complex<double>[n];
   for (int i=0; i<n; i++){
      dXdl[i] = S.dist*(Xn[i] - X[i])/delDL;
   }
   
   // Computing Fisher matrix 
   // I assume that the size of Matrix Fisher is correct: 
   
   // The order of parameters:
   // (0) M, (1) eta, (2) chi, (3) thetaS, (4) phiS, (5) log(tc), (6) psi, (7) phi0, (8) iota, (9) log(DL)
   
   Fisher(0,0) = ComputeInnerProd(n, dXm, dXm, freq, Sn);
   Fisher(0,1) = ComputeInnerProd(n, dXm, dXeta, freq, Sn);
   Fisher(0,2) = ComputeInnerProd(n, dXm, dXchi, freq, Sn);
   Fisher(0,3) = ComputeInnerProd(n, dXm, dXthS, freq, Sn);
   Fisher(0,4) = ComputeInnerProd(n, dXm, dXphS, freq, Sn);
   Fisher(0,5) = ComputeInnerProd(n, dXm, dXtc, freq, Sn);
   Fisher(0,6) = ComputeInnerProd(n, dXm, dXpsi, freq, Sn);
   Fisher(0,7) = ComputeInnerProd(n, dXm, dXphi0, freq, Sn);
   Fisher(0,8) = ComputeInnerProd(n, dXm, dXiota, freq, Sn);
   Fisher(0,9) = ComputeInnerProd(n, dXm, dXdl, freq, Sn);
   
   Fisher(1,1) = ComputeInnerProd(n, dXeta, dXeta, freq, Sn);
   Fisher(1,2) = ComputeInnerProd(n, dXeta, dXchi, freq, Sn);
   Fisher(1,3) = ComputeInnerProd(n, dXeta, dXthS, freq, Sn);
   Fisher(1,4) = ComputeInnerProd(n, dXeta, dXphS, freq, Sn);
   Fisher(1,5) = ComputeInnerProd(n, dXeta, dXtc, freq, Sn);
   Fisher(1,6) = ComputeInnerProd(n, dXeta, dXpsi, freq, Sn);
   Fisher(1,7) = ComputeInnerProd(n, dXeta, dXphi0, freq, Sn);
   Fisher(1,8) = ComputeInnerProd(n, dXeta, dXiota, freq, Sn);
   Fisher(1,9) = ComputeInnerProd(n, dXeta, dXdl, freq, Sn);
   
   Fisher(2,2) = ComputeInnerProd(n, dXchi, dXchi, freq, Sn);
   Fisher(2,3) = ComputeInnerProd(n, dXchi, dXthS, freq, Sn);
   Fisher(2,4) = ComputeInnerProd(n, dXchi, dXphS, freq, Sn);
   Fisher(2,5) = ComputeInnerProd(n, dXchi, dXtc, freq, Sn);
   Fisher(2,6) = ComputeInnerProd(n, dXchi, dXpsi, freq, Sn);
   Fisher(2,7) = ComputeInnerProd(n, dXchi, dXphi0, freq, Sn);
   Fisher(2,8) = ComputeInnerProd(n, dXchi, dXiota, freq, Sn);
   Fisher(2,9) = ComputeInnerProd(n, dXchi, dXdl, freq, Sn);
   
   Fisher(3,3) = ComputeInnerProd(n, dXthS, dXthS, freq, Sn);
   Fisher(3,4) = ComputeInnerProd(n, dXthS, dXphS, freq, Sn);
   Fisher(3,5) = ComputeInnerProd(n, dXthS, dXtc, freq, Sn);
   Fisher(3,6) = ComputeInnerProd(n, dXthS, dXpsi, freq, Sn);
   Fisher(3,7) = ComputeInnerProd(n, dXthS, dXphi0, freq, Sn);
   Fisher(3,8) = ComputeInnerProd(n, dXthS, dXiota, freq, Sn);
   Fisher(3,9) = ComputeInnerProd(n, dXthS, dXdl, freq, Sn);
   
   Fisher(4,4) = ComputeInnerProd(n, dXphS, dXphS, freq, Sn);
   Fisher(4,5) = ComputeInnerProd(n, dXphS, dXtc, freq, Sn);
   Fisher(4,6) = ComputeInnerProd(n, dXphS, dXpsi, freq, Sn);
   Fisher(4,7) = ComputeInnerProd(n, dXphS, dXphi0, freq, Sn);
   Fisher(4,8) = ComputeInnerProd(n, dXphS, dXiota, freq, Sn);
   Fisher(4,9) = ComputeInnerProd(n, dXphS, dXdl, freq, Sn);
   
   Fisher(5,5) = ComputeInnerProd(n, dXtc, dXtc, freq, Sn);
   Fisher(5,6) = ComputeInnerProd(n, dXtc, dXpsi, freq, Sn);
   Fisher(5,7) = ComputeInnerProd(n, dXtc, dXphi0, freq, Sn);
   Fisher(5,8) = ComputeInnerProd(n, dXtc, dXiota, freq, Sn);
   Fisher(5,9) = ComputeInnerProd(n, dXtc, dXdl, freq, Sn);
   
   Fisher(6,6) = ComputeInnerProd(n, dXpsi, dXpsi, freq, Sn);
   Fisher(6,7) = ComputeInnerProd(n, dXpsi, dXphi0, freq, Sn);
   Fisher(6,8) = ComputeInnerProd(n, dXpsi, dXiota, freq, Sn);
   Fisher(6,9) = ComputeInnerProd(n, dXpsi, dXdl, freq, Sn);
   
   Fisher(7,7) = ComputeInnerProd(n, dXphi0, dXphi0, freq, Sn);
   Fisher(7,8) = ComputeInnerProd(n, dXphi0, dXiota, freq, Sn);
   Fisher(7,9) = ComputeInnerProd(n, dXphi0, dXdl, freq, Sn);
   
   Fisher(8,8) = ComputeInnerProd(n, dXiota, dXiota, freq, Sn);
   Fisher(8,9) = ComputeInnerProd(n, dXiota, dXdl, freq, Sn);
   
   Fisher(9,9) = ComputeInnerProd(n, dXdl, dXdl, freq, Sn);
   
   
   for (int i=0; i<10; i++){
      for (int j=0; j<i; j++){
         Fisher(i,j) = Fisher(j,i);
      }
   }
   
    
   fout4.close();
   delete [] dXm;
   delete [] dXeta;
   delete [] dXchi;
   delete [] dXthS;
   delete [] dXphS;
   delete [] dXtc;
   delete [] dXpsi;
   delete [] dXphi0;
   delete [] dXdl;
   delete [] dXiota;
   delete [] X;
   delete [] Xn;
   delete [] Y;
   delete [] Z;
}


void ComputeFisherC::ComputeWaveXYZ(BBHTemplate S, int n, double* &freq, std::complex<double>* &X, std::complex<double>* &Y,\
         std::complex<double>* &Z)
{
   
   double tc = S.tc;
   double ThetaQ = pow(0.2*S.eta/S.Mt*tc, -0.25);
   double x = 0.25*ThetaQ*( 1. + (743./4032. + 11.*S.eta/48.)*ThetaQ );
   double fmin = pow(x, 1.5)/(LISAWP_PI*S.Mt);
   double  Dl = S.dist * (LISAWP_PC_SI/LISAWP_C_SI);
   
   
   double* tm;
   tm = new double[n];

   double cpsi = cos(2.*S.psi);
   double spsi = sin(2.*S.psi);
   std::complex<double> tmp;
   std::complex<double> img(0.0, 1.0);
   double shiftPh, fr, fmax1;
   if(tc > Tobs){
       ThetaQ = pow(0.2*S.eta/S.Mt*(tc-Tobs), -0.25);
       x = 0.25*ThetaQ*( 1. + (743./4032. + 11.*S.eta/48.)*ThetaQ );
       fmax1 = pow(x, 1.5)/(LISAWP_PI*S.Mt);      
   }
   else{
       fmax1 = fmax;
   }
   tshift = Tobs - S.tc;
   
   PhenomCwave PW_C(fmin, fmax1, df);
   std::complex<double>* Hplus = NULL;
   std::complex<double>* Hcross = NULL;
   Hplus = new std::complex<double>[n];
   Hcross = new std::complex<double>[n];
   
   int sz = PW_C.ComputeHpHc(S, n, freq, Hplus, Hcross);
   // Apply polarization angle and timeshift
   for (int i=0; i<n; i++ )
   {
      fr=i*df;
      shiftPh = tshift*LISAWP_TWOPI*fr;
      tmp = Hplus[i];
      Hplus[i] = (Hplus[i]*cpsi + Hcross[i]*spsi)*(cos(shiftPh) + img*sin(shiftPh));
      Hcross[i] = (-tmp*spsi + Hcross[i]*cpsi)*(cos(shiftPh) + img*sin(shiftPh));
   }
   
   PW_C.ComputeTime(S, freq, tm, n);
   
   // now tdi...
   
   ComputeTDIfreq tdi(L, year);
   tdi.ChooseConfiguration(conf);
   tdi.fMin = fmin;
   tdi.fMax = fmax1;
   
   
   tdi.ComputeTDIfreqXYZ(n, S.thetaS, S.phiS, tm, freq, Hplus, Hcross, X, Y, Z);
   
   delete [] Hplus;
   delete [] Hcross;
   delete [] tm;
                        
}

double ComputeFisherC::ComputeInnerProd(int sz, std::complex<double>* &x, std::complex<double>* &y, double* &freq, double* &Sn)
{
   double prod2 = 0.0;
   double df;
   for (int i=1; i<sz; i++){
      df = freq[i]-freq[i-1];
      if (freq[i] >= fminS && freq[i] <= fmaxS){
           prod2 += df*(x[i] * conj(y[i])).real()/Sn[i];
      }
   }
   prod2 *= 4.;
   return(prod2);
}

} // end of the namespace