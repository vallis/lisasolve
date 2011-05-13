/*
Copyright (C) 2011  E. Robinson, S. Babak

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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/toolbox/tools/UniformEMRIsteps.cc,v 1.2 2007/11/14 23:49:31 stba Exp $" 
*/

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include "Matrix.hh"
#include "BBHTemplate.hh"
#include "PhenomCwave.hh"
#include "Constants.hh"
#include "fftw++.h"
#include "ComputeTDI_Freq.hh"

/*** THis script computes X, Y, Z for phenomC waveform in freq. domain and transforms it to time domain */

using namespace LISAWP;

int main(){
   
   BBHTemplate H;
   double arm = 1.e9/LISAWP_C_SI; // LISA's arm in sec
   std::cout << "arm = " << arm << std::endl;
   double year = 31457280.;
   std::string spr = "    ";
   
   H.m1 = 1.e4;
   H.m2 = 1.e4;
   H.chi = 0.4;
   H.dist = 6.e9;
   H.iota = 0.87;//LISAWP_PI/3.;
   H.phi0 = 0.3;
   H.psi = 1.1;
   H.thetaS = 0.5*LISAWP_PI - 1.0;
   H.phiS = 2.0;
   H.dist = 6.e9;
   H.tc = 0.7*year;
   std::cout << "beta = " <<  H.thetaS  << std::endl;
   
   
   H.ComputeAuxParams();
   
   double Fmin = 1.e-5;
   double Fmax = 1.e0;
   double df = 1.e-6;
   double Tobs = year;
   double tc = H.tc;
   
   int n = (int) floor (Fmax / df)+1;
   double* freq;
   freq = new double[n];
   for (int i=0; i<n; i++ )
   {
      freq[i]=(double)i*df;
   }
   
   double ThetaQ = pow(0.2*H.eta/H.Mt*H.tc, -0.25);
   double x = 0.25*ThetaQ*( 1. + (743./4032. + 11.*H.eta/48.)*ThetaQ );
   double fmin = pow(x, 1.5)/(LISAWP_PI*H.Mt);
   double  Dl = H.dist * (LISAWP_PC_SI/LISAWP_C_SI);
   
   double* tm;
   tm = new double[n];

   double cpsi = cos(2.*H.psi);
   double spsi = sin(2.*H.psi);
   std::complex<double> tmp;
   std::complex<double> img(0.0, 1.0);
   double shiftPh, fr, fmax1;
   if(tc > Tobs){
       ThetaQ = pow(0.2*H.eta/H.Mt*(H.tc-Tobs), -0.25);
       x = 0.25*ThetaQ*( 1. + (743./4032. + 11.*H.eta/48.)*ThetaQ );
       fmax1 = pow(x, 1.5)/(LISAWP_PI*H.Mt);      
   }
   else{
       fmax1 = Fmax;
   }
   double tshift = Tobs - H.tc;
   
   
   PhenomCwave PW_C(fmin, fmax1, df);
   std::complex<double>* Hplus = NULL;
   std::complex<double>* Hcross = NULL;
   Hplus = new std::complex<double>[n];
   Hcross = new std::complex<double>[n];
   
   int sz = PW_C.ComputeHpHc(H, n, freq, Hplus, Hcross);
   // Apply polarization angle and timeshift
   //std::ofstream foutH("Data/TestH.dat");
   for (int i=0; i<n; i++ )
   {
      fr=i*df;
     //  foutH << fr << "    " << abs(Hplus[i]) << std::endl;
      shiftPh = tshift*LISAWP_TWOPI*fr;
      tmp = Hplus[i];
      Hplus[i] = (Hplus[i]*cpsi + Hcross[i]*spsi)*(cos(shiftPh) + img*sin(shiftPh));
      Hcross[i] = (-tmp*spsi + Hcross[i]*cpsi)*(cos(shiftPh) + img*sin(shiftPh));
      
   }
   //foutH.close();
   
   PW_C.ComputeTime(H, freq, tm, n);
   
   
   // reading the numerical Orbit
   
   double* x1;
   double* x2;
   double* x3;
   double* y1;
   double* y2;
   double* y3;
   double* z1;
   double* z2; 
   double* z3;
   double* torb;
    
   int Orsz = 17364;
   torb = new double[Orsz];
   x1 = new double[Orsz];
   y1 = new double[Orsz];
   z1 = new double[Orsz];
   x2 = new double[Orsz];
   y2 = new double[Orsz];
   z2 = new double[Orsz];
   x3 = new double[Orsz];
   y3 = new double[Orsz];
   z3 = new double[Orsz];
   
   std::ifstream finOr("Data/Orbits_HaloL1-Pos.txt");
   for (int i=0; i<Orsz; i++){
       finOr >> torb[i] >> x1[i] >> y1[i] >> z1[i] >> x2[i] >> y2[i] >> z2[i] >> x3[i] >> y3[i] >> z3[i];
   }
   finOr.close();
   
   
   // now tdi...
   
   ComputeTDIfreq tdi(arm, year);
   tdi.ChooseConfiguration("aLISA");
   tdi.fMin = fmin;
   tdi.fMax = fmax1;
   
   std::complex<double>* X;
   std::complex<double>* Y;
   std::complex<double>* Z;
   X = new  std::complex<double>[n];
   Y = new  std::complex<double>[n];
   Z = new  std::complex<double>[n];
   
   tdi.ComputeTDIfreqXYZ_NumOrb(n, H.thetaS, H.phiS, tm, freq, Hplus, Hcross, Orsz, torb, x1, y1, z1, x2, y2, z2, x3, y3, z3, X, Y, Z);
   
    for (int i=0; i<n; i++){
         if (freq[i] > fmin ){
            if (abs(X[i]) == 0.0){
                fmax1 = freq[i];
                break;
            }
         }
      }
       std::cout << "Freq. range of the signal (after): " << fmin << "   " << fmax1 << std::endl; 
   
   
   std::ofstream fout("Data/TDI_L1_F.dat");
   
   for (int i=0; i<n; i++){
      fr=(double)i*df;
      fout << std::setprecision(15) << fr << "    "  << abs(X[i]) << "   " << abs(Y[i]) << "   " << abs(Z[i]) <<  std::endl;
   }
   fout.close();
   
   
   delete [] Hplus;
   delete [] Hcross;
   delete [] tm;
   delete [] freq;
   
   delete [] torb;
   delete [] x1;
   delete [] y1;
   delete [] z1;
   delete [] x2;
   delete [] y2;
   delete [] z2;
   delete [] x3;
   delete [] y3;
   delete [] z3;
   delete [] X;
   delete [] Y;
   delete [] Z;
   
   
   return(0);
   
}