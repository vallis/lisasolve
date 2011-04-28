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
   double Dl = 6.e9; // 6 Gpc
   double arm = 5.e9/LISAWP_C_SI; // LISA's arm in sec
   std::cout << "arm = " << arm << std::endl;
   double year = 31457280.;
   std::string spr = "    ";
   
   H.m1 = 1.e6;
   H.m2 = 1.e6;
   H.chi = 0.0;
   H.dist = Dl;
   H.iota = 0.0;//LISAWP_PI/3.;
   H.phi0 = 0.0;
   H.psi = 0.0;
   H.thetaS = 0.5*LISAWP_PI - 1.0;
   H.phiS = 2.0;
   std::cout << "beta = " <<  H.thetaS  << std::endl;
   
   
   H.ComputeAuxParams();
   
   double Fmin = 1.e-4;
   double Fmax = 1.e-1;
   double df = 1.e-5;
   double dt = 15.0;
   //dt = dt/2.;
   double Tobs = 31457280.;
   //Tobs *= 0.5;
   
   Fmax = 1./(2.*dt); // assumes samling rate 15sec
   df = 1./Tobs;  // 1/year
   
   // Estimate roughly F_min based on tc:
   double tc = Tobs-5.e6; // almost 1 year
   double tshift = 1.e6;
   double ThetaQ = pow(0.2*H.eta/H.Mt*tc, -0.25);
   double x = 0.25*ThetaQ*( 1. + (743./4032. + 11.*H.eta/48.)*ThetaQ );
   Fmin = pow(x, 1.5)/(LISAWP_PI*H.Mt);
   std::cout  << "Fmin = " << Fmin << std::endl;
   std::cout << "Fmax = " << Fmax << std::endl;
   
   //!!!!!! NEED to produce warning if the sampling rate is too low !!!!!!!
   
   
   std::cout << "chirp mass = " << H.Mc << std::endl;
   Dl *= (LISAWP_PC_SI/LISAWP_C_SI); // dist in sec;
   
   PhenomCwave PW(Fmin, Fmax, df);
   
   std::complex<double>* Hplus = NULL;
   std::complex<double>* Hcross = NULL;
   
   int n = PW.ComputeHpHc(H, Hplus, Hcross);
   std::cout << "n = " << n << std::endl;
   
   double* freq;
   freq = new double[n];
   double* tm;
   tm = new double[n];
   
   for (int i=0; i<n; i++ )
   {
      freq[i]=i*df;
      tm[i] = 0.0;
   }
   //std::cout << freq[10] << "   " << tm[0] << "  " << H22[1] << std::endl;
   PW.ComputeTime(H, freq, tm, n);
   std::ofstream fout27("Data/TimeTest.dat");
   for (int i=0; i<n; i++ )
   {
      tm[i] += (5.e6-tshift);
      fout27 << std::setprecision(15) << (double)i*dt << spr << tm[i] << spr << freq[i] << std::endl;
      
   }
   fout27.close();
   // apply polarization rotation and time shift
   double cpsi = cos(2.*H.psi);
   double spsi = sin(2.*H.psi);
   std::complex<double> tmp;
   std::complex<double> img(0.0, 1.0);
   double shiftPh, fr;
   
  
   for (int i=0; i<n; i++ )
   {
      fr=i*df;
      shiftPh = tshift*LISAWP_TWOPI*fr;
      tmp = Hplus[i];
      Hplus[i] = (Hplus[i]*cpsi + Hcross[i]*spsi)*(cos(shiftPh) + img*sin(shiftPh));
      Hcross[i] = (-tmp*spsi + Hcross[i]*cpsi)*(cos(shiftPh) + img*sin(shiftPh));
   }
   
   int datSz = 2*(n-1);
   
   
 /*  std::ofstream fout26("Data/HplHcr.dat");
   double* HpT;
   double* HcT;
   HpT = new double[datSz];
   HcT = new double[datSz];
   crfft1d Backward(datSz, Hplus, HpT);
   Backward.fft(Hplus, HpT);
   Backward.fft(Hcross, HcT);
   for (int i=0; i<datSz; i++){
        fout26 << std::setprecision(15) << (double)i*dt << spr  << HpT[i] << spr << HcT[i] << std::endl;
   }
   fout26.close();
   */
   
   // Compute TDI in freq domain
   
   ComputeTDIfreq fTDI(arm, year);
    
   std::complex<double>* X;
   std::complex<double>* Y;
   std::complex<double>* Z;
   X = new  std::complex<double>[n];
   Y = new  std::complex<double>[n];
   Z = new  std::complex<double>[n];
   
   fTDI.ChooseConfiguration("aLISA");
   //fTDI.ComputeTDIfreqXYZ(n,  H.thetaS, H.phiS, tm, freq, Hplus, Hcross, X, Y, Z);
   // Or one can compute Long wavelength limit:
   fTDI.ComputeLWfreqXYZ(n,  H.thetaS, H.phiS, tm, freq, Hplus, Hcross, X, Y, Z);
   std::cout << "TDIs computed\n";
   
   // shifting the waveform...
   /*for (int i=0; i<n; i++){
      
      X[i] = X[i] * (cos(shiftPh) + img*sin(shiftPh));
   }*/
   
   // transforming in time domain
   
   double* XT;
   XT = new double[datSz];
   double* YT;
   YT = new double[datSz];
   double* ZT;
   ZT = new double[datSz];
   crfft1d Backward(datSz, X, XT);
   Backward.fft(X, XT);
   Backward.fft(Y, YT);
   Backward.fft(Z, ZT);
   std::cout << "fft computed, recording ...\n";
   
   
   std::ofstream fout("Data/TDI_T.dat");
   
   for (int i=0; i<datSz; i++){
      fout << std::setprecision(15) << (double)i*dt << "    "  << XT[i] << "   " << YT[i] << "   " << ZT[i] <<  std::endl;
   }
   fout.close();
   
   //delete [] HpT;
   //delete [] HcT;
   delete [] X;
   delete [] Y;
   delete [] Z;
   delete [] freq;
   delete [] tm;
   delete [] Hplus;
   delete [] Hcross; 

   return(0);
   
}