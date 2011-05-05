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

/**** This script computes H+ and Hx for phenomC waveform and transforms it to time domain */

using namespace LISAWP;

int main(){
   
   BBHTemplate H;
   double Dl = 6.e9; // 6 Gpc
   
   H.m1 = 1.e6;
   H.m2 = 1.e6;
   H.chi = 0.0;
   H.dist = Dl;
   H.iota = 0.0;//LISAWP_PI/3.;
   H.phi0 = 0.0;
   
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
   double tc = Tobs-5.e6; // 1 year
   double ThetaQ = pow(0.2*H.eta/H.Mt*tc, -0.25);
   double x = 0.25*ThetaQ*( 1. + (743./4032. + 11.*H.eta/48.)*ThetaQ );
   Fmin = pow(x, 1.5)/(LISAWP_PI*H.Mt);
   std::cout  << "Fmin = " << Fmin << std::endl;
   std::cout << "Fmax = " << Fmax << std::endl;
   
   
   std::cout << "chirp mass = " << H.Mc << std::endl;
   Dl *= (LISAWP_PC_SI/LISAWP_C_SI); // dist in sec;
   
   PhenomCwave PW(Fmin, Fmax, df);
   
  // std::complex<double>* H22;
   //int n = PW.ComputeH22(H, H22);
   
   std::complex<double>* Hplus;
   std::complex<double>* Hcross;

   int n = (int) floor (Fmax / df)+1;
   double* freq;
   freq = new double[n];
   for (int i=0; i<n; i++ ){
       freq[i]=(double)i*df;
   }
   
   int ck = PW.ComputeHpHc(H, n, freq, Hplus, Hcross);
   std::cout << "n = " << n << "  " << ck << std::endl;
   
   
   double* tm;
   tm = new double[n];
   
   for (int i=0; i<n; i++ )
   {
      tm[i] = 0.0;
   }
   //std::cout << freq[10] << "   " << tm[0] << "  " << H22[1] << std::endl;
   PW.ComputeTime(H, freq, tm, n);
   
   std::complex<double>* H22;    
   ck =  PW.ComputeH22(H, n, freq,  H22);
   
   //double amp1 = abs(H22[10])*pow(freq[10], 7./6.);
   
   //std::cout << "amp = " << amp1 << std::endl;
   //double norm = fact/amp1; 
   
  double fact = sqrt(5./6.)*0.25*pow(LISAWP_PI, -2./3.)*((1.+pow(cos(H.iota), 2.))/Dl) *pow(H.Mc, 5./6.);
  std::cout << "factor = " << fact << std::endl;
  std::string spr  = "   ";
  std::ofstream fout2("Data/waveF.dat");
  double shiftPh;
  std::complex<double> img(0.0, 1.0);
    for (int i=0; i<n; i++){
       fout2 << std::setprecision(15) << freq[i] << spr << abs(Hplus[i]) << spr << real(Hplus[i]) << spr << abs(H22[i]) << std::endl;
      // shiftPh = 1.e6*LISAWP_TWOPI*freq[i];
       //Hplus[i] = Hplus[i] * (cos(shiftPh) + img*sin(shiftPh));
    }       
    fout2.close();
  
  double* HpT;
  double* HcT;
  /*int datSz = 2*(n-1);
  HpT = new double[datSz];
  HcT = new double[datSz];
  std::cout << "computing FFT....." << std::endl;
  crfft1d Backward(datSz, Hplus, HpT);
  Backward.fft(Hplus, HpT);
 // Backward.fft(H22, HpT);
  std::cout << "done " << std::endl; 
  */
  
  std::cout << "computing hoft \n";
  int datSz = PW.ComputeHofT(H, HpT, HcT, 1.e6);
  std::cout << "done\n";
  std::complex<double>* Htest;
  Htest = new std::complex<double>[n];
 
   
   std::ofstream fout("Data/waveT.dat");
   
   for (int i=0; i<datSz; i++){
      //std::cout << freq[i] << spr << tm[i] << spr << abs(H22[i])*norm << std::endl;
      //fout << std::setprecision(15) << freq[i] << spr << tm[i] << spr << abs(H22[i])*norm << std::endl;
      //fout << std::setprecision(15) << freq[i] << spr << tm[i] << spr << abs(Hplus[i]) << spr << abs(Hcross[i]) <<\
               spr << fact*pow(freq[i], -7./6.) << std::endl;
               fout << std::setprecision(15) << (double)i*dt << spr  << HpT[i] << spr << HcT[i] << std::endl;
      
   }
   fout.close();
   
   std::cout << "computing Forward FFT ..... " << std::endl; 
   rcfft1d Forward(datSz, HpT, Htest);
   Forward.fftNormalized(HpT, Htest);
   std::cout << "done " << std::endl;
   std::ofstream fout3("Data/waveF_2.dat");
   for (int i=0; i<n; i++){
          fout3 << std::setprecision(15) << (double)i*df << spr << abs(Htest[i]) << spr << real(Htest[i]) << std::endl;
   }       
   fout3.close();
   
   
   
    
   
   // check
   tc = 5.*H.Mt/( 256.*H.eta * pow(LISAWP_PI*H.Mt*Fmin, 8./3.) );
   std::cout << "Tc = " << tc << std::endl;
   
   
   //delete  [] H22;
   delete  [] Hplus;   
   delete  [] Hcross;
   delete [] HpT;
   delete [] Htest;
   delete [] freq;
   delete [] tm;
   
}

