/*
Copyright (C) 2011  S. Babak

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
#include "ComputeFisher.hh"
#include "NoiseModels.hh"


/*** THis script computes inverse Fisher matrix for given setof parameters for PhenomC waveform  using X, */

using namespace LISAWP;

double ComputeInnerProd(int sz, double fmin, double fmax, std::complex<double>* &x, std::complex<double>* &y, double* &freq, double* &Sn);

int main(){
 
   BBHTemplate H;
   double arm = 5.e9;
   double year = 31457280.;
   std::string spr = "    ";
   
   int id_r, id_s;
   //#realID srcID z Mtot q betS lamS phL thL phi0 tc DL0(Mpc)  Mc  eta 
   double  z, M, q, phiS, thetaS, phiL, thetaL, phi0, tc, DL, Mc, eta;
   // lam = pi/2 - thetaS
   // beta = phiS 
   
   double Tobs = 2.*year;
   int nsources = 15360;
   double chi = 0.0;
   std::string config = "aLISA";
   
   double Fmin = 1.e-5;
   double Fmax = 1.0;
   double df = 1.e-6;
   int n = (int) floor (Fmax / df)+1;
   double* freq;
   freq = new double[n];
   for (int i=0; i<n; i++ )
   {
        freq[i]=(double)i*df;
   }
   double* tm;
   tm = new double[n];
   
   std::complex<double>* Hplus = NULL;
   std::complex<double>* Hcross = NULL;
   Hplus = new std::complex<double>[n];
   Hcross = new std::complex<double>[n];
   
   std::complex<double>* X;
   std::complex<double>* Y;
   std::complex<double>* Z;
   X = new  std::complex<double>[n];
   Y = new  std::complex<double>[n];
   Z = new  std::complex<double>[n];   
   
   std::cout << "size n = " << n << "  df =  " << df << std::endl;
   double* Sn;
   Sn = new double[n]; 
   
   
   bool galactic_bin = true;
   galactic_bin = false;
   NoiseModels NM1(galactic_bin);
   NM1.StandardLISA_X(n, freq, Sn);
   //NM1.miniLISA_C1X(n, freq, Sn);
   NM1.miniLISA_C1X(n, freq, Sn);

   std::ofstream fout23("Data/NoiseTestC3.dat");
   for (int i=0; i<n; i++){
        fout23 << std::setprecision(15) << freq[i] << spr << Sn[i] << std::endl;
   }
   fout23.close();
   exit(0);  


   // reading the noise:
   
   /*std::string noiseFile = "Data/HL1_TDIX.dat";
   int nsz = 190652;
   std::ifstream finN(noiseFile.c_str());
   if(!finN){
           std::cerr << "Cannot open the parameter file " << noiseFile << std::endl;
   	     exit(1);
   }	     
   double* ftmp;
   double* Stmp;
   ftmp = new double[nsz];
   Stmp = new double[nsz];
   for (int i=0; i<nsz; i++){
        finN >> ftmp[i] >> Stmp[i]; 
   }
   finN.close();
    
   gsl_interp_accel *acc  = gsl_interp_accel_alloc ();
   gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, nsz);
   gsl_spline_init (spline, ftmp, Stmp, nsz);
     
   for (int i=0; i<n; i++){
        Sn[i] = gsl_spline_eval(spline, freq[i], acc);
   }
     
   gsl_spline_free (spline);
   gsl_interp_accel_free (acc);
   delete [] ftmp;
   delete [] Stmp;
   
   */
   
   
   double x, ThetaQ, fmax1;
   double cpsi;
   double spsi;
   std::complex<double> tmp;
   std::complex<double> img(0.0, 1.0);
   double tshift, shiftPh, fr;
   double L = arm/LISAWP_C_SI;
   double SNR2;
   
   // reading the orbit
 /*  
   double t=0.0;
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
   for (int i=1; i<Orsz; i++){
         if (torb[i]-torb[i-1] <= 0.){
            std::cout << "negative " << i << spr << torb[i] << spr << torb[i-1] << std::endl;
         }
   }
   */
   ComputeTDIfreq tdi(L, year);
   
   
   // reading the data file
   std::ifstream fin("Data/model_HOR_MCevents_ext.dat");
  //std::ofstream fout("Data/model_HOR_MCevents_SNR_L1.dat");
   for (int i=0; i<nsources; i++){
       
      std::cout << std::endl;
      std::cout << " -------   " << i << "    ----------\n";
      std::cout << std::endl;
       
      fin >> id_r >> id_s >> z >> M >> q >> phiS >> thetaS >> phiL >> thetaL >>  phi0 >> tc >> DL >> Mc >> eta;
     // std::cout << id_r << spr << id_s << spr << z << spr << M << spr << q << spr << phiS << spr \
         << thetaS << spr <<  phiL << spr << thetaL << spr <<  phi0 << spr <<  tc << spr << DL << spr << Mc << spr << eta<< std::endl;
      //exit(0);
      
      H.M = M*(1.+z);
      H.m2 = 0.0;
      H.m1 = 0.0;
      H.q = q;
      H.chi = chi;
      H.dist = DL*1.e6;
      H.phi0 = phi0;
      H.thetaS = 0.5*LISAWP_PI - thetaS;
      H.phiS = phiS;
      H.tc = tc;
      
      H.ComputeAuxParams();
      std::cout << H.M << spr << H.q << spr <<  H.m1 << spr << H.m2 << spr << H.Mc/LISAWP_MTSUN_SI << spr << H.eta << std::endl; 
      H.iota = acos(-sin(H.thetaS)*sin(thetaL)*cos(H.phiS - phiL) - cos(H.thetaS)*cos(thetaL) );
      H.psi = atan2(cos(H.thetaS)*sin(thetaL)*cos(H.phiS - phiL) - sin(H.thetaS)*cos(thetaL), sin(thetaL)*sin(H.phiS - phiL) );
      //std::cout << H.iota << spr << H.psi << std::endl;
      //exit(0);
    
      ThetaQ = pow(0.2*H.eta/H.Mt*tc, -0.25);
      x = 0.25*ThetaQ*( 1. + (743./4032. + 11.*H.eta/48.)*ThetaQ );
      Fmin = pow(x, 1.5)/(LISAWP_PI*H.Mt);
   
      if(tc > Tobs){
          ThetaQ = pow(0.2*H.eta/H.Mt*(tc-Tobs), -0.25);
          x = 0.25*ThetaQ*( 1. + (743./4032. + 11.*H.eta/48.)*ThetaQ );
          fmax1 = pow(x, 1.5)/(LISAWP_PI*H.Mt);      
      }
      else{
          fmax1 = Fmax;
      }
      tshift = Tobs - H.tc;
      std::cout << "freq range: " << Fmin << spr << fmax1 << std::endl;
      
      PhenomCwave PW_C(Fmin, fmax1, df);

      int sz = PW_C.ComputeHpHc(H, n, freq, Hplus, Hcross);
      
      // Apply polarization angle and timeshift
      cpsi = cos(2.*H.psi);
      spsi = sin(2.*H.psi);
      for (int i=0; i<n; i++ )
      {
         fr=i*df;
         shiftPh = tshift*LISAWP_TWOPI*fr;
         tmp = Hplus[i];
         Hplus[i] = (Hplus[i]*cpsi + Hcross[i]*spsi)*(cos(shiftPh) + img*sin(shiftPh));
         Hcross[i] = (-tmp*spsi + Hcross[i]*cpsi)*(cos(shiftPh) + img*sin(shiftPh));
         //std::cout << abs(Hplus[i]) << spr;
      }
      
      std::cout << " tc = " << tc << std::endl;
      PW_C.ComputeTime(H, freq, tm, n);
      //std::cout << "h+, hx computed \n";

      // now tdi...
      //std::cout << "L = " << L << std::endl;
      
      tdi.ChooseConfiguration(config);
      tdi.fMin = Fmin;
      tdi.fMax = fmax1;

      tdi.ComputeTDIfreqXYZ_NumOrb(n, H.thetaS, H.phiS, tm, freq, Hplus, Hcross,\
                              Orsz, torb,  x1, y1, z1, x2, y2, z2, x3, y3, z3,\
                              X, Y, Z);
      //std::cout << "TDI computed \n";
      
      // no need to go beyond highest frequency, let's find it
      for (int i=0; i<n; i++){
         //std::cout << abs(Hplus[i]) << spr;
         if (freq[i] > Fmin ){
            if (abs(X[i]) == 0.0){
               fmax1 = freq[i];
               break;
            }
         }
      }
      std::cout << "freq range: " << Fmin << spr << fmax1 << std::endl;
      
      SNR2 = ComputeInnerProd(n, Fmin, fmax1, X, X, freq, Sn);   
      std::cout << "SNR^2 = " << SNR2 << std::endl;
      
   /*   tdi.ComputeLWfreqXYZ(n, H.thetaS, H.phiS, tm, freq, Hplus, Hcross, Y, Z, Z);
      
     SNR2 = ComputeInnerProd(n, Fmin, fmax1, Y, Y, freq, Sn);   */
   /*   std::cout << "SNR^2 = " << SNR2 << std::endl;
      std::ofstream fout423("Data/TestLW.dat");
      for (int i=0; i<n; i++){
         fout423 << std::setprecision(15) << freq[i] << spr << abs(X[i]) << spr << abs(Y[i]) << std::endl;
      }
      fout423.close();
   //   exit(0);   */
   
      fout << std::setprecision(15) << id_r << spr << id_s << spr << z << spr << M << spr << q << spr << phiS << spr \
            << thetaS << spr <<  phiL << spr << thetaL << spr <<  phi0 << spr <<  tc << spr << DL << spr << Mc << spr \
             << eta << spr << sqrt(SNR2) << std::endl;
   
      
   }// end of the source loop
   
   fout.close();
   
   delete [] freq;
   delete [] Hplus;
   delete [] Hcross;
   delete [] tm;
   delete [] X;
   delete [] Y;
   delete [] Z;
   delete [] Sn;
   delete [] x1;
   delete [] x2;
   delete [] x3;
   delete [] y1;
   delete [] y2;
   delete [] y3;
   delete [] z1;
   delete [] z2; 
   delete [] z3;
   delete [] torb;  
   
}


double ComputeInnerProd(int sz, double fmin, double fmax, std::complex<double>* &x, std::complex<double>* &y, double* &freq, double* &Sn)
{
   
   double prod2 = 0.0;
   double df;
   for (int i=1; i<sz; i++){
      df = freq[i]-freq[i-1];
      if (freq[i] >= fmin && freq[i] <= fmax){
           prod2 += df*(x[i] * conj(y[i])).real()/Sn[i];
      }
   }
   prod2 *= 4.;
   return(prod2);
   
}
