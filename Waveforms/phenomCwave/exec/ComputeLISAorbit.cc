#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include "Matrix.hh"
#include "BBHTemplate.hh"
#include "OrbitalMotion.hh"
#include "Constants.hh"

/*** THis script computes X, Y, Z for phenomC waveform in freq. domain and transforms it to time domain */

using namespace LISAWP;

int main(){
   
   double day = 24.0*3600.;  // day in sec. <- use it as a stepsize
   double duration = 2.0*365.; // duration in days
   int N = floor(duration);
   double arm = 2.e9/LISAWP_C_SI; // arm in sec.
   double year = 365.*day;
   std::string spr  = "    ";
   
   std::ofstream fout("Data/Orbit2years.dat");
   
   double* R;
   double** p;
   double** n;
   
   R = new double[3];
   p = new double*[3];
   n = new double*[3];

   for (int i=0; i<3; i++){
       p[i] = new double[3];
       n[i] = new double[3];
   }
   
   OrbitalMotion orbit(arm, year);
   double t = 0.0;
   for (int i=0; i<N+1; i++){
      
      orbit.EccentricLISAMotion(0.0, 0.0, t, R, p, n);
      fout << std::setprecision(15) << t << spr;
      for (int ic=0; ic<3; ic++){
         for (int ix=0; ix<3; ix++){
            fout << p[ic][ix]*LISAWP_C_SI << spr;
         }
      }
      fout << "\n";
      t += day;
   }
   fout.close();
   
   
   delete [] R;
   for (int i=0; i<3; i++){
      delete [] p[i];
      delete [] n[i];
   }
   delete p;
   delete n; 
   
}