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
#include "LinAlg.hh"
#include <gsl/gsl_linalg.h>


/*** THis script computes inverse Fisher matrix for given setof parameters for PhenomC waveform  using X, */

using namespace LISAWP;

int main(){
   
   BBHTemplate H;
   double Dl = 6.e9; // 6 Gpc
   double arm = 1.e9/LISAWP_C_SI; // LISA's arm in sec
   std::cout << "arm = " << arm << std::endl;
   double year = 31457280.;
   std::string spr = "    ";
   
   H.m1 = 1.e6;
   H.m2 = 1.e6;
   H.chi = 0.9;
   H.dist = Dl;
   H.iota = 1.2;//LISAWP_PI/3.;
   H.phi0 = 0.6;
   H.psi = 1.324;
   H.thetaS = 0.5*LISAWP_PI - 1.0;
   H.phiS = 2.0;
   H.tc = year/2.;
   std::cout << "beta = " <<  H.thetaS  << std::endl;  
   
   H.ComputeAuxParams();
   
   double Fmin = 1.e-4;
   double Fmax = 1.e-1;
   double df = 1.e-5;
   double dt = 15.0;
   //dt = dt/2.;
   double Tobs = 31457280.*2;
   //Tobs *= 0.5;
   
   Fmax = 1./(2.*dt); // assumes samling rate 15sec
   Fmax = 0.5;  //    !!!!!   check if we need to go further.   !!!!
   df = 1./Tobs;  // 1/year
   df = 1.e-6;
   
   // Compute noise:
  // std::ofstream fout("Data/CheckNoise.dat");
   int n = (int) floor (Fmax / df)+1;
   double* freq;
   freq = new double[n];
   for (int i=0; i<n; i++ )
   {
        freq[i]=(double)i*df;
   }
   std::cout << "size n = " << n << "  df =  " << df << std::endl;
   double* S_n;
   S_n = new double[n]; 
   
   bool galactic_bin = true;
   NoiseModels NM1(galactic_bin);
   //NM1.StandardLISA_X(n, freq, S_n);
   NM1.miniLISA_C2X(n, freq, S_n);
   
   double* S_n2;
   S_n2 = new double[n];
   galactic_bin = false;
   NoiseModels NM2(galactic_bin);
  // NM2.StandardLISA_X(n, freq, S_n2);
   //NM2.miniLISA_C2X(n, freq, S_n2);
   
//   for (int i=0; i<n; i++){
//      fout << std::setprecision(15) << freq[i] << spr << S_n[i] << spr << S_n2[i] << std::endl;
//   }
//   fout.close();
   
   std::string config = "aLISA";
   ComputeFisherC FishC(arm, year, config, Fmax, df, Tobs);
  
   Matrix<double> Fisher(10,10);
   FishC.ComputeRAFisher4links(H, n, freq, S_n, Fisher);
   
   std:: cout << "Fisher: \n";
   std::cout << Fisher;
   
   
   // Try that:
    int dim = 10;
   
   LinAlg la1(Fisher);
   Matrix<double> U;
   Matrix<double> S;
   Matrix<double> V;
   
   la1.SVDecomp(U, S, V);
   
   std::cout <<   " Sigma matrix is "<< std::endl;
   std::cout << S;
   
   Matrix<double> Gam1(((U*S)*V.Transpose()));
   Matrix<double> IS(dim, dim);
   IS = 0.;
   for (int i=0; i<dim; i++){
      IS(i,i) = 1./S(i,i);
//         std::cout << S(i,i)<< std::endl;
   }
   
   Matrix<double> IFisher(((V*IS)*U.Transpose()));
   //std::cout << "check: \n";
   Matrix<double> Un(Fisher*IFisher); // closer to unity
   std::cout << Un;
   std::cout << " ==========================  Inverse Fisher is =============  \n";
   std::cout << IFisher;
   
   std::cout << "sigma M/M: " << sqrt(IFisher(0,0))/H.M << std::endl;
   std::cout << "sigma eta: " << sqrt(IFisher(1,1)) << std::endl;
   std::cout << "sigma chi: " << sqrt(IFisher(2,2)) << std::endl;
   std::cout << "sigma thetaS: " << sqrt(IFisher(3,3)) << std::endl;
   std::cout << "sigma phiS: " << sqrt(IFisher(4,4)) << std::endl;
   std::cout << "sigma Tc: " << sqrt(IFisher(5,5)) << std::endl;
   std::cout << "sigma psi: " << sqrt(IFisher(6,6)) << std::endl;
   std::cout << "sigma phi0: " << sqrt(IFisher(7,7)) << std::endl;
   std::cout << "sigma iota: " << sqrt(IFisher(8,8)) << std::endl;
   std::cout << "sigma DL/DL: " << sqrt(IFisher(9,9)) << std::endl;
  
   gsl_matrix *m = gsl_matrix_alloc(dim, dim);
   gsl_vector *rhs = gsl_vector_alloc(dim);
   gsl_vector *x = gsl_vector_alloc(dim);
   
   for (int i=0; i<dim; i++){
      gsl_vector_set(rhs, i, 0.0);
	   for (int j=0; j<dim; j++){
           gsl_matrix_set(m, i, j, Fisher(i,j));
	   }	
   }
   gsl_permutation * p = gsl_permutation_alloc(dim);
   int signum;
   //gsl_vector_set(rhs, 0, 1.0);
   gsl_matrix *m1 = gsl_matrix_alloc(dim, dim);
   gsl_matrix_memcpy(m1, m);
   int info = gsl_linalg_LU_decomp (m, p, &signum);
//   gsl_linalg_LU_solve(m, p, rhs, x); // third method
   gsl_vector *residual = gsl_vector_alloc(dim);
//   info = gsl_linalg_LU_refine(m1, m, p, rhs, x, residual);
   
   for (int i=0; i<dim; i++){
  //    std::cout << i << "     "  <<gsl_vector_get(x, i) << "   " << gsl_vector_get(residual, i)<< std::endl;;   
      for (int j=0; j<dim; j++){
         gsl_vector_set(rhs, j, 0.0);
         if (i ==j )
            gsl_vector_set(rhs, i, 1.0);         
      }   
      gsl_linalg_LU_solve(m, p, rhs, x);
      info = gsl_linalg_LU_refine(m1, m, p, rhs, x, residual);
      for (int j=0; j<dim; j++){
         IFisher(j, i) = gsl_vector_get(x, j); 
      }      
   }
   
   std::cout << "refined inverse Fisher matrix \n";
   std::cout << IFisher; 
   
   std::cout << "\n check \n";
   std::cout << Fisher*IFisher << std::endl;
   
   std::cout << "sigma M/M: " << sqrt(IFisher(0,0))/H.M << std::endl;
     std::cout << "sigma eta: " << sqrt(IFisher(1,1)) << std::endl;
     std::cout << "sigma chi: " << sqrt(IFisher(2,2)) << std::endl;
     std::cout << "sigma thetaS: " << sqrt(IFisher(3,3)) << std::endl;
     std::cout << "sigma phiS: " << sqrt(IFisher(4,4)) << std::endl;
     std::cout << "sigma Tc: " << sqrt(IFisher(5,5)) << std::endl;
     std::cout << "sigma psi: " << sqrt(IFisher(6,6)) << std::endl;
     std::cout << "sigma phi0: " << sqrt(IFisher(7,7)) << std::endl;
     std::cout << "sigma iota: " << sqrt(IFisher(8,8)) << std::endl;
     std::cout << "sigma DL/DL: " << sqrt(IFisher(9,9)) << std::endl;
   
   
   gsl_permutation_free (p);
   gsl_vector_free (rhs);
   gsl_vector_free (x);
   gsl_matrix_free (m);
   gsl_matrix_free (m1);
   
   delete [] freq;
   delete [] S_n;
   delete [] S_n2;
   
}