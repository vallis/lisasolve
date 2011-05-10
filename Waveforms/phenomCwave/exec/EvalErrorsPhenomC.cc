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


/*** THis script takes as an input catalog of SMBHB observations and computes inverse Fisher matrix. It uses 
PhenomC waveform only one detector, X, and  makes cuts on the mass ration (1/20). It assumes alligned spins with L*/

using namespace LISAWP;

int main(int argc, char* argv[]){
   
    if (argc != 3){
        std::cerr << "Usage: <FileParams> <FileOut>" << std::endl;
        exit(1);
    }
    std::string spr  = "    ";
    std::cout.precision(10);
   
    BBHTemplate H;
    std::string FilePar = argv[1]; // file with parameters
    std::string FileOut = argv[2]; 
    
    // reading parameter file 
    
    //realID srcID z Mtot q betS lamS phL thL phi0 tc a1 a2 thS1 thS2 phS1 phS2 per e DL0(Mpc)  Mc  eta  thBS1 thBS2 phBS1 phBS2
    
    std::string FileSource; // file with source parameters
    int Ntot;               // total number of sources
    int realId;             // which realization we should use: [1,10]
    double q_th = 0.05;     // threshold on mass ratio, default 1/20
    double SNRth = 8.0;     // threshold on SNR for detection
    
    // Hardcodded values
    
    double Fmin = 1.e-6;
    double Fmax = 0.5;
    double df = 1.e-6;
    double year = 31457280.;
    double Tobs = 2.*year;
    
    std::string noise;
    
    
    
    
    std::ifstream fin1(FilePar.c_str());
    if(!fin1){
         std::cerr << "Cannot open the parameter file " << FilePar << std::endl;
 	      exit(1);
    }
    fin1 >> FileSource;
    fin1 >> Ntot;
    fin1 >> realId;
    fin1 >> q_th;
    fin1 >> SNRth;
    fin1 >> noise;
    
    fin1.close();
    
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
     double arm;
    
     if (noise == "LISA"){
        NM1.StandardLISA_X(n, freq, S_n);
        std::cout << "using standard LISA \n";
        arm = 5.e9/LISAWP_C_SI;
     }
     if (noise == "C1"){
        NM1.miniLISA_C1X(n, freq, S_n);
        std::cout << "using C1 configuration \n";
        arm = 1.e9/LISAWP_C_SI;
     }
     if (noise == "C2"){
        NM1.miniLISA_C2X(n, freq, S_n);
        std::cout << "using C2 config \n";
        arm = 1.e9/LISAWP_C_SI;
     }
     
    std::string config = "aLISA";
    ComputeFisherC FishC(arm, year, config, Fmax, df, Tobs);
   
    
    // reading source params
    int real_id, src_id;
    double z, M1, q, betS, lamS, phL, thL, phi0, tc, a1, a2;
    double thS1, thS2, phS1, phS2, per, e, DL0,  Mc,  eta,  thBS1, thBS2, phBS1, phBS2;
    double phiS, thetaS, M2, Mtot, chi;
    std::ifstream fin2(FileSource.c_str());
    if(!fin2){
          std::cerr << "Cannot open the source file " << FileSource << std::endl;
  	      exit(1);
    }
    double SNR2;
    int dim = 10;
    Matrix<double> Fisher(dim,dim);
    Matrix<double> IFisher(dim, dim);
    gsl_matrix *m = gsl_matrix_alloc(dim, dim);
    gsl_vector *rhs = gsl_vector_alloc(dim);
    gsl_vector *x = gsl_vector_alloc(dim);
    gsl_permutation * p = gsl_permutation_alloc(dim);
    int signum, info;
    gsl_matrix *m1 = gsl_matrix_alloc(dim, dim);
    gsl_vector *residual = gsl_vector_alloc(dim);
    
    
    std::ofstream fout(FileOut.c_str());
    for (int i=0; i< Ntot; i++)
    {
       fin2 >> real_id >> src_id >> z >> M1 >> q >> phiS >> lamS >> phL >> thL >> phi0 >> tc >> a1 >> a2;
       fin2 >> thS1 >> thS2 >> phS1 >> phS2 >> per >> e >> DL0 >>  Mc >> eta >> thBS1 >> thBS2 >> phBS1 >> phBS2;
       thetaS = 0.5*LISAWP_PI - lamS;
       DL0 *= 1.e6;
       M1 *= (1.+z);
       M2 = q*M1;
       Mtot = M1+M2;
      // std::cout <<  real_id << spr << src_id << spr << z << spr << M1 << spr << q << spr << phiS \
       << spr << thetaS << spr << phL << spr << thL << spr << phi0 << spr << tc << spr << a1 << spr << a2\
       << spr << thS1 << spr << thS2 << spr << phS1 << spr << phS2 << spr << per << spr << e << spr << DL0\
       << spr <<  Mc << spr << eta << spr << thBS1 << spr << thBS2 << spr << phBS1 << spr << phBS2 << std::endl;
       //exit(0);
       chi = M1/Mtot *a1 + M2/Mtot *a2;
       
       
       H.m1 = M1;
       H.m2 = M2;
       H.M = Mtot;
       H.q = q;
       H.chi = chi;
       H.dist = DL0;
       H.phi0 = phi0;
       H.thetaS = thetaS;
       H.phiS = phiS;
       H.tc = tc;
       
       H.iota = acos(-sin(H.thetaS)*sin(thL)*cos(H.phiS - phL) - cos(H.thetaS)*cos(thL) );
       H.psi = atan2(cos(H.thetaS)*sin(thL)*cos(H.phiS - phL) - sin(H.thetaS)*cos(thL), sin(thL)*sin(H.phiS - phL) );
       
      // std::cout << "before\n";
      // std::cout << H.M << spr << H.q << spr <<  H.m1 << spr << H.m2 << spr << H.Mc/LISAWP_MTSUN_SI << spr << H.eta << std::endl;
       H.ComputeAuxParams();
      // std::cout << "after\n";
      // std::cout << H.M << spr << H.q << spr <<  H.m1 << spr << H.m2 << spr << H.Mc/LISAWP_MTSUN_SI << spr << H.eta << std::endl;
      // std::cout << H.iota << spr << H.psi << spr << H.q << std::endl;
      // exit(0);
       std::cout << M1/Mtot *a1 << spr <<  M2/Mtot *a2 << spr << a1 << spr << a2 << std::endl; 
       std::cout << "M =  " << H.M << "  mass ratio = " << H.q << "  spin =   " << H.chi << std::endl;
       if (real_id == realId && H.q >=  q_th){
         SNR2 = FishC.ComputeRAFisher4links(H, n, freq, S_n, Fisher);
         std::cout << "SNR^2 = "<< SNR2 << std::endl;
        // exit(0);
      
         if (sqrt(SNR2) >= SNRth){
             for (int i=0; i<dim; i++){
                  gsl_vector_set(rhs, i, 0.0);
     	            for (int j=0; j<dim; j++){
                    gsl_matrix_set(m, i, j, Fisher(i,j));
     	            }	
             }
             gsl_matrix_memcpy(m1, m);
             info = gsl_linalg_LU_decomp (m, p, &signum);
      
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
             fout << std::setprecision(15) <<  real_id << spr << src_id << spr << z << spr << M1 << spr << q << spr << phiS \
              << spr << thetaS << spr << phL << spr << thL << spr << phi0 << spr << tc << spr << a1 << spr << a2\
              << spr << thS1 << spr << thS2 << spr << phS1 << spr << phS2 << spr << per << spr << e << spr << DL0\
              << spr <<  Mc << spr << eta << spr << thBS1 << spr << thBS2 << spr << phBS1 << spr << phBS2 << spr \
              << SNR2 << spr << sqrt(IFisher(0,0))/H.M << spr << sqrt(IFisher(1,1)) << spr << sqrt(IFisher(2,2)) << spr \
              << sqrt(IFisher(3,3)) << spr << sqrt(IFisher(4,4)) << spr << sqrt(IFisher(5,5)) << spr << sqrt(IFisher(6,6)) << spr\
              << sqrt(IFisher(7,7)) << spr << sqrt(IFisher(8,8)) << spr << sqrt(IFisher(9,9)) << std::endl;
         }// end of SNR-if
       }// end of the q-if
      //exit(0);
      
    }// end of the signal loop
    
    
    fin2.close();
    fout.close();
    
    gsl_permutation_free (p);
    gsl_vector_free (rhs);
    gsl_vector_free (x);
    gsl_matrix_free (m);
    gsl_matrix_free (m1);
    gsl_vector_free (residual);

    delete [] freq;
    delete [] S_n;
   
}