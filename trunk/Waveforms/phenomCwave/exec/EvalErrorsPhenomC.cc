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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


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
    
    
    
    std::string noiseFile;
    int nsz;
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
    fin1 >> noise >> noiseFile >> nsz;
    
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
     if (noise == "L1"){
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
           S_n[i] = gsl_spline_eval(spline, freq[i], acc);
        }
        
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
        delete [] ftmp;
        delete [] Stmp;
        
     }
     
    /* std::ofstream fout2254("Data/NoiseTest.dat");
     for (int i=0; i<n; i++){
        fout2254 << std::setprecision(15) << freq[i] << spr << S_n[i] << std::endl;
     }
     fout2254.close();
     exit(0); */
     
     
     // Check orbital interpolation
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
    
/*     OrbitalMotion nom(1.e9/LISAWP_C_SI, year);
    
    
    nom.InitiateSplinInterpolation(Orsz, torb, x1, y1, z1, x2, y2, z2, x3, y3, z3);
    
    std::ofstream fout2330("Data/CheckOrbitalINterp.dat");
    double** pos;
    double** nn;
    pos = new double*[3];
    nn = new double*[3];
    for (int i=0; i<3; i++){
       pos[i] = new double[3];
       nn[i] = new double[3];
    }
    double L0, L1, L2;
    while (t<torb[Orsz-1]){
       nom.NumericalData(t, pos, nn, L0, L1, L2);
       fout2330 << std::setprecision(15) << t << spr << pos[0][0] << spr << pos[0][1] << spr << pos[0][2] << spr <<\
               pos[1][0] << spr << pos[1][1] << spr << pos[1][2] << spr << pos[2][0] << spr << pos[2][1] << spr << pos[2][2] <<\
               spr << nn[0][0] << spr << nn[0][1] << spr << nn[0][2] << spr << nn[1][0] << spr << nn[1][1] << spr << nn[1][2] << spr <<\
               nn[2][0] << spr << nn[2][1] << spr << nn[2][2] << std::endl;
       t += 150.;
    }
    fout2330.close();
  */  
    
       // exit(0); 
    std::string config = "aLISA";
    config = "L1";
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
         std::cout << H.M << spr << H.eta << spr << H.chi << spr << H.thetaS << spr << H.phiS << spr <<\
            H.tc << spr << H.phi0 << spr << H.iota << spr << H.psi << spr << H.dist << std::endl;
      // std::cout << H.iota << spr << H.psi << spr << H.q << std::endl;
      // exit(0);
       std::cout << M1/Mtot *a1 << spr <<  M2/Mtot *a2 << spr << a1 << spr << a2 << std::endl; 
       //std::cout << "M =  " << H.M << "  mass ratio = " << H.q << "  spin =   " << H.chi << std::endl;
       if (real_id == realId && H.q >=  q_th){
         //SNR2 = FishC.ComputeRAFisher4links(H, n, freq, S_n, Fisher);
         SNR2 = FishC. ComputeFisher4links_NumOrb(H, n, freq, S_n, Orsz, torb,  x1, y1, z1, x2, y2, z2, x3, y3, z3, Fisher);
         std::cout << "SNR^2 = "<< SNR2 << std::endl;
         //exit(0);
      
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
     // exit(0);
      
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