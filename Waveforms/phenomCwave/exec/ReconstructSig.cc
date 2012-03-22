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
#include "Kallianpur_Striebel.h"


/*** This scripts takes the data which contains a signal and trying to reconstruct it 
 * using best known model and kallianpur-Striebel formula. */

using namespace LISAWP;

int main(int argc, char* argv[]){
   
    if (argc != 3){
        std::cerr << "Usage: <DataFile> <FileOut>" << std::endl;
        exit(1);
    }
    std::string spr  = "    ";
    std::cout.precision(10);
    std::complex<double> img(0.0, 1.0); 	
    std::string FileDat = argv[1]; // file with parameters
    std::string FileOut = argv[2]; 
 

   /* later I need to change this to read the data, but I start with 
    * using phenomC waveform generated here as a data. This is for debugging
    */

    double Fmin = 0.0035; // in units Mf for the numrel data
    double Fmax = 0.1;    // in units Mf for the numrel data
    double Tobs = 31457280.101232;
    double tc = Tobs-5.e6;  

    Fmin = 1.e-4;
    Fmax = 0.2;
    double df = 1.e-5;

    BBHTemplate H;
   
    H.M = 1.e6;
    H.q = 0.25;
    H.chi = 0.5;
    H.tc = tc;
    H.phi0 = 0.1;
    H.dist = 1.e-4;
    H.ComputeAuxParams();



    int n = (int) floor (Fmax / df)+1;
    double* freq;
    freq = new double[n];
    for (int i=0; i<n; i++ )
    {
         freq[i]=(double)i*df;
    }
    std::cout << "size n = " << n << "  df =  " << df << std::endl;

    PhenomCwave PW_C(Fmin, Fmax, df);
    std::complex<double>* h22 = NULL;
    h22 = new std::complex<double>[n];
    PW_C.ComputeH22(H, n, freq, h22);
    
    /*std::ofstream fout1446("Data/TestPhenomC.dat");
    for (int i=0; i<n; i++){
        fout1446 << std::setprecision(15) << freq[i] << spr << abs(h22[i]) << std::endl;  
    }

    fout1446.close();*/

    double tshift = Tobs - H.tc + 2054.56;
    double fr, shiftPh;
    std::cout << " tshift = " << tshift << std::endl;
    std::ofstream fout1446("Data/TestPhenomC.dat");
    for (int i=0; i<n; i++ )
    {
      fr=i*df;
      shiftPh = tshift*LISAWP_TWOPI*fr;
      //shiftPh = 0.0;
      //std::cout << shiftPh/LISAWP_TWOPI << spr << cos(shiftPh) << spr << sin(shiftPh) << std::endl; 
      h22[i] = h22[i]*(cos(shiftPh) + img*sin(shiftPh));
      fout1446 << std::setprecision(15) << freq[i] << spr << abs(h22[i]) << \
	      spr << h22[i].real() << spr << h22[i].imag() << std::endl;  
    }
    fout1446.close();
   
    KS_reconstruct ksr(h22, freq, df, Fmin, Fmax, n);
    
    ksr.SetPsd("white", 1.e-4);
    ksr.SetBaseWaveform("phenomC");
    ksr.SetIntegrationMethod("vegas");

    double* bndrs;
    bndrs = new double[6];
    bndrs[0] = 0.5e6;
    bndrs[1] = 0.2;
    bndrs[2] = 0.4;
    bndrs[3] = 1.5e6;
    bndrs[4] = 0.3;
    bndrs[5] = 0.6;
    ksr.SetBoundaries(bndrs, 3);

    ksr.ComputeNormInt();


    delete [] freq;
    delete [] h22;
    delete [] bndrs;

    return(0);

}
   

