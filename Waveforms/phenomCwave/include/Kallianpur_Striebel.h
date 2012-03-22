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


#ifndef  KSHH
#define  KSHH

#include "Constants.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include "Matrix.hh"
#include "BBHTemplate.hh"
#include "PhenomCwave.hh"
#include "Macros.hh"
#include "fftw++.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

/** Class KS_reconstruct
* it computes Bayesian estinmator for the waveform usnig the Kallianpur-Striebel 
* originally uses phenom-C waveforms as basis
* @author Stas Babak 2012
*/


namespace LISAWP{

double CompLhd(double *theta, size_t dim, void* params);	
double gg (double *k, size_t dim, void *params);


class KS_reconstruct{
   public:
     /** Constructor. the data must be in the frequency domain in the standard 
      * complex format with each point corresponding to the frequencies in 
      * the freqs array */	   
     KS_reconstruct(std::complex<double>* &input, double* freqs, double delf, \
		     double FMin, double FMax, int length);
     ~KS_reconstruct();

     /** give a noise type corresponding to the described in 
      * NoiseModels.hh */
     void SetPsd(std::string noiseType, double scale);

     /** Base waveform to be used 
      * SO far codded up only phenomC */
     void SetBaseWaveform(std::string Model);

     /** Method to be used for the N-dim integration:
      * MC, miser, vegas */
     void SetIntegrationMethod(std::string IntType);
 
     /** Here leftlower corner, upper right corner, assumes rectangular
      *  the order of parameters M, q, chi, phi0, tc, Amp */
     void SetBoundaries(double* &bndrs, int N);

     /** Computes likelihood normalization */
     double ComputeNormInt();

     /** computes estimated/reconstructed signal */
     void ComputeEstimateSgnl(std::complex<double>* &res);

   private:
     
     std::complex<double>* data;
     std::complex<double>* tmplt;
     std::string typeInt;
     double df;
     int sz;
     double fMin, fMax;
     double* Snf;
     double* freq;
     double* bndr;
     bool maxPhi0, maxTc, maxAmp;
     double ComputeInnerProduct(int n, std::complex<double>* &x, std::complex<double>* &y);
     double ComputeInnerProductV(int n, std::complex<double>* &x, std::complex<double>* &y, \
		     int &lag);
     double MaxTcPhi(int n, std::complex<double>* &x, std::complex<double>* &y, \
		     int &lag, double& phi0);
     void ComputePhenomWave(int n, BBHTemplate S, std::complex<double>* wave);  
     double ComputeL(double *theta, size_t dim, void* params);
};

} //end of the namespace

#endif


