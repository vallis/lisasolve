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


#ifndef  COMPUTEFISHERHH
#define  COMPUTEFISHERHH

#include "Constants.hh"
#include "ComputeTDI_Freq.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include "Matrix.hh"
#include "BBHTemplate.hh"
#include "ComputeFisher.hh"
#include "PhenomCwave.hh"


/** Class ComputeFsherC
* it computes Fisher matrix using supplied noise curve and 
* using RA or LW response (computed in freq. domain) and Phenom C waveform
* @author Stas Babak 2011
*/

namespace LISAWP{
class ComputeFisherC{
 
 	public:
 	   
       ComputeFisherC(double arm, double oneyear, std::string config, double Fmax, double delf, double duration);
      ~ComputeFisherC();
      void ComputeRAFisher4links(BBHTemplate S, int n, double* &freq, double* &Sn, Matrix<double> &Fisher);
      void ComputeRAFisher6links(BBHTemplate S, int n, double* &freq, double* &Sn, double** &Fisher);
      void ComputeLWFisher4links(BBHTemplate S, int n, double* &freq, double* &Sn, double** &Fisher);
      void ComputeLWFisher6links(BBHTemplate S, int n, double* &freq, double* &Sn, double** &Fisher);
      
      
   private:
      
      double L;   
      double year;
      std::string conf;
      double df;
      double fmax;
      double fmaxS;
      double fminS;
      double tshift;
      double Tobs;
      int n; // size of the frequency arrays
      bool print;
      
      void ComputeWaveXYZ(BBHTemplate S, int n, double* &freq, std::complex<double>* &X, std::complex<double>* &Y,\
         std::complex<double>* &Z);
         
      double ComputeInnerProd(int sz, std::complex<double>* &x, std::complex<double>* &y, double* &freq, double* &Sn);   
      
};
}// end of the namespace 
#endif 