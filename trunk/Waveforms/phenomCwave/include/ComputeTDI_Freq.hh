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


#ifndef  COMPUTETDIFREQHH
#define  COMPUTETDIFREQHH

#include "Constants.hh"
#include "OrbitalMotion.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>


/** Class ComputeTDIfreq
* it computes X, Y, Z TDIs for a given h+(f), hx(f), orbit and t(f).
* It assumes rigid moving triangle and adiabatic approximation 
* @author Stas Babak 2011
*/

namespace LISAWP{
class ComputeTDIfreq{
 
 	public:
 	   
       ComputeTDIfreq(double arm, double oneyear);
      ~ComputeTDIfreq();
      
      /*** Computes X,Y,Z tdis in freq domain. Here I assume that the the array of size sz are allocated
       *   I assume that t(f) is computed, which is tf[i] = t(freq[i]), and corresponding hfp (plus polar) hfc (cross polar)
       * @param sz size of the arrays
       * @param thetaS co-latitude of the source
       * @param phiS   longitude of the source
       * @param tf time at the fourier frequencies defined in freq
       * @param freq fourier frequencies
       * @param hfp SPA of hplus
       * @param hfc SPA of hcross
       * @param X, Y, Z to be filled up with corresponding TDI variables 
       ***/
      void ComputeTDIfreqXYZ(int sz, double thetaS, double phiS, double* &tf, double* &freq, \
                              std::complex<double>* &hfp, std::complex<double>* &hfc, \
                              std::complex<double>* &X, std::complex<double>* &Y, std::complex<double>* &Z);
                              
      /*** Computes long wavelength approximation, parameters are the same as for ComputeTDIfreqXYZ */                        
      void ComputeLWfreqXYZ(int sz, double thetaS, double phiS, double* &tf, double* &freq, \
                              std::complex<double>* &hfp, std::complex<double>* &hfc, \
                              std::complex<double>* &X, std::complex<double>* &Y, std::complex<double>* &Z);           
                              
      /*** This is where you should say what LISA orbit you want to use, sofar only "aLISA" is available */                        
      void ChooseConfiguration(std::string config);          
      
      /*** Generate TDI signal only between fMin and fMax ***/
      double fMin;
      double fMax;              
      
   private:
      
      double* uhat;
      double* vhat;
      double* k;
      double* kp;
      double* kn;
      double* u;
      double* v;
      double* R;
      double** p;
      double** n;
      std::string conf;
      double L;
      double year;
      
      double SINC(double x);
      
      
};   
} //end of the namespace
#endif