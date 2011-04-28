/*
Copyright (C) 2011  E. Robinson, S. Babak, F. Ohme

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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/include/Template.hh,v 1.5 2008/04/18 11:14:04 stba Exp $"
*/

#ifndef PHENOMCWAVEHH
#define PHENOMCWAVEHH

//#include <vector>
#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
//#include "Matrix.hh"
#include "Constants.hh"
#include "phenomwf.h"
#include "BBHTemplate.hh"
#include "fftw++.h"

namespace LISAWP{

  /** Template class PhenomCwave
  * Computes PhenomC waveform based on Santamaria et al: http://arxiv.org/abs/1005.3306  
  * it is a C++ wrapper class for phenomwf functions
  * @author E. Robinson, S. Babak, F. Ohme, 2011
  */

class PhenomCwave{
 
public:

  /*** Constructor 
   * @param Fmin min frequency
   * @param Fmax max frequency
   * @param delF frequency resolution
  */     
  PhenomCwave(double Fmin, double Fmax, double delF);

  /*** Computes Hplus, Hcross in radiative frame in freq. domain. Returns size of Hp, Hc
   * @param H BBHTemplate object (parameters)
   * @param Hp to be filled up with h+(f)
   * @param Hc to be filled up with hx(f)
  */
  int ComputeHpHc(BBHTemplate H, std::complex<double>* &Hp, std::complex<double>* &Hc);
  
  /*** Computes  H22, returns size of H22 */
  int ComputeH22(BBHTemplate H, std::complex<double>* &H22);
  
  /*** Computes the time[i] corresponding to the freq[i] using 3.5 PN 
   * @param H BBHtemplate object (parameters)
   * @param freq array of frequencies
   * @param tm array of time corresponding to freq (to be computed)
   * @param n size of the arrays
  */
  void ComputeTime(BBHTemplate H, double* &freq, double* &tm, int n);
  
  /*** Computes h+(t) and hx(t). Note that you need to specify shift in time to avoid folding 
     * Duation of the waveform (tc) should be less than T_obs-tShift
     * @param H BBHtemplate object (parameters)
     * @param HpT h+(t) to be filled up
     * @param HcT hx(t) to be filled up
     * @param tShift h(t) will be shifted in time by tSHift to avoid folding (iFFT of h(f) returns always waveform 
     which coalesces around t=T_obs).
  */
  int ComputeHofT(BBHTemplate H, double* &HpT, double* &HcT, double tShift);
  
private:
  
   double fMin;
   double fMax;
   double df;
  
   
};

}// END OF THE namespace

#endif
