/*
Copyright (C) 2005  S. Babak

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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/include/EMRItemplate.hh,v 1.2 2008/04/18 11:48:08 stba Exp $"
*/

#ifndef EMRITEMPLATEHH
#define EMRITEMPLATEHH

#include <vector>
#include <math.h>
#include <complex>
#include "Matrix.hh"
#include "Constants.hh"

namespace LISAWP{

  /** Template class 
  * @author S. Babak, 2006
  */

class EMRItemplate{
 
  public:

  /** Costructor 
  * @param Mtotal total mass in solar masses
  * @param eta reduced mass ratio
  * @param lambda in radian
  * @param beta in radian
  * @param Tc time to coalescence in seconds
  */

  EMRItemplate(double Mass, double mu, double spin, double lambda);
  
//  void SetInitialConditions(double nu, double e, double Phi, double gamma, double alpha);
  
  void SetPosition(double thetaS, double phiS, double thetaK, double phiK, double D);
  
  EMRItemplate(const EMRItemplate& p);
  
  EMRItemplate operator=(const EMRItemplate& p);
  
//  EMRItemplate operator=(const EMRItemplate&);
  
  double M;
  double Mt;
  double m;
  double mt;
  double a;
  double lam;
  
  double e0;
  double nu0;
  double Phi0;
  double gamma0;
  double alpha0;
//  double tStart; // start of observation 
  double t0;  // instance of time when initial conditions are defined
  double fgam0;
  double falph0;
  
  double tPl;
  double e_pl;
  double nu_pl;
  double Phi_pl;
  double alpha_pl;
  double gamma_pl;
  
  double thS;   // co-latittude !!!
  double phS;
  double thK;
  double phK;
  double stS, stK, ctS, ctK; // cos and sin of thetaS and thetaK
  double cpS, spS, cpK, spK; // cos and sin of phiS, and phiK
  
  double dist; // distance in seconds
  double Ampl; // dimensionless amplitude: mu/D
  bool SkySet;
  
  double SNR;
  double LogL;
  
  Matrix<double> XAalpha;
  Matrix<double> YEbeta;
  Matrix<double> ZTgamma;
  
};

}// end of the namespace

#endif
