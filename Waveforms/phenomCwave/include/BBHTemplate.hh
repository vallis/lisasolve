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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/include/Template.hh,v 1.5 2008/04/18 11:14:04 stba Exp $"
*/

#ifndef BBHTEMPLATEHH
#define BBHTEMPLATEHH

//#include <vector>
#include <math.h>
#include <complex>
//#include "Matrix.hh"
#include "Constants.hh"

namespace LISAWP{

  /** Template class 
  * @author S. Babak, 2011
  */

class BBHTemplate{
 
  public:

  /** Costructor 
  * @param Mtotal total mass in solar masses
  * @param eta reduced mass ratio
  * @param lambda in radian
  * @param beta in radian
  * @param Tc time to coalescence in seconds
  */

  BBHTemplate();
  
  BBHTemplate(const BBHTemplate& p);

  BBHTemplate operator=(const BBHTemplate& p);
  
  void ComputeAuxParams();  
  
  //Template(double mass1, double mass2, double lambda, double beta, double Tc);

  double M;
  double Mt;   // in seconds
  double Mc;   // in seconds
  double eta;
  double lam;
  double bet;
  double chi;
  double q;
  double tc;
  double m1;
  double m2;
  double thetaS;
  double phiS;
  double dist;
  double Ampl;
  double iota;
  double phi0;
  double psi;

  double Norm;
  double f0;
  double fEnd;
  double SNR;
  

//  private:

};
    

}//end of the namespace


#endif
