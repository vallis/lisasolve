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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/include/LISAmotion.hh,v 1.1 2007/09/18 21:50:37 stba Exp $" 
*/

#ifndef ORBITALMOTIONHH
#define ORBITALMOTIONHH

#include "Constants.hh"
#include <vector>
#include <complex>


/** Class OrbitalMotion
* it computes motion of the spaccecrafts and arms
* @author Stas Babak 2011
*/

namespace LISAWP{
	class OrbitalMotion{
	public:
		/** Constructor 
		* @param arm armlength in seconds
		* @param oneyear convention for one year in seconds
		*/
		OrbitalMotion(double arm, double oneyear);
		
		/** Standard heliocentric LISA motion, armlength is a variable*/
      void EccentricLISAMotion(double kappa0, double lambda0, double t, double* &R, double** &p, double** &n);
   	
   	// Not complete yet                         
      void GeocentricL345();
      
      // Not complete yet
      void L1ToyModel(double T_rot, double t, double* &R, double** &p, double** &n);
                                                     
		
		
	private:
   		double L;
   		double year;
         double Omega;
         double RAU;

};

} //end of the namespace

#endif