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
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


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
      
      void EccentricLISAMotion2(double kappa0, double lambda0, double t, double* &R, double** &q, double** &n);
   	
   	// Not complete yet                         
      void GeocentricL345();
      
      // Not complete yet
      void L1ToyModel(double T_rot, double t, double* &R, double** &p, double** &n);
                                                     
		void NumericalData(double t, double** &p, double** &n, double &L0, double &L1, double &L2);
		
		void InitiateSplinInterpolation(int sz, double* &torb,  double* &x1, double* &y1, double* &z1, double* &x2, double* &y2, double* &z2,\
            double* &x3, double* &y3, double* &z3);
      
      void FinalizeInterpolation();
		
	private:
   		double L;
   		double year;
         double Omega;
         double RAU;
         
         gsl_interp_accel *acc1;          
         gsl_spline *spline1;
         
         gsl_interp_accel *acc2;          
         gsl_spline *spline2;
         
         gsl_interp_accel *acc3;          
         gsl_spline *spline3;
         
         gsl_interp_accel *acc4;
         gsl_spline *spline4;
         
         gsl_interp_accel *acc5;
         gsl_spline *spline5;
         
         gsl_interp_accel *acc6;
         gsl_spline *spline6;
         
         gsl_interp_accel *acc7;
         gsl_spline *spline7;
         
         gsl_interp_accel *acc8;
         gsl_spline *spline8;
         
         gsl_interp_accel *acc9;
         gsl_spline *spline9;
         
};

} //end of the namespace

#endif
