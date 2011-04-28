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

#ifndef LISAMOTIONHH
#define LISAMOTIONHH

#include "Constants.hh"
#include <vector>
#include <complex>
#include "Matrix.hh"

/** Class LISAmotion
* it computes motion of the spaccecrafts and arms
* @author Stas Babak 2007
*/

namespace LISAWP{
	class LISAmotion{
	public:
		/** Constructor 
		* @param arm armlength in seconds
		* @param oneyear convention for one year in seconds
		*/
		LISAmotion(double arm, double oneyear);
		
		/** Computes orbital elemnets of LISA 
          * @param kappa0 initial angle of GC
          * @param lambda0 initial angle of triangle
          * @param t time at which position is evaluated
          * @param R vector-position of the GC (computed)
          * @param q matrix 3 x vector-position of the spacecraft w.r.t.
          * the GC (computed)
		  * @param n matrix 3 x link-unit_vector (computed) 
          * */
        void EccentricLISAMotion(float kappa0, float lambda0, double t, Matrix<double>& R, \
		                         Matrix<double>& q, Matrix<double>& n);
	void EccentricLISAMotion2(float kappa0, float lambda0, double t, Matrix<double>& R, \
	                         Matrix<double>& p, Matrix<double>& n);
		
		void CircularLISAMotion(float kappa0, float lambda0, double t, Matrix<double>& R, \
							                         Matrix<double>& q, Matrix<double>& n);
	private:
		double L;
		double year;
		
		
	};
	
	
}


#endif
