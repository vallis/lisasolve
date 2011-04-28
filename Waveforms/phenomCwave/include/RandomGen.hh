/*
Copyright (C) 2006  S. Babak

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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/include/RandomGen.hh,v 1.1 2007/09/18 21:50:38 stba Exp $" 
*/
#ifndef RANDOMGENHH
#define RANDOMGENHH

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"


/**** Class RandomGen: gsl wrapper
 *  computes numbers distributed according to the special law
 *  with the given range
 *  @author Stas Babak, June 2006
 */
 
namespace LISAWP{


class RandomGen{

   public:
      /*** Constructor
       * @param seed seed for the random number generator
       */
	RandomGen(int seed);
	
      /*** Destructor */
	~RandomGen();

      /*** Uniformly distribution
       * @param a lower bound
       * @param b upper bound
       */
       double  Uniform(double a, double b);
       
       double Gauss(double sigma);

   private:

     const gsl_rng_type *T;
     gsl_rng *r;



};

}

#endif

