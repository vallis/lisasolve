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

#ifndef NOISEMODELSHH
#define NOISEMODELSHH

#include "Constants.hh"
#include <vector>
#include <complex>


/** Class NoiseModels
* it computes models of the  noise for different configurations
* @author Stas Babak 2011
*/

namespace LISAWP{
	class NoiseModels{
	public:
		/** Constructor 
		* @param arm armlength in seconds
		* @param oneyear convention for one year in seconds
		*/
		NoiseModels(bool galactic_bin);
		
		/** Standard heliocentric LISA motion, armlength is a variable*/
      void StandardLISA_X(int sz, double* &freq, double* &S_n);
   	
   	// Not complete yet                         
      void miniLISA_C1X(int sz, double* &freq, double* &S_n);
      
      // Not complete yet
      void miniLISA_C2X(int sz, double* &freq, double* &S_n);
      
      void miniLISA_C4X(int sz, double* &freq, double* &S_n);
      
      void miniLISA_C5X(int sz, double* &freq, double* &S_n);
                                                     
		
		
	private:
   		bool gal;

};

} //end of the namespace

#endif