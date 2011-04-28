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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/src/RandomGen.cc,v 1.1 2007/09/18 21:50:43 stba Exp $" 
*/

///**********************************************************
///  returns uniformly distributed (with "seed") value 
///   in range "a" to "b"
//***********************************************************  

#include "RandomGen.hh"

namespace LISAWP{

RandomGen::RandomGen(int seed){
 
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

}

RandomGen::~RandomGen(){

    gsl_rng_free(r);

}

double RandomGen::Uniform(double a, double b){
    
    double u;
    u = gsl_rng_uniform(r);
    return ( a + (b-a)*u );
	
}

double RandomGen::Gauss(double sigma){
	
	double u;
	u = gsl_ran_gaussian(r, sigma);
	return(u);
}

}//end of the namespace
