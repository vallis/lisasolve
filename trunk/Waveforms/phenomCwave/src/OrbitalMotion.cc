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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/src/LISAmotion.cc,v 1.2 2007/12/07 10:36:41 stba Exp $" 
*/


#include "OrbitalMotion.hh"

namespace LISAWP{
	
OrbitalMotion::OrbitalMotion(double arm, double oneyear){
		L = arm;
		year = oneyear;
		//Omega = LISAWP_TWOPI/oneyear
		Omega = LISAWP_TWOPI*3.1709791983764586e-8;
		RAU = LISAWP_AU_SI/LISAWP_C_SI;
}


void OrbitalMotion::EccentricLISAMotion(double kappa0, double lambda0, double t, double* &R, double** &p, double** &n)
{

		double alpha = Omega*t + kappa0;
		double xi[3]; // this is beta
		xi[0] = lambda0;             //////////  MOVE THOSE CALCULATIONS OUT TO THE CONSTRUCTOR
		xi[1] = lambda0 + 2.0*LISAWP_PI/3.0;
		xi[2] = lambda0 + 4.0*LISAWP_PI/3.0;
 
		
		R[0] = RAU*cos(alpha);
		R[1] = RAU*sin(alpha);
		R[2] = 0;
		
      double q[3][3];
		
		double ecc = 0.00964838; // hardcoded
		
      double fct = 1./(2.0*sqrt(3.));
		for(int i=0; i<3; i++){
   			q[i][0] = ( sin(alpha)*cos(alpha)*sin(xi[i]) - \
   				(1.+ sin(alpha)*sin(alpha))*cos(xi[i])  )*fct;
   			q[i][1] = ( sin(alpha)*cos(alpha)*cos(xi[i]) - \
   				(1.+ cos(alpha)*cos(alpha))*sin(xi[i])  )*fct;
   		   q[i][2] = -0.5*cos(alpha - xi[i]);
   		   
            p[i][0] = R[0] + L*q[i][0];
            p[i][1] = R[1] + L*q[i][1];
            p[i][2] = R[2] + L*q[i][2];
                        
   	}
   	
		
	/*	for(int i=0; i<3; i++){
			p[i][0] = R[0] + RAU*ecc*( sin(alpha)*cos(alpha)*sin(xi[i]) - \
				(1.+ sin(alpha)*sin(alpha))*cos(xi[i])  );
			p[i][1] = R[1] + RAU*ecc*( sin(alpha)*cos(alpha)*cos(xi[i]) - \
				(1.+ cos(alpha)*cos(alpha))*sin(xi[i])  );
		   p[i][2] = -sqrt(3.)*RAU*ecc*cos(alpha - xi[i]);
		}*/
		
		for(int i=0; i<3; i++){
			 n[0][i] = (q[1][i] - q[2][i]);
		    n[1][i] = (q[2][i] - q[0][i]);
		    n[2][i] = (q[0][i] - q[1][i]);
		}
    	   	   
}

void OrbitalMotion::L1ToyModel(double T_rot, double t, double* &R, double** &p, double** &n)
{
   double xi = Omega*t;
   double eta = LISAWP_TWOPI/T_rot*3.1709791983764586e-8; // assumes that T_rot is given in years
   
   double q[3][3];
   double O2[3][3];
   double fact = L/sqrt(3.0);
   q[0][0] = 0.;
   q[0][1] = fact*0.5*sqrt(3.);
   q[0][2] = -0.5*fact;
   q[1][0] = 0.;
   q[1][1] = -fact*0.5*sqrt(3.);
   q[1][2] = -0.5*fact; 
   q[2][0] = 0.;
   q[2][1] = 0.;
   q[2][2] = fact;
   
   R[0] = 0.99*RAU*cos(xi);
	R[1] = 0.99*RAU*sin(xi);
	R[2] = 0.;
   
   
   O2[0][0] = cos(xi);
   O2[0][1] = -sin(xi) * cos(eta);
   O2[0][2] = -sin(xi) * sin(eta);
   O2[1][0] = sin(xi);
   O2[1][1] =  cos(xi) * cos(eta);
   O2[1][2] = cos(xi) * sin(eta);
   O2[2][0] = 0.0;
   O2[2][1] = -sin(eta);
   O2[2][2] = cos(eta);
   
   double qO = 0.0;
   for(int i=0; i<3; i++){
      qO = 0.0;
      for(int j=0; j<3; j++){
         for(int ii=0; ii<3; ii++){
             qO += O2[i][ii]*q[i][ii];
         }
      }
      
   }
   
   
}
	
}// end of the namespace