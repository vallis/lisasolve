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
		for(int i=0; i<3; i++){ // s/c index
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

void OrbitalMotion::InitiateSplinInterpolation(int sz, double* &torb,  double* &x1, double* &y1, double* &z1, double* &x2, double* &y2, double* &z2,\
                     double* &x3, double* &y3, double* &z3)
{
   
   acc1 = gsl_interp_accel_alloc ();
   spline1 = gsl_spline_alloc (gsl_interp_cspline, sz);
   gsl_spline_init (spline1, torb, x1, sz);
   
   acc2 = gsl_interp_accel_alloc ();
   spline2 = gsl_spline_alloc (gsl_interp_cspline, sz);
   gsl_spline_init (spline2, torb, y1, sz);
   
   acc3 = gsl_interp_accel_alloc ();
   spline3 = gsl_spline_alloc (gsl_interp_cspline, sz);
   gsl_spline_init (spline3, torb, z1, sz);
   
   acc4 = gsl_interp_accel_alloc ();
   spline4 = gsl_spline_alloc (gsl_interp_cspline, sz);
   gsl_spline_init (spline4, torb, x2, sz);
   
   acc5 = gsl_interp_accel_alloc ();
   spline5 = gsl_spline_alloc (gsl_interp_cspline, sz);
   gsl_spline_init (spline5, torb, y2, sz);
   
   acc6 = gsl_interp_accel_alloc ();
   spline6 = gsl_spline_alloc (gsl_interp_cspline, sz);
   gsl_spline_init (spline6, torb, z2, sz);
   
   acc7 = gsl_interp_accel_alloc ();
   spline7 = gsl_spline_alloc (gsl_interp_cspline, sz);
   gsl_spline_init (spline7, torb, x3, sz);
   
   acc8 = gsl_interp_accel_alloc ();
   spline8 = gsl_spline_alloc (gsl_interp_cspline, sz);
   gsl_spline_init (spline8, torb, y3, sz);
   
   acc9 = gsl_interp_accel_alloc ();
   spline9 = gsl_spline_alloc (gsl_interp_cspline, sz);
   gsl_spline_init (spline9, torb, z3, sz);
   
}

void OrbitalMotion::NumericalData( double t, double** &p, double** &n, double &L0, double &L1, double &L2)
{                     
    // Note the first spacecraft should be mather in X configuration
    //X[i] = -4.*img*sin(omL)*(-y1_32*ex2 - y231*ex1 + y123*ex2 + y3_21*ex1);
    // so it is communication between 1<->2 and 1<->3: x1, y1, z1 should be a mother
    
   
   
   p[0][0] = gsl_spline_eval (spline1, t, acc1)/LISAWP_C_SI;
   
   p[0][1] = gsl_spline_eval (spline2, t, acc2)/LISAWP_C_SI;  
   
   p[0][2] = gsl_spline_eval (spline3, t, acc3)/LISAWP_C_SI;
   
   p[1][0] = gsl_spline_eval (spline4, t, acc4)/LISAWP_C_SI;
   
   p[1][1] = gsl_spline_eval (spline5, t, acc5)/LISAWP_C_SI; 

   p[1][2] = gsl_spline_eval (spline6, t, acc6)/LISAWP_C_SI;
   
   p[2][0] = gsl_spline_eval (spline7, t, acc7)/LISAWP_C_SI;
   
   p[2][1] = gsl_spline_eval (spline8, t, acc8)/LISAWP_C_SI;
   
   p[2][2] = gsl_spline_eval (spline9, t, acc9)/LISAWP_C_SI;
   
   L0 = 0.0;
   L1 = 0.0;
   L2 = 0.0;
   for(int i=0; i<3; i++){
			 n[0][i] = (p[1][i] - p[2][i]);
		    n[1][i] = (p[2][i] - p[0][i]);
		    n[2][i] = (p[0][i] - p[1][i]);
          L0 += n[0][i]*n[0][i];
          L1 += n[1][i]*n[1][i];
          L2 += n[2][i]*n[2][i];
	}
   L0 = sqrt(L0);
   L1 = sqrt(L1);
   L2 = sqrt(L2);
	for(int i=0; i<3; i++){
			 n[0][i] = n[0][i]/L0;
		    n[1][i] = n[1][i]/L1;
		    n[2][i] = n[2][i]/L2;
	}
   
   
}

void OrbitalMotion::FinalizeInterpolation()
{
   gsl_spline_free (spline1);
   gsl_interp_accel_free (acc1);
   gsl_spline_free (spline2);
   gsl_interp_accel_free (acc2);
   gsl_spline_free (spline3);
   gsl_interp_accel_free (acc3);
   gsl_spline_free (spline4);
   gsl_interp_accel_free (acc4);
   gsl_spline_free (spline5);
   gsl_interp_accel_free (acc5);
   gsl_spline_free (spline6);
   gsl_interp_accel_free (acc6);
   gsl_spline_free (spline7);
   gsl_interp_accel_free (acc7);
   gsl_spline_free (spline8);
   gsl_interp_accel_free (acc8);
   gsl_spline_free (spline9);
   gsl_interp_accel_free (acc9);
   
}
	
}// end of the namespace