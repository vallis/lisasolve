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


#include "LISAmotion.hh"

namespace LISAWP{
	
	LISAmotion::LISAmotion(double arm, double oneyear){
		L = arm;
		year = oneyear;
	}
	void LISAmotion::EccentricLISAMotion(float kappa0, float lambda0, double t, Matrix<double>& R, \
	                         Matrix<double>& q, Matrix<double>& n){
	//	double Omega = LISAWP_TWOPI/year;
	    double Omega = LISAWP_TWOPI*3.1709791983764586e-8;
		double alpha = Omega*t + kappa0;
		Matrix<double> xi(1,3);
		xi(0) = lambda0;             //////////  MOVE THOSE CALCULATIONS OUT TO THE CONSTRUCTOR
		xi(1) = lambda0 + 2.0*LISAWP_PI/3.0;
		xi(2) = lambda0 + 4.0*LISAWP_PI/3.0;
		double RAU = LISAWP_AU_SI/LISAWP_C_SI;
		
		R(0) = RAU*cos(alpha);
		R(1) = RAU*sin(alpha);
		R(2) = 0;
		
		for(int i=0; i<3; i++){
			/*q(i,0) = ( cos(2.0*alpha - xi(i)) - 3.0*cos(xi(i)) )/(2.0*sqrt(12));
			q(i,1) = ( sin(2.0*alpha - xi(i)) - 3.0*sin(xi(i)) )/(2.0*sqrt(12));
		    q(i,2) = -0.5*cos(alpha - xi(i));*/
		    q(i,0) = ( sin(alpha)*cos(alpha)*sin(xi(i)) - \
    				(1.+ sin(alpha)*sin(alpha))*cos(xi(i))  )/(2.0*sqrt(3.));
    		 q(i,1) = ( sin(alpha)*cos(alpha)*cos(xi(i)) - \
    				(1.+ cos(alpha)*cos(alpha))*sin(xi(i))  )/(2.0*sqrt(3.));
    		 q(i,2) = -0.5*cos(alpha - xi(i)); 
		}
		for(int i=0; i<3; i++){
			n(0,i) = q(1,i) - q(2,i);
		    n(1,i) = q(2,i) - q(0,i);
		    n(2,i) = q(0,i) - q(1,i);
		}
		
	}

	void LISAmotion::EccentricLISAMotion2(float kappa0, float lambda0, double t, Matrix<double>& R, \
	                         Matrix<double>& p, Matrix<double>& n){
		//double Omega = LISAWP_TWOPI/year;
	    double Omega = LISAWP_TWOPI*3.1709791983764586e-8;
		double alpha = Omega*t + kappa0;
		Matrix<double> xi(1,3); // this is beta
		xi(0) = lambda0;             //////////  MOVE THOSE CALCULATIONS OUT TO THE CONSTRUCTOR
		xi(1) = lambda0 + 2.0*LISAWP_PI/3.0;
		xi(2) = lambda0 + 4.0*LISAWP_PI/3.0;
		double RAU = LISAWP_AU_SI/LISAWP_C_SI;
		
		R(0) = RAU*cos(alpha);
		R(1) = RAU*sin(alpha);
		R(2) = 0;
		
		double ecc = 0.00964838;
		
		for(int i=0; i<3; i++){
			p(i,0) = R(0) + RAU*ecc*( sin(alpha)*cos(alpha)*sin(xi(i)) - \
				(1.+ sin(alpha)*sin(alpha))*cos(xi(i))  );
			p(i,1) = R(1) + RAU*ecc*( sin(alpha)*cos(alpha)*cos(xi(i)) - \
				(1.+ cos(alpha)*cos(alpha))*sin(xi(i))  );
		    p(i,2) = -sqrt(3.)*RAU*ecc*cos(alpha - xi(i));
		}
		for(int i=0; i<3; i++){
			n(0,i) = (p(1,i) - p(2,i))/L;
		    	n(1,i) = (p(2,i) - p(0,i))/L;
		        n(2,i) = (p(0,i) - p(1,i))/L;
		}
		
	}

	
	void LISAmotion::CircularLISAMotion(float kappa0, float lambda0, double t, Matrix<double>& R, \
						                         Matrix<double>& q, Matrix<double>& n)
	{
			double Omega = LISAWP_TWOPI/year;
			double eta = Omega*t + kappa0;
			Matrix<double> sigma(1,3);
			double RAU = LISAWP_AU_SI/LISAWP_C_SI;

			R(0) = RAU*cos(eta);
			R(1) = RAU*sin(eta);
			R(2) = 0;				
			
			double zeta = -LISAWP_PI/6.;
			double xi = -Omega*t + lambda0;
			Matrix<double> LBrot(3,3);
			LBrot(0,0) = sin(eta)*cos(xi) - cos(eta)*sin(zeta)*sin(xi);
			LBrot(0,1) = -sin(eta)*sin(xi) - cos(eta)*sin(zeta)*cos(xi);
			LBrot(0,2) = -cos(eta)*cos(zeta);
			LBrot(1,0) = -cos(eta)*cos(xi) - sin(eta)*sin(zeta)*sin(xi);
			LBrot(1,1) = cos(eta)*sin(xi) - sin(eta)*sin(zeta)*cos(xi);
			LBrot(1,2) =  -sin(eta)*cos(zeta);
			LBrot(2,0) = cos(zeta)*sin(xi);
			LBrot(2,1) = cos(zeta)*cos(xi);
			LBrot(2,2) = -sin(zeta);
			
			sigma(0) = 1.5*LISAWP_PI;
			sigma(1) = 5.*LISAWP_PI/6.;
			sigma(2) = LISAWP_PI/6.;
			
			//In the LISA's frame:
			
			Matrix<double> qL(3,3);
			Matrix<double> nL(3,3);
			for(int i=0; i<3; i++){
				nL(i,0) = cos(sigma(i));
				nL(i,1) = sin(sigma(i));
				nL(i,2) = 0.0; 
				qL(i,0) = -cos(2.*sigma(i))/sqrt(3.);
				qL(i,1) = sin(2.*sigma(i))/sqrt(3.);
				qL(i,2) = 0.0;
			}
		/*	std::cout <<  "nL in lISA frame: " << std::endl;
			std::cout << nL;
			
			std::cout << "rotation matrix : " << std::endl;
			std::cout << LBrot;*/
			
			// In SSB frame
			n=0.0;
			q=0.0;
			for(int i=0; i<3; i++){
				for(int ii=0; ii<3; ii++){
					for(int j=0; j<3; j++){
						n(i,ii) += LBrot(ii,j)*nL(i,j);
						q(i,ii) += LBrot(ii,j)*qL(i,j);
					}	
				}
			}
			
			
/*		// Transformation to the LISA's frame
		   	double theta = 0.6981317;
			double phi = 0.52359878;
			Matrix<double> k(3);
			Matrix<double> uhat(3);
			Matrix<double> vhat(3);
			Matrix<double> kL(3);
			Matrix<double> uhatL(3);
			Matrix<double> vhatL(3);
			
			k(0) = -sin(theta)*cos(phi); 
			k(1) = -sin(theta)*sin(phi);
			k(2) = -cos(theta);

			uhat(0) = cos(theta)*cos(phi);
			uhat(1) = cos(theta)*sin(phi);
			uhat(2) = -sin(theta);

			vhat(0) = sin(phi);
			vhat(1) = -cos(phi);
			vhat(2) = 0.0;
			kL=0.0;
			uhatL=0.0;
			vhatL = 0.0;
			Matrix<double > BLrot = LBrot.Transpose();
			for(int i=0; i<3; i++){
				for(int j=0; j<3; j++){
					kL(i) += BLrot(i,j)*k(j);
					uhatL(i) += BLrot(i,j)*uhat(j);
					vhatL(i) += BLrot(i,j)*vhat(j);
				}	
			}
			
			std::cout <<"kL -> " << std::endl;
			std::cout <<  kL;
			std::cout <<"uhatL -> " << std::endl;
			std::cout <<  uhatL;
			std::cout <<"vhatL -> " << std::endl;
			std::cout <<  vhatL;
			
*/			
									
							
	}
	
}
