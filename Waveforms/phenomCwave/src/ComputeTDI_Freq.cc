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



#include "ComputeTDI_Freq.hh"

namespace LISAWP{
   
ComputeTDIfreq::ComputeTDIfreq(double arm, double oneyear){
      
      uhat = new double[3];
      vhat = new double[3];
      k = new double[3];
      kp = new double[3];
      kn = new double[3];
      u = new double[3];
      v = new double[3];
      R = new double[3];
      p = new double*[3];
      n = new double*[3];

      for (int i=0; i<3; i++){
          p[i] = new double[3];
          n[i] = new double[3];
      }
         
      L = arm;
      year = oneyear;
      fMin = 1.e-5;
      fMax = 1.0;
      
}

ComputeTDIfreq::~ComputeTDIfreq()
{
   delete [] uhat;
   delete [] vhat;
   delete [] k;
   delete [] kp;
   delete [] kn;
   delete [] u;
   delete [] v;
   delete [] R;
   
   for (int i=0; i<3; i++){
      delete [] p[i];
      delete [] n[i];
   }
   delete p;
   delete n;
}

void ComputeTDIfreq::ChooseConfiguration(std::string config)
{
   conf = config;
   
}

void ComputeTDIfreq::ComputeTDIfreqXYZ(int sz, double thetaS, double phiS, double* &tf, double* &freq, \
                        std::complex<double>* &hfp, std::complex<double>* &hfc, \
                        std::complex<double>* &X, std::complex<double>* &Y, std::complex<double>* &Z)
{

    double cS = cos(thetaS);
    double sS = sin(thetaS);
       
    double cSp = cos(phiS);
    double sSp = sin(phiS);

   // note that k = -n, where k is propagation vector and n is sky location
    k[0] = -sS*cSp; 
    k[1] = -sS*sSp;
    k[2] = -cS;

    uhat[0] = cS*cSp;
    uhat[1] = cS*sSp;
    uhat[2] = -sS;

    vhat[0] = sSp;
    vhat[1] = -cSp; 
    vhat[2] = 0.0;
    
    OrbitalMotion orbit(L, year);
    
    /*double R[3];
    double p[3][3];
    double n[3][3];*/
    
   
    
   /* if(X == NULL){
        X = new std::complex<double>[sz];
    }
    if(Y == NULL){
        Y = new std::complex<double>[sz];
    }
    if(Z == NULL){
        Z = new std::complex<double>[sz];
    }*/
    
    double tm, fr;
    double nU, nV;
    std::complex<double> fct;
    std::complex<double> Dplr;
    std::complex<double> Al1, Al2, Al3;
    std::complex<double> Al_1, Al_2, Al_3;
    std::complex<double> img(0.0, 1.0);
    std::complex<double> y123, y1_32, y231, y3_21, y2_13, y312;
    std::complex<double> ex1, ex2;
    double om2L, omL;
    double kdR;
    
   // std::ofstream fout2718("Data/CheckFpc.dat");
    for (int i=0; i<sz; i++){
       tm = tf[i];
       //std::cout << std::setprecision(15) << i << "   "  << tm << std::endl;
       fr = LISAWP_PI*freq[i]; 
       
       if (freq[i] >= fMin && freq[i] <= fMax){
          if(conf == "aLISA"){
             orbit.EccentricLISAMotion(0.0, 0.0, tm, R, p, n);
            // std::cout << tm << "   ";
            // for (int i=0; i<3; i++) std::cout << R[i] << "   ";
            // std::cout << "\n";
          }else{
             std::cerr << "Error: Unknown orbit " << conf << std::endl;
             exit(1);
          }
          kdR = 0.0; 
          for(int j =0; j<3; j++){
             kp[j] = 0.; 
             kn[j] = 0.;       
 		       nU = 0.0;
 		       nV = 0.0;        
 		       for(int ii=0; ii<3; ii++){
 		           kp[j] += k[ii]*p[j][ii];
                 kn[j] += k[ii]*n[j][ii];
 			        nU += uhat[ii]*n[j][ii];
 			        nV += vhat[ii]*n[j][ii];
		       }
		       u[j] = 0.5*(nU*nU - nV*nV);
 		       v[j] = nU*nV;
             kdR += k[j]*R[j];
          }
         /*fout2718 << std::setprecision(15) << tm << "  ";
         for (int bk=0; bk<3; bk++){
           for (int bk2=0; bk2<3; bk2++){
             fout2718 << p[bk][bk2] << "   ";
           }
         }
         fout2718 << std::endl;
       fout2718 << std::setprecision(15) << tm << "   " << u[0] << "  " << u[1] << "  " << u[2] << \
                     "   " << v[0] << "  " << v[1] << "  " << v[2] << "   "\
                     <<  kn[0] << "  " << kn[1] << "  "<< kn[2] << std::endl; */
       
          fct = -img*fr*L*(u[0]*hfp[i] + v[0]*hfc[i]);
          Al1 = fct*SINC(fr*L*(1.-kn[0]));
          Al_1 = fct*SINC(fr*L*(1.+kn[0]));
       
          fct = -img*fr*L*(u[1]*hfp[i] + v[1]*hfc[i]);
          Al2 = fct*SINC(fr*L*(1.-kn[1]));
          Al_2 = fct*SINC(fr*L*(1.+kn[1]));
       
          fct = -img*fr*L*(u[2]*hfp[i] + v[2]*hfc[i]);
          Al3 = fct*SINC(fr*L*(1.-kn[2]));
          Al_3 = fct*SINC(fr*L*(1.+kn[2]));
       
          Dplr = cos(fr*(L + kp[0] + kp[2])) - img*sin(fr*(L + kp[0] + kp[2]));
          y123 = Al2*Dplr;
          y3_21 = Al_2*Dplr;
       
          Dplr = cos(fr*(L + kp[0] + kp[1])) - img*sin(fr*(L + kp[0] + kp[1]));
          y231 = Al3*Dplr;
          y1_32 = Al_3*Dplr;
         
          Dplr = cos(fr*(L + kp[1] + kp[2])) - img*sin(fr*(L + kp[1] + kp[2]));
          y312 = Al1*Dplr;
          y2_13 = Al_1*Dplr;
       
          om2L = 4.*fr*L;
          omL = 2.*fr*L;
          ex2 = cos(om2L) - img*sin(om2L);
          ex1 = cos(omL) - img*sin(omL);
          X[i] = -4.*img*sin(omL)*(-y1_32*ex2 - y231*ex1 + y123*ex2 + y3_21*ex1);
          //X[i] = -16.*pow(fr*L, 2.)*( (u[1]-u[2])*hfp[i] + (v[1]-v[2])*hfc[i] )*(cos(fr*(-3.*L + 2.*kdR)) - img*sin(fr*(2.*kdR - 3.*L)));
          //X[i] = -16.*pow(fr*L, 2.)*( (u[1]-u[2])*hfp[i] + (v[1]-v[2])*hfc[i] )*(cos(fr*(2.*kdR)) - img*sin(fr*(2.*kdR)));
          //X[i] = 4.*img*omL*2.*img*fr*L*( (u[2]*hfp[i] + v[2]*hfc[i]) - (u[1]*hfp[i] + v[1]*hfc[i]) )*(cos(2.*fr*kdR) - img*sin(2.*fr*kdR));
       
         ///X[i] = 4.*fr*L *sin(omL)*( -(u[1]*hfp[i] + v[1]*hfc[i])*(SINC(fr*L*(1.-kn[1]))*ex1 + SINC(fr*L*(1.+kn[1])) )*\
                  (cos(fr*(3.*L + kp[0] + kp[2])) - img*sin(fr*(3.*L + kp[0] + kp[2])) ) +\
                  (u[2]*hfp[i] + v[2]*hfc[i])*(SINC(fr*L*(1.+kn[2]))*ex1 + SINC(fr*L*(1.-kn[2])) )*\
                  (cos(fr*(3.*L + kp[0] + kp[1])) - img*sin(fr*(3.*L + kp[0] + kp[1])) )  );
          Y[i] = -4.*img*sin(omL)*(-y2_13*ex2 - y312*ex1 + y231*ex2 + y1_32*ex1);
          //Y[i] = -16.*pow(fr*L, 2.)*( (u[2]-u[0])*hfp[i] + (v[2]-v[0])*hfc[i] )*(cos(2.*fr*kdR) - img*sin(2.*fr*kdR));
          Z[i] = -4.*img*sin(omL)*(-y3_21*ex2 - y123*ex1 + y312*ex2 + y2_13*ex1);
       }else{
          X[i] = 0.0;
          Y[i] = 0.0;
          Z[i] = 0.0;
       }
       
    }
    
    //fout2718.close();
   
}


void ComputeTDIfreq::ComputeLWfreqXYZ(int sz, double thetaS, double phiS, double* &tf, double* &freq, \
                        std::complex<double>* &hfp, std::complex<double>* &hfc, \
                        std::complex<double>* &X, std::complex<double>* &Y, std::complex<double>* &Z)
{
 
    double tm, fr;
    double nU, nV;
    double kdR;
    std::complex<double> img(0.0, 1.0);
    
    double cS = cos(thetaS);
    double sS = sin(thetaS);
       
    double cSp = cos(phiS);
    double sSp = sin(phiS);

   // note that k = -n, where k is propagation vector and n is sky location
    k[0] = -sS*cSp; 
    k[1] = -sS*sSp;
    k[2] = -cS;

    uhat[0] = cS*cSp;
    uhat[1] = cS*sSp;
    uhat[2] = -sS;

    vhat[0] = sSp;
    vhat[1] = -cSp; 
    vhat[2] = 0.0;
    
    std::complex<double> Dplr;
    
    
    OrbitalMotion orbit(L, year);
    for (int i=0; i<sz; i++){
        tm = tf[i];
        fr = LISAWP_PI*freq[i]; 
        if (freq[i] >= fMin && freq[i] <= fMax){
           if(conf == "aLISA"){
              orbit.EccentricLISAMotion(0.0, 0.0, tm, R, p, n);
           }else{
              std::cerr << "Error: Unknown orbit " << conf << std::endl;
              exit(1);
           }
           kdR = 0.0; 
           for(int j =0; j<3; j++){
               kp[j] = 0.; 
               kn[j] = 0.;       
  		         nU = 0.0;
  		         nV = 0.0;        
  		         for(int ii=0; ii<3; ii++){
  		            kp[j] += k[ii]*p[j][ii];
                  kn[j] += k[ii]*n[j][ii];
  			         nU += uhat[ii]*n[j][ii];
  			         nV += vhat[ii]*n[j][ii];
 		         }
 		         u[j] = 0.5*(nU*nU - nV*nV);
  		         v[j] = nU*nV;
               kdR += k[j]*R[j];
           }
           Dplr = (cos(2.*fr*kdR) - img*sin(2.*fr*kdR));
           X[i] = -16.*pow(fr*L, 2.)*( (u[1]-u[2])*hfp[i] + (v[1]-v[2])*hfc[i] )*Dplr;
           Y[i] = -16.*pow(fr*L, 2.)*( (u[2]-u[0])*hfp[i] + (v[2]-v[0])*hfc[i] )*Dplr;
           Z[i] = -16.*pow(fr*L, 2.)*( (u[0]-u[1])*hfp[i] + (v[0]-v[1])*hfc[i] )*Dplr;
        }
        else{
           X[i] = 0.0;
           Y[i] = 0.0;
           Z[i] = 0.0;
        }
    }
   
}

double ComputeTDIfreq::SINC(double x){
   if (x == 0.){
      return(1.);
   }else{
      return(sin(x)/x);
   }
   
}
   
}//end of the namespace