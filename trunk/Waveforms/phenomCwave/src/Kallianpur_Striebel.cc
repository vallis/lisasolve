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


#include "Kallianpur_Striebel.h"

namespace LISAWP{

KS_reconstruct::KS_reconstruct(std::complex<double>* &input, double* freqs, double delf, \
		     double FMin, double FMax, int length){

   sz = length;
   data = new std::complex<double>[sz];
   freq = new double[sz];
   Snf = new double[sz];
   tmplt = new std::complex<double>[sz];
   for (int i=0; i< sz; i++){
       data[i] = input[i];
       freq[i] = freqs[i];
       Snf[i] = 1.0;
   }

   fMin = FMin;
   fMax = FMax;
   df = delf;
   LISAWPAssert(freq[1]-freq[0] == df, "given df doesn't correspind to freq step in array");  

}

KS_reconstruct::~KS_reconstruct(){

  delete [] data;
  delete [] freq;
  delete [] tmplt;
  delete [] Snf;
  if (bndr == NULL)
     delete [] bndr;
}

void KS_reconstruct::SetPsd(std::string noiseType, double scale){

   if (noiseType == "white"){
      for (int i=0; i<sz; i++){
	  Snf[i] *= scale;
      }
   }else{
      std::cout << "This noise model is not implemented" << std::endl;	   
   }

}


void KS_reconstruct::SetBaseWaveform(std::string Model){

   if (Model != "phenomC"){
      std::cout << "This Model is not implemented" << std::endl;
   }

}

void KS_reconstruct::SetIntegrationMethod(std::string IntType){

    typeInt = IntType;
    std::cout << typeInt << std::endl;

    if (typeInt != "MC" && typeInt != "miser" && typeInt != "vegas"){
	 std::cout << "Unknown integration method" << std::endl;
	 exit(1);
    }
}

void KS_reconstruct::SetBoundaries(double* &bndrs, int N){

    maxPhi0= true;
    maxTc = true;
    maxAmp = true;    
    bndr = new double[2*N];

    std::cout << "N = " << N << std::endl;
    if (N > 3) maxPhi0 = false;
    if (N > 4) maxTc = false;
    if (N > 5) maxAmp = false;

    for (int i=0; i<2*N; i++){
      bndr[i] = bndrs[i];
    }
    
}

double KS_reconstruct::ComputeNormInt(){

   double* theta;
   theta = new double[6];
   theta[0] = 1.e6;
   theta[1] = 0.25;
   theta[2] = 0.5;
   theta[4] = 26457280.;
   theta[3] = 0.2;
   theta[5] = 1.0;

   void* params;
   double L = ComputeL(theta, 3, params);

   std::cout << " likelihood = " << L << std::endl;  



   //gsl_monte_function G = { &gg, 3, 0 };

   return(L);
}

double KS_reconstruct::ComputeL(double *theta, size_t dim, void* params){

   double M = theta[0];
   double q = theta[1];
   double chi = theta[2];
   double phi0 =0.0;
   double tc = 0.0;
   double  Amp = 0.0;

   if (dim > 3) phi0 = theta[3];
   if (dim > 4) tc = theta[4];
   if (dim > 5) Amp = theta[5];

   BBHTemplate S;
   
   S.M = M;
   S.q = q;
   S.chi = chi;
   S.tc = tc;
   S.phi0 = phi0;
   S.dist = 1.e-4;

   //std::cout << S.M << "  " << S.q << "  " << S.chi << "  " << S.tc << "   " <<\
	   S.phi0  << "   " << S.dist << std::endl;

   S.ComputeAuxParams();
   
   PhenomCwave phnC(fMin, fMax, df);
   phnC.ComputeH22(S, sz, freq, tmplt);

   /*std::ofstream fout1543("Data/TestKSwave.dat");
   for (int i=0; i<sz; i++ )
   {
      fout1543 << std::setprecision(15) << freq[i] << "  " << abs(tmplt[i]) << \
	      "   " << tmplt[i].real() << "  " << tmplt[i].imag() << std::endl;  
   } 
   fout1543.close();*/


   double ss = ComputeInnerProduct(sz, tmplt, tmplt);
   double sy = ComputeInnerProduct(sz, tmplt, data);
   double yy = ComputeInnerProduct(sz, data, data);
   std::cout << "ss = " << ss << "  sy = " << sy << "  yy = " << yy << std::endl;


   int lag;
   lag = 0;
   //sy =  ComputeInnerProductV(sz, tmplt, data, lag);
   //std::cout << "maximized over tc: sy = " << sy << "  lag = " << lag << std::endl; 

   
   //sy = MaxTcPhi(sz, tmplt, data, lag);
   //std::cout << "maximized over phi and tc: sy = " << sy << "  lag = " << lag << std::endl; 


   double logL = sy - 0.5*ss;
   std::cout << "logL = " << logL << std::endl;
   double dt = 1./(df*(double)sz);
   double tshift = 0.0;
   double shiftPh;
   std::complex<double> img(0.0, 1.0); 	
   if (maxAmp){
       std::cout << "maxAmp \n"; 
       logL = 0.5*sy*sy/ss;
       Amp = sy;
       std::cout << "logL = " << logL << std::endl;
   }
   if (maxTc){
      if (maxPhi0){
        std::cout << "---- maxTc and maxPhi0 ---\n"; 
        sy = MaxTcPhi(sz, tmplt, data, lag, phi0);
        std::cout << "maximized over phi and tc: sy = " << sy << "  lag = " << lag << std::endl; 
        std::cout << "shift = " << (double)lag*dt << std::endl;
        std::cout << " recovered phi0 = " << 0.5*phi0 << std::endl; //it is orbital phase -> 1/2
        logL = 0.5*sy/ss;
        std::cout << "logL = " << logL << std::endl;
        tshift = (double)lag*dt;
	//std::cout << "check -> " << fmod(tshift, LISAWP_PI) << std::endl;
        //std::cout << " recovered phi0 = " << 0.5*LISAWP_PI << "   " \
		<< phi0 -fmod(tshift, LISAWP_PI) <<\
	         "   " << phi0 +fmod(tshift, LISAWP_PI) << "\n"	<< std::endl; 

        std::cout << "------ second iteration --------\n";
        S.phi0 = 0.0;
        phnC.ComputeH22(S, sz, freq, tmplt);
        tshift = (double)lag*dt;
        for (int i=0; i<sz; i++ )
        {
           shiftPh = tshift*LISAWP_TWOPI*freq[i];
           //shiftPh = 0.0;
           tmplt[i] = tmplt[i]*(cos(shiftPh) + img*sin(shiftPh));
        }
        sy = MaxTcPhi(sz, tmplt, data, lag, phi0);
        std::cout << "maximized over phi and tc: sy = " << sy << "  lag = " << lag 
	    <<  "    out of   " <<  sz << std::endl; 
        std::cout << " recovered phi0 = " << 0.5*phi0 << std::endl; //it is orbital phase -> 1/2
        std::cout << "shift = " << (double)lag*dt << std::endl;
        logL = 0.5*sy/ss;
        std::cout << "logL = " << logL << std::endl;

 
      }else{
        std::cout << "---- maxTc  --- \n"; 
        sy =  ComputeInnerProductV(sz, tmplt, data, lag);
        std::cout << " maximized over tc: sy = " << sy << "  lag = " << lag << std::endl; 
        logL = 0.5*sy/ss;
        std::cout << "logL = " << logL << std::endl;

      }

   }

   /*std::cout << " recovered phi0 = " << 0.5*phi0 << std::endl; //it is orbital phase -> 1/2
   std::cout << "shift = " << (double)lag*dt << std::endl;
   S.phi0 = 0.0;
   phnC.ComputeH22(S, sz, freq, tmplt);
   
   tshift = (double)lag*dt;
   for (int i=0; i<sz; i++ )
   {
      shiftPh = tshift*LISAWP_TWOPI*freq[i];
      //shiftPh = 0.0;
      //std::cout << shiftPh/LISAWP_TWOPI << spr << cos(shiftPh) << spr << sin(shiftPh) << std::endl; 
      tmplt[i] = tmplt[i]*(cos(shiftPh) + img*sin(shiftPh));
  }
  std::cout << "maxTc and maxPhi0 \n"; 
  sy = MaxTcPhi(sz, tmplt, data, lag, phi0); 

  S.phi0 = phi0;
  phnC.ComputeH22(S, sz, freq, tmplt);
   
  tshift += (double)lag*dt;
   for (int i=0; i<sz; i++ )
   {
      shiftPh = tshift*LISAWP_TWOPI*freq[i];
      //shiftPh = 0.0;
      //std::cout << shiftPh/LISAWP_TWOPI << spr << cos(shiftPh) << spr << sin(shiftPh) << std::endl; 
      tmplt[i] = tmplt[i]*(cos(shiftPh) + img*sin(shiftPh));
  }
  std::cout << "maxTc and maxPhi0 \n"; 
  sy = MaxTcPhi(sz, tmplt, data, lag, phi0);
  std::cout << "------ third iteration --------\n";
  std::cout << "maximized over phi and tc: sy = " << sy << "  lag = " << lag 
	    <<  "    out of   " <<  sz << std::endl; 
  std::cout << " recovered phi0 = " << 0.5*phi0 << std::endl; //it is orbital phase -> 1/2
  std::cout << "shift = " << (double)lag*dt << std::endl;
  logL = 0.5*sy/ss;
  std::cout << "logL = " << logL << std::endl; */
 


   return(exp(logL));

}
double KS_reconstruct::ComputeInnerProduct(int n, std::complex<double>* &x,\
	       	std::complex<double>* &y){

   double prod2 = 0.0;
   // std::ofstream fout615("Data/InPrTest.dat");
   for (int i=1; i<sz; i++){
      if (freq[i] >= fMin && freq[i] <= fMax){
         prod2 += df*(x[i] * conj(y[i])).real()/Snf[i];
         /*  if (print){
              fout615 << freq[i] << "    " << prod2 << "   " <<  (x[i] * conj(y[i])).real() \
                 << "    "<< abs(x[i]) << "    " << abs(y[i]) << std::endl;
           }*/
      }
   }
//   fout615.close();
   prod2 *= 4.;
   return(prod2);

}

double KS_reconstruct::ComputeInnerProductV(int n, std::complex<double>* &x,\
	       	std::complex<double>* &y, int &lag){


   std::complex<double>* prod;
   double* rho;
   prod = new std::complex<double>[sz];
   rho = new double[sz];
 
   for (int i=1; i<sz; i++){
      prod[i] = 0.0;	   
      if (freq[i] >= fMin && freq[i] <= fMax){
         prod[i]= 2.*df*(x[i] * conj(y[i]))/Snf[i];
      }
   }
   crfft1d Backward(sz, prod, rho);
   Backward.fft(prod, rho);

   // looking for maximum
   
   double rhoMax = rho[0]*rho[0];
   lag = 0;
   for (int i=1; i<sz; i++){
        if (rho[i]*rho[i] > rhoMax){
	   rhoMax = rho[i]*rho[i];
           lag = i;
	}
   }

  // double Tobs = 1./df;
  double dt = 1./(df*(double)sz);
  std::cout << "shift = " << (double)lag*dt << std::endl; 
  delete [] prod;
  delete [] rho;
  return(rhoMax);
	
}

double KS_reconstruct::MaxTcPhi(int n, std::complex<double>* &x, std::complex<double>* &y, \
		     int &lag, double& phi0){

   std::complex<double> img(0.0, 1.0); 	
   std::complex<double>* prod1;
   //std::complex<double>* rho1;
   double* rho1;
   std::complex<double>* prod2;
   //std::complex<double>* rho2;
   double* rho2;
   
   prod1 = new std::complex<double>[sz];
   //rho1 = new std::complex<double>[sz];
   rho1 = new double[sz];
   prod2 = new std::complex<double>[sz];
   //rho2 = new std::complex<double>[sz];
   rho2 = new double[sz];
 
   for (int i=1; i<sz; i++){
      prod1[i] = 0.0;	
      prod2[i] = 0.0;   
      if (freq[i] >= fMin && freq[i] <= fMax){
         prod1[i]= 2.*df*(x[i] * conj(y[i]))/Snf[i];
         prod2[i]= 2.*df*(img*x[i] * conj(y[i]))/Snf[i];
      }
   }
   crfft1d Backward(sz, prod1, rho1);
   Backward.fft(prod1, rho1);
   Backward.fft(prod2, rho2);

   //fft1d Backward(sz, -1, prod1, rho1);
   //Backward.fft(prod1, rho1);
   //Backward.fft(prod2, rho2);

   // looking for maximum
   
   //double rhoMax = norm(rho1[0]) + norm(rho2[0]);
   double rhoMax = pow(rho1[0], 2) + pow(rho2[0], 2);
   lag = 0;
   //phi0 = -atan2(rho1[0].imag(), rho1[0].real());
   phi0 = -atan2(rho2[0], rho1[0]);
   for (int i=1; i<sz; i++){
	/*if (prod1[i].real() != 0.)
	     std::cout << prod1[i] << std::endl;*/
        //if (norm(rho1[i]) + norm(rho2[i]) > rhoMax){
	//   rhoMax = norm(rho1[i]) + norm(rho2[i]);
        if (pow(rho1[i], 2) + pow(rho2[i], 2) > rhoMax){
	   rhoMax = pow(rho1[i], 2) + pow(rho2[i], 2);
           lag = i;
           //phi0 = -atan2(rho1[i].real(), rho1[i].imag());
           phi0 = -atan2(rho2[i], rho1[i]);
	}
   }
   std::cout << " inside maximization recovered phi0 = " << phi0 << std::endl;

  // double Tobs = 1./df;
  // double dt = Tobs/(double)sz;
  delete [] prod1;
  delete [] rho1;
  delete [] prod2;
  delete [] rho2;
  return(rhoMax);
 

}
//KS_reconstruct::
//KS_reconstruct::

/*double CompLhd(double *theta, size_t dim, void* params){
  
   double M = theta[0];
   double q = theta[1];
   double chi = theta[2];
   double phi0 =0.0;
   double tc = 0.0;
   double  Amp = 0.0;

   if (dim > 3) phi0 = theta[3];
   if (dim > 4) tc = theta[4];
   if (dim > 5) Amp = theta[5];

   BBHTemplate S;
   
   S.M = M;
   S.q = q;
   S.chi = chi;
   S.tc = tc;
   S.phi0 = phi0;
   S.dist = 1.e-4;

   S.ComputeAuxParams();
   
   PhenomCwave phnC(fMin, fMax, df);
   phnC.ComputeH22(S, sz, freq, tmplt);


   double ss = ComputeInnerProduct(sz, tmplt, tmplt);
   double sy = ComputeInnerProduct(sz, tmplt, data);
   double yy = ComputeInnerProduct(sz, data, data);
   std::cout << "ss = " << ss << "  sy = " << sy << "  yy = " << yy << std::endl;


   int lag;
   lag = 0;


   double logL = sy - 0.5*ss;
   std::cout << "logL = " << logL << std::endl;
   double dt = 1./(df*(double)sz);
   double tshift = 0.0;
   double shiftPh;
   std::complex<double> img(0.0, 1.0); 	
   if (maxAmp){
       std::cout << "maxAmp \n"; 
       logL = 0.5*sy*sy/ss;
       Amp = sy;
       std::cout << "logL = " << logL << std::endl;
   }
   if (maxTc){
      if (maxPhi0){
        std::cout << "---- maxTc and maxPhi0 ---\n"; 
        sy = MaxTcPhi(sz, tmplt, data, lag, phi0);
        std::cout << "maximized over phi and tc: sy = " << sy << "  lag = " << lag << std::endl; 
        std::cout << "shift = " << (double)lag*dt << std::endl;
        std::cout << " recovered phi0 = " << 0.5*phi0 << std::endl; //it is orbital phase -> 1/2
        logL = 0.5*sy/ss;
        std::cout << "logL = " << logL << std::endl;
        tshift = (double)lag*dt;
	//std::cout << "check -> " << fmod(tshift, LISAWP_PI) << std::endl;
        //std::cout << " recovered phi0 = " << 0.5*LISAWP_PI << "   " \
		<< phi0 -fmod(tshift, LISAWP_PI) <<\
	         "   " << phi0 +fmod(tshift, LISAWP_PI) << "\n"	<< std::endl; 

        std::cout << "------ second iteration --------\n";
        S.phi0 = 0.0;
        phnC.ComputeH22(S, sz, freq, tmplt);
        tshift = (double)lag*dt;
        for (int i=0; i<sz; i++ )
        {
           shiftPh = tshift*LISAWP_TWOPI*freq[i];
           //shiftPh = 0.0;
           tmplt[i] = tmplt[i]*(cos(shiftPh) + img*sin(shiftPh));
        }
        sy = MaxTcPhi(sz, tmplt, data, lag, phi0);
        std::cout << "maximized over phi and tc: sy = " << sy << "  lag = " << lag 
	    <<  "    out of   " <<  sz << std::endl; 
        std::cout << " recovered phi0 = " << 0.5*phi0 << std::endl; //it is orbital phase -> 1/2
        std::cout << "shift = " << (double)lag*dt << std::endl;
        logL = 0.5*sy/ss;
        std::cout << "logL = " << logL << std::endl;

 
      }else{
        std::cout << "---- maxTc  --- \n"; 
        sy =  ComputeInnerProductV(sz, tmplt, data, lag);
        std::cout << " maximized over tc: sy = " << sy << "  lag = " << lag << std::endl; 
        logL = 0.5*sy/ss;
        std::cout << "logL = " << logL << std::endl;

      }

   }

   

}*/

double gg (double *k, size_t dim, void *params)
{
   double A = 1.0 / (M_PI * M_PI * M_PI);
   return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
}


} // end of the namespace
