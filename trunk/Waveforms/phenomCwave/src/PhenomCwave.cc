/*
Copyright (C) 2011  E. Robinson, S. Babak, F. Ohme

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



#include "PhenomCwave.hh"

namespace LISAWP{

PhenomCwave::PhenomCwave(double Fmin, double Fmax, double delF){
   
    fMin = Fmin;
    fMax = Fmax;
    df = delF;
   
}

int PhenomCwave::ComputeHpHc(BBHTemplate H, std::complex<double>* &Hp, std::complex<double>* &Hc){
   
   double eta = H.eta;
   double chi = H.chi;
   double Mtot = H.M;
   
   int n,i;
   n = (int) floor (fMax / df)+1;
   if (Hp == NULL){
      Hp = new std::complex<double>[n];
   }
   if (Hc == NULL){
      Hc = new std::complex<double>[n];
   }
   
  
   
   double* freq;
   freq = new double[n];
   double* amp;
   double* phs;
   amp  = new double[n];
   phs = new double[n];
   for (int i=0; i<n; i++ )
    {
      freq[i]=i*df;
      amp[i] = 0.0;
      phs[i] = 0.0;
    }
   
   phenomwf2( amp , phs , freq , n , df , fMin , fMax , eta , chi , Mtot, H.dist);

   //constructing + and x

   std::complex<double> y;
   double fact = sqrt(5./(64.*LISAWP_PI));
   std::complex<double> img(0.0, 1.0);
   double ci = cos(H.iota);
   //std::cout << "ci = " << ci << std::endl;
   //double cp = cos(2.*H.phi0);
   //double sp = sin(2.*H.phi0);
   double phase;
   double amplhp, amplhc;
   double cph, sph;
   
   for (int i=0; i<n; i++){
      phase = phs[i] + 2.*H.phi0;
      amplhp = amp[i]*fact*(1.+ci*ci); // This is |H_22Y_22| + |H_2-2 Y_2-2|, H_2-2 = H*_22  
      amplhc = -amp[i]*fact*2.*ci;         // This is -|H_22Y_22| + |H_2-2 Y_2-2|, H_2-2 = H*_22  
      cph = cos(phase);
      sph = sin(phase);
      Hp[i] = amplhp *(cph - img*sph);
      Hc[i] = img*amplhc*(cph - img*sph); // check -
      
   }
   
    

   delete [] freq;
   delete [] amp;
   delete [] phs;
   
   return(n);
   
}

int PhenomCwave::ComputeH22(BBHTemplate H, std::complex<double>* &H22){
   
   double eta = H.eta;
   double chi = H.chi;
   double Mtot = H.M;
   
  // std::cout << "eta = " << eta << "  chi = " << chi << "  Mtot = " << Mtot << std::endl;
   
   int n,i;
   n = (int) floor (fMax / df)+1;
   if (H22 == NULL){
      H22 = new std::complex<double>[n];
   }
   
   double* freq;
   freq = new double[n];
   double* hreal;
   double* himg;
   hreal  = new double[n];
   himg = new double[n];
   std::complex<double> img(0.0, 1.0);
   
   for (int i=0; i<n; i++ )
    {
      freq[i]=i*df;
      hreal[i] = 0.0;
      himg[i] = 0.0;
    }
   
   phenomwf( hreal , himg , freq , n , df , fMin , fMax , eta , chi , Mtot );
   
   for (int i=0; i<n; i++){
      H22[i] = hreal[i] - img*himg[i];
   }  
   
   delete [] freq;
   delete [] hreal;
   delete [] himg;
   
   return(n);
   
}

int PhenomCwave::ComputeHofT(BBHTemplate H, double* &HpT, double* &HcT, double tShift){
   
   /*** Here I assume that df, fmin was set up correctly
       it should be df = 1/T_{obs}, fmin is set by coalescence time
       last but not least tshift is given in sec and the time domain waveform 
       will be shifted towards the negative time. Make sure that t_c < T_{obs} -tshift ***/
   std::complex<double>* Hplus = NULL;
   std::complex<double>* Hcross = NULL;

   int n = ComputeHpHc(H, Hplus, Hcross);    
   //std::cout << "freq domain waveforms computed " << n << std::endl;
   
   int datSz = 2*(n-1);
   HpT = new double[datSz];
   HcT = new double[datSz];   
   
   // apply time shift
   std::complex<double> img(0.0, 1.0);
   double shiftPh, fr;
   for (int i=0; i<n; i++){
      fr=i*df;
      shiftPh = tShift*LISAWP_TWOPI*fr;
      Hplus[i] = Hplus[i] * (cos(shiftPh) + img*sin(shiftPh));
      Hcross[i] = Hcross[i] * (cos(shiftPh) + img*sin(shiftPh));
   }
   
   crfft1d Backward(datSz, Hplus, HpT);
   Backward.fft(Hplus, HpT);
   Backward.fft(Hcross, HcT);
   
   //std::cout << " ffts are computed\n "; 
   
   delete [] Hplus;
   delete [] Hcross;
   
   return(datSz);
       
   
}

void PhenomCwave::ComputeTime(BBHTemplate H, double* &freq, double* &tm, int n){
   
   double eta = H.eta;
   double chi = H.chi;
   double Mtot = H.Mt;
   double tf, f, om;
   double chi2 = chi*chi;
   double eta2 = eta*eta;
   double Pi = LISAWP_PI;
   
/*   if (tm = NULL){
      std::cout << "Stas\n";
      tm = new double[n];
   }*/
   
   double tprev = 0.0;
   //std::cout << fMin << " freqs   " << fMax << std::endl;
   f = fMin*Mtot;
   om = f*Pi;
   tf = -5./(256.*eta*pow(om,(8./3.))) - (5.*(743. + 924.*eta))/(64512.*eta*pow(om,2.)) \
   - (5.*(3058673. + 5472432.*eta + 4353552.*eta2 + 63504.*chi2*(-81. + 4.*eta)))/(130056192.*eta*pow(om,(4./3.))) +\ 
  (chi*(-113. + 76.*eta) + 48.*Pi)/(384.*eta*pow(om,(5./3.))) - (5.*(-1512.*chi*chi2*(-1. + 3.*eta) + \
  chi*(147101. - 137368.*eta - 17136.*eta2) + 3.*(-7729. + 1092.*eta)*Pi))/(193536.*eta*om) + \
  (5.*(18144.*chi2*chi*(19422. - 14929.*eta + 10136.*eta2) + chi*(-4074790483. + 3478848284.*eta + 2255806224.*eta2 - \ 
   404760384.*eta2*eta) + 6096384.*chi2*(-57. + 4.*eta)*Pi + \
   12.*(15419335. + 12718104.*eta - 4975824.*eta2)*Pi))/(390168576.*eta*pow(om,(1./3.))) + \
   (Pi*Pi*(10052469856691. + 206607970800.*eta2 - 462992376000.*eta*eta2 - 34927200.*chi2*(6845. - 173708.*eta + 54880.*eta2) -\ 
   1530761379840.*LISAWP_GAMMA - 62589542400.*chi*(-73. + 56.*eta)*Pi - 1001432678400.*Pi*Pi + 7700.*eta*(-3147553127. +\
   114561216.*Pi*Pi) - 1530761379840.*log(4.)) - 765380689920.*Pi*Pi*log(pow(om,(2./3.))))/(1201719214080.*eta*Pi*Pi*pow(om,(2./3.)));
   
   double tc = tf;
  
   tprev = tc-10.;
   for (int i=0; i<n; i++){
      if (freq[i] >= fMin && freq[i] <= fMax){
         //std::cout << freq[i] << std::endl;
         f = freq[i]*Mtot;
         om = f*Pi;      
         tf = -5./(256.*eta*pow(om,(8./3.))) - (5.*(743. + 924.*eta))/(64512.*eta*pow(om,2.)) \
         - (5.*(3058673. + 5472432.*eta + 4353552.*eta2 + 63504.*chi2*(-81. + 4.*eta)))/(130056192.*eta*pow(om,(4./3.))) +\ 
        (chi*(-113. + 76.*eta) + 48.*Pi)/(384.*eta*pow(om,(5./3.))) - (5.*(-1512.*chi*chi2*(-1. + 3.*eta) + \
        chi*(147101. - 137368.*eta - 17136.*eta2) + 3.*(-7729. + 1092.*eta)*Pi))/(193536.*eta*om) + \
        (5.*(18144.*chi2*chi*(19422. - 14929.*eta + 10136.*eta2) + chi*(-4074790483. + 3478848284.*eta + 2255806224.*eta2 - \ 
         404760384.*eta2*eta) + 6096384.*chi2*(-57. + 4.*eta)*Pi + \
         12.*(15419335. + 12718104.*eta - 4975824.*eta2)*Pi))/(390168576.*eta*pow(om,(1./3.))) + \
         (Pi*Pi*(10052469856691. + 206607970800.*eta2 - 462992376000.*eta*eta2 - 34927200.*chi2*(6845. - 173708.*eta + 54880.*eta2) -\ 
         1530761379840.*LISAWP_GAMMA - 62589542400.*chi*(-73. + 56.*eta)*Pi - 1001432678400.*Pi*Pi + 7700.*eta*(-3147553127. +\
         114561216.*Pi*Pi) - 1530761379840.*log(4.)) - 765380689920.*Pi*Pi*log(pow(om,(2./3.))))/(1201719214080.*eta*Pi*Pi*pow(om,(2./3.)));
         //std::cout << tc - tf << "   " << (tc-tf)*Mtot << std::endl;
         //tm[i] = (tf-tc)*Mtot;
         //tm[i] = tf*Mtot;
         
         if (tf >= tprev){
             tm[i] = (tf - tc)*Mtot;
             tprev = tf;
         }else
            tm[i] = (tprev-tc)*Mtot; 
         
      }
   }
   
}

}//end of the namespace