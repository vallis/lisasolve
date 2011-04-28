/*
Copyright (C) 2011*  S. Babak

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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/src/Template.cc,v 1.6 2008/04/18 11:14:32 stba Exp $"
*/

#include "BBHTemplate.hh"

namespace LISAWP{

BBHTemplate::BBHTemplate(){
   
   M = 0.0;
   Mt = 0.0;   // in seconds
   Mc = 0.0;   // in seconds
   eta = 0.0;
   q = 0.0;

   chi = 0.0;
   tc = 0.0;
   m1 = 0.0;
   m2 = 0.0;
   thetaS = 0.0;
   phiS = 0.0;
   dist = 0.0;
   Ampl = 0.0;
   iota = 0.0;
   phi0 = 0.0;
   psi = 0.0;

//     double Norm;
   f0 = 1.e-5;
   fEnd = 1.e-1;
   SNR = 0.0;
   
   
}

BBHTemplate::BBHTemplate(const BBHTemplate& p){

   M = p.M;
   Mt = p.Mt;   // in seconds
   Mc = p.Mc;   // in seconds
   eta = p.eta;
//   lam = ;
//   bet = ;
   chi = p.chi;
   tc = p.tc;
   m1 = p.m1;
   m2 = p.m2;
   thetaS = p.thetaS;
   phiS = p.phiS;
   dist = p.dist;
   Ampl = p.Ampl;
   iota = p.iota;
   phi0 = p.phi0;
   psi = p.psi;

//     double Norm;
   f0 = p.f0;
   fEnd = p.fEnd;
   SNR = p.SNR;
   
}

BBHTemplate BBHTemplate::operator=(const BBHTemplate& p){
   
      this->M = p.M;
      this->Mt = p.Mt;   // in seconds
      this->Mc = p.Mc;   // in seconds
      this->eta = p.eta;
   //   lam = ;
   //   bet = ;
      this->chi = p.chi;
      this->tc = p.tc;
      this->m1 = p.m1;
      this->m2 = p.m2;
      this->thetaS = p.thetaS;
      this->phiS = p.phiS;
      this->dist = p.dist;
      this->Ampl = p.Ampl;
      this->iota = p.iota;
      this->phi0 = p.phi0;
      this->psi = p.psi;

   //     double Norm;
      this->f0 = p.f0;
      this->fEnd = p.fEnd;
      this->SNR = p.SNR;
      
      return *this;
   
}  

void BBHTemplate::ComputeAuxParams(){
   
   if (m1 != 0.0 && m2 != 0.0){
      
      M = m1+m2;
      eta = m1*m2/(M*M);
      q = m2/m1;    
      
   }else{
      
      if (M != 0.0 && eta != 0.0 ){
         m1 = 0.5*M*(1.+ sqrt(1. - 4.*eta));
         m2 = 0.5*M*(1.- sqrt(1. - 4.*eta));
         q = m2/m1;
      }
      if (M != 0.0 && q != 0.0 ){
         m1 = M/(1.+q);
         m2 = q*M/(1.+q);
         eta = m1*m2/(M*M);
      }   
   }
   Mt = M*LISAWP_MTSUN_SI;
   Mc = pow(eta, 3./5.)*Mt;   
   /*
   double hpl = 1. + cos(iota)*cos(iota);
    double hcr = -2.*cos(iota);
    
    std::complex<double> img(0.,1.);
    
    double psi2 = 2.*psi;
    Fu = hpl*cos(psi2) - img*hcr*sin(psi2);
    Fv = -hpl*sin(psi2) - img*hcr*cos(psi2);
    
    uhat.resize(3);
    vhat.resize(3);
    n.resize(3);
    
    uhat[0] = cos(thetaS)*cos(phiS); 
    uhat[1] = cos(thetaS)*sin(phiS);
    uhat[2] = -sin(thetaS);
    
    vhat[0] = sin(phiS);
    vhat[1] = -cos(phiS);
    vhat[2] = 0.0;
    
    n[0] = sin(thetaS)*cos(phiS);
    n[1] = sin(thetaS)*sin(phiS);
    n[2] = cos(thetaS);
    */
   
}

}// end of the namespace

