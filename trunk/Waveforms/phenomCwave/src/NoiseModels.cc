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


#include "NoiseModels.hh"

namespace LISAWP{
	
NoiseModels::NoiseModels(bool galactic_bin){

   gal = galactic_bin;

}

void NoiseModels::StandardLISA_X(int sz, double* &freq, double* &S_n)
{
   
/*   Sacc = 2.53e-48 * fr**(-2.)
   Ssn =  1.12e-37*fr**2.
   Som = 6.33e-38*fr**2.
   S_instLISA = 16. * np.sin(2.*pi*fr*L)**2 * ( Ssn + Som + (3. + np.cos(4.*pi*fr*L)) * Sacc ) 
    (2.0 * L)**2 * (2*pi*fr)**2 * 4.0 * np.sin(om*L)**2 * (
         np.piecewise(fr,(fr >= 1.0e-4  ) & (fr < 1.0e-3  ),[lambda f: 10**-44.62 * f**-2.3, 0]) + \
         np.piecewise(fr,(fr >= 1.0e-3  ) & (fr < 10**-2.7),[lambda f: 10**-50.92 * f**-4.4, 0]) + \
         np.piecewise(fr,(fr >= 10**-2.7) & (fr < 10**-2.4),[lambda f: 10**-62.8  * f**-8.8, 0]) + \
         np.piecewise(fr,(fr >= 10**-2.4) & (fr < 10**-2.0),[lambda f: 10**-89.68 * f**-20.0,0])     )
         */

   // I presume that the memoryfor S_n was allocated elsewhere
   
   double Sacc = 2.53e-48;
   double Ssn =  1.12e-37;
   double Som = 6.33e-38;
   
   double L = 5.e9/LISAWP_C_SI;
   double omL, fr2, fr;
   double Sgal = 0.0;
   for (int i=0; i<sz; i++){
      omL = LISAWP_TWOPI*freq[i]*L;
      fr2 = pow(freq[i], 2.);
      S_n[i] = 16.*pow(sin(omL), 2.) * ((Ssn + Som)*fr2 + (3. + cos(2.*omL)) * Sacc/fr2 );    
   }
   
   if(gal){
      
      double L4 = L*L*4.;
      double om;
      for (int i=0; i<sz; i++){
         fr = freq[i];
         om = LISAWP_TWOPI*fr;
         omL = om*L;
         Sgal = 0.0;
         if ((fr >= 1.0e-4  ) && (fr < 1.0e-3 )){
            Sgal = L4 * pow(2.*om*sin(omL), 2.) * pow(10.,-44.62) * pow(fr,-2.3);
         }
         if ((fr >= 1.0e-3  ) && (fr < pow(10.,-2.7))){
            Sgal = L4 * pow(2.*om*sin(omL), 2.) * pow(10.,-50.92) * pow(fr, -4.4);
         }
         if ((fr >= pow(10.,-2.7)) & (fr < pow(10.,-2.4))){
            Sgal = L4 * pow(2.*om*sin(omL), 2.) * pow(10.,-62.8)  * pow(fr,-8.8);
         }
         if ((fr >= pow(10.,-2.4)) & (fr < 0.01)){
            Sgal = L4 * pow(2.*om*sin(omL), 2.) * pow(10.,-89.68) * pow(fr,-20.0);
         }
         S_n[i] += Sgal;
      }            
      
   }
   
/*** code by Michele 

Spm = 2.53654e-48 * f**(-2)
Sops = (1.1245e-37 * (L/5e9)**2 + 6.3253e-38) * f**2

x = 2.0 * pi * L * f

SX = 16.0 * sin(x)**2 * (2.0 * (1.0 + cos(x)**2) * Spm + Sop)
SAE = 8.0 * sin(x)**2 * (2.0 * Spm * (3.0 + 2.0*cos(x) + cos(2*x)) + Sop * (2.0 + cos(x)))
ST = 16.0 * Sop * (1.0 - cos(x)) * sin(x)**2 + 128.0 * Spm * sin(x)**2 * sin(0.5*x)**4


*/      
}

void NoiseModels::miniLISA_C1X(int sz, double* &freq, double* &S_n)
{
   
   /*# Config1

   Sacc = 8.17e-48*(1./fr + (1.8e-4/fr**2.) )**2.
   Ssn = 9.88e-37*fr**2.
   Som = 2.81e-38*fr**2.
   S_instC1 = 16. * np.sin(2.*pi*fr*L)**2 * ( Ssn + Som + (3. + np.cos(4.*pi*fr*L)) * Sacc ) */
   
   double Sacc = 8.17e-48;
   double Ssn =  9.88e-37;
   double Som = 2.81e-38;
   
   double L = 1.e9/LISAWP_C_SI;
   double omL, fr2, fr;
   double Saccf;
   double Sgal = 0.0;
   for (int i=0; i<sz; i++){
      omL = LISAWP_TWOPI*freq[i]*L;
      fr2 = pow(freq[i], 2.);
      Saccf = Sacc*pow( (1./freq[i] + (1.8e-4/fr2)), 2.);
      S_n[i] = 16.*pow(sin(omL), 2.) * ((Ssn + Som)*fr2 + (3. + cos(2.*omL)) * Saccf );    
   }
   
   if(gal){
       double L4 = L*L*4.;
       double om;
       for (int i=0; i<sz; i++){
          fr = freq[i];
          om = LISAWP_TWOPI*fr;
          omL = om*L;
          Sgal = 0.0;
          if ((fr>=4.5e-4) && (fr < 5.3e-4)){
              Sgal = L4 * pow(om*sin(omL), 2.)*0.6 *1.e-13 * pow(fr, 7.);
          }
          if ((fr>=5.3e-4) && (fr <= 5.88e-3)){
             Sgal = L4 * pow(om*sin(omL), 2.)* 0.6 *2.9714e-47 * pow(fr, -3.235);
          }
          S_n[i] += Sgal;
       }
   }
   
 /*   4.*L*L*(2.*pi*fr)**2. *3./5. * np.sin(2.*pi*L*fr)*np.sin(2.*pi*L*fr) *
       (
       	np.piecewise(fr, (fr>=4.5e-4) & (fr < 5.3e-4), [lambda f: 1.e-13 * f**7, 0]) +\
       	np.piecewise(fr, (fr>=5.3e-4) & (fr < 5.88e-3), [lambda f: 2.9714e-47 * f**(-3.235), 0]) 
       	)*/
    
   
}


void NoiseModels::miniLISA_C2X(int sz, double* &freq, double* &S_n)
{
   /*
   Sacc = 6.e-48 *(1./fr**2)
   Ssn = 6.15e-40*fr**2.
   Som = 2.81e-38*fr**2.
   S_instC2 = 16. * np.sin(2.*pi*fr*L)**2 * ( Ssn + Som + (3. + np.cos(4.*pi*fr*L)) * Sacc ) */
   
     double Sacc = 6.e-48;
     double Ssn =  6.15e-40;
     double Som = 2.81e-38;

     double L = 1.e9/LISAWP_C_SI;
     double omL, fr2, fr;
     double Sgal = 0.0;
     for (int i=0; i<sz; i++){
        omL = LISAWP_TWOPI*freq[i]*L;
        fr2 = pow(freq[i], 2.);
        S_n[i] = 16.*pow(sin(omL), 2.) * ((Ssn + Som)*fr2 + (3. + cos(2.*omL)) * Sacc/fr2 );    
     }

     if(gal){
         double L4 = L*L*4.;
         double om;
         for (int i=0; i<sz; i++){
            fr = freq[i];
            om = LISAWP_TWOPI*fr;
            omL = om*L;
            Sgal = 0.0;
            if ((fr>=4.5e-4) && (fr < 5.3e-4)){
                Sgal = L4 * pow(om*sin(omL), 2.)*0.6 *1.e-13 * pow(fr, 7.);
            }
            if ((fr>=5.3e-4) && (fr < 2.2e-3)){
               Sgal = L4 * pow(om*sin(omL), 2.)* 0.6 *2.9714e-47 * pow(fr, -3.235);
            }
            if ( (fr>=2.2e-3) && (fr < 4.e-3)){
               Sgal = L4 * pow(om*sin(omL), 2.)* 0.6 *1.517e-51 * pow(fr, -4.85);
            }
            if ((fr>=4e-3) && (fr <= 5.88e-3)){
               Sgal = L4 * pow(om*sin(omL), 2.)* 0.6 *6.706e-58 * pow(fr, -7.5);
            }
            S_n[i] += Sgal;
         }
     }
   
   /* 4.*L*L*(2.*pi*fr)**2. *3./5. * np.sin(2.*pi*L*fr)*np.sin(2.*pi*L*fr) * 
   (
   	np.piecewise(fr, (fr>=4.5e-4) & (fr < 5.3e-4), [lambda f: 1.e-13 * f**7, 0]) +\
   	np.piecewise(fr, (fr>=5.3e-4) & (fr < 2.2e-3), [lambda f: 2.9714e-47 * f**(-3.235), 0]) +\
   	np.piecewise(fr, (fr>=2.2e-3) & (fr < 4.e-3), [lambda f: 1.517e-51 * f**(-4.85), 0]) +\
   	np.piecewise(fr, (fr>=4e-3) & (fr < 5.88e-3), [lambda f: 6.706e-58 * f**(-7.5), 0])
   		)
   */
   
}

} //end of the namespace