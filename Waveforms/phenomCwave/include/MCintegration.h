/*
Copyright (C) 2012, S. Babak.  

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


#ifndef MCINTEGRATIONHH
#define MCINTEGRATIONHH

#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include "Constants.hh"
#include <gsl/gsl_rng.h>




namespace LISAWP{

 class Box{


   public:
     
     Box();
     ~Box();

     void CompCoord(int nn, int& D);
     void InitializeBox(int num);

     int n; // number of dimensions
     int* bins;     
     int* indx;
     double* xl;
     double* dx;
     double* p;
     double Fv;
     double Fv2;
     double Intg;
     double sigma;
     double prob;
     double vol; 
       
     

 };



  /** Abstract class MCint.
      And derived classes.
      @author Stas Babak, January 2012.
  */

  class MTint{

  protected:
    
    int calls;
    int Ndim;
    double* lbnd;
    double* ubnd;
    int NumBoxes;
    int* bins;
    Box* bxs;
    double** binVals;
    int maxI;

    /** pure virtual function which generates a window-vector */


    virtual void EvaluateIntegral(double &intval, double &sigval);

    virtual void ChangeBins();

    virtual void Rebin(int* &m, int I);


  public:
  
    MTint(); 

    ~MTint();
   

    virtual void SetNumBins(int* &Bins);

    virtual void SetBoundaries(int dim, double* &lowBound, double* &upperBound);

    /** Fills up supplied vector with window' values.
	@param win window-vector, \f$w_n\f$
    */

    virtual void InitializeBoxes();
    
    virtual void EvaluateBox(int i, gsl_rng *r, double* &params);

    virtual void Integrate(int seed, int numIter, int Ncall, double* & params, \
		    double &val, double &stdev, double &chi2);
   
    virtual void GetProbability(double** &coord, double* &prob);
    
    virtual double CompFunct(double* &coord, double* &pars)=0;
    int K;
    double alpha;
    bool smooth;

};




// Derive a test class
//
class TestVegas:public MTint{
   
    private:


    public:
       
       TestVegas();    
       double CompFunct(double* &coord, double* &pars); 	    

};


} // end of the namespace 


  /*template <class TYPE>
  class Hann: public Window<TYPE> {

  private:

    /// generates Hanning window

    void GenerateWindow();

  public:
   

    Hann():Window<TYPE>(){}
 
  };*/


#endif
