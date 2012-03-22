/*
Copyright (C) 2001  R. Balasubramanian

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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/src/FastFT.cc,v 1.1 2007/09/18 21:50:41 stba Exp $" 
*/

#include <cstring>
#include <cstdlib>
#include "FastFT.hh"

namespace LISAWP{

  FastFT::FastFT()
  {
    ssizeF = 0;
    ssizeR = 0;
    sfSet = 0;
    srSet = 0;
    srPlan = sfPlan = 0;
    dsizeF = 0;
    dsizeR = 0;
    dfSet = 0;
    drSet = 0;
    drPlan = 0;
    dfPlan =  0;    
    sfptr = srptr = 0;
    dfptr = drptr = 0;
  }
  
  FastFT::FastFT(const FastFT &fftWrap)
  {
    ssizeF = 0;
    ssizeR = 0;
    sfSet = 0;
    srSet = 0;
    srPlan = sfPlan = 0;
    dsizeF = 0;
    dsizeR = 0;
    dfSet = 0;
    drSet = 0;
    drPlan = 0;
    dfPlan =  0;
    sfptr = srptr = 0;
    dfptr = drptr = 0;
  }
  
  
  FastFT::~FastFT()
  {
    if(sfSet)
      fftwf_destroy_plan(sfPlan);
    if(srSet)
      fftwf_destroy_plan(srPlan);
    if(dfSet)
      fftw_destroy_plan(dfPlan);
    if(drSet)
      fftw_destroy_plan(drPlan);
    if(sfptr)
      fftwf_free(sfptr);
    if(srptr)
      fftwf_free(srptr);
    if(dfptr)
      fftw_free(dfptr);
    if(drptr)
      fftwf_free(drptr);
  }
  

  void
  FastFT::Forward(const Matrix<double> &m1, Matrix<double> &m2)
  {
    LISAWPAssert(m1.GetNRows()==1,"The matrix must be a row matrix ");
    LISAWPAssert(m2.GetNRows()==1,"The matrix must be a row matrix ");
    LISAWPAssert(m1.GetNCols() == m2.GetNCols(), "Matrix Sizes must match ");
    LISAWPAssert(m1.Size() > 0, " Matrix must be of non zero size ");
    if(((m1.Size() != dsizeF))||(!dfSet)){
      dsizeF = m1.Size();
      if(dfptr)
	fftw_free(dfptr);
      dfptr = (double *) fftw_malloc(sizeof(double)*dsizeF);
      if(dfSet){	
	fftw_destroy_plan(dfPlan);
      }
      else
	dfSet=1;
      dfPlan = fftw_plan_r2r_1d(dsizeF,dfptr,dfptr,FFTW_R2HC,FFTW_ESTIMATE);
    }   
    memcpy((char *)dfptr, (char *) m1.ConstOrigin(), dsizeF*sizeof(double));
    fftw_execute(dfPlan);
    memcpy((char *) m2.Origin(), (char *) dfptr, dsizeF*sizeof(double));
  }
  
  
  void
  FastFT::Inverse(const Matrix<double> &m1, Matrix<double> &m2)
  {    
    LISAWPAssert(m1.GetNRows()==1,"The matrix must be a row matrix ");
    LISAWPAssert(m2.GetNRows()==1,"The matrix must be a row matrix ");
    LISAWPAssert(m1.GetNCols() == m2.GetNCols(), "Matrix Sizes must match ");
    LISAWPAssert(m1.Size() > 0, " Matrix must be of non zero size ");
    if(((m1.Size() != dsizeR))||(!drSet)){
      dsizeR = m1.Size();
      if(drptr)
	fftw_free(drptr);
      drptr = (double *)fftw_malloc(sizeof(double)*dsizeR);
      if(drSet){	
	fftw_destroy_plan(drPlan);
      }
      else
	drSet=1;
      drPlan = fftw_plan_r2r_1d(dsizeR,drptr,drptr,FFTW_HC2R,FFTW_ESTIMATE);
    }   
    memcpy((char *)drptr, (char *) m1.ConstOrigin(), dsizeR*sizeof(double));
    fftw_execute(drPlan);
    memcpy((char *) m2.Origin(), (char *) drptr, dsizeR*sizeof(double));

  }


  void
  FastFT::Forward(const Matrix<float> &m1, Matrix<float> &m2)
  {
    LISAWPAssert(m1.GetNRows()==1,"The matrix must be a row matrix ");
    LISAWPAssert(m2.GetNRows()==1,"The matrix must be a row matrix ");
    LISAWPAssert(m1.GetNCols() == m2.GetNCols(), "Matrix Sizes must match ");
    LISAWPAssert(m1.Size() > 0, " Matrix must be of non zero size ");
    if(((m1.Size() != ssizeF))||(!sfSet)){
      ssizeF = m1.Size();
      if(sfptr)
	fftwf_free(sfptr);
      sfptr = (float *) fftwf_malloc(sizeof(float)*ssizeF);
      if(sfSet){	
	fftwf_destroy_plan(sfPlan);
      }
      else
	sfSet=1;
      sfPlan = fftwf_plan_r2r_1d(ssizeF,sfptr,sfptr,FFTW_R2HC,FFTW_ESTIMATE);
    }   
    memcpy((char *)sfptr, (char *) m1.ConstOrigin(), ssizeF*sizeof(float));
    fftwf_execute(sfPlan);
    memcpy((char *) m2.Origin(), (char *) sfptr, ssizeF*sizeof(float));
  }
  

  void
  FastFT::Inverse(const Matrix<float> &m1, Matrix<float> &m2)
  {    
    LISAWPAssert(m1.GetNRows()==1,"The matrix must be a row matrix ");
    LISAWPAssert(m2.GetNRows()==1,"The matrix must be a row matrix ");
    LISAWPAssert(m1.GetNCols() == m2.GetNCols(), "Matrix Sizes must match ");
    LISAWPAssert(m1.Size() > 0, " Matrix must be of non zero size ");
    if(((m1.Size() != ssizeR))||(!srSet)){
      ssizeR = m1.Size();
      if(srptr)
	fftwf_free(srptr);
      srptr = (float *)fftwf_malloc(sizeof(float)*ssizeR);
      if(srSet){	
	fftwf_destroy_plan(srPlan);
      }
      else
	srSet=1;
      srPlan = fftwf_plan_r2r_1d(ssizeR,srptr,srptr,FFTW_HC2R,FFTW_ESTIMATE);
    }   
    memcpy((char *)srptr, (char *) m1.ConstOrigin(), ssizeR*sizeof(float));
    fftwf_execute(srPlan);
    memcpy((char *) m2.Origin(), (char *) srptr, ssizeR*sizeof(float));

  }

} // namespace LISAWP





