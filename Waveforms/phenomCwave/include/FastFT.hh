/*
Copyright (C) 2001 R. Balasubramanian

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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/include/FastFT.hh,v 1.1 2007/09/18 21:50:36 stba Exp $" 
*/

#ifndef FASTFTHH
#define FASTFTHH
#include "fftw3.h"
#include "Matrix.hh"
#include "Macros.hh"
namespace LISAWP {
  /** 
   *Class FastFT: This is a  C++ Wrapper for FFTW. 
   This class uses C-functions for FFT. We use FFTW library for both precisions.

   This wrapper works only with real sequences of numbers. The Fourier Transform of the real
   sequency is stored in the same format as that of FFTW. In both the forward and inverse transforms
   the input array remains unchanged. If it is desired to have the input array and the output array the same
   then both the input matrix and output matrix can be the same. 

   Please see the test code base/test/testFastFT.cc for a test program.

   @author R. Balasubramanian, 2001.
  */
  class FastFT 
  {
  public:
    /** Constructor. */
    FastFT();
    /** Copy Constructor. */
    FastFT(const FastFT &fftWrap);
    
    /** Forward Fourier Transform
	@param mat1 reference to Input Matrix
	@param mat2 reference to output Matrix
    */
    void Forward(const Matrix<double> &mat1, Matrix<double> &mat2);

    /** Forward Fourier Transform.
	@param mat1 reference to Input Matrix
	@param mat2 reference to output Matrix (can be same as input)
    */
    void Forward(const Matrix<float> &mat1, Matrix<float>  &mat2);

    /** Inverse Fourier Transform for double matrices
	please note the input matrix will be destroyed
	@param mat1 reference to Input Matrix
	@param mat2 reference to output Matrix (can be same as input)
    */    
    void Inverse(const Matrix<double> &mat1, Matrix<double> &mat2);

    /** Inverse Fourier Transform for single precision transforms
	@param mat1 reference to Input Matrix
	@param mat2 reference to output Matrix (can be same as input)
    */
    void Inverse(const Matrix<float> &mat1, Matrix<float>  &mat2);

    /** Destructor
     */
    ~FastFT();

  private:
    bool IsPow2(unsigned num); // A private function to determine whether a number is a integral power of two.
    /// variables for single precision ffts. 
    unsigned int ssizeF; // forward size
    unsigned int ssizeR; // reverse size
    int sfSet;//flag
    int srSet;//flag
    fftwf_plan sfPlan;
    fftwf_plan srPlan;
    float *sfptr;
    float *srptr;
    // variables for double precision ffts.
    unsigned int dsizeF; 
    unsigned dsizeR;
    int dfSet;
    int drSet;
    fftw_plan dfPlan;
    fftw_plan drPlan;
    double *dfptr;
    double *drptr;
    
  };
} // namespace LISAWP
#endif










