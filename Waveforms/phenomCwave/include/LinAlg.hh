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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/include/LinAlg.hh,v 1.1 2007/09/18 21:50:37 stba Exp $" 
*/          

#ifndef LINALGHH
#define LINALGHH
#include <math.h>
#include "Macros.hh"
#include "Matrix.hh"
#include <vector>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_eigen.h"


namespace LISAWP{

  /** Collection of Linear Algebra routines. Computations are based on
   * gsl library (wraper). 
   * @author S. Babak 2005
   */

class LinAlg{
public:
   /** Constructor 
    * @param M matrix which we will work with 
    * */

    LinAlg(Matrix<double>& M);
   
    /** Destructor */
    ~LinAlg();
	
 /** Computing Singular value decomposition of M returns 
     * the vector of matrices: U, S, V:  M = U S V^{transp}
     * U orthogonal (m x n) S diagonal (n x n), V orthogonal (n x n)
     * (m >= n) I use here gsl function.
     * Matrices U,S,V must be empty 
     **/

    void SVDecomp(Matrix<double>& U, Matrix<double>& S, Matrix<double>& V);
  
 /** This function computes the eigenvalues and eigenvectors of the REAL SYMMETRIC matrix M
  * The eigenvalues are stored in the matrix Lambda and are ordered. The corresponding 
  * eigenvectors are stored in the columns of the matrix Vecs. All vectors are orhtonormal.
  * Lambda and Vecs must be empty
  */

    void Eigensys(Matrix<double>& Lambda, Matrix<double>& Vecs);

    /** This function computes LU decomposition of the Matrix (permutations)
     * returns sign of permutation. See test code for usage.
     * @param L lower part
     * @param U uppr part
     * @param perm permutation vector
     * */

    int LUDecomp(Matrix<double>& L, Matrix<double>& U, std::vector<int>& perm);

    
/** Computes determinant of the matrix */
    
    double Determinant();


  private:

    unsigned nRows;
    unsigned nCols;
    gsl_matrix *m;
    double det;

};
    
}//end of namespace

#endif
