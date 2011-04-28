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
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/src/LinAlg.cc,v 1.1 2007/09/18 21:50:42 stba Exp $" 
*/                


#include "LinAlg.hh"

namespace LISAWP{
          
LinAlg::LinAlg( Matrix<double>& M ){

   nRows = M.GetNRows();
   nCols = M.GetNCols();
   
   m = gsl_matrix_alloc(nRows, nCols);
   
   for (unsigned i=0; i<nRows; i++){
	 for (unsigned j=0; j<nCols; j++){
                gsl_matrix_set(m, i, j, M(i,j));
	   }	
   }
   
}

 LinAlg::~LinAlg(){

   gsl_matrix_free(m);

 }   
	
void LinAlg::SVDecomp(Matrix<double>& U, Matrix<double>& S, Matrix<double>& V){

	
     LISAWPAssert(nRows >= nCols, "Number of rows must be larger than columns");

           
     gsl_matrix *m1 = gsl_matrix_alloc(nRows, nCols);
     gsl_matrix *v = gsl_matrix_alloc(nCols, nCols);
     gsl_vector *s = gsl_vector_alloc(nCols);

     gsl_matrix_memcpy(m1, m);
     
     gsl_linalg_SV_decomp_jacobi(m1, v, s);
	       
     V.ExpandEmpty(nCols, nCols);
     S.ExpandEmpty(nCols, nCols);
     U.ExpandEmpty(nRows, nCols);
 
     S = 0.0;
       
     for (unsigned j=0; j<nCols; j++){
          for (unsigned i=0; i<nRows; i++){
		U(i,j) = gsl_matrix_get(m1, i, j);
		if (i < nCols){
		     V(i,j) = gsl_matrix_get(v, i, j);
		}
	  }
          S(j,j) = gsl_vector_get(s, j);
     } 
       
     gsl_matrix_free(v);
     gsl_matrix_free(m1);
     gsl_vector_free(s);

   }

 void LinAlg::Eigensys(Matrix<double>& Lambda, Matrix<double>& Vecs){

    LISAWPAssert(nRows == nCols, "Number of rows and columns must be the same");


    gsl_matrix *evec = gsl_matrix_alloc(nCols, nCols);
    gsl_vector *eval = gsl_vector_alloc(nCols);
    gsl_matrix *m1 = gsl_matrix_alloc(nCols, nCols);
    
    gsl_matrix_memcpy(m1, m);

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(nCols);

    gsl_eigen_symmv( m1, eval, evec, w);

    gsl_eigen_symmv_free(w);

    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    Lambda.ExpandEmpty(nCols, nCols);
    Vecs.ExpandEmpty(nCols, nCols);

    Lambda = 0.0;

    for (unsigned j=0; j<nCols; j++){
          for (unsigned i=0; i<nRows; i++){
		Vecs(i,j) = gsl_matrix_get(evec, i, j);
	  }
          Lambda(j,j) = gsl_vector_get(eval, j);
     } 

    
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    gsl_matrix_free(m1);
    
 }


 int LinAlg::LUDecomp(Matrix<double>& L, Matrix<double>& U, std::vector<int>& perm){

    LISAWPAssert(nRows == nCols, "Number of rows and columns must be the same");

    gsl_matrix *m1 = gsl_matrix_alloc(nCols, nCols);
    gsl_matrix_memcpy(m1, m);

    gsl_permutation *p = gsl_permutation_alloc(nCols);
    int s;
    
    gsl_linalg_LU_decomp( m1, p, &s);
    
    det  = gsl_linalg_LU_det(m1, s);

    L.ExpandEmpty(nCols, nCols);
    U.ExpandEmpty(nCols, nCols);
    L=0.0;
    U=0.0;
    
    double detCheck = 1.0;
    for (unsigned i=0; i<nRows; i++){
        for (unsigned j=0; j<nCols; j++){
	    if(i <= j){
		U(i,j) = gsl_matrix_get(m1, i, j);
	    }else{
		L(i,j) = gsl_matrix_get(m1, i, j);
	    }
	}
        L(i,i) = 1.0;
        detCheck *= U(i,i);
	perm.push_back(gsl_permutation_get(p, i));
    }
    detCheck *= s;

    LISAWPAssert((detCheck-det)<=1.e-5, "determinants do not agree");
        
    gsl_permutation_free(p);
    gsl_matrix_free(m1);

    return(s); 

 }

 double LinAlg::Determinant(){

    Matrix<double> L;
    Matrix<double> U;
    std::vector<int> p; 
    
    int s = LUDecomp(L, U, p);
    
    return(det);
 }
	
}//end of namespace

