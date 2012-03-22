/*
Copyright (C) 2012 S. Babak 

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

#include "MCintegration.h"


/**  The base class for doing multidimensional integration (vegas method) 
@author Stas Babak, Jan 2012.
*/

namespace LISAWP{


MTint::MTint(){

   lbnd = NULL;
   ubnd = NULL;
   bins = NULL;
   bxs = NULL;
   binVals = NULL;

   K = 1000;
   alpha = 1.5; // damping rate
   smooth = true;

   std::cout << "---> Use K= " << K << "  and alpha = " << alpha << std::endl;
}

MTint::~MTint(){

   if(lbnd != NULL)
      delete [] lbnd;
   if(ubnd != NULL)
      delete [] ubnd;
   if (bins != NULL)
      delete [] bins;
   if (bxs != NULL){
	   // deallocated in the destructor of boxes
   /*  for (int i=0; i< NumBoxes; i++){ 
       delete [] bxs[i].indx;       
       delete [] bxs[i].xl;
       delete [] bxs[i].dx;
       delete [] bxs[i].p;
       delete [] bxs[i].bins;
     }*/
     delete [] bxs;
   }

   if (binVals != NULL){
       for (int i=0; i<Ndim; i++){
           delete [] binVals[i];
       }
       delete [] binVals;
  }


}


void MTint::SetBoundaries(int dim, double* &lowBound, double* &upperBound)
{
   Ndim = dim;	
   if (lbnd == NULL)
         lbnd = new double[dim];
   if (ubnd == NULL)
         ubnd = new double[dim];

   for (int i =0; i<Ndim; i++){
       lbnd[i] = lowBound[i];
       ubnd[i] = upperBound[i];
   }
 	
}

void MTint::SetNumBins(int* &Bins){
  
  bins = new int[Ndim];
  NumBoxes = 1;
  maxI = Bins[0];
  for (int i= 0; i<Ndim; i++){
     bins[i] = Bins[i];
     NumBoxes *= Bins[i];
     if (bins[i] > maxI)
	maxI = bins[i];
  }

  binVals = new double*[Ndim];
  for (int i=0; i<Ndim; i++){
      binVals[i] = new double[bins[i] + 1];
  }
  
  double dx, runx;
  for (int i=0; i<Ndim; i++){
     dx = (ubnd[i] - lbnd[i])/(double)bins[i];
     runx = lbnd[i];
     for (int j=0; j<bins[i] +1; j++){
	 binVals[i][j] = runx;
	 runx += dx;
     }
  }

  /*for (int i=0; i< Ndim; i++){
     for (int j=0; j< bins[i]+1; j++){
	   std::cout << i << "   " << j << "    " << binVals[i][j] << std::endl;
     }
  }
  exit(0);*/


}


void MTint::InitializeBoxes(){

   if (bxs == NULL){
      bxs = new Box[NumBoxes];
   }

   std::cout << "  Num Boxes = " << NumBoxes << std::endl;
   double* delx;
   delx = new double[Ndim];
   double* prob;
   prob = new double[Ndim];
   for(int j = 0; j<Ndim; j++){
       delx[j] = (ubnd[j] - lbnd[j])/(double)bins[j];
       prob[j] = 1.0/((double)bins[j]*delx[j]);
   }

   for (int i=0; i< NumBoxes; i++){
       bxs[i].n = Ndim;
       bxs[i].bins = new int[Ndim];
       bxs[i].indx = new int[Ndim];       
       bxs[i].xl = new double[Ndim];
       bxs[i].dx = new double[Ndim];
       bxs[i].p = new double[Ndim];
       bxs[i].prob = 1.0;
       for (int j=0; j<Ndim; j++)
	   bxs[i].bins[j] = bins[j];    
       bxs[i].InitializeBox(i);
       bxs[i].vol = 1.0;
       //std::cout << i << "    { ";
       for (int j=0; j<Ndim; j++){
           //std::cout << i << "    " << bxs[i].bins[j];
	   bxs[i].dx[j] = delx[j];
	   bxs[i].xl[j] = (double)bxs[i].indx[j] * delx[j] + lbnd[j];
	   bxs[i].p[j] = prob[j];
	   bxs[i].prob *= prob[j];
	   bxs[i].vol *= delx[j];
	   //std::cout << bxs[i].xl[j] << ",   ";
	   //std::cout << bxs[i].p[j] << ",   ";
	   //std::cout   << bxs[i].indx[j] << ",   ";
       }
       //std::cout << " }\n";
   }
   delete [] delx;
   delete [] prob;
}


void MTint::Integrate(int seed, int numIter, int Ncall, double* &params, \
		double &val, double &stdev, double &chi2){
 
  calls = Ncall;
  gsl_rng_env_setup();
  const gsl_rng_type *T;
  gsl_rng * r;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);

  double* IntIter;
  double* SigIter;
  IntIter = new double[numIter];
  SigIter = new double[numIter];

  //std::ofstream fout1503("Data/MCdistr.dat"); 
  for(int iter = 0; iter<numIter; iter++){

     for (int k=0; k<NumBoxes; k++){	
        EvaluateBox(k, r, params);
	if (iter == 0){
           for (int i=0; i<Ndim; i++){
	      //fout1503 << bxs[k].xl[i] <<  "     "; 
	   }
	   //fout1503 << sqrt(bxs[k].Fv2) << "    " << bxs[k].prob << std::endl;
	}
	
     }
     std::cout << "\n \n" << std::endl;

     EvaluateIntegral(IntIter[iter], SigIter[iter]);

     std::cout << "iter = " << iter << "  int = " << IntIter[iter] <<\
	     "  sig = " << sqrt(SigIter[iter]) << std::endl;


     //AdjustPDF();
     ChangeBins();
     //if (iter == 0)
     //  exit(0);

  }
  //fout1503.close();
  

  // Computing the final answer

  val = 0.0;
  stdev = 0.0;

  for (int iter=0; iter<numIter; iter++){
      val += IntIter[iter]/SigIter[iter];
      //std::cout << iter << "  val =  " << val << std::endl;
      stdev += 1./SigIter[iter];
  }

//  std::cout << "stdev = " << stdev << std::endl;
  stdev = 1./stdev;

  val *= stdev;
  stdev = sqrt(stdev);

  // alternative
  
  /*for (int iter=0; iter<numIter; iter++){
      val += IntIter[iter]*pow(IntIter[iter]/SigIter[iter], 2.0);
      stdev += pow(IntIter[iter]/SigIter[iter], 2.0) ;
  }

  val = val/stdev;
  stdev = val/sqrt(stdev);*/
	
  chi2 = 0.0;
  for (int iter=0; iter<numIter; iter++){
     chi2 += pow(IntIter[iter] - val, 2.0)/SigIter[iter];
     /* for alternative */
     //chi2 +=pow((IntIter[iter] - val)*IntIter[iter]/val, 2.0)/SigIter[iter]; 
  }

  chi2 = chi2/(double)(numIter - 1.);

  gsl_rng_free(r);
  delete [] IntIter;
  delete [] SigIter;


}

void MTint::EvaluateBox(int k, gsl_rng *r, double* &params){

  double* xcoord;
  xcoord = new double[Ndim];  
  //std::cout << k <<  "  coord =  \n";

  bxs[k].Fv = 0.0;
  bxs[k].Fv2 = 0.0;

  double evf;
 /* for (int j=0; j<Ndim; j++){
      std::cout << bxs[k].xl[j] << " - " << bxs[k].xl[j] + bxs[k].dx[j] << "    ";
  }
  std::cout << std::endl;  */

  for (int j=0; j<calls; j++){
      for (int j=0; j<Ndim; j++){
         xcoord[j] = bxs[k].xl[j] +  bxs[k].dx[j] * gsl_rng_uniform(r); 
         //std::cout << xcoord[j] << ",   "; 
      }
      evf = CompFunct(xcoord, params);
      //std::cout << "  " << evf;
      bxs[k].Fv += evf;
      bxs[k].Fv2 += evf*evf;
  }
  //std::cout  << std::endl;

}

void MTint::EvaluateIntegral(double &val, double &stdev){

   double TotCalls = (double)calls*NumBoxes;
   std::cout << " tot # of calls = " << TotCalls << std::endl;  
   val = 0.0;
   stdev = 0.0;

   double check  = 0.0;
   double ttt = 1.0;
   for (int i=0; i<NumBoxes; i++){
       val += bxs[i].Fv/bxs[i].prob;
       stdev += bxs[i].Fv2/(bxs[i].prob*bxs[i].prob);
       ttt = 1.0;
       for (int I=0; I<Ndim; I++){
	   ttt *= bxs[i].dx[I];
       }
       check += bxs[i].prob*ttt;
   }
   std::cout << "check = " << check << std::endl;
   std::cout << "val = " << val << std::endl;

   val *= (1./TotCalls);
   stdev *= (1./TotCalls);

   stdev = (stdev - val*val)/(TotCalls-1.);


}

void MTint::ChangeBins(){

  double** fi = NULL;
  fi = new double*[Ndim];

  int* m = NULL;
  m = new int[maxI];

  for (int i =0; i< Ndim; i++)
      fi[i] = new double[bins[i]];

  // recomputing the probability density function
  for (int i=0; i< Ndim; i++){
     for (int j=0; j< bins[i]; j++){
	 fi[i][j] = 0.0;    
	 for (int k=0; k<NumBoxes; k++){
	     if (bxs[k].indx[i] == j){
		//fi[i][j]  += bxs[k].Fv2 * pow(bxs[k].p[i]/bxs[k].prob, 2.); 
		//falsh!!!!!
		fi[i][j]  +=  bxs[k].Fv2 * (bxs[k].p[i]/bxs[k].prob)*bxs[k].vol; 
	     }
	 }
	 fi[i][j] *= (binVals[i][j+1] - binVals[i][j]); //  * dx[I,j]
     }
  }

  // smoothening the probability
  if (smooth){
     double oldg;
     double newg, rc;
     for (int i=0; i< Ndim; i++){
        oldg = fi[i][0];
        newg = fi[i][1];
        fi[i][0] = 0.5*(oldg + newg);
        for (int j=1; j< bins[i]-1; j++){
            rc = oldg + newg;
	    oldg = newg;
	    newg = fi[i][j+1];
	    fi[i][j] = (rc + newg)/3.0;
        }
        fi[i][bins[i]-1] = 0.5*(newg + oldg);
     }
  }

  /*std::cout << " new bins: \n";
  for (int i=0; i< Ndim; i++){
     for (int j=0; j< bins[i]; j++){
	   std::cout << i << "   " << binVals[i][j] << "    " << fi[i][j] << std::endl;
     }
  }*/
   
  
  // adjust te size of bins
  
  double TotI;
  for (int i=0; i<Ndim; i++){ 
     TotI = 0.0;
     for (int j=0; j< bins[i]; j++){
	 TotI += fi[i][j];
     }	     
     for (int j=0; j< bins[i]; j++){
	 if (fi[i][j] == 0.0){
	     m[j] = 1;
	 }else{	     
	   m[j] = (int) floor(pow((double)K*(fi[i][j]/TotI - 1.0)/log(fi[i][j]/TotI), alpha));
	 }
	 //std::cout << i << "   " << j << "   " << m[j] << std::endl; 
     }

     Rebin(m, i);
     //exit(0);
  }

 /* std::cout << " new bins: \n";
  for (int i=0; i< Ndim; i++){
     for (int j=0; j< bins[i]+1; j++){
	   std::cout << i << "   " << j << "    " << binVals[i][j] << std::endl;
     }
  }*/

  // Now need to change dx and xl for each box.

  for (int k=0; k<NumBoxes; k++){	  
     bxs[k].prob = 1.0;
     bxs[k].vol = 1.0;
     for (int i=0; i<Ndim; i++){
        bxs[k].xl[i] = binVals[i][bxs[k].indx[i]];
	bxs[k].dx[i] = binVals[i][bxs[k].indx[i]+1] - bxs[k].xl[i];
	//std::cout << "dx[" << i << "]   =  " << bxs[k].dx[i] <<"   ";
	bxs[k].p[i] = 1.0/((double)bins[i]*bxs[k].dx[i]);
	bxs[k].prob *= bxs[k].p[i];
	bxs[k].vol *= bxs[k].dx[i]; 
     }
     //std::cout << std::endl;
  }

  //exit(0);

  for (int i =0; i< Ndim; i++)
      delete [] fi[i];
  delete [] fi;
  delete [] m;

}


void MTint::Rebin(int* &m, int I){

  int RebinTot = 0.0;
  for (int j=0; j<bins[I]; j++){
      RebinTot += m[j];
      //std::cout << RebinTot << "   " << m[j] << std::endl;
  }
  std::cout << I <<"-th direction \n"; 
  int NperBin = (int) floor(RebinTot/(double)bins[I]);
  std::cout << "tot num " << RebinTot << " #/bin = " << NperBin << std::endl; 

  int xind = 0; 
  double xrun = binVals[I][0];
  double dxi;

  double* tmp;
  tmp = new double[bins[I]+1];
  for (int j=1; j<bins[I]+1; j++){
      if (m[j-1] != 0){	  
         dxi = (binVals[I][j] - binVals[I][j-1])/(double)m[j-1];
      }else{
	 dxi = binVals[I][j] - binVals[I][j-1];
	 xrun += dxi; 
	 xind ++;
      }
      //std::cout << j <<  " dxi = " << dxi << "   " << m[j-1] << "   " << xrun << std::endl; 
      for (int k=0; k<m[j-1]; k++){ 
	  if(xind%NperBin == 0){
		  //std::cout << "test  " << xind << "   " << xind/NperBin << "   " << dxi \
			  << "    " << m[j-1] << "   " << xrun << std::endl; 	  
              tmp[xind/NperBin] = xrun;
	  }
	  xrun += dxi;
	  xind ++;
      }
  }
  for (int j=0; j<bins[I]; j++){
     binVals[I][j] = tmp[j];
     //std::cout << I << "  " << j << "   " << binVals[I][j] << std::endl; 
  }
  //exit(0);

  delete [] tmp;

}

void MTint::GetProbability(double** &coord, double* &prob){

    if (coord == NULL){
       std::cout << "I assume that coord: (NumBoxes x dimension) already allocated \n";
       exit(1);
    }
    if (prob == NULL){
       std::cout << "I assume that prob: (NumBoxes) already allocated \n";
       exit(1);
    }
    

    for (int k=0; k<NumBoxes; k++){
       for (int i =0; i<Ndim; i++){
	      coord[k][i]  = bxs[k].xl[i] + 0.5*bxs[k].dx[i]; 
       }
       prob[k] =  bxs[k].prob;     
    } 


}



/******      Box class   *********/

Box::Box(){

  indx = NULL;
  xl = NULL;
  dx = NULL;
  bins = NULL;
  p = NULL;



}


Box::~Box(){

   if (indx != NULL)
	delete [] indx;
   if (xl != NULL)
	delete [] xl;
   if (dx != NULL)
	delete [] dx;
   if (bins != NULL)
   	delete [] bins;
   if (p != NULL)
	delete [] p;

}


void Box::InitializeBox(int num){
   
     int D=0;
     CompCoord(num, D);     

}

void Box::CompCoord(int nn, int& D){
  
   if(indx == NULL || bins == NULL){
	std::cout << "you must allocate space for indx" << std::endl;
        exit(1);
   }
   

   int MD = 1;
   for (int i= D+1; i<n; i++){
       MD = MD * bins[i];
   }

   indx[D] = (int) floor(nn/MD);

   //std::cout << "  nn = " << nn << "  D = " << D <<  "   MD = " << MD <<  \
	  "  indx[D] = " << indx[D] << "   " << nn - indx[D]*MD  << std::endl;
   int mm; 
   if (D+1 < n){
        mm =nn - indx[D]*MD; 
	D++;
	CompCoord(mm, D);
   }	

}  

TestVegas::TestVegas(){

   std::cout << "Poehali\n";
}

double TestVegas::CompFunct(double* &coord, double* &pars){

    //double A = 1.0 / (LISAWP_PI * LISAWP_PI * LISAWP_PI);
    //return A / (1.0 - cos (coord[0]) * cos (coord[1]) * cos (coord[2])); 
   
    // sinc(x)*sinc(y)
    /*if (coord[0] == 0.0){
	if (coord[1] == 0.0){
	   return(1.0);
	}else{
	   return(sin(coord[1])/coord[1]);
	}
    }else{
	if (coord[1] == 0.0){
	   return(sin(coord[0])/coord[0]);
	}else
           return( sin(coord[0])*sin(coord[1])/(coord[0]*coord[1]) );
    }*/

    /*if (coord[0] == 0.0 || coord[1] == 0.0){
	return(1.0);
    }else{
        return( sin(coord[0]*coord[1])/(coord[0]*coord[1]) );
    } */

    //return( exp(coord[0]*cos(coord[1])) );
    //return( sin(coord[0])*cos(coord[1]) );
    //
    // torus

    double psi, theta; 
    double x = coord[0];
    double y = coord[1];

    
    if (x != 0.0)
       theta = atan2(y,x);
    else
       theta = 0.5*LISAWP_PI;

    double cpsi;
    //std::cout << theta + LISAWP_TWOPI  << std::endl;
    //exit(0); 
    if (x == 0.0 && y== 0.0)
	   return(0.0);
    if (x != 0.0){
 	 cpsi = 0.5*(x/cos(theta) -4.0);
	 if (fabs(cpsi) > 1.0)
	      return(0.0);
	 else
	     psi = acos(cpsi);
    }else{
	 cpsi = 0.5*(y/sin(theta) -4.0);
         if (fabs(cpsi) > 1.0)
	      return(0.0);
	 else
	     psi = acos(cpsi);
    }
    //std::cout  << "psi = " << psi << std::endl;
    return(2.0*sin(psi));
        

    //double A = 1.0 / (LISAWP_PI * LISAWP_PI * LISAWP_PI);
   //return A / (1.0 - cos (coord[0]) * cos (coord[1]) * cos (coord[2]));

}


}// end of the namespace


