
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include "Constants.hh"
#include "MCintegration.h"
#include "GSL_vegas.h"



using namespace LISAWP;

double g (double *k, size_t dim, void *params)
{
       double A = 1.0 / (LISAWP_PI * LISAWP_PI * LISAWP_PI);
       return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
}

int main(){

   /* The answewr for integration 
    * of sinc(x)*sinc(y) in the interval, {-4pi, -4pi} {4pi, 4pi}
    * is   8.906180496 
    *
    * The integral sin(x*y)/(x*y) in the same limits is
    * 35.43242500 
    * */

   
    std::string spr = "  ";
    TestVegas TV;
   
   //int dim = 3;
   int dim = 2;
   int* NumBins;
   NumBins = new int[dim];
   NumBins[0] = 250;
   NumBins[1] = 250;
   //NumBins[2] = 50;

   double* lx;
   double* ux;
   lx = new double[dim];
   ux = new double[dim];

   double tiny = 1.e-10;
   /*lx[0] = tiny;
   lx[1] = tiny;
   lx[2] = tiny;

   ux[0] = LISAWP_PI-tiny;
   ux[1] = LISAWP_PI-tiny;
   ux[2] = LISAWP_PI-tiny;
*/

   lx[0] = -4.*LISAWP_PI;
   lx[1] = -4.*LISAWP_PI;
   //lx[0] = -2.;

   ux[0] = 4.*LISAWP_PI;
   ux[1] = 4.*LISAWP_PI;
   //ux[0] = 2.;
   

   TV.SetBoundaries(dim, lx, ux);

   TV.SetNumBins(NumBins);

   TV.InitializeBoxes();
 
   TV.smooth = false;

   int seed = 0;
   int numIter = 5;
   int Ncalls = 10;
   double *params;
   double val, stdev, chi2;
   
   double th = atan(-2./2.);
   std::cout << "stas " << th << std::endl; 
   std::cout << 0.5*(-2./cos(th) -4.0) << spr <<  0.5*(2.0/sin(th) -4.0) << std::endl;  
   double* co;
   co = new double[dim];
   th = 6.0;
   double psi = 6.0;

   co[0] = 4.*cos(th) + 2.*cos(psi)*cos(th);
   co[1] = 4.*sin(th) + 2.*cos(psi)*sin(th);
   std::cout << co[0] << spr << co[1] << spr <<  TV.CompFunct(co, params) \
	   << spr << 2.*sin(psi) << spr << cos(-1.14159) << spr << cos(th) << std::endl;
   delete [] co; 
   //exit(0);

   TV.Integrate(seed, numIter, Ncalls, params, val, stdev, chi2);

   double** coord = NULL;
   double* prob = NULL;

   int numBoxes = 1;
   for (int i=0; i<dim; i++)
	numBoxes *= NumBins[i];
   coord = new double*[numBoxes];
   prob = new double[numBoxes];

   for (int i =0; i< numBoxes; i++){
       coord[i] = new double[dim];
   }

   TV.GetProbability(coord, prob);

   std::ofstream fout1543("Data/MCprob.dat");
   for (int i = 0; i< numBoxes; i++){
       for (int j=0; j<dim; j++){	   
          fout1543 << coord[i][j] << spr;
       }
       fout1543 << prob[i] << spr << TV.CompFunct(coord[i], params) << std::endl;
   }
   fout1543.close();


   std::cout << "Int = " << val << std::endl;
   std::cout << "sigma = " << stdev << std::endl;
   std::cout << "chi^2 = " << chi2 << std::endl; 

   for (int i =0; i< numBoxes; i++){
       delete [] coord[i];
   }
   delete [] coord;
   delete [] prob;
 

   exit(0);
    double res, err;
     
    double xl[3] = { 0, 0, 0 };
    double xu[3] = { M_PI, M_PI, M_PI };
     
    const gsl_rng_type *T;
    gsl_rng *r;
     
    gsl_monte_function G = { &g, 3, 0 };
     
    //size_t calls = 500000;
    size_t calls = 625000;
     
    gsl_rng_env_setup ();
     
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);


    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
    

    double exact = 1.3932039296856768591842462603255; 
    std::cout << "Exact = " << exact << std::endl;

    GSL_vegas_integrate (&G, xl, xu, 3, calls, r, s,
                                    &res, &err);
    std::cout << "vegas warm-up  " << res << "   " <<  err << std::endl;
     
    printf ("converging...\n");
    do
      {
        GSL_vegas_integrate (&G, xl, xu, 3, calls, r, s,
                                        &res, &err);
         std::cout << "result = " << res << " sigma =  " << err <<\
                     "   chisq/dof = " <<  gsl_monte_vegas_chisq (s) << std::endl;
       }
    while (fabs (gsl_monte_vegas_chisq (s) - 0.3) > 0.1);
    std::cout << "result = " << res << " sigma =  " << err <<\
                     "   chisq/dof = " <<  gsl_monte_vegas_chisq (s) << std::endl;

    gsl_rng_free (r);

   gsl_monte_vegas_free (s);
   delete [] NumBins;
   delete [] lx;
   delete [] ux;

   return(0);

}
