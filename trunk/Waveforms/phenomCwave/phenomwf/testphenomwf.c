#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phenomwf.h"

int main( int argc, char *argv[] )
{

  if ( argc == 9 ) 
  {

    double eta = atof(argv[1]);
    double chi = atof(argv[2]);
    double Mtot = atof(argv[3]);
    double fMin = atof(argv[4]);
    double fMax = atof(argv[5]);
    double df = atof(argv[6]);
    double dL = atof(argv[7]);
    FILE *file = fopen( argv[8], "w" );

    int n,i;
    /* will start at f=0 and zero out below fMin */
    n = floor (fMax / df);

    double *wfreal=NULL;
    double *wfimag=NULL;
    double *f=NULL;

    wfreal = malloc(sizeof(double)*n);
    wfimag = malloc(sizeof(double)*n);
    f = malloc(sizeof(double)*n);

    /* set up frequencies */
    for ( i=0; i<n; i++ )
    {
      f[i]=i*df;
    }
    /* call waveform generation */
    phenomwf( wfreal , wfimag , f , n , df , fMin , fMax , eta , chi , Mtot , dL);

    /* write to file */
    if ( file == 0 )
    {
        fprintf(stdout, "Error: could not open file\n" );
        return(1);
    }
    else 
    {
    /* write  */
      fprintf(file, "# freq (Hz)\treal part\timaginary part\n");
      for ( i=0; i<n; i++)
      {
        fprintf(file, "%e\t%e\t%e\n",f[i],wfreal[i],wfimag[i]);
      }
      fclose( file );
    }

  /* free memory */
  free(wfreal);
  free(wfimag);
  free(f);

  }
  else
  {
    fprintf(stdout,"Error: wrong number of arguments.\nUsage: testphenomwf eta chi Mtot fMin fMax df dL outputfile\n");
    return(1);
  }

  return(0);
}
