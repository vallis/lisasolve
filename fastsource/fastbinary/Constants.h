/****************************************************************************/
/*                                                                          */
/* TITLE: Preprocessor Macros For The Galactic Background                   */
/*                                                                          */
/* ABSTRACT: This header file contains a number of constants used           */
/* throughout the making of the galactic backgound realizations.            */
/*                                                                          */
/****************************************************************************/



          /* --------------  MATHEMATICAL CONSTANTS  -------------- */

 /* Set the value of pi */
#define pi 3.141592653589793

#define sq2pi 2.506628274631

 /* Square root of 3 */
#define sq3 1.73205080757
          /* ----------------  NATURAL CONSTANTS  ----------------- */

 /* Speed of light (m/s) */
#define clight 299792458.

 /* Number of seconds in a sidereal year */
#define year 3.15581498e7

 /* Newton's gravitational constant (mks) */
#define G 6.67259e-11

 /* Astronomical unit (meters) */
#define AU 1.49597870660e11

 /* Number of meters in a parsec */
#define pc 3.0856775807e16

 /* Number of meters in a kiloparsec */
#define kpc 3.0856775807e19

 /* Distance bewteen the Sun and the center of the galaxy (kpc) */
#define Rsun 8.5

 /* Mass of the Sun (kg) */
#define Msun 1.9889e30


          /* ----------------  DETECTOR CONSTANTS  ---------------- */

 /* Observation time (seconds) */
//#define T 3.15581498e7 //(1.0*year)
//#define T 6.31162996e7 //(2.0*year)
#define T 6.2914560e7
//#define T pow(2.0,22.0)*15.0  //~ 2 years (2^22 points at 15 s cadence)

 /* Number of data points */
#define NFFT 4194304

#define dt 15.0

 /* Initial azimuthal position of the guiding center */
#define kappa 0.0

 /* Initial orientation of the LISA constellation */
#define lambda 0.0

 /* Orbital radius of the guiding center */
#define Rgc (1.0*AU)

 /* Mean arm length of the LISA detector (meters) */
#define L 5.0e9

 /* Photon shot noise power */
#define Sps 4.0e-22
 
 /* Acceleration noise power */
#define Sacc 9.0e-30

 /* Transfer frequency */ 
#define fstar 0.00954269032

 /* LISA orbital eccentricity */
#define ec 0.009648370435

 /* LISA modulation frequency */
#define fm 3.168753575e-8

#define ec 0.009648370435

#define sq3 1.73205080757

#define GEOM 4.92569043916e-6

