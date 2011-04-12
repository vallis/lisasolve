/* constants from Constants.h */

// mathematical
const double pi = 3.141592653589793;
const double sq3 = 1.73205080757;

// physical
const double clight = 299792458;      // m/s
const double AU = 1.49597870660e11;   // m

// LISA
const double L = 5.0e9;             // armlength
const double fstar = 0.00954269032; // transfer frequency
const double fm = 3.168753575e-8;   // modulation frequency
const double kappa = 0.0;           // initial azimuthal position of the guiding center
const double lambda = 0.0;          // initial orientation of the LISA constellation
const double ec = 0.009648370435;   // eccentricity

const double Sps = 4.0e-22;         // photon shot noise power (Neil's value)
const double Sacc = 9.0e-30;        // acceleration noise power (Neil's value)

/* #include <complex>
   typedef std::complex<double> cdouble; */

extern "C" {
    #include <fftw3.h>
    
    /* double AEnoise(double f); */
}

double AEnoise(double f);

class FastResponse {
private:
    long N, M;
    double T, dt;
    
    double *u,*v,*k;            // Gravitational Wave basis vectors
    double *kdotx, **kdotr;     // Dot products
    double *xi, *fonfs;     // Distance, gravitational wave frequency & ratio of f and transfer frequency f*
    double **eplus, **ecross, **dplus, **dcross;    // Polarization basis tensors, convenient quantities    
    double *x, *y, *z;                              // Spacecraft position and separation vectors
    double **xv, **yv, **zv;
    double *r12, *r13, *r21, *r23, *r31, *r32;
    double **TR, **TI;                              // Time varying quantities (Re & Im) broken up into convenient segments
    double *data12, *data13, *data21, *data23, *data31, *data32;    // Fourier coefficients before FFT and after convolution:
                                                                    // Time series of slowly evolving terms at each vertex   
    double *a12, *a13, *a21, *a23, *a31, *a32;                      // Fourier coefficients of slowly evolving terms (numerical)
    double *b;                                                      // Fourier coefficients of rapidly evolving terms (analytical)
    double *an, *bn;                                                // MV: Fourier transforms for convolve_fft
    double *c12, *c13, *c21, *c23, *c31, *c32;                      // Fourier coefficients of entire response (convolution)

    double *ReA, *ImA, *ReB, *ImB, *ReC, *ImC;    

    double *X, *Y, *Z;
    
    fftw_complex *in, *out;
    fftw_plan plan_forward, plan_backward;

    void spacecraft(double t);
    
    void convolve(double *a, double *b, double *cn, int method);

    void XYZ(double f0, long q, double *XLS, double *XSL, double *YLS, double *YSL, double *ZLS, double *ZSL);

public:
    FastResponse(long Nreq,double Treq,double dtreq);
    ~FastResponse();
    
    void Response(double f0,double fdot,double theta,double phi,double A,double iota,double psi,double phio,
                  double *XLS,long XLSlen,double *XSL,long XSLlen,
                  double *YLS,long YLSlen,double *YSL,long YSLlen,
                  double *ZLS,long ZLSlen,double *ZSL,long ZSLlen,
                  int method);
};
