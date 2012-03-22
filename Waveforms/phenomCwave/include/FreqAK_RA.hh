#ifndef FREQAKRAHH
#define FREQAKRAHH

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "Macros.hh"
#include "Constants.hh"
#include "Matrix.hh"
#include "EMRItemplate.hh"
#include "OrbitalMotion.hh"
#include <vector>


/*** Class FreqAK_RA
 *
 * Generator of analytic kludge waveforms described in 
 * Barack & Cutler gr-qc/0310125. it generates waveform in the freq. domain using rigid adiabatic approximation
 *
 * @author  Jon Gair, Stas Babak 2010
 */

namespace LISAWP{
  
class FreqAK_RA{
   
  public:
   
   /*** Constructor
    * @param timestep cadence for the time domain 
    * @param maxDur maximum duration of the observation, stop waveform if plunge didn't occure
    * @param dtPhase cadence (internal) for evolving orbit and LISA motion, (usually 2048 sec).
    */
    FreqAK_RA(double timestep, double maxDur, double dtPhase, double L); 
    
    /***  destructor */
    ~FreqAK_RA();
    
    /*** Evolves orbit and responce (RA) with large time step dtPhase
     * Note! that integration starting from any point is not yet implemented it is hardcoded
     * that template starts at 0 at up to maxDuration long.... 
     * @param  S  EMRItemplate object with parameters etc.
     * @param t0  initial moment of observations
     * @param dur duration of the template
     * @param tStart moment of time at which parameters of the template are specified
     */
     
    void PhaseEv_RA(EMRItemplate& S, double t0, double dur, double tStart);
   
    /** The same as above but for long-wavelength limit response */ 
    void PhaseEv_LW(EMRItemplate& S, double t0, double dur, double tStart);
    
    /** Constructing waveform in time domain using LW response */
    void ConstructWaveTime_LW(EMRItemplate& S, Matrix<double>& X, Matrix<double>& Y, Matrix<double>& Z);
    
    /** Constructing waveform in time domain using RA response */
    void ConstructWaveTime_RA(EMRItemplate& S, Matrix<double>& X, Matrix<double>& Y, Matrix<double>& Z);
     
    /**   Constructing waveform in time domain using LW response */
    void ConstructWaveTime_LW(EMRItemplate& S, Matrix<double>& X, Matrix<double>& Z_Y);

   /** Constructing waveform in freq domain using RA response */
    void ConstructWaveFreq_RA(EMRItemplate& S, Matrix<double>& Xf, Matrix<double>& Yf, Matrix<double>& Zf);

   /** Constructing waveform in freq domain using LW response */
    void ConstructWaveFreq_LW(EMRItemplate& S, Matrix<double>& Xf, Matrix<double>& Z_Yf);
    
  private:
     
     double* tm_ph;
     double* phiev;
     double* alpev;
     double* gamev;
     double* eccev;
     double* nuev;
     double* gamdotev;
     double* alpdotev;
     double* phiddotev;
     double* gamddotev;
     double* alpddotev;   
     double* Ampev;
     
     double* FplusI;
     double* FplusII;
     double* FcrosI;
     double* FcrosII;
     
     std::complex<double>** Xp;
     std::complex<double>** Xc;
     std::complex<double>** Yp;
     std::complex<double>** Yc;
     std::complex<double>** Zp;
     std::complex<double>** Zc;
     
     int Nph;     // number of points for phase and response evaluation
     int Nps;    // number of points for waveform evaluation
   
     double dt_w;
     double dt_ph;
     double df;    // freq resolution df = 1/Tobs
     double Tobs;  // duration of observation in sec
     int imax;  // maximum index of the orbital evolution
     double arm; // armlength
     double year;
     
     double ArcTan(double up, double down);
     
     double Sinc(double y);
     
     
};
} // end of the namespace

#endif