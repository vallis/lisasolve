import math
import numpy as N

from FrequencyArray import FrequencyArray

def AET(X,Y,Z):
    return (Z - X)/math.sqrt(2.0), (X - 2.0*Y + Z)/math.sqrt(6.0), (X + Y + Z)/math.sqrt(3.0)


class TDIf(object):
    """TDI triple-observable object. Can be initialized with X,Y,Z or A,E,T, but
       will keep A,E,T internally."""
    
    def __init__(self,aet=None,xyz=None):
        if aet != None:
            self.Af, self.Ef, self.Tf = aet
        elif xyz != None:
            self.Af, self.Ef, self.Tf = AET(*xyz)
            self.Xf = xyz[0]                        # keep also X (temporary?)
        
        self._Sae, self._St, self._Sx = None, None, None
    
    
    def pad(self,leftpad=1,rightpad=1):
        self.Af = self.Af.pad(leftpad,rightpad)
        self.Ef = self.Ef.pad(leftpad,rightpad)
        self.Tf = self.Tf.pad(leftpad,rightpad)
        
        self.Xf = self.Xf.pad(leftpad,rightpad)
        
        self._Sae, self._St, self._Sx = None, None, None
        
        return self
    
    
    # TO DO: memoize the properties?    
    @property
    def Sae(self):
        if self._Sae == None: self._Sae = noisepsd_AE(self.Af)
        return self._Sae
    
    @property
    def St(self):
        if self._St  == None: self._St  = noisepsd_T(self.Tf)
        return self._St
    
    @property
    def Sx(self):
        if self._Sx  == None: self._Sx  = noisepsd_X(self.Xf)
        return self._Sx
    
    @property
    def kmin(self):
        return self.Af.kmin
    
    def __len__(self):
        return len(self.Af)
    
    def __add__(self,other):
        return TDIf(aet=(self.Af + other.Af,self.Ef + other.Ef,self.Tf + other.Tf))
    
    def __sub__(self,other):
        return TDIf(aet=(self.Af - other.Af,self.Ef - other.Ef,self.Tf - other.Tf))
    
    def __div__(self,other):
        if isinstance(other,TDIf):
            return TDIf(aet=(self.Af/other.Af,self.Ef/other.Ef,self.Tf/other.Tf))
        else:
            return TDIf(aet=(self.Af/other,self.Ef/other,self.Tf/other))
    
    def __iadd__(self,other):
        self.Af += other.Af; self.Ef += other.Ef; self.Tf += other.Tf
    
    def __isub__(self,other):
        self.Af -= other.Af; self.Ef -= other.Ef; self.Tf -= other.Tf
    
    def normsq(self,noisepsd = None,extranoise = [0,0,0]):
        if noisepsd == None:
            return (4.0 / self.Af.df) * ( N.sum(N.abs(self.Af)**2 / (self.Sae + extranoise[0])) + 
                                          N.sum(N.abs(self.Ef)**2 / (self.Sae + extranoise[1])) +
                                          N.sum(N.abs(self.Tf)**2 / (self.St  + extranoise[2])) )
        else:
            return (4.0 / self.Af.df) * ( N.sum(N.abs(self.Af)**2 / noisepsd[0]) + 
                                          N.sum(N.abs(self.Ef)**2 / noisepsd[1]) +
                                          N.sum(N.abs(self.Tf)**2 / noisepsd[2]) )
    
    def normsqx(self,noisepsd = None):
        if noisepsd == None:
            return (4.0 / self.Xf.df) * ( N.sum(N.abs(self.Xf)**2 / (self.Sx)) )
        else:
            return (4.0 / self.Xf.df) * ( N.sum(N.abs(self.Xf)**2 / noisepsd) )
    
    def dotprod(self,other):
        return (4.0 / self.Af.df) * N.real( N.sum(N.conj(self.Af) * other.Af / self.Sae) +
                                            N.sum(N.conj(self.Ef) * other.Ef / self.Sae) +
                                            N.sum(N.conj(self.Tf) * other.Tf / self.St ) )    
    
    def logL(self,other):
        return -0.5 * (4.0 / self.Af.df) * ( N.sum(N.abs(self.Af.rsub(other.Af))**2 / self.Sae) +
                                             N.sum(N.abs(self.Ef.rsub(other.Ef))**2 / self.Sae) +
                                             N.sum(N.abs(self.Tf.rsub(other.Tf))**2 / self.St ) )
    



def SNR(state):
    signal = state.model.datamodel(state)
    
    return math.sqrt(signal.normsq())

    
def fisher(center,scale):
    x0 = center.asarray(); d = len(x0)
    dx = center.fromarray(scale)    
    
    delta = lambda i: N.identity(d)[i,:]    # is there a good numpy builtin for the *vector* delta_ij?
    gettdi = lambda s: center.model.datamodel(center.model.state(s))
    
    derivs = [ ( gettdi(x0 + delta(i)*dx[i]) - gettdi(x0 - delta(i)*dx[i]) ) / (2.0*dx[i]) for i in range(d) ]
        
    prods = N.zeros((d,d),'d')
    for i in range(d):
        for j in range(i,d):
            prods[i,j] = prods[j,i] = derivs[i].dotprod(derivs[j])
    
    return prods


defaultmodel = 'lisareq'

def lisanoises(f,noisemodel=None):
    if noisemodel == None: noisemodel = defaultmodel
    
    if   noisemodel == 'mldc':
        Spm = 2.5e-48 * (1.0 + (f/1.0e-4)**-2) * f**(-2)
        Sop = 1.8e-37 * f**2
    elif noisemodel == 'mldc-nominal':
        Spm = 2.53654e-48 * (1.0 + (f/1.0e-4)**-2) * f**(-2)    # 3e-15 m/s^2/sqrt(Hz)
        Sop = 1.75703e-37 * f**2                                # 2e-11 m/sqrt(Hz)
    elif noisemodel == 'lisareq':
        Spm = 2.53654e-48 * (1.0 + (f/1.0e-4)**-1) * (1.0 + (f/0.008)**4) * f**(-2)
        Sop = 1.42319e-37 * (1.0 + (f/0.002)**-4) * f**2
    
    return Spm, Sop

def lisanoise(f,noisemodel=None):
    if noisemodel == 'cutler':
        return (20.0/3.0)*(9.18e-52 * f**-4 + 1.59e-41 + 9.18e-38 * f**2)
    
    # otherwise use the LISA requirements model
    c = 299792458; L = 16.6782 * c  # note: it's really 16.67 in the requirements
    
    Sa  = 3e-15 * N.sqrt(1.0 + (f/1.0e-4)**-1) * N.sqrt(1.0 + (f/0.008)**4)
    Sac = Sa * 2.0 / (2.0 * math.pi * f)**2
    
    So  = 18e-12 * N.sqrt(1 + (f/0.002)**-4)
    
    ft = c / (2.0 * L)
    T2 = 1.0 + (f/(0.41 * ft))**2
    
    return (20.0/3.0) * T2 * (Sac**2 + So**2) / L**2

def simplesnr(f,h,i=None,years=1,noisemodel=None):
    if i == None:
        h0 = h * math.sqrt(16.0/5.0)    # rms average over inclinations
    else:
        h0 = h * math.sqrt((1 + math.cos(i)**2)**2 + (2*math.cos(i))**2)        
    
    return h0 * math.sqrt(years * 365.25*24*3600) / math.sqrt(lisanoise(f,noisemodel))


includewd = False

# from lisatools makeTDIsignal-synthlisa2.py
def noisepsd_X(frequencydata):
    f = frequencydata.f
    
    L  = 16.6782
    x = 2.0 * math.pi * L * f
    
    Spm, Sop = lisanoises(f)
    
    Sx = 16.0 * N.sin(x)**2 * (2.0 * (1.0 + N.cos(x)**2) * Spm + Sop)
    # Sxy = -4.0 * N.sin(2*x) * N.sin(x) * (Sop + 4.0*Spm)
    # Sa = Sx - Sxy
        
    if includewd:
        Sx += (2.0 * L)**2 * (2*math.pi*f)**2 * 4.0 * N.sin(x)**2 * (
                N.piecewise(f,(f >= 1.0e-4  ) & (f < 1.0e-3  ),[lambda f: 10**-44.62 * f**-2.3, 0]) + \
                N.piecewise(f,(f >= 1.0e-3  ) & (f < 10**-2.7),[lambda f: 10**-50.92 * f**-4.4, 0]) + \
                N.piecewise(f,(f >= 10**-2.7) & (f < 10**-2.4),[lambda f: 10**-62.8  * f**-8.8, 0]) + \
                N.piecewise(f,(f >= 10**-2.4) & (f < 10**-2.0),[lambda f: 10**-89.68 * f**-20.0,0])     )
    
    return FrequencyArray(Sx,kmin=frequencydata.kmin,df=frequencydata.df)

def noisepsd_AE(frequencydata):    
    f = frequencydata.f
    
    L  = 16.6782
    x = 2.0 * math.pi * L * f
    
    Spm, Sop = lisanoises(f)
    
    Sa = 8.0 * N.sin(x)**2 * (2.0 * Spm * (3.0 + 2.0*N.cos(x) + N.cos(2*x)) +
                              Sop * (2.0 + N.cos(x)))
    
    if includewd:
        Swd = (2.0 * L)**2 * (2*math.pi*f)**2 * 4.0 * N.sin(x)**2 * (
                N.piecewise(f,(f >= 1.0e-4  ) & (f < 1.0e-3  ),[lambda f: 10**-44.62 * f**-2.3, 0]) + \
                N.piecewise(f,(f >= 1.0e-3  ) & (f < 10**-2.7),[lambda f: 10**-50.92 * f**-4.4, 0]) + \
                N.piecewise(f,(f >= 10**-2.7) & (f < 10**-2.4),[lambda f: 10**-62.8  * f**-8.8, 0]) + \
                N.piecewise(f,(f >= 10**-2.4) & (f < 10**-2.0),[lambda f: 10**-89.68 * f**-20.0,0])     )
    
        # from lisatools makeTDIsignal-synthlisa2.py, adjusted for our definition of A
        Sa += 1.5 * Swd
    
    return FrequencyArray(Sa,kmin=frequencydata.kmin,df=frequencydata.df)

# TO DO: currently not including WD background here... probably OK
def noisepsd_T(frequencydata):
    f = frequencydata.f
    
    L  = 16.6782
    x = 2.0 * math.pi * L * f
    
    Spm, Sop = lisanoises(f)
    
    St = 16.0 * Sop * (1.0 - N.cos(x)) * N.sin(x)**2 + 128.0 * Spm * N.sin(x)**2 * N.sin(0.5*x)**4
    
    return FrequencyArray(St,kmin=frequencydata.kmin,df=frequencydata.df)


# === note on FFT normalization: it's always fun, isn't it? ===
#
# numpy.fft.fft and numpy.fft.ifft are inverses of each other, but
# numpy.fft.fft(X)/sqrt(N) has the same sum-squared as X
# numpy.fft.rfft is a subset of numpy.fft.fft, and numpy.fft.irfft is its inverse
#
# now Parseval says
#   int |h(f)|^2 df = int h(t)^2 dt
# and discretizing
#   sum |h(f_i)|^2 / T = sum h(t_i)^2 T / N
# hence
#   sum |h(f_i)|^2 (N/T^2) = sum |fft(X)_i|^2 / N
# whence
#   h(f_i) = fft(X)_i / (N df) -> int |h(f)|^2 df = sum |h(f_i)|^2 df
#                                                 = sum |fft(X)_i|^2 / (N^2 df)
#
# we're using a compact representation of the frequency-domain signal
# with no negative frequencies, and only 2*N - 1 components; if we were really
# integrating over all frequencies, we should take half the DC and Nyquist contributions,
# but that won't be the case here (and the trapezoidal rule would take care of it anyway)
#
# last question is what is returned by FastBinary.cpp. Since
# irfft(ret(f_i)) * N seems to be the correct time-domain signal,
# then ret(f_i) must be h(f_i) * df, and therefore the SNR integral is
#
# 4 Re int |h(f)|^2 / Sn(f) df = 4 sum |ret(f_i) / df|^2 / Sn(f_i) df = (4/df) sum |ret(f_i)|^2 / Sn(f_i) 
#
# now how's h(f_i) related to S_h(f)? I have
#
# int_0^\infty S_h(f) = (1/T) \int_0^T |h(t)|^2 dt = (2/T) \int_0^\infty |h(f)|^2 df
# so S_h(f_i) = (2/T) |h(f_i)|^2 = (2/T) |ret(f_i)|^2 / df^2 = (2/df) |ret(f_i)^2|
