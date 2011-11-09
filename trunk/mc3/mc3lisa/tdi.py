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
    
    def __mul__(self,other):
        if isinstance(other,TDIf):
            return TDIf(aet=(self.Af*other.Af,self.Ef*other.Ef,self.Tf*other.Tf))
        else:
            return TDIf(aet=(self.Af*other,self.Ef*other,self.Tf*other))
    
    def __rmul__(self,other):
        return TDIf(aet=(self.Af*other,self.Ef*other,self.Tf*other))
    
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
    
    def cprod(self,other):
        return (4.0 / self.Af.df) * ( N.sum(N.conj(self.Af) * other.Af / self.Sae)
                                    + N.sum(N.conj(self.Ef) * other.Ef / self.Sae)
                                    + N.sum(N.conj(self.Tf) * other.Tf / self.St) )
    
    def dotprod(self,other):
        return (4.0 / self.Af.df) * N.real( N.sum(N.conj(self.Af) * other.Af / self.Sae) +
                                            N.sum(N.conj(self.Ef) * other.Ef / self.Sae) +
                                            N.sum(N.conj(self.Tf) * other.Tf / self.St ) )
    
    def logL(self,other):
        return -0.5 * (4.0 / self.Af.df) * ( N.sum(N.abs(self.Af.rsub(other.Af))**2 / self.Sae) +
                                             N.sum(N.abs(self.Ef.rsub(other.Ef))**2 / self.Sae) +
                                             N.sum(N.abs(self.Tf.rsub(other.Tf))**2 / self.St ) )
    


def SNR(center,gettdi):
    return math.sqrt(gettdi(center).normsq())


def fisher(center,scale,gettdi):    
    x0, dx = N.asarray(center), N.asarray(scale)
    d = len(x0)
    
    # is there a better numpy builtin for the *vector* delta_ij?
    delta = lambda i: N.identity(d)[i,:]
    
    derivs = [ ( gettdi(x0 + delta(i)*dx[i]) - gettdi(x0 - delta(i)*dx[i]) ) / (2.0*dx[i]) for i in range(d) ]
    
    prods = N.zeros((d,d),'d')
    for i in range(d):
        for j in range(i,d):
            prods[i,j] = prods[j,i] = derivs[i].dotprod(derivs[j])
    
    return prods



class model(object):
    c = 299792458
    
    defaultnoise = 'lisareq'
    noisemodel = defaultnoise
    
    defaultL = 16.6782
    lisaL = defaultL
    
    defaultD = 0.4
    lisaD = defaultD
    
    defaultP = 1.0
    lisaP = defaultP
    
    defaultWD = None
    lisaWD = defaultWD
    
    @staticmethod
    def setL(Lm):
        # set L in meters (same as in fastsource/fastbinary)
        model.lisaL = Lm / model.c
    
    @staticmethod
    def setmodel(mymodel,myarm=None):
        model.noisemodel = model.defaultnoise
        model.lisaL  = model.defaultL
        model.lisaD  = model.defaultD
        model.lisaP  = model.defaultP
        model.lisaWD = model.defaultWD
        
        if myarm != None:
            model.setL(myarm) # will be overridden by models that assume arm
        
        if mymodel in ['lisa-classic','default']:
            pass
        elif mymodel == 'CLISA1_P005c_LPF':
            model.noisemodel = 'newlpf'
            model.lisaL  = 1e9 / model.c
            model.lisaP  = 0.05
        elif mymodel == '10LISA1_P2_DRS':
            model.noisemodel = 'newdrs-wrong'
            model.lisaL  = 1e9 / model.c
            model.lisaP  = 2
        elif mymodel == '10LISA1_P07_D25_DRS_4L':
            model.noisemodel = 'newdrs'
            model.lisaL  = 1e9 / model.c
            model.lisaP  = 0.7
            model.lisaD  = 0.25
        elif mymodel == '10LISA1_P2_D25_DRS_4L':
            model.noisemodel = 'newdrs'
            model.lisaL  = 1e9 / model.c
            model.lisaP  = 2
            model.lisaD  = 0.25
        elif mymodel == '10LISA1_P07_D25_RDRS_4L':
            model.noisemodel = 'reddrs'
            model.lisaL  = 1e9 / model.c
            model.lisaP  = 0.7
            model.lisaD  = 0.25
        elif mymodel in ['mldc','mldc-nominal','lisareq','toy','newlpf','newdrs','reddrs','lpf','wind','ax50']:
            model.noisemodel = mymodel
        else:
            raise NotImplementedError(mymodel)
    


# currently only the lpf model allows for the part of optical noise that does not change with armlength
def lisanoises(f,noisemodel=None):
    if noisemodel == None: noisemodel = model.noisemodel
    
    if   noisemodel == 'mldc':
        Spm = 2.5e-48 * (1.0 + (f/1.0e-4)**-2) * f**(-2)
        Sop = 1.8e-37 * (model.lisaL/model.defaultL)**2 * f**2
    elif noisemodel == 'mldc-nominal':
        Spm = 2.53654e-48 * (1.0 + (f/1.0e-4)**-2) * f**(-2)        # 3e-15 m/s^2/sqrt(Hz), divided by (2 pi c) and squared
        Sop = 1.75703e-37 * (model.lisaL/model.defaultL)**2 * f**2  # 20 pm/sqrt(Hz), mult. by (2 pi / c) and squared; all optical noise scales as shot noise
    elif noisemodel == 'lisareq':
        Spm = 2.53654e-48 * (1.0 + (f/1.0e-4)**-1) * (1.0 + (f/0.008)**4) * f**(-2)
        Sop = 1.42319e-37 * (model.lisaL/model.defaultL)**2 * (1.0 + (f/0.002)**-4) * f**2
    elif noisemodel == 'toy':
        Spm = 2.53654e-48 * f**(-2)                                 # 3e-15 m/s^2/sqrt(Hz), no reddening
        
        Sops = 1.1245e-37 * (model.lisaL/model.defaultL)**2 * (model.defaultD/model.lisaD)**4 * (model.defaultP/model.lisaP)   # 16 pm/sqrt(Hz)
        Sopo = 6.3253e-38                                                                                                      # 12 pm/sqrt(Hz)
        Sop = (Sops + Sopo) * f**2
    elif noisemodel == 'newlpf':  # lisalight, to be used with lisaL = 1Gm, lisaP = 0.05
        Spm = 8.17047e-48 * (1.0 + (f/1.8e-4)**-1)**2 * f**(-2)     # 5.3e-15 m/s^2/sqrt(Hz)
        
        Sops = 6.15e-38 * (model.lisaL/model.defaultL)**2 * (model.defaultD/model.lisaD)**4 * (model.defaultP/model.lisaP)      # 11.83 pm/sqrt(Hz)
        Sopo = 2.81e-38                                                                                                         # 8 pm/sqrt(Hz)
        Sop = (Sops + Sopo) * f**2
    elif noisemodel == 'newdrs-wrong':  # lisalight, to be used with lisaL = 1Gm, lisaP = 2
        Spm = 6.00314e-48 * f**(-2)                                 # 4.6e-15 m/s^2/sqrt(Hz)
        
        Sops = 3.07e-38 * (model.lisaL/model.defaultL)**2 * (model.defaultD/model.lisaD)**4 * (model.defaultP/model.lisaP)      # 8.36 pm/sqrt(Hz) - was scaled wrong with P
        Sopo = 2.81e-38                                                                                                         # 8 pm/sqrt(Hz)
        Sop = (Sops + Sopo) * f**2
    elif noisemodel == 'newdrs':  # lisalight, to be used with lisaL = 1Gm, lisaP = 2
        Spm = 6.00314e-48 * f**(-2)                                 # 4.6e-15 m/s^2/sqrt(Hz)
        
        Sops = 6.15e-38 * (model.lisaL/model.defaultL)**2 * (model.defaultD/model.lisaD)**4 * (model.defaultP/model.lisaP)      # 11.83 pm/sqrt(Hz)
        Sopo = 2.81e-38                                                                                                         # 8 pm/sqrt(Hz)
        Sop = (Sops + Sopo) * f**2        
    elif noisemodel == 'reddrs':  # used for lisalight C6
        Spm = 6.0e-48 * (1 + (1e-4/f)) * f**(-2)    # 4.61e-15 m/s^2/sqrt(Hz)
        
        Sops = 6.17e-38 * (model.lisaL/model.defaultL)**2 * (model.defaultD/model.lisaD)**4 * (model.defaultP/model.lisaP)
        Sopo = 2.76e-38
        Sop = (Sops + Sopo) * f**2        
    elif noisemodel == 'lpf':
        # LPF CBE curve from ? via Oliver; the coefficient in front is [10^-14.09 * c / (2 pi)]^2
        Spm = 1.86208e-47 * (1.0 + (f/10**-3.58822)**-1.79173) * (1.0 + (f/10**-2.21652)**3.74838) * f**(-2)
        
        # see LISA-variant Mathematica notebook; include only shot-noise correction due to armlength
        # constants are 7.7 and 5.15 pm (formerly 9.25) squared and multiplied by (2 pi / c)^2
        # notice the standard curve includes 18 pm of which 35% is margin
        Sop = (1.16502e-38 + 2.60435e-38*(model.lisaL/model.defaultL)**2) * f**2    
    elif noisemodel == 'wind':
        Spm = 1.76e-50 * f**-0.75 * f**(-2)
        Sop = 1.42319e-37 * (model.lisaL/model.defaultL)**2 * (1.0 + (f/0.002)**-4) * f**2
    elif noisemodel == 'ax50':
        Spm = 50 * 2.53654e-48 * (1.0 + (f/1.0e-4)**-1) * (1.0 + (f/0.008)**4) * f**(-2)
        Sop = 1.42319e-37 * (model.lisaL/model.defaultL)**2 * (1.0 + (f/0.002)**-4) * f**2
    else:
        raise NotImplementedError(noisemodel)
    
    return Spm, Sop


def phinneyswitch(Sinst,Sgwdb,switch):
    return N.minimum(Sinst*switch,Sinst + Sgwdb)


class phinneybackground(object):
    def __init__(self,Sh=1.4e-44,dNdf=2e-3,koverT=1.5,Sh_exp=-7.0/3.0,dNdf_exp=-11.0/3.0,dNdf_func=N.exp,dNdf_switch=phinneyswitch):
        self.Sh,   self.Sh_exp   = Sh,   Sh_exp
        self.dNdf, self.dNdf_exp = dNdf, dNdf_exp
        self.koverT = koverT / (365.25 * 24 * 3600)
        self.dNdf_func, self.dNdf_switch = dNdf_func, dNdf_switch
    
    def __call__(self,f,Sinst = None):
        Sgwdb = self.Sh * f**self.Sh_exp
        dNdf  = self.dNdf * f**self.dNdf_exp
        
        if Sinst != None:
            return self.dNdf_switch(Sinst,Sgwdb,self.dNdf_func(self.koverT*dNdf))
        else:
            return Sgwdb
    


# note: not updated for shorter armlengths, LPF noise, better optical noise model
def lisanoise(f,noisemodel=None,includewd=None):
    if noisemodel == None: noisemodel = model.noisemodel
    if includewd == None: includewd = model.lisaWD
        
    if noisemodel == 'cutler':
        # compare to Eq. (25) of Barack and Cutler Phys.Rev. D 70, 122002 (2004)
        # their Sh, defined by <n n> = 3/40 Sh, is 6.12e-51 f**-4,
        # so the <n n> noise is enhanced (as "seen" by signals) because of signal averaging
        #
        # that's the same as defining Sh by <n n> = 1/2 Sh as 9.18e-52 f**-4 (as we do)
        # and then enhancing it by 20/3 because of signal averaging
        #
        # either way the S_h used in averaged SNR expressions is 6.12e-51
        # (and I hope I never have to think about this again)
        
        Sh = (20.0/3.0)*(9.18e-52 * f**-4 + 1.59e-41 + 9.18e-38 * f**2)
        
        if includewd == True:
            pb = phinneybackground()
            return pb(f,Sh)
        elif includewd == None:
            return Sh
        else:
            raise NotImplementedError
    else:
        optscale = (model.lisaL/model.defaultL)**2 * (model.defaultD/model.lisaD)**4 * (model.defaultP/model.lisaP)
        
        if noisemodel == 'lisareq':
            Sa = 3e-15 * N.sqrt(1.0 + (f/1.0e-4)**-1) * N.sqrt(1.0 + (f/0.008)**4)
            So = 18e-12 * optscale * N.sqrt(1 + (f/0.002)**-4)
        elif noisemodel == 'lpf':
            Sa = 10**-14.09 * N.sqrt((1.0 + (f/10**-3.58822)**-1.79173) * (1.0 + (f/10**-2.21652)**3.74838))
            So = N.sqrt((7.7e-12)**2  * optscale + (5.15e-12)**2)
        elif noisemodel == 'toy':
            Sa = 3e-15
            So = N.sqrt((1.6e-11)**2  * optscale + (1.2e-11)**2)
        elif noisemodel == 'newlpf':
            Sa = 5.3e-15 * (1.0 + (f/1.8e-4)**-1)
            So = N.sqrt((1.18e-11)**2 * optscale + (8.0e-12)**2)
        elif noisemodel == 'newdrs-wrong':
            Sa = 4.6e-15
            So = N.sqrt((8.36e-12)**2 * optscale + (8.0e-12)**2)
        elif noisemodel == 'newdrs':
            Sa = 4.6e-15
            So = N.sqrt((1.18e-11)**2 * optscale + (8.0e-12)**2)
        elif noisemodel == 'wind':
            Sa = 2.5e-16 * f**-0.75
            So = 18e-12 * optscale * N.sqrt(1 + (f/0.002)**-4)
        elif noisemodel == 'ax50':
            Sa = 50 * 3e-15 * N.sqrt(1.0 + (f/1.0e-4)**-1) * N.sqrt(1.0 + (f/0.008)**4)
            So = 18e-12 * optscale * N.sqrt(1 + (f/0.002)**-4)        
        else:
            raise NotImplementedError(noisemodel)
        
        Sac = Sa * 2.0 / (2.0 * math.pi * f)**2
        
        c = 299792458; L = model.lisaL * c
        
        # how to fit WD noise in here? Let's reason for X first...
        # TDI-X WD noise could be added to Sop (as defined in lisanoises) after dividing by 16 sin(x)^2
        # now Sop = So^2 * (2*pi/c)^2 f^2; hence we can add Swd/(16 sin(x)^2) / (2*pi*f/c)^2 to So
        # note that confusion noise is L dependent because of subtraction... but it probably still
        # makes sense to factor the L and add L^2 * Swd / (16 sin(x)^2 x^2/c^2) with x = 2 pi L f
        
        # transfer function
        ft = 0.5 / model.lisaL
        T2 = 1.0 + (f/(0.41 * ft))**2
                
        if includewd == None:
            Swd = 0
        elif includewd == 'cutler':
            pb = phinneybackground()
            return pb(f,(20.0/3.0) * T2 * (Sac**2 + So**2) / L**2)
        elif isinstance(includewd,phinneybackground):
            return includewd(f,(20.0/3.0) * T2 * (Sac**2 + So**2) / L**2)
        else:
            x = 2.0 * math.pi * model.lisaL * f            
            Swd = makewdnoise(f,includewd,obs='X') * L**2 / (16.0 * N.sin(x)**2 * x**2)
        
        return (20.0/3.0) * T2 * (Sac**2 + So**2 + Swd) / L**2


def simplesnr(f,h,i=None,years=1,noisemodel=None,includewd=None):
    if i == None:
        h0 = h * math.sqrt(16.0/5.0)    # rms average over inclinations
    else:
        h0 = h * math.sqrt((1 + math.cos(i)**2)**2 + (2*math.cos(i))**2)        
    
    return h0 * math.sqrt(years * 365.25*24*3600) / math.sqrt(lisanoise(f,noisemodel,includewd))


wdnoise = {}

# fit between 1e-4 and 5e-3 for X, between 1e-4 and 4e-4 for AET; all SNR = 5
wdnoise['tau2']   = (('rat42',[-1.2503, -13.3508, -94.1852, -296.6416, -313.8596, 4.9418, 6.1323]),
                     ('rat42',[-1.2599, -13.8309, -97.7703, -311.5419, -336.4092, 5.0691, 6.4637]))
wdnoise['opt']    = (('rat42',[-1.0865, -11.2113, -83.9764, -271.5378, -287.9153, 4.8456, 5.8931]),
                     ('rat42',[-1.0781, -11.3477  -85.3638, -279.6701, -301.9440, 4.9496, 6.1504]))
wdnoise['pess']   = (('rat42',[-1.2649, -13.5895, -95.5196, -301.0872, -319.7566, 4.9740, 6.2117]),
                     ('rat42',[-1.2813, -14.1556, -99.5091, -316.7877, -342.7881, 5.1004, 6.5392]))
wdnoise['hybrid'] = (('poly4',[-2.4460, -33.4121,-171.5341, -390.7209, -373.5341]),
                     ('poly4',[-2.7569, -38.0938,-197.8030, -455.9119, -433.8260]))

def makewdnoise(f,wdstyle,obs='X'):
    if wdstyle == 'mldc':
        x = 2.0 * math.pi * model.lisaL * f
        t = 4 * x**2 * N.sin(x)**2 * (1.0 if obs == 'X' else 1.5)
        
        return t * ( N.piecewise(f,(f >= 1.0e-4  ) & (f < 1.0e-3  ),[lambda f: 10**-44.62 * f**-2.3, 0]) + \
                     N.piecewise(f,(f >= 1.0e-3  ) & (f < 10**-2.7),[lambda f: 10**-50.92 * f**-4.4, 0]) + \
                     N.piecewise(f,(f >= 10**-2.7) & (f < 10**-2.4),[lambda f: 10**-62.8  * f**-8.8, 0]) + \
                     N.piecewise(f,(f >= 10**-2.4) & (f < 10**-2.0),[lambda f: 10**-89.68 * f**-20.0,0])     )
    elif wdstyle in wdnoise:
        mod, p = wdnoise[wdstyle]
        p = p[0] if obs == 'X' else p[1] # assume AE if not X
        y = N.log10(f)
    
        if mod == 'rat42':
            return 10.0**( (p[0]*y**4+p[1]*y**3+p[2]*y**2+p[3]*y+p[4])/(y**2+p[5]*y+p[6]) )
        elif mod == 'poly4':
            return 10.0**( p[0]*y**4+p[1]*y**3+p[2]*y**2+p[3]*y+p[4] )
        else:
            raise NotImplementedError
    else:
        if '.txt' in wdstyle:
            conf = N.loadtxt(wdstyle)
            conf[N.isnan(conf[:,1]),1] = 0
        
            return N.interp(f,conf[:,0],conf[:,1])
        else:
            raise NotImplementedError


# from lisatools makeTDIsignal-synthlisa2.py
def noisepsd_X(frequencydata,includewd=None):
    if includewd == None: includewd = model.lisaWD
    
    f = frequencydata.f if isinstance(frequencydata,FrequencyArray) else frequencydata
    x = 2.0 * math.pi * model.lisaL * f
    
    Spm, Sop = lisanoises(f)
    
    Sx = 16.0 * N.sin(x)**2 * (2.0 * (1.0 + N.cos(x)**2) * Spm + Sop)
    # Sxy = -4.0 * N.sin(2*x) * N.sin(x) * (Sop + 4.0*Spm)
    # Sa = Sx - Sxy
    
    if includewd != None:
        Sx += makewdnoise(f,includewd,'X')
    
    if isinstance(frequencydata,FrequencyArray):
        return FrequencyArray(Sx,kmin=frequencydata.kmin,df=frequencydata.df)
    else:
        return Sx


def noisepsd_AE(frequencydata,includewd=None):
    if includewd == None: includewd = model.lisaWD
    
    f = frequencydata.f if isinstance(frequencydata,FrequencyArray) else f
    x = 2.0 * math.pi * model.lisaL * f
    
    Spm, Sop = lisanoises(f)
    
    Sa = 8.0 * N.sin(x)**2 * (2.0 * Spm * (3.0 + 2.0*N.cos(x) + N.cos(2*x)) +
                              Sop * (2.0 + N.cos(x)))
    
    if includewd != None:
        Sa += makewdnoise(f,includewd,'AE')
    
    if isinstance(frequencydata,FrequencyArray):
        return FrequencyArray(Sa,kmin=frequencydata.kmin,df=frequencydata.df)
    else:
        return Sa

# TO DO: currently not including WD background here... probably OK
def noisepsd_T(frequencydata):
    f = frequencydata.f if isinstance(frequencydata,FrequencyArray) else f
    x = 2.0 * math.pi * model.lisaL * f
    
    Spm, Sop = lisanoises(f)
    
    St = 16.0 * Sop * (1.0 - N.cos(x)) * N.sin(x)**2 + 128.0 * Spm * N.sin(x)**2 * N.sin(0.5*x)**4
    
    if isinstance(frequencydata,FrequencyArray):
        return FrequencyArray(St,kmin=frequencydata.kmin,df=frequencydata.df)
    else:
        return St


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
