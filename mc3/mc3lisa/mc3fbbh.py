from __future__ import division
import math
import numpy as N
import FrequencyArray

from functools import partial

class memoize(object):
    """cache the return value of an instance method -- DOES NOT WORK WITH NUMPY ARRAY ARGUMENTS
    [based on http://code.activestate.com/recipes/577452-a-memoize-decorator-for-instance-methods]"""
    
    def __init__(self,func):
        self.func = func
    
    def __get__(self,obj,objtype=None):
        if obj is None:
            return self.func            # if called on the class, return only the function
        else:
            return partial(self,obj)    # if called on the instance, return a partially bound function
    
    def __call__(self,*args,**kwargs):
        obj = args[0]
        
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        
        key = (self.func,args[1:],frozenset(kwargs.items()))
        
        try:
            res = cache[key]
        except AttributeError:
            res = cache[key] = self.func(*args, **kw)
        
        return res
    


# let's see... there are two use cases:
#
# - for LISA, generally we will have a signal of fixed length (1 year), unless the
#   beginning is below a limiting frequency, in which case it's OK to truncate
# - for LIGO, we generally want a signal starting at a minimum frequency,
#   whatever the length

# - premature optimization is the root of all evil!


class fbbh(object):
    def __init__(self,state,minf=1e-5,deltat=None,samples=None,fiducialf=None,debug=False):
        self.state, self.debug = state, debug
        
        self.deltat          = self.checknyquist(deltat)
        self.minf, self.maxt = self.checkduration(minf)        
        self.samples         = self.checksamples(samples,self.maxt,self.deltat)
        
        self.df = 1.0 / (self.samples * self.deltat)
        self.kmin = int(self.minf/self.df)
        
        if self.debug: print "frequency resolution: %s, minimum discrete frequency: %s" % (self.df,self.df*self.kmin)
        
        self.fidf, self.fidk = self.checkfiducialf(fiducialf,self.df,self.kmin)
        
        self.f = self.df * N.arange(self.kmin,self.samples/2+1)      # go up to Nyquist frequency
    
    def checkfiducialf(self,fiducialf,df,kmin):
        if fiducialf:
            if fiducialf > self.maxf():
                print "fiducial frequency too high... proceeding without it"
                return None, None
            else:
                fidk = int(fiducialf/df - kmin)
                fidf = self.df * (kmin + fidk)
            
                if self.debug: print "fiducial frequency %s, adjusted to %s" % (fiducialf,fidf)
                
                return fidf, fidk
        else:
            return None, None
    
    # TO DO: need more breathing space for Nyquist and duration?
    
    def checknyquist(self,deltat):
        maxf = self.maxf()
        
        if deltat is None:
            deltat = 2**math.floor(math.log(0.5/maxf,2))       # 2**n multiples or submultiples of 1 s
            if self.debug: print "computed deltat:", deltat
        elif 0.5 / deltat < maxf:
            deltat /= 2**math.ceil(math.log(maxf/(0.5/deltat),2))
            if self.debug: print "adjusted deltat:", deltat
        
        if self.debug: print "max f: %s, Nyquist: %s" % (maxf,0.5/deltat)
        
        return deltat
    
    def checkduration(self,minf):
        # compute the duration of the waveform from fmin and maxf
        maxf = self.maxf()
        maxt = self.timing(maxf) - self.timing(minf)
        
        if self.debug: print "computed duration: %s" % maxt
        
        # if the requested duration is shorter than we have, adjust the initial frequency
        # note that setting state['T'] to None disables this logic
        if getattr(self.state,'T',None) and maxt > self.state.T:
            maxt = self.state.T
            
            fs = N.linspace(minf,maxf,1000)     # a little rough, but should be OK
            ts = self.timing(fs)
            
            minf = N.interp(-maxt,ts - ts[-1],fs)
            
            if self.debug:
                print ("requested duration: %s, adjusted minf: %s, check: %s" %
                       (maxt,minf,self.timing(maxf) - self.timing(minf)))
        
        return minf, maxt
    
    def checksamples(self,samples,maxt,deltat):
        if samples is None:
            samples = 2**math.ceil(math.log(maxt/deltat,2))
            if self.debug: print "computed samples: %s" % samples
        # elif samples * deltat < maxt:
        #     samples *= 2**math.ceil(math.log(maxt/(samples*deltat),2))
        #     if self.debug: print "adjusted samples: %s" % samples
        
        return samples
    
    
    Msun = 4.9257e-6
    Mpcs = 1.02927e14
    
    gamma = 0.577215664901
    PNlambda = -1987/3080
    PNtheta  = -11831/9240
    
    def PNalpha(self):
        eta, spinbeta, spinsigma = self.state.eta, self.state.spinbeta, self.state.spinsigma
        
        alpha = [1,
                 0,
                 (20/9) * (743/336 + 11/4*eta),
                 -16 * math.pi + 4 * spinbeta,
                 10 * (3058673/1016064 + 5429/1008*eta + 617/144*eta**2 - spinsigma),
                 0,                                                                 # omit constant term at 2.5PN
                 ( (11583231236531/4694215680 - 640/3*math.pi**2 - 6848/21*fbbh.gamma) +
                   (-15335597827/3048192 + 2255/12*math.pi**2 - 1760/3*fbbh.PNtheta + 12320/9*fbbh.PNlambda) * eta +
                   76055/1728*eta**2 - 127825/1296*eta**3 - 6848/21*math.log(4) ),
                 math.pi * (77096675/254016 + 378515/1512*eta - 74045/756*eta**2)
                 ]
        
        alphalog = [0,0,0,0,0,
                    math.pi * (38645/252 - 65/9*eta*3), # 2.5PN
                    -6848/21,                           # 3PN
                    0]
        
        return alpha, alphalog
    
    def PNbeta(self):
        alpha, alphalog = self.PNalpha()
        
        beta = ([1] +
                [((5-k)/5) * alpha[k] for k in xrange(1,5)] +
                [(-1/5) * alphalog[5], (-1/5) * (alpha[6] + alphalog[6]), (-2/5) * alpha[7]])
        betalog = [0,0,0,0,0,0,(-1/5) * alphalog[6],0]
        
        return beta, betalog
    
    def maxf(self):
        # waveform-ending logic
        try:
            return self.state.fend
        except AttributeError:
            Mc, eta = self.state.Mc, self.state.eta
            M = Mc * eta**(-3/5)            
            return 4400 / M                 # default to Schwarzschild
    
    
    def v(self,f):
        Mc, eta = self.state.Mc, self.state.eta
        
        M = Mc * eta**(-3/5)
        v = (math.pi * M * fbbh.Msun * f)**(1/3)        
        
        vs = [1]
        for i in xrange(8):
            vs.append(vs[-1] * v)
        logv = N.log(v)
        
        return vs,logv
    
    
    def phasing(self,f):
        Mc, eta, phi0, t0 = self.state.Mc, self.state.eta, self.state.phi0, self.state.t0
        pnorder = getattr(self.state,'pnorder',7)   # default is 7 for 3.5PN
        
        vs, logv = self.v(f)
        alpha, alphalog = self.PNalpha()
        
        phasing = 0
        for i in xrange(pnorder + 1):
            phasing += (alpha[i] + alphalog[i] * logv) * vs[i]
        
        phasing *= (3/(128*eta) / vs[5])
        
        # alternative-gravity correction
        # Cornish et al. define it as an additional phasing term of beta u^b, where u = pi Mc f = eta^(3/5) v^3
        if hasattr(self.state,'b'):
            # logu = (3/5) * math.log(eta) + 3 * logv
            # phasing += self.state.beta * N.exp(self.state.b * logu)
            phasing += self.state.beta * (math.pi * Mc * fbbh.Msun * f)**self.state.b
        
        phasing -= phi0                     # following Maggiore 5.275, but omitting constant -pi/4
        phasing += (2 * math.pi * t0 * f)
        
        # fiducialf adjustment
        if self.fidf:
            phasing -= (phasing[self.fidk] + 2 * math.pi * t0 * self.fidf)
        
        return phasing
    
    def amplitude(self,f):
        # Mc in Msun, d in Mpc
        Mc, d = self.state.Mc, self.state.d
        
        # checked against Eq. (5.274) of Maggiore
        A = math.sqrt(5/24) * math.pi**(-2/3) * (Mc * fbbh.Msun)**(5/6) / (d * fbbh.Mpcs)
        
        # debugging Mc derivative
        # A = getattr(self.state,'A',A); self.state.A = A
        
        return A * f**(-7/6) * (f <= self.maxf()) * (f >= self.minf)
    
    def timing(self,f):
        Mc, eta = self.state.Mc, self.state.eta
        pnorder = getattr(self.state,'pnorder',7)   # default is 7 for 3.5PN
        
        M = Mc * eta**(-3/5)
        
        vs, logv = self.v(f)
        beta, betalog = self.PNbeta()
        
        timing = 0        
        for i in xrange(pnorder + 1):
            timing += (beta[i] + betalog[i] * logv) * vs[i]
                
        timing *= (-(5*M*fbbh.Msun)/(256*eta) / vs[8])
        
        # alternative-gravity correction - TO DO: check derivation with change in definition of b
        # if hasattr(self.state,'b'):
        #     logu = (3/5) * math.log(eta) + 3 * logv
        #     timing += self.state.beta * ((5 - self.state.b)/5) * N.exp(self.state.b * logu)
        
        # find approximate ending of waveform, set timing to zero there
        # maxf = self.maxf()
        # endpos = N.nonzero(f <= maxf)[0][-1]
        # timing -= timing[endpos]
        
        # TO DO: fiducial f stuff
        
        return timing
    
    
    def hpc(self):
        inc = self.state.inc
        
        # checked against Eq. (5.274) of Maggiore        
        ap, ac = 0.5 * (1 + math.cos(inc)**2), math.cos(inc)
        # add (2 pi f tc) to phasing [or -1j * (2 pi f tc) to exponent] for time of coalescence
        # the sign of the exponential is consistent with numpy's fft, h(f) = int exp(-2 pi i f t) h(t) dt
        # per Maggiore's 5.275, add pi/2 to hp phasing to get hc
        hcp = self.amplitude(self.f) * N.exp(-1j * self.phasing(self.f))        
        
        hp = FrequencyArray.FrequencyArray(ap * hcp,      kmin=self.kmin,df=self.df)
        hc = FrequencyArray.FrequencyArray(ac * -1j * hcp,kmin=self.kmin,df=self.df)
        
        return hp,hc
    
    def hpcderiv(self,par):
        if par not in ['Mc','eta','phi0','t0','d','cosi','inc','b','beta']:
            raise NotImplementedError
        
        Mc, eta, d, inc = self.state.Mc, self.state.eta, self.state.d, self.state.inc
        spinbeta, spinsigma = self.state.spinbeta, self.state.spinsigma
        
        pnorder = getattr(self.state,'pnorder',7)   # default is 7 for 3.5PN
        
        if par == 'Mc':
            alpha = [-5/128,
                     0,
                     (-5*(743 + 924*eta))/32256,
                     -spinbeta/16 + math.pi/4,
                     (-5*(3058673 + 5472432*eta + 4353552*eta**2 - 1016064*spinsigma))/65028096,
                     (-5*(-7729 + 1092*eta)*math.pi)/32256,
                     ( 10052469856691/600859607040 + (76055*eta**2)/221184 - (127825*eta**3)/165888 -
                       (107*fbbh.gamma)/42 - (5*math.pi**2)/3 + eta*(-15335597827/390168576 +
                       (385*fbbh.PNlambda)/36 + (2255*math.pi**2)/1536 - (55*fbbh.PNtheta)/12) - (107*math.log(4))/42 ),
                     (-5*(-15419335 - 12718104*eta + 4975824*eta**2)*math.pi)/16257024
                     ]
            
            alphalog = [0,0,0,0,0,0,-107/42,0]
            
            vs, logv = self.v(self.f)
            
            dphasing = 0
            for i in xrange(pnorder + 1):
                dphasing += (alpha[i] + alphalog[i] * logv) * vs[i]
                        
            # since the factor is complex, the in-place operator would not promote the array
            dphasing = dphasing * (-1j/(eta * Mc) / vs[5])
            
            if hasattr(self.state,'b'):
                # logu = (3/5) * math.log(eta) + 3 * logv
                # dphasing += -1j/Mc * self.state.beta * self.state.b * N.exp(self.state.b * logu)
                
                u = math.pi * Mc * fbbh.Msun * self.f
                dphasing += -1j/Mc * self.state.beta * self.state.b * u**self.state.b
            
            if self.fidf:
                dphasing -= dphasing[self.fidk]
            
            damplitude = 5/(6 * Mc)
            # damplitude = 0
        elif par == 'eta':
            alpha = [0,
                     0,
                     (-743 + 1386*eta)/16128,
                     -(9/160)*spinbeta + (9*math.pi)/40,
                     (-3058673 + 1368108*eta + 6530328*eta**2 + 1016064*spinsigma)/5419008,
                     ((-7729 + 1092*eta)*math.pi)/10752,
                     ( -11328104339891/166905446400 + (15211*eta**2)/18432 -
                       (25565*eta**3)/6144 + (321*fbbh.gamma)/35 + 6*math.pi**2 +
            	       eta*(15335597827/650280960 - (77*fbbh.PNlambda)/12 -
            	       (451*math.pi**2)/512 + (11*fbbh.PNtheta)/4) + (321*math.log(4))/35 ),
            	     -((15419335 + 3633744*eta + 2132496*eta**2)*math.pi)/1548288
            	     ]
            
            alphalog = [0,0,0,0,0,(-38645*math.pi)/10752,321/35,0]
            
            vs, logv = self.v(self.f)
            
            dphasing = 0
            for i in xrange(pnorder + 1):
                dphasing += (alpha[i] + alphalog[i] * logv) * vs[i]
                        
            dphasing = dphasing * (-1j/(eta**2) / vs[5])
            
            if self.fidf:
                dphasing -= dphasing[self.fidk]
            
            # since the alt-gravity phasing term is a function of pi Mc f, there's no dependence on eta
            # and no additional dphasing term
            
            damplitude = 0
        elif par == 'phi0':
            dphasing = 1j
            damplitude = 0
        elif par == 't0':
            if self.fidf:
                dphasing = -2j * math.pi * (self.f - self.fidf)
            else:
                dphasing = -2j * math.pi * self.f    
            
            damplitude = 0
        elif par == 'd':
            dphasing = 0
            damplitude = -1/d
        elif par == 'b':
            u = math.pi * Mc * fbbh.Msun * self.f
            dphasing = -1j * self.state.beta * N.log(u) * u**self.state.b
            
            if self.fidf:
                dphasing -= dphasing[self.fidk]
            
            damplitude = 0
        elif par == 'beta':
            u = math.pi * Mc * fbbh.Msun * self.f
            dphasing = -1j * u**self.state.b
            
            if self.fidf:
                dphasing -= dphasing[self.fidk]
            
            damplitude = 0
        if par == 'cosi':
            ap, ac = math.cos(inc), 1
            dphasing, damplitude = 0, 1
        elif par == 'inc':
            ap, ac = -math.sin(inc) * math.cos(inc), -math.sin(inc)
            dphasing, damplitude = 0, 1
        else:
            ap, ac = 0.5 * (1 + math.cos(inc)**2), math.cos(inc)
        
        amp = self.amplitude(self.f)
        
        hcp = (amp * dphasing + amp * damplitude) * N.exp(-1j * self.phasing(self.f))
        
        hp = FrequencyArray.FrequencyArray(ap * hcp,      kmin=self.kmin,df=self.df)
        hc = FrequencyArray.FrequencyArray(ac * -1j * hcp,kmin=self.kmin,df=self.df)
        
        return hp,hc
    
