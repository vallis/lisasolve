import mc3
import searchfuncs

import BBH

import math
import numpy

# how can we make it easier to search over subspaces?
# - use default parameters?
class bbhmodel(mc3.model):
    def signal(self,state):
        bbh = BBH.BlackHoleBinary()
        
        # these can later be collected if desired
        state['m2'] = state['m1'] / state['q']
        state['beta'] = math.asin(state['sinb'])
        
        bbh.Mass1 = state['m1']
        bbh.Mass2 = state['m2']
        bbh.EclipticLatitude  = state['beta']
        bbh.EclipticLongitude = state['lambda']
                
        bbh.CoalescenceTime = state['ctime']
        
        # these parameters are recovered by the Fstat, so we can fix them to conventional values
        bbh.Distance, bbh.Inclination, bbh.Polarization, bbh.InitialAngularOrbitalPhase = 1.0e9, 0.0, 0.0, 0.0
        
        # this is standard
        bbh.TaperApplied = 7.0
        
        return bbh


class lisadata(mc3.data):
    def __init__(self,filename,wd=False,sr=30000,fmin=1.0e-5,synthlisatime=10000):
        [self.N,self.dt,sx,sy,sz] = searchfuncs.readfile(filename)
        
        (sxt, syt, szt) = map(numpy.fft.rfft,(sx,sy,sz))
        self.sIf = sxt; self.sIIf = (syt - szt) / math.sqrt(3.0)
        
        if wd:
            self.lisapsd = searchfuncs.new_tdiXmultsinwd_psd(self.N,self.dt)
        else:
            self.lisapsd = searchfuncs.new_tdiXmultsin_psd(self.N,self.dt)
        
        self.sr = sr
        self.fmin = fmin
        self.synthlisatime = synthlisatime
        
    def logL(self,bbh,state=None):
        F2t = searchfuncs.Ftnumpy(self.sIf,self.sIIf,bbh,
                                  self.lisapsd,self.N,self.dt,self.fmin,
                                  synthlisatime=self.synthlisatime)
        
        # look for maximum on either side of CoalescenceTime
        imax, F2tmax = searchfuncs.shift_to_max(F2t,self.sr)
        tcbest = bbh.CoalescenceTime + self.dt * (imax - self.sr)
        
        # re-evaluate the Fstat at that time
        bbh.CoalescenceTime = tcbest
        F2t = searchfuncs.Ftnumpy(self.sIf,self.sIIf,bbh,
                                  self.lisapsd,self.N,self.dt,self.fmin,
                                  synthlisatime=self.synthlisatime)
        
        # interpolate maximum
        imax, xmax, F2tmax, F2intmax = searchfuncs.shift_to_max_then_interpolate(F2t,self.sr)
        tcbest = bbh.CoalescenceTime + self.dt * (imax + xmax - self.sr)
        
        # this can be later collected if desired. It would be nice to get also the Fstat parameters
        if state != None:
            state['tc'] = tcbest
        
        return 0.5*F2intmax
    
    
