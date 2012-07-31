from __future__ import division
import math
import numpy as N
import FrequencyArray

# see LIGO T010110-00-Z
# but mostly from Nissanke et al ApJ 725, 496 (2010) -- TO DO: check delays

response = {}
position = {}   # meters, geocentric

response['LLO']   = N.matrix([[0.41128090024,   0.14020970464,   0.24729430676],
                              [0.14020970464,  -0.10900560021,  -0.18161569536],
                              [0.24729430676,  -0.18161569536,  -0.30227550864]])
position['LLO']   = N.matrix([-7.42760440000e+04, -5.49628372100e+06, 3.22425701800e+06])
        
response['LHO']   = N.matrix([[-0.39261409640,  -0.07761300355,  -0.24738860130],
                              [-0.07761300355,   0.31952440739,   0.22799809277],
                              [-0.2473886013 ,   0.22799809277,   0.07309029996]])
position['LHO']   = N.matrix([-2.16141490000e+06, -3.83469520000e+06, 4.60035020000e+06])

response['VIRGO'] = N.matrix([[ 0.24387399852,  -0.09908380359,  -0.23257620633],
                              [-0.09908380359,  -0.44782578945 ,  0.18783310056],
                              [-0.23257620633,   0.18783310056 ,  0.20395180583]])
position['VIRGO'] = N.matrix([4.54637409863e+06, 8.42989697467e+05, 4.37857696275e+06])

response['AIGO']  = N.matrix([[0.3804453308496,   0.248416253670,   0.096425872565],
                              [0.2484162536700,  -0.016170849085,  -0.200358164090],
                              [0.0964258725647,  -0.200358164090,  -0.364274481765]])
position['AIGO']  = N.matrix([-4.072708092382e+06, 8.46246235124e+06, -5.69421822793e+06])

response['LCGT']  = N.matrix([[ 0.140195381798,  -0.198339093258,   0.321822699094 ],
                              [-0.198339093258,  -0.348250870274,   0.1223989559938],
                              [ 0.321822699094,   0.122398955992,   0.208055488476 ]])
position['LCGT']  = N.matrix([-3.19388097472e+06, 2.94619660227e+06, 3.18561477478e+06])


speedoflight = 299792458  # m/s

# returns F+, Fx and delays
def FpcD(detector,skyposition):    
    theta, phi, psi = skyposition
        
    ex = N.matrix([ math.sin(phi) * math.cos(psi) - math.sin(psi) * math.cos(phi) * math.cos(theta),
                   -math.cos(phi) * math.cos(psi) - math.sin(psi) * math.sin(phi) * math.cos(theta),
                    math.sin(psi) * math.sin(theta)])
    
    ey = N.matrix([-math.sin(phi) * math.sin(psi) - math.cos(psi) * math.cos(phi) * math.cos(theta),
                    math.cos(phi) * math.sin(psi) - math.cos(psi) * math.sin(phi) * math.cos(theta),
                    math.cos(psi) * math.sin(theta)])
    
    R = response[detector]
    
    Fp = float(ex * R * ex.T - ey * R * ey.T)
    Fc = float(ex * R * ey.T + ey * R * ex.T)
    
    k = N.matrix([math.cos(theta)*math.cos(phi),math.cos(theta)*math.sin(phi),math.sin(theta)])    # vector to the source
    D = -float(k * position[detector].T) / speedoflight
    
    return Fp, Fc, D


class networkf(object):
    def __init__(self,signal=None,dets=None,data=None,deriv=None,fisher=False):
        if data:
            self.data, self.dets = data, data.keys()
        elif signal and dets:
            self.dets = dets
            self.data = {}
            
            skyposition = (signal.state.theta, signal.state.phi, signal.state.psi)
            
            if deriv:
                hp, hc = signal.hpcderiv(deriv,fisher=fisher)
            else:
                hp, hc = signal.hpc(fisher=fisher)
            
            f = signal.f
            
            for det in dets:
                # TO DO: check the sign here vs the definition of the Fourier transform
                #        also remember that FpcD returns a delay
                Fp, Fc, D = FpcD(det,skyposition)
                self.data[det] = N.exp(2*math.pi*-1j*f*D) * (Fp * hp + Fc * hc) * hp.df     # note that we store h(f_i) * df
                                                                                            # (same as for Galactic binaries)
        self._Sh = None

    def __add__(self,other):
        # in python 2.7, can do {k: v for (k,v) in ...}
        return networkf(data=dict((det,self.data[det] + other.data[det]) for det in self.dets))
    
    def __sub__(self,other):
        # in python 2.7, can do {k: v for (k,v) in ...}
        return networkf(data=dict((det,self.data[det] - other.data[det]) for det in self.dets))
    
    def __mul__(self,other):
        return networkf(data=dict((det,self.data[det]*other) for det in self.dets))
    
    def __rmul__(self,other):
        return networkf(data=dict((det,self.data[det]*other) for det in self.dets))
    
    def __div__(self,other):
        return networkf(data=dict((det,self.data[det]/other) for det in self.dets))
    
    def Sh(self,det):
        if self._Sh is None:
            self._Sh = dict((det,noisepsd(self.data[det])) for det in self.dets)
        return self._Sh[det]
    
    def normsq(self):
        norm = 4.0 / self.data[self.dets[0]].df
        return norm * sum(N.sum(N.abs(self.data[det])**2 / self.Sh(det)) for det in self.dets)
    
    def dotprod(self,other):
        norm = 4.0 / self.data[self.dets[0]].df
        return norm * N.real( sum(N.sum(N.conj(self.data[det]) * other.data[det] / self.Sh(det)) for det in self.dets) )
    
    def logL(self,other):
        norm = 4.0 / self.data[self.dets[0]].df
        return -0.5 * norm * sum(N.sum(N.abs(self.data[det].rsub(other.data[det]))**2 / self.Sh(det)) for det in self.dets)
    

# for the moment use advanced LIGO for all (from CSYP)
def noisepsd(frequencydata):
    f = frequencydata.f if isinstance(frequencydata,FrequencyArray.FrequencyArray) else f
    
    x = f / 215.0
    
    return 1e-49 * (x**-4.14 - 5.0*x**-2 + 111.0*(2.0 - 2.0*x**2 + x**4)/(2.0 + x**2))
