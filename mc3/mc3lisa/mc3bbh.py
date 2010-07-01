import numpy as N
import fastbbhc
import FrequencyArray
import mc3, tdi

class bbh(mc3.model):
    deltat, duration = 15.0, 2**21
    minf = 1e-5
    withsynthlisa, synthlisatime = True, 100000
    
    # we would only need this to call BBH directly, but let's write it down anyway...
    translate = {'m1': 'Mass1', 'm2': 'Mass2',
                 'lat': 'EclipticLatitude', 'lon': 'EclipticLongitude',
                 'ctime': 'CoalescenceTime', 'd': 'Distance',
                 'inc': 'Inclination', 'pol': 'Polarization', 'phi': 'InitialAngularOrbitalPhase'}
    
    @staticmethod
    def datamodel(state):
        binary = fastbbhc.fastBBH(state['m1'],state['m2'],state['ctime'],
                                  state['phi'],state['inc'],
                                  state['d'],7.0,
                                  state['lat'],state['lon'],-state['pol'],bbh.deltat)
        
        X = binary.Xm(bbh.duration,withsynthlisa=bbh.withsynthlisa,synthlisatime=bbh.synthlisatime)
        Y = binary.Ym(bbh.duration,withsynthlisa=bbh.withsynthlisa,synthlisatime=bbh.synthlisatime)
        Z = binary.Zm(bbh.duration,withsynthlisa=bbh.withsynthlisa,synthlisatime=bbh.synthlisatime)
        
        df = 1.0 / (bbh.duration * bbh.deltat)
        kmin = int(bbh.minf/df + 1.0)
        
        # consistently with mc3gwdb, we should return h(f_i) * df;
        # on the other hand, h(f_i) = rfft(X)_i / (N df) [for positive frequencies], so:
        
        Xf,Yf,Zf = map(lambda a: FrequencyArray.FrequencyArray(a[kmin:]/len(X),kmin=kmin,df=df),
                       map(N.fft.rfft,(X,Y,Z)))
        
        # to integrate correctly, we take half the Nyquist contribution
        Xf[-1] *= 0.5; Yf[-1] *= 0.5; Zf[-1] *= 0.5
        
        return tdi.TDIf(xyz=(Xf,Yf,Zf))
    
