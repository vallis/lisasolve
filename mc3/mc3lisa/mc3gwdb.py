import numpy as N
import FastBinary
import mc3, tdi

class gwdb(mc3.model):
    translate = {'f': 'Frequency', 'fdot': 'FrequencyDerivative',
                 'lat': 'EclipticLatitude', 'lon': 'EclipticLongitude',
                 'A': 'Amplitude',
                 'inc': 'Inclination', 'psi': 'Polarization', 'phi': 'InitialPhase'}
    
    @staticmethod
    def datamodel(state):
        initdict = dict((gwdb.translate[par],state[par]) for par in gwdb.translate if par in state.keys())
        binary = FastBinary.FastGalacticBinary(init = initdict)
        
        return tdi.TDIf(xyz=binary.fourier())
    


# will use the mc3.multimodel.__init__ by default; that's OK
# TO DO: certainly a synthetic TDI object would help here...
class multigwdb(mc3.multimodel):
    @staticmethod
    def datamodel(multistate):
        return reduce( operator.add,
                       (gwdb.datamodel(multistate.singlestate(i)) for i in range(multistate.model.dim)) )
    
    @property
    def singlemodel(self):
        return gwdb
    
