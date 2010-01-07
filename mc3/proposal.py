import random, math

from state import state

# need to keep the random state somehow
# perhaps it can be done by the base class

class proposal(object):
    pass

class Gaussian(proposal):
    def __init__(self):
        pass
    
    def propose(self,current):
        proposed = state(current.model)
        
        norm = math.sqrt(len(current.model.parameters))
        for par in current.model.parameters:
            proposed[par] = current[par] + random.normalvariate(0,1) / norm
            
            if hasattr(par,'renormalize'):
                proposed[par] = par.renormalize(proposed[par])
        
        return proposed
