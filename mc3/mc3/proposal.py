import math

import numpy

from state import state

# TO DO: need to keep the random state somehow
#        perhaps it can be done by the base class

class proposal(object):
    def __init__(self,func=None):
        self.func = func
    
    def propose(self,current):
        return self.func(current)
        

class Gaussian(proposal):
    def __init__(self,scale=None):
        self.scale = scale
        self.norm  = None
    
    def propose(self,current):
        x = current.asarray()
        
        # the first time we run, initialize the proposal normalization matrix
        if self.norm == None:
            if self.scale == None:
                # by default, do no scaling
                self.scale = numpy.ones(len(x),'d')
            else:
                # otherwise use the vector or matrix passed at initialization
                self.scale = current.fromarray(self.scale)
            
            # if we're carrying a vector, use it as a diagonal scaling matrix
            # if we're carrying a matrix, treat it as a covariance matrix,
            #    and use it square root (Cholesky) to produce correlated Gaussian displacements 
            self.norm = (numpy.diag(self.scale) if len(self.scale.shape) == 1
                                                else numpy.linalg.cholesky(self.scale)) / math.sqrt(len(x))
        
        y = x + numpy.dot(self.norm,numpy.random.randn(len(x)))
        
        return current.model.state(y)
    

class fromprior(proposal):
    def propose(self,current):
        return current.model.state([p.random() for p in current.model.parameters])
    
