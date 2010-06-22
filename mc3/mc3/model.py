import operator, random, math

import numpy as N

from state import state, multistate

class parameter(str):
    """Instances of this class represent parameters in the model.
    The only basic properties of a parameter are its name, and its prior
    distribution. The default prior for a parameter is uniform and unlimited
    (hence it's an improper prior)."""
    
    # because str is immutable, it is *created* (rather than initialized)
    # with a value; thus to subclass str we need to use __new__ to set the
    # value of the string. After that, __init__ will be called with *args
    # and **kwargs 
    def __new__(cls,parname,*args,**kwargs):
        return str.__new__(cls,parname)
    
    def __init__(self,*args,**kwargs):
        pass
    
    def prior(self,value):
        """The most conservative prior possible, to be overridden
        by parameter subclasses with more defined behavior (such as uniform)
        or by the user directly."""
        return 1.0
    
    def renorm(self,value):
        """Called to remap a parameter into its standard range; to be overridden."""
        return value
    


# note 1: since parameter is a subclass of str, it will hash like a string,
#         and therefore it can be used as key in a dictionary, interchangeably
#         with the parameter name

# note 2: if prior or renorm are replaced at runtime, as in
#         parinstance.prior = lambda x: f(x)
#         the bound instance method will be replaced by a function, which means
#         that it won't have access to "self". That's probably OK. Otherwise one would do
#         parinstance.method = new.instancemethod(lambda self, x: f(x), parinstance, parinstance.__class__)


class uniform(parameter):
    """Represents a uniformly distributed parameter, which may have a prescribed
    range, and may be periodic across that range."""
    
    def __init__(self,parname,xmin=None,xmax=None,periodic=False):
        parameter.__init__(self,parname)
        
        self.xmin, self.xmax = xmin, xmax
        self.periodic = periodic
    
    def prior(self,value):
        """If range lower or upper limits are set, prohibits parameter values
        below or above them, respectively. If both limits are set, returns the
        proper uniform prior for in-range values; otherwise returns unity."""
        
        if self.xmin != None and value < self.xmin:
            return 0
        elif self.xmax != None and value > self.xmax:
            return 0
        elif self.xmin != None and self.xmax != None:
            return 1.0 / (self.xmax - self.xmin)
        else:
            return 1.0
    
    def renorm(self,value):
        """Remaps a uniformly distributed, periodic parameter into its standard
        range."""
        
        if self.periodic == True and self.xmin != None and self.xmax != None:
            if value < self.xmin or value > self.xmax:
                # thanks to the behavior of %, this works even if (value - xmin) < 0
                return self.xmin + (value - self.xmin) % (self.xmax - self.xmin)
        
        return value
    
    def random(self):
        return random.uniform(self.xmin,self.xmax)
    


class azimuthangle(uniform):
    """Represents a longitude-like parameter, with range [0,2pi]."""
    
    def __init__(self,parname):
        uniform.__init__(self,parname,0,2.0*math.pi,periodic=True)
    


class polarangle(uniform):
    """Represents a colatitude-like parameter, with range [0,pi]."""
    
    def __init__(self,parname):
        # we cannot make it periodic and define renorm for it, since that would involve
        # changing the longitude at the same time
        uniform.__init__(self,parname,0,math.pi,periodic=False)
    
    def prior(self,value):
        if 0 < value < math.pi:
            return 0.5 * math.sin(value)
        else:
            return 0
    
    def random(self):
        return math.acos(random.uniform(-1.0,1.0))
    

class latitudeangle(uniform):
    """Represents a latitude-like parameter, with range [0,pi]."""
    
    def __init__(self,parname):
        # same story here as for polarangle
        uniform.__init__(self,parname,-0.5*math.pi,0.5*math.pi,periodic=False)
    
    def prior(self,value):
        if -0.5*math.pi < value < 0.5*math.pi:
            return 0.5 * math.cos(value)
        else:
            return 0
    
    def random(self):
        return math.asin(random.uniform(-1.0,1.0))
    


class computed(parameter):
    """Represents a parameter that is computed from the current state, and possibly
    from the corresponding data model."""
    
    def __init__(self,parname,computefunc):
        parameter.__init__(self,parname)
        
        self.compute = computefunc
    

class default(computed):
    """Represents a parameter that does not take part in the search, and has a default value."""
    
    def __init__(self,parname,value):
        computed.__init__(self,parname,lambda state: value)
    


class model(object):
    """Instances of this class represent a model of the data, as characterized
    by a set of parameters."""
    
    def __init__(self,searchpars,otherpars=[],prior=None,logL=None):
        # be generous, and normalize parameters that are passed simply as strings...
        self.parameters = [p if isinstance(p,parameter) else parameter(p) for p in searchpars]
        self.d = len(self.parameters)
        
        # otherpars may be passed as a list of computed/default objects, or as a shorthand
        # dictionary of callable and noncallable objects
        
        if isinstance(otherpars,dict):
            computedvars = [self.choose(*p) for p in otherpars.iteritems()]
        else:
            # normalize computed and default parameters that are passed as strings,
            # when a function or variable of the same name exists
            computedvars = [p if isinstance(p,parameter) else self.choose(p,globals()[p]) for p in otherpars]
        
        # it's actually better to hold these as a dictionary
        self.computed = dict((str(p),p.compute) for p in computedvars)
        
        # treat the prior as a standard computed parameter
        self.computed['prior'] = prior if prior else self.prior
        
        # if we're passed a model logL, also treat it as a computed parameter
        if logL:
            self.computed['logL'] = logL
        else:
            self.computed['logL'] = lambda s: 1.0
    
    # is it a computed parameter or a default value?
    def choose(self,x,y):
        return computed(x,y) if callable(y) else default(x,y)
    
    def state(self,init=None):
        """Returns a state for this model. Same as calling state(thismodel,init)."""
        
        return state(self,init)
    
    @staticmethod
    def datamodel(state):
        """Returns the 'value' of the model (to be compared with the data)
        for the parameters specified by 'state'. Defaults to a vector of
        the parameter values."""
        
        return tuple(state[par] for par in state.model.parameters)
    
    def prior(self,state):
        """Composes the prior for a model state by multiplying the priors
        for the individual parameters."""
        
        return reduce(operator.mul,(par.prior(state[par]) for par in self.parameters))
    

class multimodel(model):
    # need to force default parameters into numpy arrays, and hope they have the right dimension
    def __init__(self,dim,searchpars,otherpars=[],prior=None,logL=None):
        self.dim = dim
        model.__init__(self,searchpars,otherpars,prior,logL)
        
        self._singlemodel = self.singlemodel(searchpars,otherpars,prior,logL)
    
    def choose(self,x,y):
        return computed(x,y) if callable(y) else default(x,N.array(y if isinstance(y,(tuple,list)) else [y]*self.dim))
    
    def state(self,init=None):
        return multistate(self,self.dim,init)
    
    def singlestate(self,init=None):
        return state(self._singlemodel,init)
    
    def prior(self,state):
        return reduce(operator.mul,(reduce(operator.mul,(par.prior(singlepar)
                                                         for singlepar in state[par]))
                                    for par in self.parameters))
    
