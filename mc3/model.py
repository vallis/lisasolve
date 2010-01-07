from state import state

# this is OK, because all string subclasses will hash like a string,
# so state can remain a dict
# however subclassing immutable builtins is a little tricky...
class parameter(str):
    def __new__(cls,value,*args,**kwargs):
        return str.__new__(cls,value)
    

class uniform(parameter):
    def __init__(self,parname,xmin=None,xmax=None,periodic=False):
        parameter.__init__(parname)
        
        self.xmin, self.xmax = xmin, xmax
        self.periodic = periodic
    
    def prior(self,value):
        if self.xmin != None and value < self.xmin:
            return 0
        elif self.xmax != None and value > self.xmax:
            return 0
        elif self.xmin != None and self.xmax != None:
            return 1.0 / (self.xmax - self.xmin)
        else:
            return 1.0
    
    def renormalize(self,value):
        if self.periodic == True and self.xmin != None and self.xmax != None:
            if value < self.xmin or value > self.xmax:
                return self.xmin + (value - self.xmin) % (self.xmax - self.xmin)
        
        return value
    

class model(object):
    def __init__(self,parameters,init,computed=[],**kwargs):
        self.parameters = parameters
        self.computed = computed
        
        self.initstate = state(self)
        for i,par in enumerate(parameters):
            self.initstate[par] = init[i]
        
        self.extra = kwargs
    
    def signal(self,state):
        return tuple(state[par] for par in self.parameters)
    
    def prior(self,value):
        pi = 1.0
        
        for par in self.parameters:
            if hasattr(par,'prior'):
                pi = pi * par.prior(value[par])
        
        return pi
    
