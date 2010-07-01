import math, operator

import numpy as N

# would be fun to use metaprogramming to have model-specific states
class state(dict):
    def __init__(self,model,init=None):
        """Initialize a state for the given model; accepts a tuple, list,
        dict, or another state."""
        
        self.__dict__['model'] = model
        
        if isinstance(init,dict):
            for p in self.model.parameters:
                self[p] = p.renorm(init[p])
        elif isinstance(init,(tuple,list,N.ndarray)):
            for p,v in zip(self.model.parameters,init):
                self[p] = p.renorm(v)
    
    def __getitem__(self,attr):
        try:
            return dict.__getitem__(self,attr)
        except KeyError:
            # if we don't have it, maybe we know how to compute it,
            # and we'll cache it we're were at it!
            self[attr] = self.model.computed[attr](self)
            return dict.__getitem__(self,attr)
    
    # TO DO: what kind of access is needed by "for par in state"?
    def keys(self):
        return self.model.parameters + self.model.computed.keys()
    
    # the primary interface for using state is with class-like
    # attribute accessors, but we actually maintain the items
    # in a dict, for ease of transferring them across MPI
    def __getattr__(self,attr):
        try:
            return self.__dict__[attr]
        except KeyError:
            try:
                return self.__getitem__(attr)
            except KeyError:
                raise AttributeError, "state cannot access attribute %s" % attr
    
    def __setattr__(self,attr,value):
        self.__setitem__(attr,value)
    
    def asarray(self):
        """Return numpy vector of parameters; assumes they are all scalar."""
        return N.fromiter((self[p] for p in self.model.parameters),'d')
    
    def fromarray(self,expr):
        """Verifies that an array (passed as a tuple, list, dict, or numpy array)
        has the required number of parameter slots, and returns it as a numpy array.
        Covers both 1-D and 2-D cases; assumes all parameters are scalar."""
        
        if isinstance(expr,(tuple,list,N.ndarray)):
            array = N.array(expr,'d')
        # I suppose we could also support iterators
        elif isinstance(expr,dict):
            array = N.fromiter((expr[p] for p in self.model.parameters),'d')
        else:
            raise TypeError
        
        if array.shape[0] != self.model.d or (len(array.shape) == 2 and array.shape[1] != self.model.d):
            raise ValueError, "state.checkarray: wrong-sized array or vector for this model"
        
        return array
    
    def datamodel(self):
        return self.model.datamodel(state)
    
    # the following two methods support the sample storage mechanism of chain
    def select(self,pars):
        return [self[par] for par in pars]
    
    def index(self,tostore,item):
        return tostore.index(item)
    

# for the moment, we'll assume all the parameters are multiple, and that the dimensionality is fixed
# we're also probably assuming that all parameters are floats
class multistate(state):
    # the multistate should be initialized with ((a[0],b[0],c[0]),(a[1],b[1],c[1]),...)
    #                                      or   ({'a': a[0],'b': b[0],'c': c[0]},{'a': a[1],'b': b[1],'c': c[1]},...)
    #                                      or   (numpy[a[0],b[0],c[0]],numpy[a[1],b[1],c[1]])
    #                                      or   numpy[a[0],b[0],c[0],a[1],b[1],c[1],...]
    def __init__(self,model,dim,init=None):
        self.__dict__['model'] = model
        self.__dict__['dim']   = dim
        
        if isinstance(init,(tuple,list)) and len(init) == self.dim:
            for i,p in enumerate(self.model.parameters):
                self[p] = N.array([p.renorm(single[p]) if isinstance(single,dict) else p.renorm(single[i]) for single in init])
        elif isinstance(init,N.ndarray):
            for i,p in enumerate(self.model.parameters):
                self[p] = N.vectorize(p.renorm)(init[i::self.model.d])
    
    def singlestate(self,i):
        return singlestate(self,i)
    
    # support returning a single parameter list
    def asarray(self):
        ret = N.zeros(self.dim * self.model.d,'d')
        
        for i,p in enumerate(self.model.parameters):
            ret[i::self.model.d] = getattr(self,p)
        
        return ret
    
    def fromarray(self,expr):
        # for 1D, accept a list/tuple of the accepted single-state formats 
        if isinstance(expr,(tuple,list)) and len(expr) == self.dim:
            return [state.fromarray(self,e) for e in expr]
        elif isinstance(expr,N.ndarray):
            if expr.shape[0] != self.model.d * self.dim or (len(expr.shape) == 2 and expr.shape[1] != self.model.d * self.dim):
                raise ValueError, "state.checkarray: wrong-sized array or vector for this model"
        
        return expr
    
    def select(self,pars):
        return reduce(operator.add,([self[par][i] for par in pars] for i in range(self.dim)))
    
    def index(self,tostore,item):
        return [tostore.index(item) + i*len(tostore) for i in range(self.dim)]
    

# this is a really a proxy for a section of a multistate;
# still somewhat experimental, since it's not clear if .model should be multistate.model
# or the "single" model
class singlestate(state):
    def __init__(self,multistate,i):
        self.__dict__['multistate'] = multistate
        self.__dict__['i'] = i
        
        self.__dict__['model'] = multistate.model
    
    def asarray(self):
        return N.array([getattr(self.multistate,p)[self.i] for p in self.multistate.model.parameters])
    
    def __getitem__(self,attr):
        return self.multistate[attr][self.i]
    
    def __setitem__(self,attr,value):
        self.multistate[attr][self.i] = value
    
    def __getattr__(self,attr):
        return getattr(self.multistate,attr)[self.i]
    
    def __setattr__(self,attr,value):
        getattr(self.multistate,attr)[self.i] = value
    
