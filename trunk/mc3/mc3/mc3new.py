import math, random, itertools, collections, operator
import countdown

import numpy as N, matplotlib.pyplot as P

# new definition of Parameter
# - every parameter is now a different class; its instances hold values (and are derived types of float, etc.)
# - this allows range checking, etc.
# - the parameter name is simply the class name


# newest-style parameters. Classes are parameters, instances are values. Probably overkill?

class Parameter(object):
    @property
    def name(self):
        return self.__class__.__name__
    
    def prior(self):    # by default, anything is possible
        return 1.0
    

def ParameterModel(name,baseclass,**kwargs):
    return type(name,(baseclass,Parameter),kwargs)

def UniformReal(name,range):
    def prior(self):
        if range[0] < self < range[1]:
            return 1.0 / (range[1] - range[0])
        else:
            return 0
    
    def random(self):
        return random.uniform(range[0],range[1])    # TO DO: standard seeding of randoms
    
    return ParameterModel(name,float,prior=prior,random=random)


# --- newer-style parameters. Instances are parameters, values are kept elsewhere.

# the old inheritance from str was not a bad idea afterall...
class Parameter(object):
    def __new__(cls,*args,**kwargs):
        if 'type' in kwargs:
            return super(Parameter,kwargs['type']).__new__(kwargs['type'])
        elif 'range' in kwargs:
            return super(Parameter,UniformReal).__new__(UniformReal)            
        else:
            return object.__new__(cls)
    
    def __init__(self,name):
        self.name = name
    
    def __str__(self):
        return self.name
    
    def prior(self,value):  # by default, anything is possible
        return 1.0
    

class UniformReal(Parameter):
    def __init__(self,name,range,periodic=False,**kwargs):
        Parameter.__init__(self,name)
        
        self.xmin, self.xmax, self.dx = range[0], range[1], range[1] - range[0]
        self.periodic = periodic
    
    def prior(self,value):
        if self.xmin < value < self.xmax:
            return 1.0 / self.dx
        else:
            return 0
    
    def random(self):
        return random.uniform(self.xmin,self.xmax)    # TO DO: standard seeding of randoms
    
    def normalize(self,value):
        if value < self.xmin or value > self.xmax and self.periodic == True:
            # thanks to the behavior of %, this works even if (value - xmin) < 0
            return self.xmin + (value - self.xmin) % self.dx
        else:
            return value
    

class PolarAngle(UniformReal):
    def __init__(self,name,**kwargs):
        UniformReal.__init__(self,name,range=(0,math.pi),**kwargs)
        
    def prior(self,value):
        if 0 < value < math.pi:
            return 0.5*math.sin(value)
        else:
            return 0
    
    def random(self):
        return math.acos(random.uniform(-1,1))  # TO DO: standard seeding of randoms
    

class CopolarAngle(UniformReal):
    def __init__(self,name,**kwargs):
        UniformReal.__init__(self,name,range=(-0.5*math.pi,0.5*math.pi),**kwargs)
    
    def prior(self,value):
        if -0.5*math.pi < value < 0.5*math.pi:
            return 0.5*math.cos(value)
        else:
            return 0
    
    def random(self):
        return math.asin(random.uniform(-1,1))  # TO DO: standard seeding of randoms
    

class Angle(UniformReal):
    def __init__(self,name,**kwargs):
        UniformReal.__init__(self,name,range=(0,2*math.pi),periodic=True,**kwargs)
    


# --- new-style states. Models are classes, instances are individual states.

class MemoFunc(object):
    # when used to decorate a method, it will perform
    # as a memoized property (with the cache kept in the instance under the same name)
    def __init__(self,name,func):
        self.name, self.func = name, func
    
    # TO DO: recognize when we're being called from the class, return function
    def __get__(self,instance,owner):
        cached = self.func(instance)
        setattr(instance,self.name,cached)
        return cached
    

class State(object):
    def __init__(self,init={},**kwargs):
        # initialize state from a dictionary or sequence (using the order in self.parameters)
        # will leave attributes undefined if not given, in which case defaults will apply
        
        if init == 'random':
            for par in self.parameters:
                setattr(self,str(par),par.random())
        elif isinstance(init,dict):
            self.__dict__.update(init)
        else:
            for par,val in zip(self.parameters,init):
                setattr(self,str(par),val)              # TO DO: renormalization of parameters?
                                                        #        currently done in proposal
        
        self.__dict__.update(kwargs)
    
    # TO DO: does not copy all!
    def copy(self):
        return self.__class__(self.parlist())
    
    def parlist(self,pars=None):
        return [getattr(self,str(par))                for par in (self.parameters if pars is None else pars)]
    
    def pardict(self,pars=None):
        return dict((str(par),getattr(self,str(par))) for par in (self.parameters if pars is None else pars))
    
    @property
    def prior(self):
        # since Python does not have "mul" in analogy to "sum", we have to make it like this
        return reduce(operator.mul,[par.prior(getattr(self,str(par))) for par in self.parameters])
    
    def normalize(self):
        for par in self.parameters:
            setattr(self,str(par),par.normalize(getattr(self,str(par))))
        
        return self
    
    # TO DO: should we enforce instance attributes as write-once? This would protect the memoization...
    

def Model(name,pars,defs={},**kwargs):
    # a factory for State classes, which teaches them:
    #   - their name
    #   - a list of parameters (required), will be coerced to Parameter class
    #   - a dictionary of default or computed (when callable) attribute values;
    #     the latter will be installed as memoized properties
    #   - additional default or computed values given as keyword arguments (for syntactic sugar)
    
    namespace = {'parameters': [par if isinstance(par,Parameter) else Parameter(par) for par in pars], 'n': len(pars)}
    # namespace = {'parameters': [par if issubclass(par,Parameter) else ParameterModel(par,float) for par in pars]}
    
    for par,val in itertools.chain(defs.items(),kwargs.items()):
        namespace[par] = MemoFunc(par,val) if isinstance(val,collections.Callable) else val
    
    # if we don't worry about the caching, could do
    # import types
    # namespace[par] = types.MethodType(val,None,State)
    # although it may be a problem that we're binding to the base class
    # maybe we should add these _after_ creating the new class below
    # but not really, because we really want to make a property...
    # so we should just use property(val)?
    
    return type(name,(State,),namespace)
    


# TO DO: we're assuming all models have the same parameters! otherwise we cannot de/serialize to/from lists
class StateArray(object):
    def __init__(self,init={}):
        if init == 'random':
            self.states = [model('random') for model in self.models]
        if len(init) == self.m:
            self.states = [model(single) for model,single in zip(self.models,init)]
        elif len(init) == self.n:
            self.states = [model(init[i*self.ncomp:(i+1)*self.ncomp]) for (i,model) in enumerate(self.models)]
        else:
            self.states = [model(init) for model in self.models]
    
    def __getitem__(self,index):
        return self.states[index]
    
    def __getattr__(self,par):
        return [getattr(state,par) for state in self.states]
    
    # the parlist is serialized, e.g. [s[0].a, s[0].b, s[1].a, s[1].b, ...]
    def parlist(self,pars=None):
        return reduce(operator.add,[state.parlist(pars) for state in self.states])
    
    def pardict(self,pars=None):
        return [state.pardict(pars) for state in self.states]
    
    @property
    def prior(self):
        return reduce(operator.mul,[state.prior for state in self.states])
    
    def logL(self):
        return reduce(operator.mul,[state.logL for state in self.states])
    
    def normalize(self):
        for state in self.states:
            state.normalize()
        
        return self
    

# TO DO: in future may allow for Array parameters (i.e., searchable) currently allowing only defs
def ModelArray(name,models,defs={},**kwargs):
    # TO DO: assuming all models are the same; idiomatic [model] * m initialization is possible!
    namespace = {'models': models, 'parameters': models[0].parameters,
                 'm': len(models), 'ncomp': models[0].n, 'n': len(models)*models[0].n}
    
    for par,val in itertools.chain(defs.items(),kwargs.items()):
        namespace[par] = MemoFunc(par,val) if isinstance(val,collections.Callable) else val
    
    return type(name,(StateArray,),namespace)


# --- new-style proposers. Functions. Will need temperature?

def RandomPropose(state):
    return type(state)('random')

class GaussianPropose(object):
    def __init__(self,scale=1,covmat=None):
        if covmat:
            self.norm = N.linalg.cholesky(N.asarray(covmat)) / math.sqrt(len(covmat))
        else:
            self.norm = N.asarray(scale)            
    
    def __call__(self,current):
        proposed = N.asarray(current.parlist()) + N.dot(self.norm,N.random.randn(current.n))
        
        return current.__class__(proposed).normalize()
    

class MultiPropose(object):
    # take a list of 2-tuples (rule,prob)
    def __init__(self,*args):
        self.rules, self.probs = [], []
        
        cum = 0
        for rule,prob in args:
            self.rules.append(rule)
            
            cum = cum + prob
            self.probs.append(cum)
        
        for i in range(len(self.probs)):
            self.probs[i] /= cum
    
    def __call__(self,current):
        x = N.random.random()
        
        for i,prob in enumerate(self.probs):
            if x < prob:
                return self.rules[i](current)
    

# --- new-style steppers. Functions.

def AlwaysStep(current,proposed):
    return proposed

def MetropolisStep(current,proposed):
    if proposed.prior == 0:
        return current
    else:
        logrho = ( (math.log(proposed.prior) + proposed.logL) -
                   (math.log(current.prior)  + current.logL ) )     # / temperature
        logx = math.log(random.random())
        
        return proposed if logrho > logx else current


# --- new-style chain

class Chain(object):
    def __init__(self,model,proposer=RandomPropose,stepper=MetropolisStep,init='random',store=None): 
        self.model, self.proposer, self.stepper = model, proposer, stepper        
        self.current = init if isinstance(init,self.model) else self.model(init)
        self.store = store
        
        self.N, self.accepted = 0, 0
        
        self.samples = [self.current.parlist(self.store)]
        self.cols = _NamedColumns(self.samples,self.model,self.store)
    
    def run(self,iterations):
        # TO DO: can we use countdown with "with"? or even better, somehow embed it in the for loop?
        C = countdown.countdown(iterations,100)
        for iter in range(iterations):
            C.status()
        
            # make a proposal for next sample
            proposed = self.proposer(self.current)
            
            # apply the stepping criterion
            self.current = self.stepper(self.current,proposed)
            if self.current is proposed:
                self.accepted += 1
            
            # store the new sample
            self.samples.append(self.current.parlist(self.store))
            self.N += 1
        C.end()
    
    def plot(self,xbins=50,window=100,p={},fignum=1,title='MCMC parameter plots'):
        P.figure(fignum)
        
        subplots = len(self.cols.names)
        for i,par in enumerate(self.cols.names):
            x = getattr(self.cols,par)
            
            # histogram
            P.subplot(subplots,3,3*i+1); P.hist(x,bins=xbins,normed=True); P.xlabel(par)
            if par in p:
                z = N.linspace(N.min(x),N.max(x),100)
                P.hold(True); P.plot(z,N.vectorize(p[par])(z),'r'); P.hold(False)
            
            # trajectory
            P.subplot(subplots,3,3*i+2); P.plot(x); P.xlabel(par)
            
            # correlation -- TO DO: try to find window that yields decent correlation times
            xm = x - N.mean(x)
            c = N.correlate(xm,xm,'full')[len(x)-1:len(x)+window]; c = c / c[0]
            P.subplot(subplots,3,3*i+3); P.plot(c); P.xlabel(par)        
        
        P.suptitle(title)
    
    def stat(self,skip=1):
        print "%d samples, %d moves accepted (%d%%)" % (self.N,self.accepted,100*self.accepted/self.N)
        
        # TO DO: reformat once we have storage backends
        parray = N.array(self.samples[::skip],'d')
        
        means  = N.mean(parray,axis=0)
        covmat = N.cov(parray.T)
        errors = N.sqrt(N.diagonal(covmat))
        
        maxlen = max(len(par) for par in self.cols.names)
        for i,par in enumerate(self.cols.names):
            print "%s: %g (%g)" % (par.ljust(maxlen),means[i],errors[i])
        
        return means, errors, covmat
    


# proxy object to access a column in a list of lists
class _NamedColumns(object):
    def __init__(self,samples,model,store):
        self.samples = samples
        self.model = model
        self.indices = dict((str(par),i) for i,par in enumerate(store if store is not None else self.model.parameters))
        
        # TO DO: homogenize this with numpy or pytables usage for column names
        if issubclass(self.model,StateArray):
            self.names = reduce(operator.add,(['%s[%d]' % (str(par),i) for par in self.model.parameters] for i in range(self.model.m)))
        else:
            self.names = [str(par) for par in self.model.parameters]
    
    def __getattr__(self,par):
        if issubclass(self.model,StateArray):
            if '[' in par:
                p, i = par.split('['); i = int(i[:-1])
                return [sample[i*self.model.ncomp + self.indices[p]] for sample in self.samples]
            else:
                return [[sample[i*self.model.ncomp + self.indices[par]] for sample in self.samples] for i in range(self.model.m)]
        else:
            return [sample[self.indices[par]] for sample in self.samples]
    


# in storing array values, need to use some kind of parameter proxy that looks like a[0], a[1], etc., for easy plotting and stat-ing
# or maybe multistate offers access to parameters that returns tuples!

# offer various storage backends: State, numpy, tables
# collect important stuff in mc3.core, then extensions
