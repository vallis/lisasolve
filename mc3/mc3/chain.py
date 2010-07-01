import sys, time

from state import state
from stepper import Metropolis

# TO DO: I suppose chain doesn't need to know about "data", since we could
# initialize model with a logL that uses the data... however there may be
# consequences for parallel operation... 

class chain(object):
    def __init__(self,model,proposal,stepper=Metropolis,data=None,init=None,store=None,temperature=1.0,progress=False):
        # these are always required
        self.model    = model
        self.proposal = proposal
        self.stepper  = stepper
        
        # data may not be present
        self.data = data
        
        # init may be passed as a state or as a dict/list/tuple...
        self.state      = init if isinstance(init,state) else self.model.state(init)
        
        self.state.logL = self.data.logL(self.state) if self.data else self.model.computed['logL'](self.state)
        
        # by default, we'll store only the search parameters
        self.tostore = store if store else self.model.parameters
        
        self.temperature = temperature
        
        self.progress = progress
        
        self.samples  = []
        
    def iterate(self,current,temperature):
        newstate = self.proposal.propose(current)
        newstate.logL = self.data.logL(newstate) if self.data else self.model.computed['logL'](newstate)
        
        return self.stepper.step(current,newstate,temperature)
    
    def store(self,state):
        return state.select(self.tostore)
    
    def run(self,iterations):
        t0, tlast = [time.time()] * 2
        
        for iter in range(iterations):
            self.state = self.iterate(self.state,self.temperature)  # advance the chain
            self.samples.append(self.store(self.state))             # store desired parameters
            
            t = time.time()
            if self.progress and iter > 0 and (t - tlast) > 1:
                tlast = t
                print '\r%d/%d complete (%d%%), %ds ETA...          ' % (iter+1,iterations,
                                                                         int(100.0 * iter / iterations),
                                                                         int(1.0 * (iterations - iter) * (t - t0) / iter)),
                sys.stdout.flush()
        
        if self.progress:
            t = time.time()
            print "\rCompleted %d in %ds, %.2f logL/s." % (iterations,t - t0,iterations / (t - t0))
    
    def __getitem__(self,item):
        try:
            index = self.state.index(self.tostore,item)
        except ValueError:
            raise KeyError
        
        if isinstance(index,(tuple,list)):
            return [[sample[i] for i in index] for sample in self.samples]
        else:
            return [sample[index] for sample in self.samples]
    
    def __getattr__(self,attr):
        try:
            return self.__dict__[attr]
        except KeyError:
            try:
                return self.__getitem__(attr)
            except KeyError:
               raise AttributeError
    
    def __len__(self):
        return len(self.samples)
    
