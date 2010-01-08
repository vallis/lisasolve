import math, copy

class chain(object):
    def __init__(self,model,proposal,stepper,data,temperature=1.0):
        self.model    = model
        self.proposal = proposal
        self.stepper  = stepper
        self.data     = data
        
        self.temperature = temperature
        
        self.state = None
        self.samples = []
    
    def run(self,iterations):
        if len(self.samples) == 0:
            self.state = self.model.initstate
            self.state.logL = self.data.logL(self.model.signal(self.state),self.state)
                
        for iter in range(iterations):
            newstate = self.proposal.propose(self.state)        
            newstate.logL = self.data.logL(self.model.signal(newstate),newstate)
            # should we cache the computation of the signal? It's called twice...
            
            self.state = self.stepper.step(self.state,newstate,self.temperature)
            
            for par in self.model.computed:
                if par not in self.state and hasattr(par,'compute'):
                    self.state[par] = par.compute(self.model.signal(self.state))
            
            pars =  [self.state[par] for par in self.model.saved]
            pars += [self.state.prior, self.state.logL]
            
            self.samples.append(pars)
        
    def getpar(self,par):
        # see if we have stored the parameter...
        try:
            ind = (self.model.saved + ['prior','logL']).index(par)
            return [sample[ind] for sample in self.samples]
        except:
            return []
    
