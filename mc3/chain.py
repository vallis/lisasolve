import math

class chain(object):
    def __init__(self,model,proposal,stepper,data):
        self.model    = model
        self.proposal = proposal
        self.stepper  = stepper
        self.data     = data
        
        self.state = None
        
        self.samples = []
    
    def run(self,iterations):
        if len(self.samples) == 0:
            self.state = self.model.initstate            
            self.state.logL = self.data.logL(self.model.signal(self.state),self.state)
        
        for iter in range(iterations):
            newstate = self.proposal.propose(self.state)        
            newstate.logL = self.data.logL(self.model.signal(newstate),newstate)
            
            self.state = self.stepper.step(self.state,newstate)
            
            pars =  [self.state[par] for par in self.model.parameters + self.model.computed]
            pars += [self.state.prior, self.state.logL, self.state.p]
            
            self.samples.append(pars)
    
    def getpar(self,par):
        if par in self.model.parameters:
            ind = self.model.parameters.index(par)
        elif par in self.model.computed:
            ind = len(self.model.parameters) + self.model.computed.index(par)
        else:
            return []
        
        return [sample[ind] for sample in self.samples]
