class chain(object):
    def __init__(self,model,proposal,stepper,data,temperature=1.0):
        self.model    = model
        self.proposal = proposal
        self.stepper  = stepper
        self.data     = data
        
        self.temperature = temperature
        
        self.state = None
        self.samples = []
    
    def init(self):
        self.state = self.model.initstate            
        self.state.logL = self.data.logL(self.model.signal(self.state),self.state)
    
    def iterate(self,current,temperature):
        newstate = self.proposal.propose(current)        
        newstate.logL = self.data.logL(self.model.signal(newstate),newstate)
        
        return self.stepper.step(current,newstate,temperature)

    def compute(self,state):
        for par in self.model.computed:
            if par not in state and hasattr(par,'compute'):
                state[par] = par.compute(self.model.signal(state))
    
    def store(self,state):
        pars =  [state[par] for par in self.model.saved]
        pars += [state.prior, state.logL]
        
        return pars
    
    def run(self,iterations):
        if len(self.samples) == 0:
            self.init()             # initialize from model.initstate
        
        for iter in range(iterations):
            self.state = self.iterate(self.state,self.temperature)  # advance the chain
            self.compute(self.state)                                # compute derived parameters
            self.samples.append(self.store(self.state))             # store desired parameters
    
    def getpar(self,par):
        # see if we have stored the parameter...
        try:
            ind = (self.model.saved + ['prior','logL']).index(par)
            return [sample[ind] for sample in self.samples]
        except:
            return []
    
