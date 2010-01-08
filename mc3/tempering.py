import math, random
import mc3

class tempering(mc3.chain):
    def __init__(self,model,proposal,stepper,data,temperature,alpha):
        mc3.chain.__init__(self,model,proposal,stepper,data,temperature)
        
        self.N = len(temperature)
        self.alpha = alpha
        
        if self.N < 2 or (not 0 < self.alpha < 1.0):
            raise ValueError
    
    def run(self,iterations):
        if len(self.samples) == 0:
            self.init()
            self.states = [self.state] * self.N
        
        attempted, swapped = 0, 0
        
        for iter in range(iterations):
            if random.random() < self.alpha:
                attempted += 1
                i = random.randrange(0,self.N-1)
                
                logp1 = math.log(self.states[i  ].prior) + self.states[i  ].logL
                logp2 = math.log(self.states[i+1].prior) + self.states[i+1].logL
                
                logrho = (logp1 - logp2) / (self.temperature[i+1] - self.temperature[i])
                
                if logrho > math.log(random.random()):
                    self.states[i], self.states[i+1] = self.states[i+1], self.states[i]
                    swapped += 1
            else:
                for i in range(self.N):
                    self.states[i] = self.iterate(self.states[i],self.temperature[i])
            
            for state in self.states:
                self.compute(state)
            
            self.samples.append(sum([self.store(state) for state in self.states],[]))
        
        print "Swaps: %d/%d/%d" % (swapped,attempted,iterations)
    
    def getpar(self,par,chain=0):
        try:
            allsaved = self.model.saved + ['prior','logL']
            ind = chain * len(allsaved) + allsaved.index(par)
            return [sample[ind] for sample in self.samples]
        except:
            return []
    
