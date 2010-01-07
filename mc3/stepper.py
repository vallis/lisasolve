import random

# need to keep the random state somehow
# perhaps it can be done by the base class

class stepper(object):
    pass

class Metropolis(stepper):
    def step(self,current,proposed):
        rho = proposed.p / current.p
        
        if rho > 1 or rho > random.random():
            return proposed
        else:
            return current
