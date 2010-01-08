import random
import math

# need to keep the random state somehow
# perhaps it can be done by the base class

class stepper(object):
    pass

class Metropolis(stepper):
    def step(self,current,proposed,temperature=1.0):
        if proposed.prior == 0.0:
            return current
        
        logrho = ( (math.log(proposed.prior) + proposed.logL) -
                   (math.log(current.prior)  + current.logL ) ) / temperature
        
        if logrho > math.log(random.random()):
            return proposed
        else:
            return current
