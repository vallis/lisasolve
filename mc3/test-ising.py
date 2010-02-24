import sys,random, math
import numpy
import mc3

if len(sys.argv) > 1:
    T = float(sys.argv[1])
else:
    T = 2.0 / math.log(1.0 + math.sqrt(2.0))    # critical temperature

# gridsize

d = 32

# model

mymodel = mc3.model(parameters=[mc3.parameter('spins')],
                    init=[numpy.ones([d,d])],
                    computed=[mc3.computed('m', lambda X: numpy.sum(X[0]) / d**2)],
                    saved=['m'])

# likelihood

mydata = mc3.data()
mydata.logL = lambda X: numpy.sum(X[0][:-1,:]*X[0][1:,:]) + numpy.sum(X[0][:,:-1]*X[0][:,1:])

# proposal

def spinflip(current):
    newspins = current['spins'].copy()
    newspins[random.randrange(0,len(newspins)),random.randrange(0,len(newspins))] *= -1.0
    
    proposed = mc3.state(current.model)
    proposed['spins'] = newspins
    
    return proposed

myprop = mc3.proposal()
myprop.propose = spinflip

# MCMC

mychain = mc3.chain(model=mymodel,
                    proposal=myprop,
                    stepper=mc3.Metropolis(),
                    data=mydata,
                    temperature=T)

mychain.run(iterations=50000)

import pylab

m = numpy.array(mychain.getpar('m'))

pylab.figure(1); pylab.hist(m,bins=50)
pylab.figure(2); pylab.plot(m)
pylab.show()
