import mc3

mymodel = mc3.model(parameters=[mc3.uniform('x',-2,2,periodic=True),
                                mc3.uniform('y',-4,4)],
                    init=[0,0])

mydata = mc3.data()
mydata.logL = lambda X: -0.5 * ((X[0] - 1)**2 + X[1]**2)

mychain = mc3.chain(model=mymodel,
                    proposal=mc3.Gaussian(),
                    stepper=mc3.Metropolis(),
                    data=mydata)

mychain.run(iterations=30000)

import pylab

pylab.figure(1); pylab.hist(mychain.getpar('x'),bins=50)
pylab.figure(2); pylab.hist(mychain.getpar('y'),bins=50)
pylab.show()

# what's next?
# - get the BH example to work
# - try a parallel tempering example

# MPI runs - main initializes data and model
#          - then nodes > 0 call a single function that repeatedly evaluates signal/logL
#            blocking on receives of a state variable; and terminate when they receive None
#          - this can be used for "wasteful" parallelism (multiple proposals)
#            as well as parallel tempering

# if init is not given, then model can try to make it, using parameter calls

# need to save the state of the random number generator... probably in state
# (and we need discipline in using always the same one, and always from the main)

# checkpointing: after each step, pickle the state structure (for atomicity: do new file, move onto old)
# write out sample line (binary or ascii), appending to file. It may simplify things to keep the
# samples in a numpy array. (or lisaXML?) Include resume option in chain call.

# need some management of computed variables
