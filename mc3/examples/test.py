import math
import matplotlib.pyplot as P
import numpy as N

import mc3

# TO DO: use the generic function for normal distributions (it's in plot)

# TO DO: see if we can fix the overall title of multiple plots
# TO DO: show the estimated mean and variance (with est. errors), also after culling

# MODEL #1: create a simple bivariate Gaussian model with x in [0,4], y in [-5,5]

parameters_1 = [mc3.uniform('x',0,4),mc3.uniform('y',-5,5)]

def logL_1(s):
    return -0.5 * ((s.x - 2.0)**2 + (s.y/2.0)**2)   # note this is unnormalized, but that's OK for MCMC

model_1 = mc3.model(searchpars=parameters_1,logL=logL_1)

# CHAIN #1: propose following a simple unnormalized Gaussian random walk,
#           accept with the original Metropolis rule, run for 30,000 steps;
#           the resulting histograms follow the likelihood correctly,
#           and the autocorrelation plots suggest that we should use every
#           other 20th and 60th sample (respectively) to estimate moments

chain_1 = mc3.chain(model=model_1,
                    init=[0,0],
                    proposal=mc3.Gaussian(),
                    stepper=mc3.Metropolis(),
                    progress=True)

chain_1.run(iterations=30000)

# get statistics...
mc3.statchain(chain_1)

# ...and plot everything, and also the theoretical posteriors
# (the lambda-form is just a short-hand for defining functions in-line)
mc3.plotchain(chain_1,'p(x,y) = exp[-(x-2)^2/2 - (y/2)^2/2]',fignum=1,
              overplot=[mc3.normaldist(2.0,1.0),mc3.normaldist(0,2.0**2)])

# MODEL #2: extend the range of y to [-100,100], and its variance to 50;
#           now the chain explores y very poorly, and the autocorrelation
#           is very high; we expect ~ (2*50)^2 steps are necessary to walk across the range

parameters_2 = [mc3.uniform('x',0,4),mc3.uniform('y',-100,100)]
logL_2  = lambda s: -0.5 * ((s.x - 2.0)**2 + (s.y/50.0)**2)
model_2 = mc3.model(searchpars=parameters_2,logL=logL_2)
chain_2 = mc3.chain(model=model_2,init=[0,0],proposal=mc3.Gaussian(),stepper=mc3.Metropolis(),progress=True)
chain_2.run(iterations=30000)

mc3.statchain(chain_2)
mc3.plotchain(chain_2,'p(x,y) = exp[-(x-2)^2/2 - (y/50)^2/2]',fignum=2,corrwindow=[100,1000],
              overplot=[mc3.normaldist(2.0,1.0),mc3.normaldist(0,50.0**2)])


# CHAIN #3: scaling the proposal improves things considerably

chain_3 = mc3.chain(model=model_2,init=[0,0],proposal=mc3.Gaussian([1,50]),stepper=mc3.Metropolis(),progress=True)
chain_3.run(iterations=30000)

mc3.statchain(chain_3)
mc3.plotchain(chain_3,'p(x,y) = exp[-(x-2)^2/2 - (y/50)^2/2] (scaled prop.)',fignum=3,
              overplot=[mc3.normaldist(2.0,1.0),mc3.normaldist(0,50.0**2)])

# CHAIN #4: let's try out a correlated distribution
#           random-walking over the original variables ain't great
#           note the covariance matrix is [[101/4,-99/4],[-99/4,101/4]]

cxx = 101.0/4; cxy = -99.0/4

parameters_4 = [mc3.uniform('x',-20,20),mc3.uniform('y',-20,20)]
logL_4  = lambda s: -0.5 * ((s.x + s.y)**2 + ((s.x - s.y)/10.0)**2)
model_4 = mc3.model(searchpars=parameters_4,logL=logL_4,
                    otherpars=[mc3.computed('xpy', lambda s: s.x + s.y),    # to watch what x+y and x-y are doing,
                               mc3.computed('xmy', lambda s: s.x - s.y)])   # we teach model to compute them...
chain_4 = mc3.chain(model=model_4,init=[0,0],proposal=mc3.Gaussian(),stepper=mc3.Metropolis(),
                    store=['x','y','xpy','xmy'],                            # and we ask chain to store them
                    progress=True)
chain_4.run(iterations=30000)

mc3.statchain(chain_4)
mc3.plotchain(chain_4,'p(x,y) = exp[-(x+y)^2/2 - ((x-y)/10)^2/2]',fignum=4,otherpars=['xpy','xmy'],
              overplot=[mc3.normaldist(0,cxx),mc3.normaldist(0,cxx),
                        mc3.normaldist(0,1.0),mc3.normaldist(0,100.0)])

# CHAIN #5: we can do better by using a correlated proposal
# TO DO: estimate covariance from previous?

chain_5 = mc3.chain(model=model_4,init=[0,0],proposal=mc3.Gaussian([[cxx,cxy],[cxy,cxx]]),stepper=mc3.Metropolis(),
                    store=['x','y','xpy','xmy'],                            # and we ask chain to store them
                    progress=True)
chain_5.run(iterations=30000)

mc3.statchain(chain_5)
mc3.plotchain(chain_5,'p(x,y) = exp[-(x+y)^2/2 - ((x-y)/10)^2/2] (corr. prop.)',fignum=5,otherpars=['xpy','xmy'],
              overplot=[mc3.normaldist(0,cxx),mc3.normaldist(0,cxx),
              mc3.normaldist(0,1.0),mc3.normaldist(0,100.0)])

P.show()
