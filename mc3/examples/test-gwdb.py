import sys, math

import matplotlib.pyplot as P
import numpy as N

import mc3
import mc3lisa

import FastBinary
import FrequencyArray

# DATA: the target source is a monochromatic Galactic binary

#       first we define the source parameters in a Python dictionary
truepars = {'f': 0.001,
            'lat': 0.25*math.pi,
            'lon': 0.25*math.pi,
            'A': 1e-22,
            'inc': 0.45*math.pi,
            'psi': 0.35*math.pi,
            'phi': 1.75*math.pi}

#       get the frequency-domain representation of the signal, and pad the data on both sides
data = mc3lisa.gwdb.datamodel(truepars)
data.pad()

#       then the log-likelihood of an MCMC state is simply
def datalogL(state):
    return data.logL(mc3lisa.gwdb.datamodel(state))

# CHAIN 1: search only over lat and lon, using arbitrary Gaussian scaling for the steps

#          these are parameters to search on
skypos = [mc3.latitudeangle('lat'),mc3.azimuthangle('lon')]

#          and these are parameters that we fix, using the mc3 "default" parameter
others = [mc3.default(par,truepars[par]) for par in ('f','A','inc','psi','phi')]

model_1 = mc3lisa.gwdb(searchpars=skypos,
                       otherpars=others,
                       logL=datalogL)

chain_1 = mc3.chain(model=model_1,
                    init=[truepars['lat'],truepars['lon']], # start the chain at the true values
                    proposal=mc3.Gaussian([0.005,0.005]),   # propose uncorrelated Gaussian jumps with arbitrary amplitude...
                    stepper=mc3.Metropolis(),               # ...and use the simple Metropolis rule
                    progress=True)                          # show progress while we run

print "\n=== CHAIN 1: search over lat and lon ==="

chain_1.run(iterations=10000)

mc3.statchain(chain_1)

#          before plotting, let's also see what the SNR and Fisher-matrix errors look like
#          note that the Fisher matrix does not know about the prior on lat
truestate = model_1.state(truepars)
snr, covmat = mc3lisa.SNR(truestate), N.linalg.inv(mc3lisa.fisher(truestate,[0.005,0.005])) # arbitrary scaling for finite differences

print 'SNR: %f, Fisher D(lat) = %f, D(lon) = %f' % (snr,math.sqrt(covmat[0,0]),math.sqrt(covmat[1,1]))
print 'Fisher covariance:'; print covmat

mc3.plotchain(chain_1,'lat-lon search, arbitrary steps',fignum=1,
              overplot=[mc3.normaldist(truepars['lat'],covmat[0,0]),    # a simple way to plot a normal distribution
                        mc3.normaldist(truepars['lon'],covmat[1,1])])   # of given mean and variance

# CHAIN 2: not bad, but we're not mixing too well---the correlation time is high.
#          Let's rescale the proposal using the Fisher matrix

chain_2 = mc3.chain(model=model_1,
                    init=[truepars['lat'],truepars['lon']],
                    proposal=mc3.Gaussian(covmat),
                    stepper=mc3.Metropolis(),
                    progress=True)

print "\n=== CHAIN 2: search over lat and lon with Fisher scaling ==="

chain_2.run(iterations=10000)

mc3.statchain(chain_2)
mc3.statchain(chain_2,report=False,skip=20)     # see if things change by culling samples

mc3.plotchain(chain_2,'lat-lon search, Fisher steps',fignum=2,
              overplot=[mc3.normaldist(truepars['lat'],covmat[0,0]),
                        mc3.normaldist(truepars['lon'],covmat[1,1])])

# CHAIN 3: now try a seven-dimensional search (all parameters);
#          again we'll use the Fisher matrix to scale the proposal
#          match with Fisher isn't bad, perhaps the f-derivative could be done more accurately
#          correlation time is long, but acceptance is below 50% already

allpars = [mc3.uniform('f',0.9e-3,1.1e-3),                       # not really uniform of course
           mc3.latitudeangle('lat'),mc3.azimuthangle('lon'),
           mc3.uniform('A',0.5e-22,1.5e-22),                     # ditto
           mc3.polarangle('inc'),mc3.uniform('psi',0,math.pi,periodic=True),mc3.azimuthangle('phi')]

model_3 = mc3lisa.gwdb(searchpars = allpars,
                       logL = datalogL)

truestate = model_3.state(truepars)
snr, covmat = mc3lisa.SNR(truestate), N.linalg.inv(mc3lisa.fisher(truestate,[1e-8,0.005,0.005,1e-23,0.005,0.005,0.005]))

chain_3 = mc3.chain(model=model_3,
                    init=truepars,
                    proposal=mc3.Gaussian(covmat),
                    stepper=mc3.Metropolis(),
                    progress=True)

print "\n=== CHAIN 3: search over all seven parameters with Fisher scaling ==="

chain_3.run(iterations=10000)

N.set_printoptions(precision=3)     # reduce output precision
mc3.statchain(chain_3)

print 'SNR: %f, Fisher D = %s (f,lat,lon,A,inc,psi,phi)' % (snr,N.sqrt(N.diagonal(covmat)))
print 'Fisher covariance:'; print covmat

mc3.plotchain(chain_3,'7-dim search, Fisher steps',fignum=3,
              overplot=[mc3.normaldist(truepars[p],covmat[i,i]) for (i,p) in enumerate(allpars)])   # do all distributions

P.show()
