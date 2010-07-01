import sys, math

import matplotlib.pyplot as P
import numpy as N

import mc3
import mc3lisa

# DATA: the target source is a "challenge 1B" MBH binary

#       first we define the source parameters in a Python dictionary
truepars = {'m1': 1204622.0708,
            'm2': 401333.588406,
            'lat': 0.156694915526,
            'lon': 5.99832784615,
            'ctime': 13636993.4823,
            'd': 44784411151.1,
            'inc': 1.1500606131,
            'pol': 2.60961000859,
            'phi': 3.67499142949}

#       get the frequency-domain representation of the signal;
#       mc3lisa.bbh uses Michele's "fastbbhc" rigid adiabatic approximation
data = mc3lisa.bbh.datamodel(truepars)

#       then the log-likelihood of an MCMC state is simply
def datalogL(state):
    return data.logL(mc3lisa.bbh.datamodel(state))

# CHAIN 1: search only over lat and lon, using two-dimensional Fisher scaling

#          these are parameters to search on
skypos = [mc3.latitudeangle('lat'),mc3.azimuthangle('lon')]

#          and these are parameters that we fix, using the mc3 "default" parameter
others = [mc3.default(par,truepars[par]) for par in ('m1','m2','ctime','d','inc','pol','phi')]

model_1 = mc3lisa.bbh(searchpars=skypos,
                      otherpars=others,
                      logL=datalogL)

#          get SNR and Fisher
truestate = model_1.state(truepars)
snr, covmat = mc3lisa.SNR(truestate), N.linalg.inv(mc3lisa.fisher(truestate,[0.005,0.005])) # arbitrary scaling for finite differences

print 'SNR: %f, Fisher D(lat) = %f, D(lon) = %f' % (snr,math.sqrt(covmat[0,0]),math.sqrt(covmat[1,1]))
print 'Fisher covariance:'; print covmat

chain_1 = mc3.chain(model=model_1,
                    init=[truepars['lat'],truepars['lon']], # start the chain at the true values
                    proposal=mc3.Gaussian(covmat),          # propose Fisher-correlated Gaussian jumps
                    stepper=mc3.Metropolis(),               # ...and use the simple Metropolis rule
                    progress=True)                          # show progress while we run

print "\n=== CHAIN 1: search over lat and lon ==="

chain_1.run(iterations=10)

mc3.statchain(chain_1)

mc3.plotchain(chain_1,'lat-lon search, Fisher steps',fignum=1,
              overplot=[mc3.normaldist(truepars['lat'],covmat[0,0]),    # a simple way to plot a normal distribution
                        mc3.normaldist(truepars['lon'],covmat[1,1])])   # of given mean and variance

N.array(chain_1.samples,'d').tofile('/Users/vallis/Dropbox/bbh2-samples.bin')

# CHAIN 2: search over all parameters with Fisher scaling

#          take a 10% uniform window for non-angular parameters
allpars = [mc3.uniform('m1',truepars['m1']*0.95,truepars['m1']*1.05),
           mc3.uniform('m2',truepars['m2']*0.95,truepars['m2']*1.05),
           mc3.latitudeangle('lat'),mc3.azimuthangle('lon'),
           mc3.uniform('ctime',truepars['ctime']*0.95,truepars['ctime']*1.05),
           mc3.uniform('d',truepars['d']*0.95,truepars['d']*1.05),
           mc3.polarangle('inc'),mc3.uniform('pol',0,math.pi,periodic=True),mc3.azimuthangle('phi')]

model_2 = mc3lisa.bbh(searchpars = allpars,
                      logL = datalogL)

truestate = model_2.state(truepars)
deltas = [truestate['m1']*1e-6,truestate['m2']*1e-6,0.005,0.005,10.0,truestate['d']*1e-3,0.005,0.005,0.005]
snr, covmat = mc3lisa.SNR(truestate), N.linalg.inv(mc3lisa.fisher(truestate,deltas))

chain_2 = mc3.chain(model=model_2,
                    init=truepars,
                    proposal=mc3.Gaussian(covmat),
                    stepper=mc3.Metropolis(),
                    progress=True)

print "\n=== CHAIN 2: search over all nine parameters with Fisher scaling ==="

chain_2.run(iterations=10)

N.set_printoptions(precision=3)     # reduce output precision
mc3.statchain(chain_2)

print 'SNR: %f, Fisher D = %s (m1,m2,lat,lon,ctime,d,inc,pol,phi)' % (snr,N.sqrt(N.diagonal(covmat)))
print 'Fisher covariance:'; print covmat

mc3.plotchain(chain_2,'9-dim search, Fisher steps',fignum=3,
              overplot=[mc3.normaldist(truepars[p],covmat[i,i]) for (i,p) in enumerate(allpars)])   # do all distributions

N.array(chain_2.samples,'d').tofile('/Users/vallis/Dropbox/bbh9-samples.bin')

P.show()
