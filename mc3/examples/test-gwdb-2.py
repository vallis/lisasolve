import sys, math

import matplotlib.pyplot as P
import numpy as N

import mc3
import mc3lisa

import FastBinary
import FrequencyArray

# DATA: two target sources, separated by a few frequency bins

year = 3.15581498e7

src1 = {'f': 0.001,
        'lat': 0.25*math.pi,'lon': 0.25*math.pi,
        'A': 1e-22,
        'inc': 0.45*math.pi,'psi': 0.35*math.pi,'phi': 1.75*math.pi}

src2 = {'f': 0.001 + 3.5/year,
        'lat': -0.4*math.pi,'lon': 1.65*math.pi,
        'A': 7e-23,
        'inc': 0.15*math.pi,'psi': 0.65*math.pi,'phi': 0.25*math.pi}

#       sum the two sources in the frequency domain, and pad the resulting array   
data = mc3lisa.gwdb.datamodel(src1) + mc3lisa.gwdb.datamodel(src2)
data.pad()

#       the log likelihood is obtained by comparison with the (multi)model
def datalogL(state):
    return data.logL(mc3lisa.multigwdb.datamodel(state))

# CHAIN 1: 7-dimensional search

#          setup parameters as in the single-source case
allpars = [mc3.uniform('f',0.9e-3,1.1e-3),                       # not really uniform of course
           mc3.latitudeangle('lat'),mc3.azimuthangle('lon'),
           mc3.uniform('A',0.5e-22,1.5e-22),                     # ditto
           mc3.polarangle('inc'),mc3.uniform('psi',0,math.pi,periodic=True),mc3.azimuthangle('phi')]

#          create a "multimodel" with dimension two
model_1 = mc3lisa.multigwdb(dim=2,searchpars=allpars,logL=datalogL)

#           get the Fisher matrix for the joint model... use arbitrary scaling for the finite differences
truestate = model_1.state([src1,src2])

steps = [1e-8,0.005,0.005,1e-23,0.005,0.005,0.005]
covmat = N.linalg.inv(mc3lisa.fisher(truestate,steps * 2))

#           OK, let's run the chain

chain_1 = mc3.chain(model=model_1,
                    init=truestate,
                    proposal=mc3.Gaussian(covmat),
                    stepper=mc3.Metropolis(),
                    progress=True)

chain_1.run(iterations=10000)

#           get statistics...
N.set_printoptions(precision=3)
mc3.statchain(chain_1,cov=False,err=True)

#           plot distributions independently...
mc3.plotchain(chain_1,index=0,fignum=1,overplot=[mc3.normaldist(src1[p],covmat[i,i])     for (i,p) in enumerate(allpars)])
mc3.plotchain(chain_1,index=1,fignum=2,overplot=[mc3.normaldist(src2[p],covmat[i+7,i+7]) for (i,p) in enumerate(allpars)])

#           and again together
mc3.plotchain(chain_1,index=0,fignum=3)
mc3.plotchain(chain_1,index=1,fignum=3)

#           let's check the SNRs
SNR1 = mc3lisa.SNR(model_1.singlestate(src1))
SNR2 = mc3lisa.SNR(model_1.singlestate(src2))

print "SNR1 = ", SNR1, "SNR2 = ", SNR2

#           let's compare the double- and single-Fisher expected errors
covmat1 = N.linalg.inv(mc3lisa.fisher(model_1.singlestate(src1),steps))
covmat2 = N.linalg.inv(mc3lisa.fisher(model_1.singlestate(src2),steps))

print "Double-Fisher errors (src 1):"; print N.sqrt(N.diagonal(covmat)[0:7])
print "Single-Fisher errors (src 1):"; print N.sqrt(N.diagonal(covmat1))
print "Double-Fisher errors (src 2):"; print N.sqrt(N.diagonal(covmat)[7:14])
print "Single-Fisher errors (src 2):"; print N.sqrt(N.diagonal(covmat2))

P.show()
