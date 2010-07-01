import sys,random, math

import numpy as N
import matplotlib.pyplot as P

import mc3

# MODEL #1: (gridsize x gridsize) Ising model
#           there's only one parameter, "spins", which holds the state of the entire grid
#           as a numpy array

gridsize = 16
T = 4.5

# computing the total energy for up-down-left-right nearest-neighbor pairs
# for homogeneous ferromagnetic interaction
def energy(s):
    return -(N.sum(s.spins[:-1,:]*s.spins[1:,:]) + N.sum(s.spins[:,:-1]*s.spins[:,1:]) +    # bulk
             N.sum(s.spins[-1,:]*s.spins[0,:])   + N.sum(s.spins[:,-1]*s.spins[:,0]))       # boundary

# computing the average magnetization
def magnetization(s):
    return N.sum(s.spins) / gridsize**2

def absmag(s):
    return abs(magnetization(s))

model_1 = mc3.model(searchpars = [mc3.parameter('spins')],
                    otherpars  = [mc3.computed('m',magnetization),  # we teach the model how to compute m
                                  mc3.computed('am',absmag)],
                    logL = lambda s: -energy(s) / T)                # remember p \propto exp{-E/kt}

# CHAIN #1: single spin-flip proposal; does not mix very well...

def spinflip(current):
    newspins = current['spins'].copy()
    newspins[random.randrange(0,gridsize),random.randrange(0,gridsize)] *= -1.0
    
    proposed = current.model.state({'spins': newspins})
    
    return proposed

chain_1 = mc3.chain(model = model_1,
                    init  = {'spins': N.sign(N.random.randn(gridsize,gridsize))},
                    proposal = mc3.proposal(spinflip),
                    stepper = mc3.Metropolis(),
                    store = ['m','am'],                             # we'd like to store only the magnetization
                    progress = True)

chain_1.run(iterations=10000)

mc3.statchain(chain_1)
mc3.plotchain(chain_1,otherpars=['m','am'],fignum=1,xbins=20)

# CHAIN #2: Wolff flipping procedure

def neighbors(pos):
    return [((pos[0] + D[0]) % gridsize,
             (pos[1] + D[1]) % gridsize) for D in ((1,0),(-1,0),(0,1),(0,-1))]

# Wolff cluster-update algorithm (from Sethna 2006)
def wolff(current):
    spins = current['spins'].copy()
    
    begin = (random.randrange(0,gridsize),random.randrange(0,gridsize))
    spinbegin = spins[begin]
    
    toflip = [begin]
    
    # Wolff detailed-balance probability
    p = 1.0 - math.exp(-2.0/T)
    
    while toflip:
        if spins[toflip[0]] == spinbegin:
            spins[toflip[0]] *= -1
            
            for pos in neighbors(toflip[0]):
                if spins[pos] == spinbegin and random.random() < p:
                    toflip.append(pos)
        
        del toflip[0]
    
    return current.model.state({'spins': spins})

chain_2 = mc3.chain(model = model_1,
                    init  = {'spins': N.sign(N.random.randn(gridsize,gridsize))},
                    proposal = mc3.proposal(wolff),
                    stepper = mc3.always(), # the Wolff procedure does not require the Metropolis arbitration
                    store = ['m','am'],
                    progress = True)

chain_2.run(iterations=10000)

mc3.statchain(chain_2)
mc3.plotchain(chain_2,otherpars=['m','am'],fignum=2,xbins=20)

P.show()

