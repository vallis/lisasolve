import math

import mc3
import mc3lisa

mydata  = mc3lisa.lisadata('/Users/vallis/Desktop/mldc/challenge1B.2.1-training-frequency.xml')

mymodel = mc3lisa.bbhmodel([mc3.uniform('m1',1.0e6,5.0e6),
                            mc3.uniform('q',1.0,4.0),                           # q = m1/m2
                            mc3.uniform('sinb',-1.0,1.0),                       # sin(beta) in [-1,1]
                            mc3.uniform('lambda',0,2.0*math.pi,periodic=True)],
                            [3.0e6,2.5,0,0],                                    # init guesses
                            ['m2','beta','tc'],                                 # computed parameters to collect
                            ctime=15552000.0)                                   # tc_guess = about six months

print mymodel.initstate
print mydata.logL(mymodel.signal(mymodel.initstate),mymodel.initstate)
print mymodel.initstate
