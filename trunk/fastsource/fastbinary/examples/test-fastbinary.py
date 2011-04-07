import math
import mc3lisa
import FastBinary

fastbin = FastBinary.FastGalacticBinary()

# -- binary one

# parameters: f [Hz], fdot [Hz/s], lat [-pi/2,pi/2], lon [0,2pi], A [strain], inc [0,pi], psi [0,pi], phi [0,2pi]
par1 = [1e-3, 0, 0.25*math.pi, 0.25*math.pi, 1e-22, 0.45*math.pi, 0.35*math.pi, 1.75*math.pi]
# T and dt are the default values; best to have T/dt = power of 2
tdi1 = mc3lisa.tdi.TDIf(xyz=fastbin.onefourier(vector=par1,T=6.2914560e7,dt=15))

# -- binary two

year = 3.15581498e7
par2 = [1e-3 + 3.5/year, 0, -0.4*math.pi, 1.65*math.pi, 7e-23, 0.15*math.pi, 0.65*math.pi, 0.25*math.pi]
tdi2 = mc3lisa.tdi.TDIf(xyz=fastbin.onefourier(vector=par2))

# -- A-E-T SNRs and match

print "SNR(h1) = ", math.sqrt(tdi1.normsq())
print "SNR(h2) = ", math.sqrt(tdi2.normsq())
print "(h1,h2) = ", math.sqrt(tdi1.dotprod(tdi2))

try:
    import matplotlib.pyplot as P
    
    print "Plotting |A| for h1 and h2..."
    
    P.plot(tdi1.Af.f,abs(tdi1.Af)); P.hold(True)
    P.plot(tdi2.Af.f,abs(tdi2.Af)); P.hold(False)
    P.show()
except:
    pass
