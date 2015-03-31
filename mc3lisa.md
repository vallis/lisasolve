# LISA signal processing: mc3/mc3lisa.tdi #

## Installation ##

  * Prerequisites: numpy, mc3 (but see below).
  * Get the lisasolve svn from lisasolve.googlecode.com. Or update to the latest if you have it already.
  * Go into lisasolve/mc3 and run `python setup.py install` (or `python setup.py install --prefix=...` if you've been using a user Python library directory). This will install mc3lisa as well as its prerequisite mc3.

## Usage ##

  * The basic object is the `TDIf` class, which contains a triad of TDI-observable datasets, each described as `FrequencyArray` frequency-domain objects. Internally, the triad is stored as A, E, T; it can be initialized as `TDIf(aet=(A,E,T))` or as `TDIf(xyz=(X,Y,Z))`. In the latter case the X array is also stored internally.

  * Currently `TDIf` objects can be added and subtracted, which will thread these operations over A, E, and T. In addition, the following operations are supported (where `t1` and `t2` are both `TDIf` objects:
```
t1.normsq()      # the noise-weighted inner product of t1 with itself;
                 # LISA noise PSDs are computed automatically over the required frequency range

t1.normsq(noisepsd=(Sa,Se,St))  # noise-weighted inner product with the given PSDs;
                                # these need to be arrays with the same length and offset as t1

t1.normsqx(...)  # same as normsq, but will only use the X channel

t1.dotprod(t2)   # noise-weighted inner product (t1,t2)

t1.logL(t2)      # noise-weighted likelihood -(t1-t2,t1-t2)/2 
```

  * Regarding noise models: the default setting is to use the current LISA-science-requirements prescriptions for the LISA secondary noises. This can be changed by setting `mc3lisa.tdi.defaultmodel` to `mldc` (the MLDC as-implemented values) or `mldc-nominal` (the slightly different MLDC documented values). `mc3lisa.tdi` provides the PSD functions `noisepsd_X(f)` and `noisepsd_AE(f)` and `noisepsd_T(f)`, which are used internally by `TDIf` and take an array of frequencies as argument.