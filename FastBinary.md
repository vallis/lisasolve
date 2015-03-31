# Fast waveform generators: fastsource/FastBinary #

## Installation ##

  * Prerequisites: numpy, lisatools (only lisatools/lisaXML/io-python, really), FFTW-3.
  * Get the lisasolve svn from lisasolve.googlecode.com. Or update to the latest if you have it already.
  * Go into lisasolve/common, and run `python setup.py install` (or `python setup.py install --prefix=...` if you've been using a user Python library directory).
  * Do the same in lisasolve/fastsource/fastbinary. If FFTW is not accessible in a system-wide location, pass its root directory to setup as --with-fftw=XXX (i.e., if the FFTW-3 files are in XXX/lib and XXX/include, pass XXX).

## Usage ##

  * Quickstart: see the [test-fastbinary.py](http://code.google.com/p/lisasolve/source/browse/trunk/fastsource/fastbinary/examples/test-fastbinary.py) example script.

  * To make a single binary, create a new `FastBinary.FastGalacticBinary` object:
```
>>> fastbin = FastBinary.FastGalacticBinary()
```

  * Then assign its standard parameters: `['Frequency', 'FrequencyDerivative', 'EclipticLatitude', 'EclipticLongitude', 'Amplitude', 'Inclination', 'Polarization', 'InitialPhase']`. Frequency and its derivative in Hz, Hz/s; angles in radians; amplitude as strain at SSB; amplitude/inclination related by the MLDC convention.
```
>>> fastbin.Frequency = 0.001
>>> fastbin.FrequencyDerivative = 0
>>> ...
```

  * You can also initialize the binary when you create it, by passing a dictionary with all the above parameters:
```
>>> pardict = {'Frequency': 0.001, 'FrequencyDerivative': 0, ...}
>>> fastbin = FastBinary.FastGalacticBinary(init=pardict)
```

  * Then the method `(X,Y,Z) = fastbin.onefourier(T=...,dt=...)` will return a 3-tuple consisting of the X, Y, and Z FFTs for the binary. The arguments T and dt, if not given, default to two "MLDC years" (62914560 s) and 15 s respectively. The FFTs are returned as lisasolve/common/FrequencyArray objects.

  * To describe what these objects are, let's take X for definiteness. X is a complex numpy array with bins spaced by `X.df`, beginning at frequency `X.kmin * X.df`, and ending at frequency `(X.kmin + len(X) - 1)*X.df`. The length and kmin of the array are computed by onefourier to enclose all the "interesting" part of the binary signal. Negative frequencies are not given since the underlying time-domain signal is real. If you sum two such FrequencyArrays for different binaries, the result will be an extended array that contains both.

  * If you want to do many binaries at once, don't bother to initialize fastbin, and call the method `(X,Y,Z) = fastbin.fourier(table=mytable,T=...,dt=...)` where mytable is a numpy array where each row contains the 8 parameters described above, in that order. In fact, mytable could be any Python "sequence" of "sequences", such as a list of lists, but numpy arrays will be fastest. In this case the FrequencyArray objects that you get back are based at kmin = 0 and extend up to the Nyquist frequency 0.5/dt.

  * If you want the time-domain X, Y, Z signals instead, call `(X,Y,Z) = fastbin.TDI(T=...,dt=...)` instead, which internally calls "fourier" and then does an inverse FFT.