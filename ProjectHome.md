A collaborative, open-source repository of algorithms, tools, libraries, and documents to analyze LISA datasets, detect gravitational-wave sources, and estimate their parameters.

Not much here yet... currently the repository contains:
  * **mc3**, Michele Vallisneri's LISA-specific framework for MCMC simulations (work-in-progress, with some unresolved dependencies)
  * a suite of **common** Python utilities by Michele, including a numpy array subclass to model Fourier-domain signals (`FrequencyArray`), a utility class to "count down" iterations with ETA (`countdown`), and simple MPI master/slave and workflow management tools (`masterslave` and `workflow`)
  * **fastbinary**, a Python/C++ module to generate Galactic binary waveforms in the frequency domain using Michele's variant of the rapid MLDC code by Cornish and Lyttenberg
  * **pymultinest**, a Python bridge (by Michele) for Cambridge University's [MultiNest](http://ccpforge.cse.rl.ac.uk/gf/project/multinest/)

Material will continue to accrete here throughout 2011...