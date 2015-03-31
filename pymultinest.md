# Python multinest bridge: pymultinest #

## Installation ##

  * Prerequisites: slightly messy unfortunately; here I describe my latest OS X Snow Leopard install.
    * A working, gcc-compatible Fortran installation (such as gfortran 4.2.1 from http://r.research.att.com/tools).
    * You may have to reinstall numpy once the Fortran is in (for the sake of f2py).
    * A Fortran-aware MPI installation (e.g., Open MPI from http://www.open-mpi.org); the Open MPI shipped with OS X is not.
    * mpi4py from http://mpi4py.scipy.org.
  * The current svn version of the bridge includes the core MultiNest 2.8 files (from http://ccpforge.cse.rl.ac.uk/gf/project/multinest/), plus patched versions of `nested.F90` (which becomes `pynested.F90`), `Makefile`, and the f2py wrapper `nested.pyf`.
  * Then run `make pymultinest.so` in `lisasolve/multinest`. Some changes will be needed in the makefile for other OSes and architectures; please e-mail me and I'll try to help you (and augment the SVN)

## Usage ##

  * Look at `example_eggbox.py`. After initializing MPI and importing `pymultinest`, we need to define a `LogLike(cube, ndim)` function that returns a log likelihood where:
    * `cube` will be a Python list containing the `ndim` search parameters, rescaled to [0,1];
    * the function needs to set them (in place) to their physical values (implementing priors with nonlinear transformations if needed);
    * note that `len(cube)` may be larger than `ndim`; this allows returning "nonphysical" parameters that are not used in the search, but reported in output.
  * We then call the MultiNest executable using `pymultinest.nested.nestRun`. Parameters are explained in the MultiNest `readme.txt` and in `example_eggbox.py`. Remember to have an empty directory `chains` on hand.
  * With Open MPI, the parallel call would be something like `mpiexec -np 4 python example_eggbox.py`.

## Compilation problems ##

  * If the library complains about MPI libraries upon `import`ing, f2py may be doing the final library-compilation step with the regular Fortran compiler (e.g., gfortran) instead of the MPI version. One way to fix this is to modify the f2py line in the Makefile, by adding the linking directives that you can get from "mpif90 -link-info". For instance, on my Linux system it would be `-L/usr/lib -lmpichf90 -lmpichf90 -lmpich -lopa -lmpl -lrt -lcr -lpthread`.

## Compiling on Ubuntu 11/12 ##

  * On a fresh install (+ updates), install gfortran (sudo apt-get install gfortran), mpich2, subversion, python-dev, python-numpy, liblapack (liblapack-dev on Ubuntu 12). Download lisasolve as described in http://code.google.com/p/lisasolve/source/checkout.
  * May need to add -fPIC to the C and C++ compiler flags in the Makefile.

Michele, 5/22/2012