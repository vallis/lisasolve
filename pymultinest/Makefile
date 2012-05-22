# === OS X configuration, using GCC compilers ===
# == using gcc and gfortran; on Linux, probably need to add -fPIC to CC and CXX ==
CC  = mpicc
CXX = mpicxx
FC  = mpif90 -DMPI
CFLAGS   += -I. -O2
CXXFLAGS +=  -I. -O2
FFLAGS   += -O3 -ffree-form -ffree-line-length-none
# warning: f2py needs LDFLAGS to be unset
AR = ar r
# == linker options for any executables that use multinest ==
LIBS = -lmultinest -lpthread -framework vecLib
# == f2py options ==
COMPILER = gnu95
LINKLIB  = --link-lapack_opt

# === Python install location ===
PYTHON_INCLUDE = $(shell python -c "import distutils.sysconfig; print distutils.sysconfig.get_config_var('CONFINCLUDEPY')")
PYTHON_LIBNAME = python$(shell python -c "import distutils.sysconfig; print distutils.sysconfig.get_config_var('VERSION')")
PYTHON_LIBDIR  = $(shell python -c "import distutils.sysconfig; print distutils.sysconfig.get_python_lib()")

# === multinest files ===
LIBOBJ = utils.o utils1.o kmeans_clstr.o xmeans_clstr.o posterior.o priors.o nested.o
PYOBJ  = utils.o utils1.o kmeans_clstr.o xmeans_clstr.o posterior.o priors.o pynested.o

# === rules ===
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^ 

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^ 

pymultinest.so: $(PYOBJ)
	f2py --fcompiler=$(COMPILER) --f90exec=$(FC) --f90flags="$(FFLAGS)" $(LINKLIB) --opt=-O2 -c nested.pyf $^
	cp pymultinest.so $(PYTHON_LIBDIR)

multinest.a: $(LIBOBJ)
	$(AR) $@ $^

all: pymultinest.so multinest.a

clean: 
	-rm *.a *.so *.o *.mod
