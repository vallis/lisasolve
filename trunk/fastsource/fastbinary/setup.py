#!/usr/bin/env python

modulename  = 'FastBinary'
version     = '$Id: $'
description = 'A MLDC plugin to create Galactic-binary waveforms'
author      = 'Michele Vallisneri, based on work by Neil Cornish and others'
email       = 'vallis@vallis.org'
url         = 'http://lisasolve.googlecode.com'

csourcefiles = ['FastBinary.cpp','AEnoise.c','arrays.c']
sourcefiles  = csourcefiles + ['FastBinary.i']
headers      = ['FastBinary.h','arrays.h']

# please don't change anything below without contacting vallis@vallis.org

# note: setup.cfg specifies that we want SWIG to create a C extension

import os
import sys

prefix = ''; fftw_prefix = ''; argv_replace = []
for arg in sys.argv:
    if arg.startswith('--prefix='):
        prefix = arg.split('=', 1)[1]
        argv_replace.append(arg)
    elif arg.startswith('--with-fftw='):
        fftw_prefix = arg.split('=', 1)[1]
    else: 
        argv_replace.append(arg)
sys.argv = argv_replace

if fftw_prefix == '':
    if prefix == '':
        print >> sys.stderr, "You must specify FFTW location --with-fftw=<fftw_path>"
        sys.exit(1)
    else:
        fftw_prefix = prefix

from distutils.core import setup, Extension
from distutils.command.build import build
from distutils.command.install_lib import install_lib
from distutils.spawn import spawn
from distutils.file_util import copy_file
from distutils.util import get_platform

# process arguments...

argv_replace = []
make_clib = False

for arg in sys.argv:
    if arg == '--make-clib':
        make_clib = True
    else:
        argv_replace.append(arg)

sys.argv = argv_replace

# build the Python modules last so that the include the SWIG interface

class swig_build(build):
    # we assume build_py is the first in the sequence
    sub_commands = build.sub_commands[1:] + [build.sub_commands[0]]

# build also a C library

class qm_install_lib(install_lib):
    def run(self):
        install_lib.run(self)

        if self.distribution.has_c_libraries():
            build_clib = self.get_finalized_command('build_clib')
            libs = build_clib.get_library_names()

            clib_dir = build_clib.build_clib

            for lib in libs:
                clib = 'lib' + lib + '.a'

                src_file = os.path.join(clib_dir, clib)
                dest_file = os.path.join(self.install_dir, clib)

                copy_file(src_file, dest_file)

                if sys.platform[:6] == "darwin":
                    spawn(['ranlib'] + [dest_file])

# get the numpy installation dir so we know where to find the header
# (an interesting approach)

from numpy import __path__ as numpypath
numpyinclude = numpypath[0] + '/core/include'

# now run the setup

if make_clib == True:
    clibrary = [(modulename,{'sources': csourcefiles,
                             'depends': headers,
                             'include_dirs': [numpyinclude]})]

    # this hack needed on OS X for universal binaries since ar fails if the .a is already present...

    try:
        os.remove('build/temp.' + get_platform() + '-%s.%s' % sys.version_info[0:2] + '/lib' + modulename + '.a')
    except:
        pass
else:
    clibrary = []

setup(name = modulename,
      version = version,
      description = description,
      author = author,
      author_email = email,
      url = url,

      py_modules = [modulename],
      
      ext_modules = [Extension('_' + modulename,
                               sourcefiles,
                               include_dirs = [numpyinclude,fftw_prefix + '/include'],
                               library_dirs = [fftw_prefix + '/lib'],
                               libraries=['fftw3'],
                               depends = headers)],

      libraries = clibrary,      
                               
      cmdclass = {'build': swig_build, 'install_lib' : qm_install_lib}
      )
