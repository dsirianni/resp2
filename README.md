# RESP2
[![Build Status](https://travis-ci.org/rmcgibbo/resp2.svg?branch=master)](https://travis-ci.org/rmcgibbo/resp2)

_Restrained electrostatic potential (RESP) fitting for PSI4_

Overview
--------

This is a partially-complete plugin for [PSI4](https://github.com/psi4/psi4public) that implements restrained electrostatic potential fitting (RESP) [1] to compute atomic partial charges, potentially for use in MD.

This plugin still probably needs to be verified against a known-good implementation of this method, such as the code in AmberTools.

RESP is really a pretty simple method, and PSI4 has a _very_ nice plugin architecture for adding in new functionality that
can make use of the core PSI4 objects (the Wavefunction, etc), without needing to recompile all of PSI4 as you change your plugin.

The basic structure of the calculation is:

1. Generate a bunch of points in R^3 around the molecule.
   * This requires computing points with a ~uniform density on the VdW surface -- see ``src/vdwsurface.cc``
2. Calculate the QM ESP at these points from the wavefunction.
3. Calculate the MM ESP at these points as a function of the atomic partial charges.
4. Minimize some norm of the difference between these two sets of charges w.r.t the atomic partial charges.

For the minimzer, this code  uses [nlopt](http://ab-initio.mit.edu/wiki/index.php/NLopt) (cpp interface), so you'll
need the headers / libraries for that.

Installation
------------
1. Install psi4 using conda
   * `conda install -c psi4 psi4`
   * This requires having the Anaconda or Miniconda Python distribution installed.
   * See http://www.psicode.org/psi4manual/master/conda.html#quick-installation
   * You'll need psi4 >= 0.3.534.
2. You'll need to install the nlopt library and its development headers
   * On Debian-based distros, install the `libnlopt-dev` package
   * On RPM-based distros, I think it's `NLopt-devel`.
   * It's also very easy to compile from source, and in that case you'll get a static
     library by default which makes everything easier.
3. You need a C++11-capable compiler and boost.
   * If you're using a cluster with an old linux, just run
     `conda install gcc boost`
4. Checkout this repository and run `./configure; make` in this directory
5. Run `psi4` in this directory to run the example input file.

CMake Installation
------------------
1. Install psi4 using conda
   * `conda install psi4 -c psi4/label/dev -c psi4`
   * This requires having the Anaconda or Miniconda Python distribution installed.
   * See http://www.psicode.org/psi4manual/master/conda.html#how-to-install-a-psi4-binary-into-an-ana-miniconda-distribution
   * You'll need psi4 >= 1.1
2. You'll need to install the nlopt library and its development headers
   * `conda install nlopt -c conda-forge`
   * If CMake doesn't find it, pass `cmake ... -DNLOPT_PREFIX=${CONDA_PREFIX}`
3. Also, install Boost 
   * `conda install boost`
3. You need a C++11-capable compiler
   * psi4 plugin machinery should take care of this, as the g++ distributed w/psi4 works just fine
4. Checkout this repository
   * Run `psi4 --plugin-compile` in this directory.
   * Run the resulting line, adding any other options you like
   * Run `make`
5. Run `python input.py` in this directory to run the example PsiAPI input file.

Example
-------
See the [build log](https://travis-ci.org/rmcgibbo/resp2) on Travis-CI.

References
----------
[1] Bayly, Cieplak, Cornell, and Kollman (1993) http://pubs.acs.org/doi/abs/10.1021/j100142a004

License
-------
To the maximum extent possible, everything in this repository that is my (=rmcgibbo) own work
is released under CC0 (see LICENSE file).
