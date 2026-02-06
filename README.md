# C++ Wrapper of fastMTT

## Content

This package contains the following components:
* `fastmtt_cpp.cpp` - C++ wrapper of fastMTT algorithm
* `kinfit_3pr.cpp` - C++ wrapper of kinematic fit for tau(X)+tau(a1) decays
* `kinfit_3pr_3pr.cpp` - C++ wrapper of kinematic fit for tau(a1)+tau(a1) decays
* `functions.h` - collection of functions used in fastmtt and kinfit
* `KinFit.cpp` - C++ code used for compilation of 
* `compile_kinfit.bash` - compilation script (creates shared library)
* `test_kinfit.py` - testing script
* `pybind11` - folder of pybind11 package (C++ binding to python)
* `Plot.py` - plotting macro

## Getting code from git

```
cd $CMSSW_BASE/src
git clone https://github.com/raspereza/DiTauMass.git
cd $CMSSW_BASE/
scramv1 b -j 4
```

## Compiling shared library

```
cd $CMSSW_BASE/src/
./compile_kinfit.bash
```

## Running test

Example:
```
./check_kinfit.py --channel tt 
```
The script will create RooT file named `kinfit_3pr_tt.root` which you can inspect.

