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
git clone --recurse-submodules https://github.com/raspereza/DiTauMass.git
cd $CMSSW_BASE/
scramv1 b -j 4
```

## Compiling shared library

```
cd $CMSSW_BASE/src/
./compile_kinfit.bash
```

## Running test

### Evaluating FastMTT:
Example for tt channel:
```
python3 FastMTT/examples/batch_processing.py /eos/cms/store/group/phys_tau/lrussell/forAliaksei/CPSignalStudies/Run3_2022EE/tt/GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay/nominal/merged.root --entry-stop 200000
```

this will create file merged.root in the working directory.

Then you have to run the test file and give the adress to the file. For example:
```
./test_kinfit.py merged.root --channel tt
```
The script will create RooT file named `kinfit_3pr_tt.root` which you can inspect.

