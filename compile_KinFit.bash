#!/bin/bash
c++ -O3 -ftemplate-depth=5000 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) KinFit.cpp -o KinFit$(python3-config --extension-suffix)
mv KinFit$(python3-config --extension-suffix) ${CMSSW_BASE}/lib/el9_amd64_gcc12
