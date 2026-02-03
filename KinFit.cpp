#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "functions.h"

#include "kinfit_3pr.cpp"
#include "kinfit_3pr_3pr.cpp"
#include "fastmtt_cpp.cpp"

PYBIND11_MODULE(KinFit, m) {
  m.def("fastmtt_cpp", fastmtt_cpp, "FastMTT_cpp");
  m.def("kinfit_3pr", kinfit_3pr, "KinFit_3pr");
  m.def("kinfit_3pr_3pr",kinfit_3pr_3pr, "KinFit_3pr_3pr");
}

