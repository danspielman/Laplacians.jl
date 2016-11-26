#=

This sets the matlab path correctly.

=#

#include("$(Pkg.dir("Laplacians"))/matlab/matlabSolvers.jl")

using MATLAB
matdir = "$(Pkg.dir("Laplacians"))/matlab"
mat"""
  addpath($matdir)
"""
  

