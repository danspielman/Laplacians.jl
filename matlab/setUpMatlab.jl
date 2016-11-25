#=

This sets the matlab path correctly.

=#

using MATLAB
matdir = "$(Pkg.dir("Laplacians"))/matlab"
mat"""
  addpath($matdir)
"""
  
