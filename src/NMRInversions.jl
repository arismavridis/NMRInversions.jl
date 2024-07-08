module NMRInversions

using DelimitedFiles
using LinearAlgebra
using SparseArrays
using Statistics
using Optimization, OptimizationOptimJL
using NativeFileDialog
using PolygonOps
using GLMakie
using RipQP
using QuadraticModels
using JuMP
using HiGHS

#These are for multiple dispatch purposes
struct inversion1D end
IR = inversion1D()
CPMG = inversion1D()
PFG = inversion1D()
export IR
export CPMG
export PFG

struct inversion2D end
IRCPMG = inversion2D()
export IRCPMG

export inversion1D
export inversion2D


# Include the files 
include("regularization.jl")
include("inversions_1D.jl")
include("inversions_2D.jl")
include("gui.jl")
include("svds.jl")
include("data_import.jl")


export invert
export svdcompress
export import_1D
export import_spinsolve
export select_peaks

end