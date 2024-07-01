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

# Include the files 

include("regularization.jl")
include("svds.jl")
include("inversions_1D.jl")
include("inversions_2D.jl")
include("data_import.jl")

#These are for multiple dispatch purposes
struct inversion1D end
struct inversion2D end

T1 = inversion1D()
T2 = inversion1D()
D = inversion1D()
T1T2 = inversion2D()

export T1T2map

end

