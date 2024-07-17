module NMRInversions

# Import dependencies
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using Statistics
using NativeFileDialog
using PolygonOps
using GLMakie
import JuMP
import HiGHS
import Optimization, OptimizationOptimJL


## The following are custom types for multiple dispatch purposes
# Supported 1D experiment types 
struct inversion1D end
export inversion1D

IR = inversion1D()
export IR
CPMG = inversion1D()
export CPMG
PFG = inversion1D()
export PFG

# Supported 2D experiment types
struct inversion2D end
export inversion2D

IRCPMG = inversion2D()
export IRCPMG

# Supported solvers (via extensions)
abstract type tikhonov_solver end
export tikhonov_solver

struct brd_solver <: tikhonov_solver end
brd = brd_solver()
export brd

struct ip_solver <: tikhonov_solver end
IP = ip_solver()
export IP

struct ripqp_solver <: tikhonov_solver end
ripqp = ripqp_solver()
export ripqp



## Include the package files 
include("regularization.jl")
include("inversions_1D.jl")
include("inversions_2D.jl")
include("gui.jl")
include("svds.jl")
include("data_import.jl")
include("tikhonov_brd.jl")
include("tikhonov_jump.jl")


function foo(a::NMRInversions.brd_solver=brd)
    display("type is brd_solver")
end


# Export useful functions
export invert
export svdcompress
export import_1D
export import_spinsolve
export select_peaks

end