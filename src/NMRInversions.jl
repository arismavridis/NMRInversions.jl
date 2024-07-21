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
import Ipopt
import Optimization, OptimizationOptimJL


## The following are custom types for multiple dispatch purposes

# Pulse sequences
customtypes = Dict(
    :IR => :inversion1D,
    :CPMG => :inversion1D,
    :PFG => :inversion1D,
    :IRCPMG => :inversion2D
)
struct inversion1D end
struct inversion2D end
export inversion1D, inversion2D

for (v, t) in customtypes
    eval(
        quote
            $v = $t()
            export $v
        end
    )
end


# Supported solvers (via extensions)
abstract type regularization_solver end
export regularization_solver

reg_types = Dict(
    :brd => :brd_solver,
    :ripqp => :ripqp_solver,
    :IP => :ip_solver,
    :ipoptL1 => :jump_L1_solver
    :pdhgm => :pdhgm_solver
)

for (v, t) in reg_types
    eval(
        quote
            struct $t <: regularization_solver end
            $v = $t()
            export $v
        end
    )
end


## Include the package files 
include("misc.jl")
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