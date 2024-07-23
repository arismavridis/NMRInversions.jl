module NMRInversions

# Import dependencies
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using NativeFileDialog
using PolygonOps
using GLMakie
import JuMP
import HiGHS
import Ipopt
import Optimization, OptimizationOptimJL

"""
to do list:
- Move the makie gui to extension
- rename svd function to create_kernel, and make it more user friendly
- add L curve method 
- change multiple dispatch, from ::inversion1D to ::Type{inversion1D}, and use the stuct as input instead of defining variable with struct name
abstract type A end
struct a <: A end
foo(a::Type{<:A})

"""

## The following are custom types for multiple dispatch purposes

# Pulse sequences
customtypes = Dict(
    :IR => :inversion1D,
    :CPMG => :inversion1D,
    :PFG => :inversion1D,
    :IRCPMG => :inversion2D
)

abstract type inversion1D end
abstract type inversion2D end
export inversion1D, inversion2D

for (a, A) in customtypes
    eval(
        quote
            struct $a <: $A end
            export $a
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
    :ipoptL1 => :jump_L1_solver,
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
include("inversions_io.jl")
include("kernels.jl")
include("inversions_1D.jl")
include("inversions_2D.jl")
include("gui.jl")
include("tikhonov_brd.jl")
include("tikhonov_jump.jl")


function foo(a::NMRInversions.brd_solver=brd)
    display("type is brd_solver")
end


# Export useful functions
export invert
export create_kernel
export import_1D
export import_spinsolve
export select_peaks

end