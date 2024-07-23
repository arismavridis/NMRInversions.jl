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
- add L curve method 
- add SVD to 1D kernels, implement GCV for 1D inversion

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

# Supported solvers 
abstract type regularization_solver end
export regularization_solver

for x in [:song, :ripqp, :pdhgm] 
    eval(
        quote
            struct $x <: regularization_solver end
            export $x
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
include("tikhonov_song.jl")
include("tikhonov_jump.jl")



# Export useful functions
export invert
export create_kernel
export import_1D
export import_spinsolve
export select_peaks

end