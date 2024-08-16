module NMRInversions

# Import dependencies
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using NativeFileDialog
using PolygonOps
import Optimization, OptimizationOptimJL

"""
to do list:
- add L curve method 
- add gcv for reci method
"""

## The following are custom types for multiple dispatch purposes

# Pulse sequences
abstract type inversion1D end
abstract type inversion2D end
struct IR <: inversion1D end
struct CPMG <: inversion1D end
struct PFG <: inversion1D end
struct IRCPMG <: inversion2D end
struct PFGCPMG <: inversion2D end
export inversion1D, inversion2D, IR, CPMG, PFG, IRCPMG, PFGCPMG

# Supported solvers 
abstract type regularization_solver end
struct song <: regularization_solver end
struct ripqp <: regularization_solver end
struct pdhgm <: regularization_solver end
export regularization_solver, song, ripqp, pdhgm

# Supported methods to determine regularization α parameter
abstract type smoothing_optimizer end
struct gcv <: smoothing_optimizer end
struct brd <: smoothing_optimizer end
struct lcurve <: smoothing_optimizer end
export smoothing_optimizer, gcv, brd, lcurve



## Include the package files 
include("inversions_io.jl")
include("kernels.jl")
include("misc.jl")
include("inversions_1D.jl")
include("inversions_2D.jl")
include("tikhonov_song.jl")
include("L1_regularization.jl")


# Export useful functions
export invert
export create_kernel
export import_1D
export import_spinsolve, import_geospec

# The following functions are modified from extension files
function select_peaks end
function pubfig end
export select_peaks, pubfig

end