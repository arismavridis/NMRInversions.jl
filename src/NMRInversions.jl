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
- fix L curve method 
- add gcv for reci method
- consider RegularizedLeastSquares.jl
- fix readresults and writeresults
- add 1D gui
"""

## The following are custom types for multiple dispatch purposes

# Pulse sequences
abstract type pulse_sequence1D end
abstract type pulse_sequence2D end
struct IR <: pulse_sequence1D end
struct CPMG <: pulse_sequence1D end
struct PFG <: pulse_sequence1D end
struct IRCPMG <: pulse_sequence2D end
struct PFGCPMG <: pulse_sequence2D end
export pulse_sequence1D, pulse_sequence2D, IR, CPMG, PFG, IRCPMG, PFGCPMG

# Supported solvers 
abstract type regularization_solver end 
struct brd <: regularization_solver end
struct brd <: regularization_solver end
struct ripqp <: regularization_solver end
struct pdhgm <: regularization_solver end
struct optim_nnls <: regularization_solver end
export regularization_solver, brd, ripqp, pdhgm, optim_nnls
export regularization_solver, brd, ripqp, pdhgm, optim_nnls

# Supported methods to determine regularization α parameter
abstract type smoothing_optimizer end
struct gcv <: smoothing_optimizer end
struct lcurve <: smoothing_optimizer end
export smoothing_optimizer, gcv, lcurve
export smoothing_optimizer, gcv, lcurve


## Include the package files 
include("inversions_io.jl")
include("kernels.jl")
include("finding_alpha.jl")
include("optim_regularizations.jl")
include("L1_regularization.jl")
include("inversion_functions.jl")
include("misc.jl")

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