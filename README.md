# NMRInversions.jl
<p align="center">
    <img width=450 src="./logo/logo.svg"/>
</p>

This package can be used to easily perform numerical inversions for 1D and 2D NMR relaxation and diffusion measurements.

The default method is Tikhonov regularization, where the optimal smoothing term is determinded through generalized cross-validation, as described in [Mitchell et al 2012](https://doi.org/10.1016/j.pnmrs.2011.07.002).

A GUI is implemented as an extension through the use of the [GLMakie](https://github.com/MakieOrg/Makie.jl) package, which enables interactive visualisation for 2D NMR experiments.

For more details, please refer to the documentation.

If you have any problems or questions, please feel free to [submit an issue](https://github.com/arismavridis/NMRInversions.jl/issues).
