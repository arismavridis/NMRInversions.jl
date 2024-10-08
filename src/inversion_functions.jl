"""
# Inversion for 1D pulse sequences:
    invert(seq, x, y ; lims, alpha, order, solver)

This function will build a kernel and use it to perform an inversion using the algorithm of your choice.
The output is an `inv_out_1D` structure.

 Necessary (positional) arguments:
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PGSE)
- `x` is the experiment x axis (time or b factor etc.)
- `y` is the experiment y axis (intensity of the NMR signal)

 Optional (keyword) arguments:
- `lims=(a,b,c)` will set the "limits" of the output X, so that it starts from 10^a, ends in 10^b and consists of c logarithmically spaced values (default is (-5, 1, 128)). Alternatiively, a vector of values can be used directly (e.g. `lims=[10^-5, ... ,10]` ).
- `alpha` determines the smoothing term. Use a real number for a fixed alpha.  No selection will lead to automatically determining alpha through the defeault method, which is `gcv`.
- `order` is the order of the regularization. Default is 0.
- `solver` is the algorithm used to do the inversion math. Default is `brd`.

"""
function invert(seq::Type{<:pulse_sequence1D}, x::AbstractArray, y::Vector;
    lims=(-5, 1, 128), alpha=gcv, order=0, solver=brd)

    if typeof(lims) == Tuple{Int,Int,Int}
        X = exp10.(range(lims...))
    elseif typeof(lims) == AbstractVector
        X = lims
    else
        error("lims keyword argument must be a tuple or a vector")
    end

    α = 1 #placeholder, will be replaced below 

    if isa(alpha, Real)

        α = alpha
        ker_struct = create_kernel(seq, x, X, y)

        f, r = solve_regularization(ker_struct.K, ker_struct.g, α, solver, order)

    elseif alpha == gcv

        ker_struct = create_kernel(seq, x, X, y)
        f, r, α = solve_gcv(ker_struct, solver, order)

    else
        error("alpha must be a real number or gcv")

    end

    x_fit = collect(range(0, x[end] + 0.1 * x[end], 128))
    y_fit = create_kernel(seq, x_fit, X) * f

    isreal(y) ? SNR = NaN : SNR = calc_snr(y)

    weighted_average =  f'X / sum(f)

    return inv_out_1D(seq, x, y, x_fit, y_fit, X, f, r, SNR, α, weighted_average)

end

"""
    invert(data::input1D ; kwargs...)

Instead of the positional arguments `seq`, `x` and `y`,
you can use a single `input1D` structure, which contains the same information. 
Especially useful if you're using the output of one 
of the import functions (look documentation tutorial section).
"""
function invert(data::input1D; kwargs...)

    return invert(data.seq, data.x, data.y; kwargs...)

end




"""
# Inversion for 2D pulse sequences:
    invert(seq, x_direct, x_indirect, X_direct, X_indirect, Data)

This function will build a kernel and use it to perform an inversion using the algorithm of your choice.
The output is an `inv_out_2D` structure.

 Necessary (positional) arguments:
- `seq` is the 2D pulse sequence (e.g. IRCPMG)
- `x_direct` is the direct dimension acquisition parameter (e.g. the times when you aquire CPMG echoes).
- `x_indirect` is the indirect dimension acquisition parameter (e.g. all the delay times τ in your IR sequence).
- `Data` is the 2D data matrix of complex data.


 Optional (keyword) arguments:
- `rdir` determines the output "range" of the inversion in the direct dimension (e.g. T₂ times in IRCPMG)
- `rindir` determines the output "range" of the inversion in the indirect dimension (e.g. T₁ times in IRCPMG)
 In both cases above, you can use a tuple specifying the limits of the range, or a vector of values, same as the `lims` argument in the 1D inversion.

- `alpha` determines the smoothing term. Use a real number for a fixed alpha.  No selection will lead to automatically determining alpha through the defeault method, which is `gcv`.
- `order` is the order of the regularization. Default is 0.
- `solver` is the algorithm used to do the inversion math. Default is `brd`.

"""
function invert(
    seq::Type{<:pulse_sequence2D}, x_direct::AbstractVector, x_indirect::AbstractVector, Data::AbstractMatrix;
    alpha=gcv, rdir=(-5, 1, 100), rindir=(-5, 1, 100),
    solver=brd, order=0)

    if typeof(rdir) == Tuple{Int,Int,Int}
        X_direct = exp10.(range(rdir...))
    elseif typeof(rdir) == AbstractVector
        X_direct = lims
    else
        error("rdir keyword argument must be a tuple or a vector")
    end

    if typeof(rindir) == Tuple{Int,Int,Int}
        X_indirect = exp10.(range(rindir...))
    elseif typeof(rindir) == AbstractVector
        X_indirect = lims
    else
        error("rindir keyword argument must be a tuple or a vector")
    end


    ker_struct = create_kernel(seq, x_direct, x_indirect, X_direct, X_indirect, Data)

    α = 1 #placeholder, will be replaced below

    if isa(alpha, Real)
        α = alpha
        f, r = solve_regularization(ker_struct.K, ker_struct.g, α, solver, order)

    elseif alpha == gcv
        f, r, α = solve_gcv(ker_struct, solver, order)

    else
        error("alpha must be a real number or gcv")

    end

    F = reshape(f, length(X_direct), length(X_indirect))

    results = inv_out_2D(seq, X_direct, X_indirect, F, r, calc_snr(Data), α, ones(size(F)), [])

    return results

end


"""
    invert(data::input2D ; kwargs...)

Instead of the positional arguments `seq`, `x_direct` ,
`x_indirect` and `Data`, you can use a single `input2D` 
structure, which contains the same information. 
Especially useful if you're using the output of one of the 
import functions (look documentation tutorial section).
"""
function invert(input::input2D; kwargs...)

    invert(input.seq, input.x_direct, input.x_indirect, input.data; kwargs...)

end
