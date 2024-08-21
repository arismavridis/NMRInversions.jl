
function invert(data::input1D; kwargs...)

    return invert(data.exptype, data.x, data.y; kwargs...)

end


function invert(exptype::Type{<:inversion1D}, x::AbstractArray, y::Vector;
    lims=(-5, 1, 128), α=1, order=0, solver=song,
    savedata=false, makeplot=false)

    if typeof(lims) == Tuple{Int,Int,Int}
        X = exp10.(range(lims...))
    elseif typeof(lims) == AbstractVector
        X = lims
    else
        error("lims must be a tuple or a vector")
    end


    if isa(α, Real)

        ker_struct = create_kernel(exptype, x, X, y)

        f, r = solve_regularization(ker_struct.K, ker_struct.g, α, solver, order)

    elseif α == gcv

        ker_struct = create_kernel(exptype, x, X, y)
        f, r, α = solve_gcv(ker_struct, solver, order)

    else
        error("α must be a real number or gcv")

    end

    if savedata == true
        open("delim_file.txt", "w") do io
            writedlm(io, ["inv_x" "inv_y"])
            writedlm(io, [T_range f])
        end
    end

    x_fit = collect(range(0, x[end] + 0.1 * x[end], 258))
    y_fit = create_kernel(exptype, x_fit, X) * f

    isreal(y) ? SNR = NaN : SNR = calc_snr(y)

    return invres1D(exptype, X, f, r, SNR, α, x_fit, y_fit)

end

