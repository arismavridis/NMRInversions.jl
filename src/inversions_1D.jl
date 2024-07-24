struct invres1D
    ps::String
    x::Vector
    f::Vector
    r::Vector
    SNR::Real
    alpha::Real
end

function invert(exptype::Type{<:inversion1D}, file::String=pick_file(pwd()); kwargs...)

    csvmatrix = readdlm(file, ',')
    x = csvmatrix[:, 1]

    if size(csvmatrix, 2) == 2 # If we got only real data
        y = csvmatrix[:, 2]

    elseif size(csvmatrix, 2) == 3  # If we got both real and imaginary data

        if exptype in [IR, IRCPMG]
            positive_start = false
        else
            positive_start = true
        end

        yre, yim, ϕ = autophase(csvmatrix[:, 2], csvmatrix[:, 3], maxre1=positive_start)
        display("Data phase corrected by $(round(ϕ, digits=3)) radians.")

        y = complex.(yre, yim)
    else
        error("The 1D experiment data file must have two or three columns.")
    end

    invert(exptype, x, y; kwargs...)

end



function invert(exptype::Type{<:inversion1D}, x::AbstractArray, y::Vector;
    lims=(-5, 1, 128), α=1, order=0, solver=song,
    savedata=false, makeplot=false)

    if isa(α, Real)

        K = create_kernel(exptype, x, exp10.(range(lims...)))

        f, r = solve_regularization(K, y, α, solver, order)

    elseif α == gcv
        ker_struct = create_kernel(exptype, x, exp10.(range(lims...)), y)

        f, r = solve_gcv(ker_struct, solver, order)

    end

    if savedata == true
        open("delim_file.txt", "w") do io
            writedlm(io, ["inv_x" "inv_y"])
            writedlm(io, [T_range f])
        end
    end

    if makeplot == true
        plotrange = collect(range(0, x[end] + 0.1 * x[end], 258))
        fit = create_kernel(exptype, plotrange, exp10.(range(lims...))) * f

    end

    return f, r

end


## Conversion to standard form, copied from RegularizationTools
# L = Γ(size(A, 2), order)
# L⁺ = pinv(L)
# p, n = size(L)
# Iₙ = Matrix{Float64}(I, n, n)
# Q, R = qr(L')
# K0 = Q[:, p+1:end]
# H, T = qr(A * K0)
# H0 = H[:, 1:n-p]
# K₀T⁻¹H₀ᵀ = K0 * T^-1 * H0'
# L⁺ₐ = (Iₙ - K₀T⁻¹H₀ᵀ * A) * L⁺
# Ā = A * L⁺ₐ
# b̄ = b - A * K₀T⁻¹H₀ᵀ * b
# Solve
# f, r = solve_regularization(Ā, b̄, α, solver=solver)
# Go back to general form
# f = L⁺ * f + K₀T⁻¹H₀ᵀ * (b - A * L⁺ * f)