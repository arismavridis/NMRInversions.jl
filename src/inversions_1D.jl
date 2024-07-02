struct invres1D
    ps::String
    x::Vector
    f::Vector
    r::Vector
    SNR::Real
    alpha::Real
end

function invert(exptype::inversion1D, file::String=pick_file(pwd()); kwargs...)

    csvmatrix = readdlm(file, ',')
    x = csvmatrix[:, 1]
    y = csvmatrix[:, 2]

    invert(exptype, x, y; kwargs...)

end

function invert(exptype::inversion1D, x, y;
    lims=(-5, 1, 128),
    α=0, order=0, solver=:brd,
    savedata=false, saveplot=false)

    if exptype == IR 
        kernel_eq = (t, T) -> 1 - 2 * exp(-t / T)

    elseif exptype == CPMG 
        kernel_eq = (t, T) -> exp(-t / T)

    elseif exptype == PFG
        kernel_eq = (t, D) -> exp(-t / D)

    end

    T_range = exp10.(range(lims...))
    A = kernel_eq.(x, T_range')

    display("Regularization of order $(order).")

    b = y
    # Conversion to standard form, copied from RegularizationTools
    L = Γ(size(A, 2), order)
    L⁺ = pinv(L)
    p, n = size(L)
    Iₙ = Matrix{Float64}(I, n, n)
    Q, R = qr(L')
    K0 = Q[:, p+1:end]
    H, T = qr(A * K0)
    H0 = H[:, 1:n-p]
    K₀T⁻¹H₀ᵀ = K0 * T^-1 * H0'
    L⁺ₐ = (Iₙ - K₀T⁻¹H₀ᵀ * A) * L⁺
    Ā = A * L⁺ₐ
    b̄ = b - A * K₀T⁻¹H₀ᵀ * b

    # Solve
    f, r = solve_tikhonov(Ā, b̄, α, solver=solver)

    # Go back to general form
    f = L⁺ * f + K₀T⁻¹H₀ᵀ * (b - A * L⁺ * f)


    if savedata == true
        open("delim_file.txt", "w") do io
            writedlm(io, [T_range f])
        end
    end

    plotrange = range(0, x[end] + 0.1 * x[end], 258)
    fit = sum(f .* kernel_eq.(plotrange', T_range), dims=1)'

    av = sum(f .* T_range) / sum(f)

    display(plot(log10.(T_range), f))

    return T_range, f

end