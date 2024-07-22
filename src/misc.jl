function create_kernel(exptype::inversion1D, x::Vector, X::Vector=exp10.(range(-5, 1, 100)))
    if exptype == IR
        kernel_eq = (t, T) -> 1 - 2 * exp(-t / T)
    elseif exptype == CPMG
        kernel_eq = (t, T) -> exp(-t / T)
    elseif exptype == PFG
        kernel_eq = (t, D) -> exp(-t / D)
    end

    return kernel_eq.(x, X')
end

function create_kernel(exptype::inversion2D, 
    x_dir::Vector, x_indir::Vector,
    X_dir::Vector=exp10.(range(-5, 1, 100)), X_indir::Vector=exp10.(range(-5, 1, 100)))

    if exptype == IRCPMG

    end

end

function gcv_score(α, r, s, x) # where r is the residuals of the solution and x=Ṽ₀'f

    ñ = length(s)
    c = r ./ α
    σ² = s .^ 2
    m̂ = sum(σ² ./ (σ² .+ α))
    dm̂ = sum(σ² ./ ((σ² .+ α) .^ 2))
    t̂ = sum(x .^ 2 ./ (σ² .+ α))
    φₙ = (α .^ 2 * c' * c .* ñ) ./ ((ñ - m̂) .^ 2)  # GCV score to be minimized
    αₙ = (α .^ 2 * c' * c * dm̂) / (t̂ * (ñ - m̂))  # α update value, test this one next

    return φₙ, αₙ
end

function Γ(m::Int, order::Int)
    # Eilers, P. H. C. (2003). Analytical Chemistry, 75(14), 3631–3636. (Supporting Information)
    # Taken from RegularizationTools.jl
    if order == 0
        return Array{Float64}(LinearAlgebra.I, (m, m))
    end
    return diff(Γ(m, order - 1), dims=1)
end

